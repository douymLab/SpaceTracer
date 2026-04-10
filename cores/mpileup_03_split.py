#!/usr/bin/env python3
"""
Split Mpileup Step - Split mpileup file when the file is so large 
(this will not work when you use the one-step "run" function, but when you run separately)
"""

import os
from pathlib import Path
import sqlite3
from typing import Dict
import json

from collections import defaultdict
from SpaceTracer.utils.logger import get_logger
from SpaceTracer.utils.utils import check_dir
model_name=__name__
logger = get_logger(model_name)

class SplitMpileupStep():
    """
    split filtered_mpileup file
    
    1. check filtered_mpileup file size
    2. if file length >= threshold, split file
    3. split by chrom[chrM will be specific handle]
    
    Parameters:
        input_file: the input mpileup file
        output_dir: the output dir
        split_manifest: the manifest file recoding split information
        genome_details: the genome information
        chrom_chunk_size: the chunk size for auto chromosome (default:5000)
        chrM_chunk_size: the chunk size for mitochondrial (default:100)
    """
    def __init__(self, input_file: str, output_dir: str, db_path:str,genome_details: Dict, config: dict):
        self.input_file = input_file
        self.output_dir =output_dir
        self.db_path=db_path
        self.genome_details=genome_details
        self.config=config


    def get_step_config(self) -> Dict:
        return self.config.get('steps', {}).get('mpileup', {})
    
    def _run(self, chrom_chunk_size, chrM_chunk_size, read_len, max_cost):
        """Split mpileup and write chunk index db directly."""
        try:
            split_finish = self._split_by_chromosome_and_chunk(
                chrom_chunk_size=chrom_chunk_size,
                chrM_chunk_size=chrM_chunk_size,
                read_len=read_len,
                max_cost=max_cost
            )

            return True

        except Exception:
            raise
        
    def _init_chunk_index_db(self, db_path: str):
        if os.path.exists(db_path):
            os.remove(db_path)

        conn = sqlite3.connect(db_path)
        cur = conn.cursor()

        cur.execute("""
        CREATE TABLE chunks (
            chunk_id    TEXT PRIMARY KEY,
            chrom       TEXT NOT NULL,
            chrom_type  TEXT,
            chunk_idx   INTEGER NOT NULL,
            chunk_file  TEXT NOT NULL,
            start_pos   INTEGER,
            end_pos     INTEGER,
            records     INTEGER,
            span_bp     INTEGER,
            max_depth   INTEGER,
            mean_depth  REAL,
            cost        REAL
        )
        """)

        cur.execute("""
        CREATE TABLE chroms (
            chrom          TEXT PRIMARY KEY,
            chrom_type     TEXT,
            total_records  INTEGER,
            num_chunks     INTEGER,
            chunk_size     INTEGER,
            total_cost     REAL
        )
        """)

        cur.execute("CREATE INDEX idx_chunks_chrom ON chunks(chrom)")
        cur.execute("CREATE INDEX idx_chunks_chrom_idx ON chunks(chrom, chunk_idx)")

        conn.commit()
        return conn

    def _split_by_chromosome_and_chunk(self, chrom_chunk_size, chrM_chunk_size, read_len, max_cost=10000):
        mpileup_file = self.input_file
        genome_info = self.genome_details
        db_path=self.db_path

        split_dir = os.path.join(self.output_dir, "split_for_candidates")
        check_dir(split_dir)

        conn = self._init_chunk_index_db(db_path)
        cur = conn.cursor()

        chrom_config = genome_info["chromosomes"]
        autosomes = set(chrom_config["autosomes"])
        sex_chromosomes = set(chrom_config["sex_chromosomes"])
        mitochondrial = set(chrom_config["mitochondrial"])
        contigs = set(chrom_config["contigs"])

        handlers = {}
        # split_files = []

        def infer_chrom_type_and_chunk_size(chrom: str):
            if chrom in mitochondrial:
                return "mitochondrial", int(chrM_chunk_size)
            elif chrom in sex_chromosomes:
                return "sex_chromosome", int(chrom_chunk_size)
            elif chrom in autosomes:
                return "autosome", int(chrom_chunk_size)
            elif chrom in contigs:
                return "contig", int(chrom_chunk_size)
            else:
                logger.warning(f"Chromosome {chrom} not found in genome_info, fallback to unknown")
                return "unknown", int(chrom_chunk_size)

        def get_handler(chrom: str):
            if chrom not in handlers:
                chrom_type, chunk_size = infer_chrom_type_and_chunk_size(chrom)
                chunk_idx = 0
                file_path = os.path.join(split_dir, f"{chrom}_chunk{chunk_idx:04d}.mpileup")

                handlers[chrom] = {
                    "chrom": chrom,
                    "chrom_type": chrom_type,
                    "chunk_size": chunk_size,
                    "chunk_idx": chunk_idx,
                    "file_obj": open(file_path, "w"),
                    "current_file": str(file_path),

                    # current chunk runtime state
                    "count": 0,
                    "chunk_start_pos": None,
                    "last_pos": None,
                    "max_depth": 0,
                    "depth_sum": 0,

                    # chrom-level summary
                    "total_records": 0,
                    "num_chunks": 0,
                    "total_cost": 0.0,
                }

                # split_files.append(str(file_path))

            return handlers[chrom]

        def reset_chunk_state(handler):
            handler["count"] = 0
            handler["chunk_start_pos"] = None
            handler["last_pos"] = None
            handler["max_depth"] = 0
            handler["depth_sum"] = 0

        def open_next_chunk_file(handler):
            handler["chunk_idx"] += 1
            file_path = os.path.join(
                split_dir,
                f"{handler['chrom']}_chunk{handler['chunk_idx']:04d}.mpileup"
            )
            handler["file_obj"] = open(file_path, "w")
            handler["current_file"] = str(file_path)
            # split_files.append(str(file_path))

        def finalize_chunk(handler):
            if handler["count"] == 0:
                return

            start_pos = handler["chunk_start_pos"]
            end_pos = handler["last_pos"]
            span_bp = end_pos - start_pos + 1 if start_pos is not None and end_pos is not None else 0
            mean_depth = handler["depth_sum"] / handler["count"] if handler["count"] > 0 else 0.0
            cost = handler["max_depth"] * (span_bp / read_len) if read_len > 0 else 0.0
            chunk_id = f"{handler['chrom']}_chunk{handler['chunk_idx']:04d}"

            cur.execute("""
            INSERT INTO chunks (
                chunk_id, chrom, chrom_type, chunk_idx, chunk_file,
                start_pos, end_pos, records, span_bp, max_depth, mean_depth, cost
            ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            """, (
                chunk_id,
                handler["chrom"],
                handler["chrom_type"],
                handler["chunk_idx"],
                handler["current_file"],
                start_pos,
                end_pos,
                handler["count"],
                span_bp,
                handler["max_depth"],
                mean_depth,
                cost,
            ))

            handler["total_records"] += handler["count"]
            handler["num_chunks"] += 1
            handler["total_cost"] += cost

            reset_chunk_state(handler)

        with open(mpileup_file, "r") as f:
            for line_num, line in enumerate(f, 1):
                if line.startswith("#"):
                    continue

                fields = line.rstrip("\n").split("\t")
                if len(fields) < 9:
                    continue

                chrom = fields[0]
                pos = int(fields[1])
                depth = int(fields[5])

                handler = get_handler(chrom)

                if handler["chunk_start_pos"] is None:
                    next_count = 1
                    next_span = 1
                    next_max_depth = depth
                else:
                    next_count = handler["count"] + 1
                    next_span = pos - handler["chunk_start_pos"] + 1
                    next_max_depth = max(handler["max_depth"], depth)

                next_cost = next_max_depth * (next_span / read_len) if read_len > 0 else 0.0

                # need to split before writing current line into current chunk
                if handler["count"] > 0 and (
                    next_count > handler["chunk_size"] or next_cost > max_cost
                ):
                    finalize_chunk(handler)
                    handler["file_obj"].close()
                    open_next_chunk_file(handler)

                    # current line becomes the first line of the new chunk
                    next_count = 1
                    next_max_depth = depth

                handler["file_obj"].write(line)

                if handler["chunk_start_pos"] is None:
                    handler["chunk_start_pos"] = pos

                handler["count"] = next_count
                handler["last_pos"] = pos
                handler["max_depth"] = next_max_depth
                handler["depth_sum"] += depth

        # finalize remaining chunks and chromosome summaries
        for chrom, handler in handlers.items():
            if handler["count"] > 0:
                finalize_chunk(handler)

            if handler["file_obj"]:
                handler["file_obj"].close()

            cur.execute("""
            INSERT INTO chroms (
                chrom, chrom_type, total_records, num_chunks, chunk_size, total_cost
            ) VALUES (?, ?, ?, ?, ?, ?)
            """, (
                handler["chrom"],
                handler["chrom_type"],
                handler["total_records"],
                handler["num_chunks"],
                handler["chunk_size"],
                handler["total_cost"],
            ))

        conn.commit()
        conn.close()

        return  True


# def save_manifest(manifest: Dict, output_manifest: Path):
#     """save manifest into JSON file"""
#     with open(output_manifest, 'w') as f:
#         json.dump(manifest, f, indent=2)
    
#     logger.info(f"Manifest saved to: {output_manifest}")
    

def finalize_chunk(handler,read_len):
    if handler['count'] == 0 or handler['chunk_start_pos'] is None or handler['last_pos'] is None:
        return

    span_bp = handler['last_pos'] - handler['chunk_start_pos'] + 1
    mean_depth = handler['depth_sum'] / handler['count'] if handler['count'] > 0 else 0.0
    cost = mean_depth * max(1.0, span_bp / read_len)

    chunk_meta = {
        'chunk_idx': handler['chunk_idx'],
        'file': handler['files'][-1],
        'records': handler['count'],
        'start_pos': handler['chunk_start_pos'],
        'end_pos': handler['last_pos'],
        'span_bp': span_bp,
        'max_depth': handler['max_depth'],
        'mean_depth': mean_depth,
        'cost': cost,
    }
    handler['chunks'].append(chunk_meta)

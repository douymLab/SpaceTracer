#!/usr/bin/env python3
"""
Split Mpileup Step - Split mpileup file when the file is so large
(this will not work when you use the one-step "run" function, but when you run separately)
"""

from pathlib import Path
from typing import Dict
import csv

from SpaceTracer.utils.logger import get_logger

model_name = __name__
logger = get_logger(model_name)


class SplitMpileupStep:
    """
    Split filtered mpileup file.

    1. check filtered_mpileup file size
    2. if file length >= threshold, split file
    3. split by chromosome; chrM handled with specific chunk size

    Parameters:
        input_file: the input mpileup file
        output_dir: the output dir
        manifest_path: the output TSV manifest file
        genome_details: the genome information
        config: pipeline config
    """

    def __init__(
        self,
        input_file: str,
        output_dir: str,
        manifest_path: str,
        genome_details: Dict,
        config: dict
    ):
        self.input_file = input_file
        self.output_dir = output_dir
        self.manifest_path = manifest_path
        self.genome_details = genome_details
        self.config = config

    def get_step_config(self) -> Dict:
        return self.config.get("steps", {}).get("mpileup", {})

    def _run(self, chrom_chunk_size, chrM_chunk_size, read_len, max_cost):
        """
        Split mpileup and write chunk manifest TSV directly.
        """
        try:
            self._split_by_chromosome_and_chunk(
                chrom_chunk_size=chrom_chunk_size,
                chrM_chunk_size=chrM_chunk_size,
                read_len=read_len,
                max_cost=max_cost
            )
            return True
        except Exception:
            raise

    def _split_by_chromosome_and_chunk(
        self,
        chrom_chunk_size,
        chrM_chunk_size,
        read_len,
        max_cost=10000
    ):
        """
        Scan original mpileup once and write UTF-8 TSV manifest.

        Chunking rules:
          1) chromosome switch => finalize current chunk
          2) next_count > chunk_size => finalize current chunk
          3) next_cost > max_cost => finalize current chunk

        cost = depth_sum / read_len
        """

        mpileup_file = str(Path(self.input_file).resolve())
        manifest_path = str(Path(self.manifest_path).resolve())
        genome_info = self.genome_details

        if not Path(mpileup_file).exists():
            raise FileNotFoundError(f"Input mpileup file not found: {mpileup_file}")

        chrom_config = genome_info.get("chromosomes", {})
        autosomes = set(chrom_config.get("autosomes", []))
        sex_chromosomes = set(chrom_config.get("sex_chromosomes", []))
        mitochondrial = set(chrom_config.get("mitochondrial", []))
        contigs = set(chrom_config.get("contigs", []))

        manifest_fields = [
            "chunk_id",
            "chrom",
            "chrom_type",
            "chunk_idx",
            "source_file",
            "start_offset",
            "end_offset",
            "start_pos",
            "end_pos",
            "records",
            "span_bp",
            "max_depth",
            "mean_depth",
            "cost",
        ]

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
                logger.warning(
                    f"Chromosome {chrom} not found in genome_info, fallback to unknown"
                )
                return "unknown", int(chrom_chunk_size)

        chrom_summaries = {}

        def ensure_chrom_summary(chrom: str):
            if chrom not in chrom_summaries:
                chrom_type, chunk_size = infer_chrom_type_and_chunk_size(chrom)
                chrom_summaries[chrom] = {
                    "chrom": chrom,
                    "chrom_type": chrom_type,
                    "chunk_size": chunk_size,
                    "total_records": 0,
                    "num_chunks": 0,
                    "total_cost": 0.0,
                }
            return chrom_summaries[chrom]

        def new_chunk_state(chrom: str, pos: int, offset: int, depth: int):
            chrom_summary = ensure_chrom_summary(chrom)
            chunk_idx = chrom_summary["num_chunks"]
            return {
                "chrom": chrom,
                "chrom_type": chrom_summary["chrom_type"],
                "chunk_size": chrom_summary["chunk_size"],
                "chunk_idx": chunk_idx,
                "start_offset": offset,
                "start_pos": pos,
                "last_pos": pos,
                "count": 1,
                "max_depth": depth,
                "depth_sum": depth,
            }

        def finalize_chunk(writer, chunk_state, end_offset_exclusive: int):
            if chunk_state is None or chunk_state["count"] == 0:
                return

            start_pos = chunk_state["start_pos"]
            end_pos = chunk_state["last_pos"]
            records = chunk_state["count"]
            span_bp = end_pos - start_pos + 1 if start_pos is not None and end_pos is not None else 0
            mean_depth = chunk_state["depth_sum"] / records if records > 0 else 0.0
            cost = chunk_state["depth_sum"] / read_len if read_len > 0 else 0.0
            chunk_id = f"{chunk_state['chrom']}_chunk{chunk_state['chunk_idx']:04d}"

            row = {
                "chunk_id": chunk_id,
                "chrom": chunk_state["chrom"],
                "chrom_type": chunk_state["chrom_type"],
                "chunk_idx": int(chunk_state["chunk_idx"]),
                "source_file": mpileup_file,
                "start_offset": int(chunk_state["start_offset"]),
                "end_offset": int(end_offset_exclusive),
                "start_pos": int(start_pos),
                "end_pos": int(end_pos),
                "records": int(records),
                "span_bp": int(span_bp),
                "max_depth": int(chunk_state["max_depth"]),
                "mean_depth": round(mean_depth, 6),
                "cost": round(cost, 6),
            }
            writer.writerow(row)

            chrom_summary = ensure_chrom_summary(chunk_state["chrom"])
            chrom_summary["total_records"] += records
            chrom_summary["num_chunks"] += 1
            chrom_summary["total_cost"] += cost

        current_chunk = None
        current_chrom = None

        Path(manifest_path).parent.mkdir(parents=True, exist_ok=True)

        logger.info(f"Start splitting mpileup: {mpileup_file}")
        logger.info(f"Manifest output: {manifest_path}")
        logger.info(
            f"Chunk params: chrom_chunk_size={chrom_chunk_size}, "
            f"chrM_chunk_size={chrM_chunk_size}, read_len={read_len}, max_cost={max_cost}"
        )

        # manifest 强制 UTF-8 TSV
        with open(manifest_path, "w", encoding="utf-8", newline="") as mf:
            writer = csv.DictWriter(mf, fieldnames=manifest_fields, delimiter="\t")
            writer.writeheader()

            # mpileup 用二进制读取，保证 tell()/offset 精确
            with open(mpileup_file, "rb") as f:
                while True:
                    line_start = f.tell()
                    line = f.readline()
                    if not line:
                        break

                    if line.startswith(b"#"):
                        continue

                    try:
                        text = line.decode("utf-8").rstrip("\n")
                    except UnicodeDecodeError:
                        logger.warning(
                            f"Skip undecodable mpileup line at byte offset {line_start}"
                        )
                        continue

                    fields = text.split("\t")
                    if len(fields) < 9:
                        continue

                    chrom = fields[0]

                    try:
                        pos = int(fields[1])
                        depth = int(fields[5])
                    except ValueError:
                        continue

                    # 染色体切换：先收尾旧 chunk，再开启新 chunk
                    if current_chrom is not None and chrom != current_chrom:
                        finalize_chunk(writer, current_chunk, line_start)
                        current_chunk = None
                        current_chrom = None

                    if current_chunk is None:
                        current_chunk = new_chunk_state(chrom, pos, line_start, depth)
                        current_chrom = chrom
                        continue

                    next_count = current_chunk["count"] + 1
                    next_depth_sum = current_chunk["depth_sum"] + depth
                    next_cost = next_depth_sum / read_len if read_len > 0 else 0.0

                    if (
                        next_count > current_chunk["chunk_size"]
                        or next_cost > max_cost
                    ):
                        finalize_chunk(writer, current_chunk, line_start)
                        current_chunk = new_chunk_state(chrom, pos, line_start, depth)
                        current_chrom = chrom
                        continue

                    current_chunk["count"] += 1
                    current_chunk["depth_sum"] += depth
                    current_chunk["max_depth"] = max(current_chunk["max_depth"], depth)
                    current_chunk["last_pos"] = pos
                    current_chrom = chrom

                file_end = f.tell()
                finalize_chunk(writer, current_chunk, file_end)

        total_chunks = sum(summary["num_chunks"] for summary in chrom_summaries.values())
        total_records = sum(summary["total_records"] for summary in chrom_summaries.values())

        logger.info(
            f"Split finished. Total chromosomes={len(chrom_summaries)}, "
            f"total_chunks={total_chunks}, total_records={total_records}"
        )

        for chrom, summary in chrom_summaries.items():
            logger.info(
                f"[{chrom}] type={summary['chrom_type']} "
                f"chunks={summary['num_chunks']} "
                f"records={summary['total_records']} "
                f"total_cost={summary['total_cost']:.4f}"
            )

        return True

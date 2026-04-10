#!/usr/bin/env python3
"""
UMI combine
"""

from functools import partial
import json
import time
import pyarrow as pa
import pyarrow.parquet as pq
import multiprocessing
import os
from typing import List, Tuple, Dict, Union
from concurrent.futures import ProcessPoolExecutor, as_completed
import csv

import pandas as pd
from SpaceTracer.steps.base import BaseStep
from SpaceTracer.utils.logger import get_logger
from SpaceTracer.utils.read_files import load_chunk_files_from_db, read_bam

from SpaceTracer.cores.UMI_combine import process_single_region_for_umi_combine
from SpaceTracer.utils.utils import build_region_tasks_for_UMI_combine

model_name=__name__
logger = get_logger(model_name)

# ── Column definitions ────────────────────────────────────────────────────────

# Matches the barcode_list structure built inside the worker function:
# [chrom, pos, ".", ref, alt, barcode, consensus_read_count_str,
#  qual_A, qual_T, qual_C, qual_G]
COLUMNS_MAIN = [
    "#chrom", "pos", "strand", "ref", "alt",
    "barcode", "consensus_read_count", "qA", "qT", "qC", "qG",
]

# error_allele file: one row per allele observation at each position
COLUMNS_ERROR_ALLELE = ["#chrom", "pos", "ref", "alt", "strand"]

# ── Step class ────────────────────────────────────────────────────────────────
class UMICombineStep(BaseStep):
    def get_inputs(self, context: Dict) -> Dict[str, str]:
        """input"""
        inputs = {
            'in_filter_bam': context.get('in_filter_bam'),
            'reference': self.config.get('genome_fasta'),
            'filter_mpileup_file':context.get('filter_mpileup_file'),
            'db_path': context.get('db_path','')
        }

        return inputs
    
    def get_outputs(self,context: Dict) -> Dict[str, str]:
        """output"""
        # return {
        #     'spot_count_file': os.path.join(self.step_dir, 'spot.count.parquet'),
        #     # 'spot_count_parquet': os.path.join(self.step_dir, 'spot.count.parquet'),
        #     'error_count_file': os.path.join(self.step_dir, 'error.count.parquet')
        # }
        return {
            'spot_count_file': os.path.join(self.step_dir, 'spot.count.list.csv'),
            # 'spot_count_parquet': os.path.join(self.step_dir, 'spot.count.parquet'),
            'error_count_file': os.path.join(self.step_dir, 'error.count.list.csv')
        }
    
    def _run(self,context: Dict):
        inputs=self.get_inputs(context)
        bam_file=inputs["in_filter_bam"]
        mpileup_file=inputs["filter_mpileup_file"]
        db_path=inputs['db_path']

        
        seq_type=self.config.get('sequence_type')
        thread=self.threads

        outputs=self.get_outputs(context)
        spot_count_file=outputs['spot_count_file']
        # spot_count_parquet=outputs['spot_count_parquet']
        error_count_file=outputs['error_count_file']

        # df = pd.read_csv(mpileup_file, sep="\t")
        
        max_region_size=50000
        max_variants_per_region=5000
        chunksize=5
        max_workers=self.threads
        thread_per_chunk=2
        main_flush_rows=100000
        err_flush_rows=5000
        compression="snappy"

        chunk_output_dir = os.path.join(self.step_dir, "chunk_outputs")
        final_output_dir = os.path.join(self.step_dir, "final_outputs")

        os.makedirs(chunk_output_dir, exist_ok=True)
        os.makedirs(final_output_dir, exist_ok=True)
        t0 = time.time()

        spot_files, error_files = _run_umi_combine_parallel(
            db_path=db_path,
            bam_file=bam_file,
            seq_type=seq_type,
            output_dir=chunk_output_dir,
            max_workers=max_workers,
            thread_per_chunk=thread_per_chunk,
            main_flush_rows=main_flush_rows,
            err_flush_rows=err_flush_rows,
            max_region_size=max_region_size,
            max_variants_per_region=max_variants_per_region,
            chunksize=chunksize,
            compression=compression
        )
        logger.info(f"[parallel umi_combine] {time.time()-t0:.2f}s")
        final_spot = spot_count_file
        final_error = error_count_file
        t1 = time.time()
        save_parquet_file_list(spot_files, final_spot)
        # merge_parquet_files(spot_files, final_spot, compression=compression)
        
        logger.info(f"[final spot file merge] {time.time()-t1:.2f}s")
        t2 = time.time()
        save_parquet_file_list(error_files, final_error)
        # merge_parquet_files(error_files, final_error, compression=compression)
        logger.info(f"[final error file merge] {time.time()-t2:.2f}s")

      

def _build_region_tasks_from_one_chunk_file_for_UMI_combine(chunk_file,max_region_size,max_variants_per_region):
    df=pd.read_csv(chunk_file,sep="\t",header=None,names=["#chrom", "pos", "type", "ref", "alt1", "total_depth", "ref_depth", "alt1_depth", "vaf"],comment='#')
    region_tasks=build_region_tasks_for_UMI_combine(df,max_region_size,max_variants_per_region)
    return region_tasks


def _umi_combine_to_parquet_buffered(
        bam_file,
        seq_type,
        region_tasks,
        thread,
        spot_count_parquet,
        error_count_parquet,
        main_flush_rows,
        err_flush_rows,
        chunksize,
        compression
    ):  
    ## write 
    main_writer = None
    err_writer = None

    main_buffer = []
    err_buffer = []

    try:
        with multiprocessing.Pool(
            thread,
            initializer=_init_worker,
            initargs=(bam_file,)
        ) as pool:

            partial_func = partial(_worker_wrapper, seq_type)

            for result in pool.imap_unordered(partial_func, region_tasks, chunksize=chunksize):
                if result is None:
                    continue

                for spot_count_list, error_list in result:
                    if spot_count_list:
                        main_buffer.extend(spot_count_list)

                    if error_list:
                        err_buffer.append(error_list)

                if len(main_buffer) >= main_flush_rows:
                    main_writer = _flush_rows_to_parquet(
                        rows=main_buffer,
                        columns=COLUMNS_MAIN,
                        writer=main_writer,
                        out_file=spot_count_parquet,
                        compression=compression
                    )
                    main_buffer.clear()

                if len(err_buffer) >= err_flush_rows:
                    err_writer = _flush_rows_to_parquet(
                        rows=err_buffer,
                        columns=COLUMNS_ERROR_ALLELE,
                        writer=err_writer,
                        out_file=error_count_parquet,
                        compression=compression
                    )
                    err_buffer.clear()

        if main_buffer:
            main_writer = _flush_rows_to_parquet(
                rows=main_buffer,
                columns=COLUMNS_MAIN,
                writer=main_writer,
                out_file=spot_count_parquet,
                compression=compression
            )

        if err_buffer:
            err_writer = _flush_rows_to_parquet(
                rows=err_buffer,
                columns=COLUMNS_ERROR_ALLELE,
                writer=err_writer,
                out_file=error_count_parquet,
                compression=compression
            )

    finally:
        if main_writer is not None:
            main_writer.close()
        if err_writer is not None:
            err_writer.close()


    # bam_handle = read_bam(bam_file)

    # for region_task in region_tasks:
    #     result = process_single_region_for_umi_combine(
    #         bam_handle, region_task, seq_type
    #     )

    #     if result is None:
    #         continue

    #     for spot_count_list, error_list in result:
    #         if spot_count_list:
    #             main_buffer.extend(spot_count_list)
    #         if error_list:
    #             err_buffer.append(error_list)

    #     if len(main_buffer) >= main_flush_rows:
    #         main_writer = _flush_rows_to_parquet(
    #             rows=main_buffer,
    #             columns=COLUMNS_MAIN,
    #             writer=main_writer,
    #             out_file=spot_count_parquet,
    #             compression=compression
    #         )
    #         main_buffer.clear()

    #     if len(err_buffer) >= err_flush_rows:
    #         err_writer = _flush_rows_to_parquet(
    #             rows=err_buffer,
    #             columns=COLUMNS_ERROR_ALLELE,
    #             writer=err_writer,
    #             out_file=error_count_parquet,
    #             compression=compression
    #         )
    #         err_buffer.clear()


def _run_umi_combine_one_chunk(
    chunk_file,
    bam_file,
    seq_type,
    thread_per_chunk,
    output_dir,
    max_region_size,
    max_variants_per_region,
    main_flush_rows,
    err_flush_rows,
    chunksize,
    compression
):
    os.makedirs(output_dir, exist_ok=True)
    chunk_base = os.path.basename(chunk_file)
    chunk_prefix = os.path.splitext(chunk_base)[0]

    spot_count_parquet = os.path.join(output_dir, f"{chunk_prefix}.spot.parquet")
    error_count_parquet = os.path.join(output_dir, f"{chunk_prefix}.error.parquet")
    
    t0 = time.time()
    # logger.info(f"chunk start: {chunk_file}")

    region_tasks = _build_region_tasks_from_one_chunk_file_for_UMI_combine(chunk_file,max_region_size,max_variants_per_region)

    t1 = time.time()
    # logger.info(f"chunk={chunk_file} build_region_tasks done, n={len(region_tasks)}, elapsed={t1-t0:.2f}s")

    if not region_tasks:
        return {
            "chunk_file": chunk_file,
            "spot_parquet": None,
            "error_parquet": None
        }

    _umi_combine_to_parquet_buffered(
        bam_file=bam_file,
        seq_type=seq_type,
        region_tasks=region_tasks,
        thread=thread_per_chunk,
        spot_count_parquet=spot_count_parquet,
        error_count_parquet=error_count_parquet,
        main_flush_rows=main_flush_rows,
        err_flush_rows=err_flush_rows,
        chunksize=chunksize,
        compression=compression
    )

    t2 = time.time()
    # logger.info(f"chunk={chunk_file} processing done, elapsed={t2-t1:.2f}s total={t2-t0:.2f}s")

    return {
        "chunk_file": chunk_file,
        "spot_parquet": spot_count_parquet,
        "error_parquet": error_count_parquet,
    }


def _run_umi_combine_parallel(
    db_path,
    bam_file,
    seq_type,
    output_dir,
    max_workers=4,
    thread_per_chunk=1,
    main_flush_rows=100000,
    err_flush_rows=50000,
    max_region_size=20000,
    max_variants_per_region=100,
    chunksize=20,
    compression="snappy"
):
    chunk_files = load_chunk_files_from_db(db_path)

    if not chunk_files:
        print("[warn] no chunk files found in db")
        return [], []

    # print(f"[info] total chunk files = {len(chunk_files)}")

    all_spot_parquets = []
    all_error_parquets = []

    with ProcessPoolExecutor(max_workers=max_workers) as exe:
        futures = []

        for chunk_file in chunk_files:
            futures.append(
                exe.submit(
                    _run_umi_combine_one_chunk,
                    chunk_file,
                    bam_file,
                    seq_type,
                    thread_per_chunk,
                    output_dir,
                    max_region_size,
                    max_variants_per_region,
                    main_flush_rows,
                    err_flush_rows,
                    chunksize,
                    compression
                )
            )

        for f in as_completed(futures):
            result = f.result()

            if result["spot_parquet"] is not None:
                all_spot_parquets.append(result["spot_parquet"])

            if result["error_parquet"] is not None:
                all_error_parquets.append(result["error_parquet"])

    return all_spot_parquets, all_error_parquets

    
    
# ── Module-level worker (must NOT be inside the class) ───────────────────────
_bam_handle = None
def _init_worker(bam_file):
    global _bam_handle
    _bam_handle = read_bam(bam_file)
    

def _worker_wrapper(seq_type, region_info):
    return process_single_region_for_umi_combine(_bam_handle,region_info,seq_type)

            
def _flush_rows_to_parquet(rows, columns, writer, out_file, compression="snappy"):

    if not rows:
        return writer

    df = pd.DataFrame(rows, columns=columns)
    df = df.fillna("NA")

    table = pa.Table.from_pandas(df, preserve_index=False)

    if writer is None:
        writer = pq.ParquetWriter(
            out_file,
            table.schema,
            compression=compression
        )

    writer.write_table(table)
    return writer
    
def save_parquet_file_list(chunk_files: List[Union[str, Tuple[str, str]]], output_file: str):
    """
    Args:
        chunk_files:
            支持两种格式：
            1. ["/path/chr10_chunk0000.spot.parquet", ...]
            2. [("chr10_chunk0000", "/path/chr10_chunk0000.spot.parquet"), ...]
        output_file:
            输出 manifest tsv 文件
    """
    valid_rows = []

    for item in chunk_files:
        chunk = None
        file = None

        # 情况1：tuple/list -> (chunk, file)
        if isinstance(item, (tuple, list)) and len(item) == 2:
            chunk, file = item

        # 情况2：str -> 从文件名解析 chunk
        elif isinstance(item, str):
            file = item
            basename = os.path.basename(file)

            # 例如: chr10_chunk0000.spot.parquet -> chr10_chunk0000
            if basename.endswith(".spot.parquet"):
                chunk = basename[:-len(".spot.parquet")]
            elif basename.endswith(".parquet"):
                chunk = basename[:-len(".parquet")]
            else:
                chunk = os.path.splitext(basename)[0]

        else:
            continue

        if not file or not isinstance(file, str) or file.strip() == "":
            continue
        if not os.path.exists(file):
            continue
        if not chunk:
            continue

        valid_rows.append((chunk, file))

    out_dir = os.path.dirname(output_file)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    with open(output_file, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f, delimiter="	")
        writer.writerow(["chunk", "parquet_file"])
        writer.writerows(valid_rows)

    if not valid_rows:
        print(f"[manifest] no valid parquet files, wrote empty manifest: {output_file}")


def merge_parquet_files(parquet_files, output_file, compression="snappy"):
    parquet_files = [f for f in parquet_files if f is not None]

    if not parquet_files:
        print(f"[merge] no parquet files to merge for {output_file}")
        return

    writer = None
    try:
        for file in parquet_files:
            if not file or (isinstance(file, str) and file.strip() == ""):
                # logger.warning(f"[{model_name}]-[merge] skip empty file path: {file}")
                continue

            if not os.path.exists(file):
                # logger.warning(f"[{model_name}]-[merge] file not found, skip: {file}")
                continue

            table = pq.read_table(file)

            if writer is None:
                writer = pq.ParquetWriter(
                    output_file,
                    table.schema,
                    compression=compression
                )

            writer.write_table(table)

    finally:
        if writer is not None:
            writer.close()
# def _process_single_region_for_umi_combine(bam_handle,region_info,seq_type):
#     """
#     sites:
#         [(pos, ref, alt, check_mosaic, check_error), ...]
#     """

#     chrom,start,end,sites = region_info
#     reads=bam_handle.fetch(chrom,start,end)

#     results=[]
#     for site in sites:
#         pos, ref, _, check_mosaic, check_error=site

#         try:
#             mosaic_spot_list, error_list = combine_UMI_spot_for_both_mosaic_and_error(
#                 reads,
#                 check_mosaic,   # "" when this position is not candidate_mosaic
#                 check_error,    # "" when this position is not candidate_error
#                 seq_type,
#                 (chrom,pos,ref),
#             )
#             # new_list  : list of rows (multiple barcodes per position)
#             # error_allele : single char, e.g. "A"  — or "" if not applicable
#             # strand       : single char, e.g. "+"  — or "" if not applicable
#             results.append([ mosaic_spot_list, error_list])
        
#         except Exception as exc:
#             identifier=(chrom,pos,ref)
#             logger.warning(f"[UmiCombine] worker failed for {identifier}: {exc}")
#             results.append([ None, None])

#     return results


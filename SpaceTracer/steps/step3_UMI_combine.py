#!/usr/bin/env python3
"""
UMI combine
"""

import io
import time
import os
import gc
import csv
from functools import partial
from typing import List, Tuple, Dict, Union

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

from SpaceTracer.steps.base import BaseStep
from SpaceTracer.utils.logger import get_logger
from SpaceTracer.utils.parallel import compute_events_threshold, parallel_map
from SpaceTracer.utils.read_files import (
    load_manifest_files_for_chunk_tasks_only_umi_combine,
    read_bam,
    read_chunk_bytes_by_offset,
)
from SpaceTracer.cores.UMI_combine import process_single_region_for_umi_combine_auto
from SpaceTracer.utils.utils import build_region_tasks_for_UMI_combine

model_name = __name__
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

# the configuration for per chunck
small_cfg = {
    "task_type": "small",
    "thread_per_chunk": 1,
    "main_flush_rows": 20000,
    "err_flush_rows": 1000,
    "max_region_size": 20000,
    "max_variants_per_region": 30,
}

medium_cfg = {
    "task_type": "medium",
    "thread_per_chunk": 1,
    "main_flush_rows": 10000,
    "err_flush_rows": 1000,
    "max_region_size": 5000,
    "max_variants_per_region": 10,
}

huge_cfg = {
    "task_type": "huge",
    "thread_per_chunk": 1,
    "main_flush_rows": 5000,
    "err_flush_rows": 1000,
    "max_region_size": 1000,
    "max_variants_per_region": 3,
}

mito_cfg = {
    "task_type": "mito",
    "thread_per_chunk": 1,
    "main_flush_rows": 5000,
    "err_flush_rows": 1000,
    "max_region_size": 150,
    "max_variants_per_region": 3,
}

def _pick_runtime_cfg_for_umi_chunk(task, small_cfg, medium_cfg, huge_cfg, mito_cfg=None):
    chrom = task.get("chrom", "")
    is_mito = chrom in ("chrM", "MT", "M")

    max_depth = task.get("max_depth", 0) or 0
    mean_depth = task.get("mean_depth", 0) or 0
    cost = task.get("cost", 0) or 0

    if is_mito and mito_cfg is not None:
        return mito_cfg, "mito"

    if max_depth >= 100000 or mean_depth >= 50000 or cost >= 150000:
        return huge_cfg, "huge"

    if max_depth >= 20000 or mean_depth >= 5000 or cost >= 20000:
        return medium_cfg, "medium"

    return small_cfg, "small"


# ── Step class ────────────────────────────────────────────────────────────────
class UMICombineStep(BaseStep):
    def get_inputs(self, context: Dict) -> Dict[str, str]:
        """input"""
        inputs = {
            "in_filter_bam": context.get("in_filter_bam"),
            "reference": self.config.get("genome_fasta"),
            "filter_mpileup_file": context.get("filter_mpileup_file"),
            "manifest_path": context.get("manifest_path", ""),
        }
        return inputs

    def get_outputs(self, context: Dict) -> Dict[str, str]:
        """output"""
        return {
            "spot_count_file": os.path.join(self.step_dir, "spot.count.list.csv"),
            "error_count_file": os.path.join(self.step_dir, "error.count.list.csv"),
        }

    def _run(self, context: Dict):
        inputs = self.get_inputs(context)
        bam_file = inputs["in_filter_bam"]
        # mpileup_file = inputs["filter_mpileup_file"]
        manifest_path = inputs["manifest_path"]

        seq_type = self.config.get("sequence_type")
        # bin_size = int(self.config.get("bin_size"))
        bin_size = 1
        thread = self.threads
        max_mem = self.memory
        tmp_root= self.step_dir

        outputs = self.get_outputs(context)
        spot_count_file = outputs["spot_count_file"]
        error_count_file = outputs["error_count_file"]

        chunk_output_dir = os.path.join(self.step_dir, "chunk_outputs")
        # final_output_dir = os.path.join(self.step_dir, "final_outputs")
        if seq_type=="visium":
            events_mem_fraction=1
            events_threshold=compute_events_threshold(
                    memory_limit_bytes=max_mem,
                    n_workers=thread,
                    bytes_per_event = 300,
                    events_mem_fraction = events_mem_fraction,
            )
        else:
            events_mem_fraction=0.7
            events_threshold=compute_events_threshold(
                memory_limit_bytes=max_mem,
                n_workers=thread,
                bytes_per_event = 300,
                events_mem_fraction = events_mem_fraction,
                min_threshold=1000000,
                max_threshold=3000000
            )

        os.makedirs(chunk_output_dir, exist_ok=True)
        # os.makedirs(final_output_dir, exist_ok=True)

        print("************events_threshold",events_threshold)
        t0 = time.time()
        spot_files, error_files = _run_umi_combine_parallel(
            manifest_path=manifest_path,
            bam_file=bam_file,
            seq_type=seq_type,
            bin_size=bin_size,
            output_dir=chunk_output_dir,
            max_workers=thread,
            max_mem=max_mem,
            small_cfg=small_cfg,
            medium_cfg=medium_cfg,
            huge_cfg=huge_cfg,
            mito_cfg=mito_cfg,
            compression="snappy",
            progress_interval=0.05,
            logger=logger,
            max_in_flight=thread,
            tmp_root=tmp_root,
            events_threshold=events_threshold
        )
        logger.info(f"[parallel umi_combine] {time.time() - t0:.2f}s")

        final_spot = spot_count_file
        final_error = error_count_file

        t1 = time.time()
        save_parquet_file_list(spot_files, final_spot)
        logger.info(f"[final spot file merge] {time.time() - t1:.2f}s")

        t2 = time.time()
        save_parquet_file_list(error_files, final_error)
        logger.info(f"[final error file merge] {time.time() - t2:.2f}s")


# ── Task building ─────────────────────────────────────────────────────────────
def _build_region_tasks_from_one_chunk_file_for_UMI_combine(
    chunk_task,
    max_region_size,
    max_variants_per_region,
):
    raw = read_chunk_bytes_by_offset(
        source_file=chunk_task["source_file"],
        start_offset=chunk_task["start_offset"],
        end_offset=chunk_task["end_offset"],
    )

    if not raw:
        return []

    df = pd.read_csv(
        io.BytesIO(raw),
        sep="\t",
        header=None,
        names=[
            "#chrom", "pos", "type", "ref", "alt1",
            "total_depth", "ref_depth", "alt1_depth", "vaf",
        ],
        comment="#",
        dtype={
            "#chrom": "string",
            "pos": "int64",
            "type": "string",
            "ref": "string",
            "alt1": "string",
            "total_depth": "int64",
            "ref_depth": "int64",
            "alt1_depth": "int64",
            "vaf": "float64",
        },
    )

    region_tasks = build_region_tasks_for_UMI_combine(
        df,
        max_region_size,
        max_variants_per_region,
    )

    del df
    return region_tasks


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
            compression=compression,
        )

    writer.write_table(table)

    del df
    del table
    return writer


def _umi_combine_to_parquet_buffered(
    bam_file,
    seq_type,
    bin_size,
    region_tasks,
    spot_count_parquet,
    error_count_parquet,
    main_flush_rows,
    err_flush_rows,
    compression,
    events_threshold,
    tmp_root="/tmp",
    chunk_id=None,
):
    main_writer = None
    err_writer = None
    main_buffer = []
    err_buffer = []
    bam_handle = None

    try:
        bam_handle = read_bam(bam_file)

        for region_idx, region_info in enumerate(region_tasks, 1):
            # if region_idx == 1 or region_idx % 10 == 0:
                # log_worker_mem(
                #         logger,
                #         "before_region",
                #         extra=(
                #             f"chunk_id={chunk_id} region_idx={region_idx}/{len(region_tasks)} "
                #             f"main_buf={len(main_buffer)} err_buf={len(err_buffer)} "
                #             f"{obj_size_text(main_buffer, 'main_buffer')} {obj_size_text(err_buffer, 'err_buffer')}"
                #         ),
                # )
            # result = process_single_region_for_umi_combine(
            #     bam_handle,
            #     region_info,
            #     seq_type,
            #     bin_size,
            # )

            result = process_single_region_for_umi_combine_auto(
                bam_handle=bam_handle,
                region_info=region_info,
                seq_type=seq_type,
                bin_size=bin_size,
                tmp_root=tmp_root,          # 可选
                threshold=3,
                weigh=0.5,
                events_threshold=events_threshold, 
                cell_dict={},
                debug_log=False
            )
            if result is None:
                continue

            for spot_count_list, error_list in result:
                if spot_count_list:
                    main_buffer.extend(spot_count_list)

                if error_list:
                    err_buffer.append(error_list) 

            if len(main_buffer) >= main_flush_rows:
                #***************
                # log_worker_mem(
                #     logger,
                #     "before_flush_main",
                #     extra=(
                #         f"chunk_id={chunk_id} rows={len(main_buffer)} "
                #         f"{obj_size_text(main_buffer, 'main_buffer')}"
                #     ),
                # )

                main_writer = _flush_rows_to_parquet(
                    rows=main_buffer,
                    columns=COLUMNS_MAIN,
                    writer=main_writer,
                    out_file=spot_count_parquet,
                    compression=compression,
                )
                main_buffer = []
                #*****************
                # log_worker_mem(
                #     logger,
                #     "after_flush_main",
                #     extra=f"chunk_id={chunk_id} rows={len(main_buffer)}",
                # )

            if len(err_buffer) >= err_flush_rows:
                #*****************
                # log_worker_mem(
                #     logger,
                #     "before_flush_err",
                #     extra=(
                #         f"chunk_id={chunk_id} rows={len(err_buffer)} "
                #         f"{obj_size_text(err_buffer, 'err_buffer')}"
                #     ),
                # )
                err_writer = _flush_rows_to_parquet(
                    rows=err_buffer,
                    columns=COLUMNS_ERROR_ALLELE,
                    writer=err_writer,
                    out_file=error_count_parquet,
                    compression=compression,
                )
                err_buffer=[]
                #*****************
                # log_worker_mem(
                #     logger,
                #     "after_flush_err",
                #     extra=f"chunk_id={chunk_id} rows={len(err_buffer)}",
                # )

        if main_buffer:
            #*****************
            # log_worker_mem(
            #     logger,
            #     "before_final_flush_main",
            #     extra=f"chunk_id={chunk_id} rows={len(main_buffer)} {obj_size_text(main_buffer, 'main_buffer')}",
            # )
            main_writer = _flush_rows_to_parquet(
                rows=main_buffer,
                columns=COLUMNS_MAIN,
                writer=main_writer,
                out_file=spot_count_parquet,
                compression=compression,
            )
            main_buffer.clear()

        if err_buffer:
            #***********
            # log_worker_mem(
            #     logger,
            #     "before_final_flush_err",
            #     extra=f"chunk_id={chunk_id} rows={len(err_buffer)} {obj_size_text(err_buffer, 'err_buffer')}",
            # )
            err_writer = _flush_rows_to_parquet(
                rows=err_buffer,
                columns=COLUMNS_ERROR_ALLELE,
                writer=err_writer,
                out_file=error_count_parquet,
                compression=compression,
            )
            err_buffer.clear()

    finally:
        # log_worker_mem(logger, "before_cleanup", extra=f"chunk_id={chunk_id}")
        try:
            if bam_handle is not None:
                bam_handle.close()
        except Exception:
            pass

        if main_writer is not None:
            main_writer.close()
        if err_writer is not None:
            err_writer.close()

        del main_buffer
        del err_buffer
        gc.collect()


def _run_umi_combine_one_chunk(
    chunk_task,
    bam_file,
    seq_type,
    bin_size,
    output_dir,
    max_region_size,
    max_variants_per_region,
    main_flush_rows,
    err_flush_rows,
    compression,
    tmp_root,
    events_threshold
):
    os.makedirs(output_dir, exist_ok=True)

    chunk_id = chunk_task["chunk_id"]
    chunk_prefix = chunk_id

    spot_count_parquet = os.path.join(output_dir, f"{chunk_prefix}.spot.parquet")
    error_count_parquet = os.path.join(output_dir, f"{chunk_prefix}.error.parquet")
    
    # **********
    # log_worker_mem( 
    #     logger,
    #     "chunk_start",
    #     extra=(
    #         f"chunk_id={chunk_id} chrom={chunk_task.get('chrom')} "
    #         f"records={chunk_task.get('records')} span_bp={chunk_task.get('span_bp')} "
    #         f"max_depth={chunk_task.get('max_depth')} mean_depth={chunk_task.get('mean_depth')}"
    #     ),
    # )

    region_tasks = _build_region_tasks_from_one_chunk_file_for_UMI_combine(
        chunk_task,
        max_region_size,
        max_variants_per_region,
    )

    # **********
    # log_worker_mem(
    #     logger,
    #     "after_build_region_tasks",
    #     extra=f"chunk_id={chunk_id} n_region_tasks={len(region_tasks)} {obj_size_text(region_tasks, 'region_tasks')}",
    # )

    if len(region_tasks) == 0:
        return {
            "chunk_id": chunk_id,
            "spot_parquet": None,
            "error_parquet": None,
        }

    _umi_combine_to_parquet_buffered(
        bam_file=bam_file,
        seq_type=seq_type,
        bin_size=bin_size,
        region_tasks=region_tasks,
        spot_count_parquet=spot_count_parquet,
        error_count_parquet=error_count_parquet,
        main_flush_rows=main_flush_rows,
        err_flush_rows=err_flush_rows,
        compression=compression,
        tmp_root=tmp_root,
        events_threshold=events_threshold
    )

    # ***************
    # log_worker_mem(logger, "after_chunk_compute", extra=f"chunk_id={chunk_id}")

    del region_tasks
    gc.collect()

    # ***************
    # log_worker_mem(logger, "after_chunk_gc", extra=f"chunk_id={chunk_id}")

    return {
        "chunk_id": chunk_id,
        "spot_parquet": spot_count_parquet,
        "error_parquet": error_count_parquet,
    }


def _build_umi_combine_tasks(manifest_path, small_cfg,medium_cfg,huge_cfg, mito_cfg):
    tasks = load_manifest_files_for_chunk_tasks_only_umi_combine(manifest_path)
    if not tasks:
        return []

    merged = []
    for task in tasks:
        chrom = task.get("chrom", "")
        is_mito = chrom in ("chrM", "MT", "M")
        runtime_cfg, risk_type = _pick_runtime_cfg_for_umi_chunk(
            task=task,
            small_cfg=small_cfg,
            medium_cfg=medium_cfg,
            huge_cfg=huge_cfg,
            mito_cfg=mito_cfg,
        )
        
        new_task = dict(task)
        new_task["is_mito"] = is_mito
        new_task["risk_type"] = risk_type
        new_task["runtime_cfg"] = runtime_cfg
        merged.append(new_task)

    # sort by cost, large->small
    merged.sort(key=lambda x: (-x["cost"], x["chrom"], x["chunk_idx"]))

    # max n_heavy 
    n_heavy = min(20, len(merged))
    heavy = merged[:n_heavy]

    # sort by sort, small->large
    light = sorted(
        merged[n_heavy:],
        key=lambda x: (x["cost"], x["chrom"], x["chunk_idx"])
    )

    small_gap = 30

    ordered = []
    i = 0
    j = 0

    while i < len(heavy) or j < len(light):
        if i < len(heavy):
            ordered.append(heavy[i])
            i += 1

        for _ in range(small_gap):
            if j < len(light):
                ordered.append(light[j])
                j += 1

    return ordered


def _run_single_umi_combine_task(
    task,
    bam_file,
    seq_type,
    bin_size,
    output_dir,
    events_threshold,
    compression="snappy",
    tmp_root="/tmp",
):
    cfg = task["runtime_cfg"]

    return _run_umi_combine_one_chunk(
        chunk_task=task,
        bam_file=bam_file,
        seq_type=seq_type,
        bin_size=bin_size,
        output_dir=output_dir,
        max_region_size=cfg["max_region_size"],
        max_variants_per_region=cfg["max_variants_per_region"],
        main_flush_rows=cfg["main_flush_rows"],
        err_flush_rows=cfg["err_flush_rows"],
        compression=compression,
        tmp_root=tmp_root,
        events_threshold=events_threshold
    )


def _run_umi_combine_parallel(
    manifest_path,
    bam_file,
    seq_type,
    bin_size,
    output_dir,
    max_workers,
    max_mem,
    small_cfg,
    medium_cfg,
    huge_cfg,
    mito_cfg,
    events_threshold,
    compression="snappy",
    progress_interval=0.05,
    logger=None,
    max_in_flight=None,
    tmp_root="/tmp",
):
    tasks = _build_umi_combine_tasks(
        manifest_path=manifest_path,
        small_cfg=small_cfg,
        medium_cfg=medium_cfg,
        huge_cfg=huge_cfg,
        mito_cfg=mito_cfg,
    )

    if not tasks:
        if logger:
            logger.warning("[umi_combine] no tasks found")
        return [], []

    if logger:
        mito_n = sum(1 for t in tasks if t["is_mito"])
        logger.info(
            f"[umi_combine] unified queue: total={len(tasks)}, mito={mito_n}, "
            f"normal={len(tasks) - mito_n}, max_workers={max_workers}, "
            f"max_in_flight={max_in_flight if max_in_flight is not None else 'default'}"
        )

    worker_fn = partial(
        _run_single_umi_combine_task,
        bam_file=bam_file,
        seq_type=seq_type,
        bin_size=bin_size,
        output_dir=output_dir,
        compression=compression,
        tmp_root=tmp_root,
        events_threshold=events_threshold
    )

    results = parallel_map(
        items=tasks,
        worker_fn=worker_fn,
        max_workers=max_workers,
        memory_limit_bytes=max_mem,
        desc="umi_combine",
        raise_on_error=True,
        backend="process",
        progress_interval=progress_interval,
        logger=logger,
        worker_takes_tuple=False,
        max_in_flight=max_in_flight if max_in_flight is not None else  max_workers,
        # top_n_children=1, #************
        # debug_children_every_print=True, #************
    )

    all_spot_parquets = []
    all_error_parquets = []

    for r in results:
        if not r:
            continue
        if r.get("spot_parquet") is not None:
            all_spot_parquets.append(r["spot_parquet"])
        if r.get("error_parquet") is not None:
            all_error_parquets.append(r["error_parquet"])

    return all_spot_parquets, all_error_parquets


# ── Output helpers ────────────────────────────────────────────────────────────
def save_parquet_file_list(chunk_files: List[Union[str, Tuple[str, str]]], output_file: str):
    valid_rows = []

    for item in chunk_files:
        chunk = None
        file = None

        if isinstance(item, (tuple, list)) and len(item) == 2:
            chunk, file = item
        elif isinstance(item, str):
            file = item
            basename = os.path.basename(file)

            if basename.endswith(".spot.parquet"):
                chunk = basename[:-len(".spot.parquet")]
            elif basename.endswith(".error.parquet"):
                chunk = basename[:-len(".error.parquet")]
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
        writer = csv.writer(f, delimiter="\t")
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
                continue

            if not os.path.exists(file):
                continue

            table = pq.read_table(file)

            if writer is None:
                writer = pq.ParquetWriter(
                    output_file,
                    table.schema,
                    compression=compression,
                )

            writer.write_table(table)

            del table

    finally:
        if writer is not None:
            writer.close()

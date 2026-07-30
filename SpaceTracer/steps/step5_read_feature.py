import gc
import os
import time
from functools import partial
from typing import Dict, List

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

from SpaceTracer.cores.read_feature_01_process import handel_bam_file_for_region
from SpaceTracer.cores.read_feature_02_extract import readLevelFeatures
from SpaceTracer.steps.base import BaseStep
from SpaceTracer.utils.get_read_level_feature import detect_read_length
from SpaceTracer.utils.logger import get_logger
from SpaceTracer.utils.parallel import parallel_map
from SpaceTracer.utils.read_files import read_bam
from SpaceTracer.utils.utils import (
    barcode_cell_mapping,
    get_regions,
    load_manifest_tsv,
    save_manifest_tsv,
)

model_name = __name__
logger = get_logger(model_name)

VARIANT_KEY_COLUMNS = ["#chrom", "pos", "ref", "alt"]

small_cfg = {
    "task_type": "small",
    "flush_rows": 500,
    "max_region_size": 20000,
    "max_variants_per_region": 100,
}

medium_cfg = {
    "task_type": "medium",
    "flush_rows": 200,
    "max_region_size": 5000,
    "max_variants_per_region": 30,
}

huge_cfg = {
    "task_type": "huge",
    "flush_rows": 100,
    "max_region_size": 1000,
    "max_variants_per_region": 5,
}

mito_cfg = {
    "task_type": "mito",
    "flush_rows": 100,
    "max_region_size": 150,
    "max_variants_per_region": 3,
}


def _load_mutation_identifiers(mutation_list_file: str) -> List[str]:
    identifiers = []
    with open(mutation_list_file, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            line = line.strip()
            if line:
                identifiers.append(line)
    return identifiers


def _chunk_has_mito_variants(identifiers: List[str]) -> bool:
    mito_chroms = {"chrM", "MT", "M", "chrMT"}
    for identifier in identifiers:
        chrom = identifier.split("\t", 1)[0]
        if chrom in mito_chroms:
            return True
    return False


def _pick_runtime_cfg_for_read_chunk(
    n_variants: int,
    is_mito: bool,
    small_cfg,
    medium_cfg,
    huge_cfg,
    mito_cfg=None,
):
    if is_mito and mito_cfg is not None:
        return mito_cfg, "mito"

    if n_variants >= 5000:
        return huge_cfg, "huge"

    if n_variants >= 1000:
        return medium_cfg, "medium"

    return small_cfg, "small"


def _estimate_chunk_variant_count(row: Dict[str, str]) -> int:
    mutation_list_file = row.get("ind_geno_filter_mutation_list", "")
    if not mutation_list_file or not os.path.exists(mutation_list_file):
        return 0

    count = 0
    with open(mutation_list_file, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.strip():
                count += 1
    return count


def _build_read_feature_tasks(manifest_rows, small_cfg, medium_cfg, huge_cfg, mito_cfg):
    merged = []
    for row in manifest_rows:
        identifiers = _load_mutation_identifiers(row["ind_geno_filter_mutation_list"])
        n_variants = len(identifiers)
        is_mito = _chunk_has_mito_variants(identifiers)
        runtime_cfg, risk_type = _pick_runtime_cfg_for_read_chunk(
            n_variants=n_variants,
            is_mito=is_mito,
            small_cfg=small_cfg,
            medium_cfg=medium_cfg,
            huge_cfg=huge_cfg,
            mito_cfg=mito_cfg,
        )

        new_task = dict(row)
        new_task["n_variants"] = n_variants
        new_task["is_mito"] = is_mito
        new_task["risk_type"] = risk_type
        new_task["runtime_cfg"] = runtime_cfg
        new_task["identifiers"] = identifiers
        merged.append(new_task)

    merged.sort(key=lambda x: (-x["n_variants"], x.get("chunk", "")))

    n_heavy = min(20, len(merged))
    heavy = merged[:n_heavy]
    light = sorted(
        merged[n_heavy:],
        key=lambda x: (x["n_variants"], x.get("chunk", "")),
    )

    ordered = []
    i = 0
    j = 0
    small_gap = 30

    while i < len(heavy) or j < len(light):
        if i < len(heavy):
            ordered.append(heavy[i])
            i += 1

        for _ in range(small_gap):
            if j < len(light):
                ordered.append(light[j])
                j += 1

    return ordered


def _process_single_region(
    region_info,
    bam_handle,
    bam_file,
    run_type,
    bins,
    cell_dict,
    readLen,
    downsample,
    target_depth,
    seed,
):
    var_result_dict = handel_bam_file_for_region(
        bam_file=bam_file,
        region_dict=region_info,
        run_type=run_type,
        bins=bins,
        cell_dict=cell_dict,
        readLen=readLen,
        downsample=downsample,
        target_depth=target_depth,
        seed=seed,
        bam_handle=bam_handle,
    )

    region_features = []
    for identifier, read_info_dict in var_result_dict.items():
        feature_dict = readLevelFeatures.from_read_info_to_dict(
            identifier=identifier,
            read_info_dict=read_info_dict,
        )
        if feature_dict is not None:
            region_features.append(feature_dict)

    return region_features


def _normalize_feature_rows(rows: List[Dict]) -> pd.DataFrame:
    df = pd.DataFrame(rows)
    if df.empty:
        return df
    if "chrom" in df.columns:
        df = df.rename(columns={"chrom": "#chrom"})
    return df


def _format_feature_value_for_parquet(value):
    if value is None:
        return pd.NA
    if isinstance(value, str):
        if value == "" or value == "NA":
            return pd.NA
        return value
    if isinstance(value, float) and pd.isna(value):
        return pd.NA
    return str(value)


def _sanitize_feature_df_for_parquet(df: pd.DataFrame) -> pd.DataFrame:
    """
    read_level feature columns may mix floats and the literal string 'NA'
    (e.g. from refine_mean). Coerce value columns to nullable string so
    buffered parquet writes keep a stable schema across flushes.
    """
    if df.empty:
        return df

    sanitized = df.copy()
    for col in sanitized.columns:
        sanitized[col] = sanitized[col].map(_format_feature_value_for_parquet)
        sanitized[col] = sanitized[col].astype("string")
    return sanitized


def _write_empty_chunk_output(output_file: str):
    empty_df = pd.DataFrame(columns=VARIANT_KEY_COLUMNS)
    empty_df.to_csv(output_file, sep="\t", index=False)

    parquet_file = str(output_file).replace(".txt", ".parquet")
    empty_indexed = empty_df.set_index(VARIANT_KEY_COLUMNS)
    empty_indexed.to_parquet(
        parquet_file, index=True, engine="pyarrow", compression="snappy"
    )


def _flush_feature_rows_to_parquet(rows, writer, parquet_schema, out_file, compression="snappy"):
    if not rows:
        return writer, parquet_schema

    df = _normalize_feature_rows(rows)
    df = df.set_index(VARIANT_KEY_COLUMNS)
    df = _sanitize_feature_df_for_parquet(df)
    table = pa.Table.from_pandas(df, preserve_index=True)

    if writer is None:
        parquet_schema = table.schema
        writer = pq.ParquetWriter(out_file, parquet_schema, compression=compression)
    else:
        table = table.cast(parquet_schema)

    writer.write_table(table)
    del df
    del table
    return writer, parquet_schema


def _append_feature_rows_to_txt(rows, txt_handle, header_written: bool) -> bool:
    if not rows:
        return header_written

    df = _normalize_feature_rows(rows)
    df.to_csv(txt_handle, sep="\t", index=False, header=not header_written)
    return True


def _read_feature_to_outputs_buffered(
    bam_file,
    region_dict_list,
    run_type,
    bins,
    cell_dict,
    readLen,
    downsample,
    target_depth,
    seed,
    output_file,
    flush_rows,
    compression="snappy",
) -> int:
    parquet_file = str(output_file).replace(".txt", ".parquet")
    parquet_writer = None
    parquet_schema = None
    feature_buffer = []
    bam_handle = None
    total_features = 0
    header_written = False

    try:
        bam_handle = read_bam(bam_file)

        with open(output_file, "w", encoding="utf-8") as txt_handle:
            for region_info in region_dict_list:
                region_features = _process_single_region(
                    region_info=region_info,
                    bam_handle=bam_handle,
                    bam_file=bam_file,
                    run_type=run_type,
                    bins=bins,
                    cell_dict=cell_dict,
                    readLen=readLen,
                    downsample=downsample,
                    target_depth=target_depth,
                    seed=seed,
                )

                if not region_features:
                    continue

                feature_buffer.extend(region_features)
                total_features += len(region_features)

                if len(feature_buffer) >= flush_rows:
                    header_written = _append_feature_rows_to_txt(
                        feature_buffer, txt_handle, header_written
                    )
                    parquet_writer, parquet_schema = _flush_feature_rows_to_parquet(
                        rows=feature_buffer,
                        writer=parquet_writer,
                        parquet_schema=parquet_schema,
                        out_file=parquet_file,
                        compression=compression,
                    )
                    feature_buffer = []

            if feature_buffer:
                header_written = _append_feature_rows_to_txt(
                    feature_buffer, txt_handle, header_written
                )
                parquet_writer, parquet_schema = _flush_feature_rows_to_parquet(
                    rows=feature_buffer,
                    writer=parquet_writer,
                    parquet_schema=parquet_schema,
                    out_file=parquet_file,
                    compression=compression,
                )
                feature_buffer.clear()

    finally:
        try:
            if bam_handle is not None:
                bam_handle.close()
        except Exception:
            pass

        if parquet_writer is not None:
            parquet_writer.close()

        del feature_buffer
        gc.collect()

    if total_features == 0:
        _write_empty_chunk_output(output_file)
        return 0

    return total_features


def _run_read_feature_one_chunk(
    task,
    bam_file,
    seq_type,
    bins,
    cell_dict,
    readLen,
    downsample,
    target_depth,
    seed,
    output_file,
    compression="snappy",
):
    chunk = task["chunk"]
    mutation_list_file = task["ind_geno_filter_mutation_list"]
    cfg = task["runtime_cfg"]

    if not mutation_list_file or not os.path.exists(mutation_list_file):
        raise FileNotFoundError(
            f"mutation list file not found for chunk={chunk}: {mutation_list_file}"
        )

    parquet_file = str(output_file).replace(".txt", ".parquet")
    identifiers = task.get("identifiers")
    if identifiers is None:
        identifiers = _load_mutation_identifiers(mutation_list_file)

    if len(identifiers) == 0:
        _write_empty_chunk_output(output_file)
        return {
            "chunk": chunk,
            "ind_geno_filter_mutation_list": mutation_list_file,
            "read_feature_txt": output_file,
            "read_feature_parquet": parquet_file,
        }

    region_dict_list = get_regions(
        identifiers,
        cfg["max_region_size"],
        cfg["max_variants_per_region"],
    )

    if not region_dict_list:
        _write_empty_chunk_output(output_file)
        return {
            "chunk": chunk,
            "ind_geno_filter_mutation_list": mutation_list_file,
            "read_feature_txt": output_file,
            "read_feature_parquet": parquet_file,
        }

    _read_feature_to_outputs_buffered(
        bam_file=bam_file,
        region_dict_list=region_dict_list,
        run_type=seq_type,
        bins=bins,
        cell_dict=cell_dict,
        readLen=readLen,
        downsample=downsample,
        target_depth=target_depth,
        seed=seed,
        output_file=output_file,
        flush_rows=cfg["flush_rows"],
        compression=compression,
    )

    del region_dict_list
    gc.collect()

    return {
        "chunk": chunk,
        "ind_geno_filter_mutation_list": mutation_list_file,
        "read_feature_txt": output_file,
        "read_feature_parquet": parquet_file,
    }


def _run_single_read_feature_task(
    task,
    bam_file,
    seq_type,
    bins,
    cell_dict,
    readLen,
    downsample,
    target_depth,
    seed,
    step_dir,
    sample,
    compression="snappy",
):
    chunk = task["chunk"]
    chunk_dir = os.path.join(step_dir, sample, chunk)
    os.makedirs(chunk_dir, exist_ok=True)
    output_file = os.path.join(chunk_dir, "read_feature.txt")

    return _run_read_feature_one_chunk(
        task=task,
        bam_file=bam_file,
        seq_type=seq_type,
        bins=bins,
        cell_dict=cell_dict,
        readLen=readLen,
        downsample=downsample,
        target_depth=target_depth,
        seed=seed,
        output_file=output_file,
        compression=compression,
    )


def _run_read_feature_parallel(
    valid_rows,
    bam_file,
    seq_type,
    bins,
    cell_dict,
    readLen,
    downsample,
    target_depth,
    seed,
    step_dir,
    sample,
    max_workers,
    max_mem,
    small_cfg,
    medium_cfg,
    huge_cfg,
    mito_cfg,
    compression="snappy",
    progress_interval=0.05,
    max_in_flight=None,
):
    tasks = _build_read_feature_tasks(
        valid_rows,
        small_cfg=small_cfg,
        medium_cfg=medium_cfg,
        huge_cfg=huge_cfg,
        mito_cfg=mito_cfg,
    )

    if not tasks:
        logger.warning("[read_feature] no tasks found")
        return []

    mito_n = sum(1 for t in tasks if t["is_mito"])
    logger.info(
        f"[read_feature] unified queue: total={len(tasks)}, mito={mito_n}, "
        f"normal={len(tasks) - mito_n}, max_workers={max_workers}, "
        f"max_in_flight={max_in_flight if max_in_flight is not None else max_workers}"
    )

    worker_fn = partial(
        _run_single_read_feature_task,
        bam_file=bam_file,
        seq_type=seq_type,
        bins=bins,
        cell_dict=cell_dict,
        readLen=readLen,
        downsample=downsample,
        target_depth=target_depth,
        seed=seed,
        step_dir=step_dir,
        sample=sample,
        compression=compression,
    )

    return parallel_map(
        items=tasks,
        worker_fn=worker_fn,
        max_workers=max_workers,
        memory_limit_bytes=max_mem,
        desc="read_feature",
        raise_on_error=True,
        backend="process",
        progress_interval=progress_interval,
        logger=logger,
        max_in_flight=max_in_flight if max_in_flight is not None else max_workers,
    )


class ReadFeatureStep(BaseStep):
    def get_inputs(self, context):
        return {
            "raw_bam": context.get("bam_file"),
            "genotype_results": (
                context.get("genotype_results", "")
                or context.get("ind_geno_filter_mutation_list", "")
            ),
        }

    def get_outputs(self, context):
        sample = context.get("sample", "sample")
        out_dir = os.path.join(self.step_dir, sample)
        os.makedirs(out_dir, exist_ok=True)

        return {
            "read_feature_results": os.path.join(out_dir, "read_feature_chunk_manifest.tsv")
        }

    def get_step_config(self):
        return self.config.get("steps", {}).get("read_feature", {})

    def _run(self, context):
        t0 = time.time()

        inputs = self.get_inputs(context)
        bam_file = inputs["raw_bam"]
        chunk_manifest = inputs["genotype_results"]

        if not bam_file or not os.path.exists(bam_file):
            raise FileNotFoundError(f"bam file not found: {bam_file}")
        if not chunk_manifest or not os.path.exists(chunk_manifest):
            raise FileNotFoundError(f"chunk manifest not found: {chunk_manifest}")

        output_file = self.get_outputs(context)["read_feature_results"]

        step_config = self.get_step_config()
        cell_info_file = step_config.get("cell_info", "")
        downsample = step_config["downsample"]
        target_depth = step_config["downsample_target_depth"]
        seed = step_config["seed"]

        seq_type = self.config.get("sequence_type")
        bins = self.config.get("bin_size")
        sample = context.get("sample", "sample")

        readLen = detect_read_length(bam_file)

        if cell_info_file:
            cell_dict = barcode_cell_mapping(cell_info_file)
        else:
            cell_dict = {}

        rows = load_manifest_tsv(chunk_manifest)
        if not rows:
            raise ValueError(f"No chunk records found in chunk manifest: {chunk_manifest}")

        valid_rows = []
        for row in rows:
            if not row:
                continue
            if not row.get("chunk"):
                continue
            if not row.get("ind_geno_filter_mutation_list"):
                continue
            valid_rows.append(row)

        if not valid_rows:
            raise ValueError(f"No valid chunk rows found in manifest: {chunk_manifest}")

        max_workers = self.config.get("runtime", {}).get("max_parallel", self.threads)

        chunk_results = _run_read_feature_parallel(
            valid_rows=valid_rows,
            bam_file=bam_file,
            seq_type=seq_type,
            bins=bins,
            cell_dict=cell_dict,
            readLen=readLen,
            downsample=downsample,
            target_depth=target_depth,
            seed=seed,
            step_dir=self.step_dir,
            sample=sample,
            max_workers=max_workers,
            max_mem=self.memory,
            small_cfg=small_cfg,
            medium_cfg=medium_cfg,
            huge_cfg=huge_cfg,
            mito_cfg=mito_cfg,
            max_in_flight=max_workers,
        )

        save_manifest_tsv(chunk_results, output_file)
        logger.info(f"[parallel read_feature] {time.time() - t0:.2f}s")

        return {
            "read_feature_results": output_file
        }

import subprocess
from typing import Dict
import pandas as pd
import os
import time
from datetime import datetime
from pandas.errors import EmptyDataError

from SpaceTracer.cores.read_feature_01_process import handel_bam_file_for_region
from SpaceTracer.cores.read_feature_02_extract import readLevelFeatures
from SpaceTracer.steps.base import BaseStep
from SpaceTracer.utils.get_read_level_feature import detect_read_length
from SpaceTracer.utils.parallel import parallel_map
from SpaceTracer.utils.utils import (
    get_regions,
    barcode_cell_mapping,
    load_manifest_tsv,
    save_manifest_tsv
)
from SpaceTracer.utils.logger import get_logger

model_name = __name__
logger = get_logger(model_name)


def _now():
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def _rss_mb():
    try:
        import psutil
        return psutil.Process(os.getpid()).memory_info().rss / 1024 / 1024
    except Exception:
        return None


def _fmt_rss():
    rss = _rss_mb()
    if rss is None:
        return "NA"
    return f"{rss:.1f}MB"


class SimpleTimer:
    def __init__(self, name: str):
        self.name = name
        self.t0 = None

    def __enter__(self):
        self.t0 = time.perf_counter()
        logger.info(f"[PROFILE_START][{_now()}] {self.name} | rss={_fmt_rss()}")
        return self

    def __exit__(self, exc_type, exc, tb):
        cost = time.perf_counter() - self.t0
        status = "OK" if exc_type is None else f"ERROR={exc_type.__name__}"
        logger.info(
            f"[PROFILE_END][{_now()}] {self.name} | "
            f"cost={cost:.2f}s | rss={_fmt_rss()} | status={status}"
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

    def _get_chunk_output_path(self, context: Dict, chunk: str) -> str:
        sample = context.get("sample", "sample")
        chunk_dir = os.path.join(self.step_dir, sample, chunk)
        os.makedirs(chunk_dir, exist_ok=True)
        return os.path.join(chunk_dir, "read_feature.txt")

    def _write_empty_chunk_output(self, output_file: str):
        t_total = time.perf_counter()

        colnames = ['#chrom', 'pos', 'ref', 'alt']

        t0 = time.perf_counter()
        empty_df = pd.DataFrame(columns=colnames)
        df_cost = time.perf_counter() - t0

        t0 = time.perf_counter()
        empty_df.to_csv(output_file, sep='\t', index=False)
        write_txt_cost = time.perf_counter() - t0

        parquet_file = str(output_file).replace('.txt', '.parquet')

        t0 = time.perf_counter()
        empty_df = empty_df.set_index(['#chrom', 'pos', 'ref', 'alt'])
        set_index_cost = time.perf_counter() - t0

        t0 = time.perf_counter()
        empty_df.to_parquet(parquet_file, index=True, engine='pyarrow', compression='snappy')
        write_parquet_cost = time.perf_counter() - t0

        total_cost = time.perf_counter() - t_total

    def _run_one_chunk(
        self,
        row: Dict[str, str],
        bam_file: str,
        seq_type,
        bins,
        cell_dict: Dict,
        readLen: int,
        downsample,
        target_depth,
        seed,
        max_region_size,
        max_variants_per_region,
        context: Dict,
    ) -> Dict[str, str]:
        t_chunk_total = time.perf_counter()

        chunk = row["chunk"]
        mutation_list_file = row["ind_geno_filter_mutation_list"]

        if not mutation_list_file or not os.path.exists(mutation_list_file):
            raise FileNotFoundError(f"mutation list file not found for chunk={chunk}: {mutation_list_file}")

        output_file = self._get_chunk_output_path(context, chunk)
        parquet_file = str(output_file).replace('.txt', '.parquet')

        t0 = time.perf_counter()
        identifiers = []
        with open(mutation_list_file, "r", encoding="utf-8", errors="ignore") as f:
            for line in f:
                line = line.strip()
                if line:
                    identifiers.append(line)
        load_mutation_cost = time.perf_counter() - t0

        if len(identifiers) == 0:
            self._write_empty_chunk_output(output_file)

            total_cost = time.perf_counter() - t_chunk_total

            return {
                "chunk": chunk,
                "ind_geno_filter_mutation_list": mutation_list_file,
                "read_feature_txt": output_file,
                "read_feature_parquet": parquet_file,
            }

        t0 = time.perf_counter()
        region_dict_list = get_regions(identifiers, max_region_size, max_variants_per_region)
        get_regions_cost = time.perf_counter() - t0

        if not region_dict_list:
            self._write_empty_chunk_output(output_file)

            total_cost = time.perf_counter() - t_chunk_total

            return {
                "chunk": chunk,
                "ind_geno_filter_mutation_list": mutation_list_file,
                "read_feature_txt": output_file,
                "read_feature_parquet": parquet_file,
            }

        t0 = time.perf_counter()
        self._process_and_save(
            bam_file=bam_file,
            region_dict_list=region_dict_list,
            run_type=seq_type,
            bins=bins,
            cell_dict=cell_dict,
            readLen=readLen,
            downsample=downsample,
            target_depth=target_depth,
            seed=seed,
            n_processes=1,
            output_file=output_file,
        )
        process_save_cost = time.perf_counter() - t0

        t0 = time.perf_counter()
        if not os.path.exists(output_file) or os.path.getsize(output_file) == 0:
            self._write_empty_chunk_output(output_file)
        output_check_cost = time.perf_counter() - t0

        total_cost = time.perf_counter() - t_chunk_total


        return {
            "chunk": chunk,
            "ind_geno_filter_mutation_list": mutation_list_file,
            "read_feature_txt": output_file,
            "read_feature_parquet": parquet_file,
        }

    def _run(self, context):
        t_run_total = time.perf_counter()

        t0 = time.perf_counter()
        inputs = self.get_inputs(context)
        bam_file = inputs["raw_bam"]
        chunk_manifest = inputs["genotype_results"]
        get_inputs_cost = time.perf_counter() - t0

        t0 = time.perf_counter()
        if not bam_file or not os.path.exists(bam_file):
            raise FileNotFoundError(f"bam file not found: {bam_file}")
        if not chunk_manifest or not os.path.exists(chunk_manifest):
            raise FileNotFoundError(f"chunk manifest not found: {chunk_manifest}")
        check_inputs_cost = time.perf_counter() - t0

        t0 = time.perf_counter()
        output_file = self.get_outputs(context)["read_feature_results"]
        get_outputs_cost = time.perf_counter() - t0

        t0 = time.perf_counter()
        step_config = self.get_step_config()
        cell_info_file = step_config.get("cell_info", "")
        downsample = step_config["downsample"]
        target_depth = step_config["downsample_target_depth"]
        max_region_size = step_config["max_region_size"]
        max_variants_per_region = step_config["max_variants_per_region"]
        seed = step_config["seed"]

        seq_type = self.config.get("sequence_type")
        bins = self.config.get("bin_size")
        load_config_cost = time.perf_counter() - t0

        t0 = time.perf_counter()
        readLen = detect_read_length(bam_file)
        detect_read_length_cost = time.perf_counter() - t0

        t0 = time.perf_counter()
        if cell_info_file:
            cell_dict = barcode_cell_mapping(cell_info_file)
        else:
            cell_dict = {}
        cell_mapping_cost = time.perf_counter() - t0

        t0 = time.perf_counter()
        rows = load_manifest_tsv(chunk_manifest)
        load_manifest_cost = time.perf_counter() - t0

        if not rows:
            raise ValueError(f"No chunk records found in chunk manifest: {chunk_manifest}")

        t0 = time.perf_counter()
        valid_rows = []
        for row in rows:
            if not row:
                continue
            if not row.get("chunk"):
                continue
            if not row.get("ind_geno_filter_mutation_list"):
                continue
            valid_rows.append(row)
        filter_rows_cost = time.perf_counter() - t0

        if not valid_rows:
            raise ValueError(f"No valid chunk rows found in manifest: {chunk_manifest}")

        max_workers = self.config.get("runtime", {}).get("max_parallel", self.threads)

        def worker(row: Dict[str, str]) -> Dict[str, str]:
            return self._run_one_chunk(
                row=row,
                bam_file=bam_file,
                seq_type=seq_type,
                bins=bins,
                cell_dict=cell_dict,
                readLen=readLen,
                downsample=downsample,
                target_depth=target_depth,
                seed=seed,
                max_region_size=max_region_size,
                max_variants_per_region=max_variants_per_region,
                context=context,
            )

        t0 = time.perf_counter()
        chunk_results = parallel_map(
            valid_rows,
            worker_fn=worker,
            max_workers=max_workers,
            desc=f"read_feature",
            raise_on_error=True,
        )
        parallel_map_cost = time.perf_counter() - t0

        t0 = time.perf_counter()
        save_manifest_tsv(chunk_results, output_file)
        save_manifest_cost = time.perf_counter() - t0

        total_cost = time.perf_counter() - t_run_total

        return {
            "read_feature_results": output_file
        }

    @staticmethod
    def _process_single_region(
        region_info,
        bam_file,
        run_type,
        bins,
        cell_dict,
        readLen,
        downsample,
        target_depth,
        seed
    ):
        t0 = time.perf_counter()
        var_result_dict = handel_bam_file_for_region(
            bam_file=bam_file,
            region_dict=region_info,
            run_type=run_type,
            bins=bins,
            cell_dict=cell_dict,
            readLen=readLen,
            downsample=downsample,
            target_depth=target_depth,
            seed=seed
        )
        bam_process_cost = time.perf_counter() - t0

        t0 = time.perf_counter()
        region_features = []
        for identifier, read_info_dict in var_result_dict.items():
            feature_dict = readLevelFeatures.from_read_info_to_dict(
                identifier=identifier,
                read_info_dict=read_info_dict
            )
            if feature_dict is not None:
                region_features.append(feature_dict)
        feature_extract_cost = time.perf_counter() - t0

        return region_features

    def _process_and_save(
        self,
        bam_file,
        region_dict_list,
        run_type,
        bins,
        cell_dict,
        readLen,
        downsample,
        target_depth,
        seed,
        n_processes,
        output_file
    ):
        from functools import partial

        t_total = time.perf_counter()

        t0 = time.perf_counter()
        process_func = partial(
            self._process_single_region,
            bam_file=bam_file,
            run_type=run_type,
            bins=bins,
            cell_dict=cell_dict,
            readLen=readLen,
            downsample=downsample,
            target_depth=target_depth,
            seed=seed
        )
        make_partial_cost = time.perf_counter() - t0

        tasks = region_dict_list
        total_features = 0
        header_written = False

        process_regions_cost = 0.0
        build_df_cost = 0.0
        rename_cost = 0.0
        write_txt_cost = 0.0

        n_regions = 0
        n_nonempty_regions = 0
        n_empty_regions = 0

        t0 = time.perf_counter()
        with open(output_file, 'w', encoding='utf-8') as f:
            open_output_cost = time.perf_counter() - t0

            for region_info in tasks:
                n_regions += 1

                t_region = time.perf_counter()
                region_features = process_func(region_info)
                process_regions_cost += time.perf_counter() - t_region

                if not region_features:
                    n_empty_regions += 1
                    continue

                n_nonempty_regions += 1

                t_df = time.perf_counter()
                df = pd.DataFrame(region_features)
                build_df_cost += time.perf_counter() - t_df

                if df.empty:
                    continue

                t_rename = time.perf_counter()
                if 'chrom' in df.columns:
                    df = df.rename(columns={'chrom': '#chrom'})
                rename_cost += time.perf_counter() - t_rename

                t_write = time.perf_counter()
                if not header_written:
                    df.to_csv(f, sep='\t', index=False, header=True)
                    header_written = True
                else:
                    df.to_csv(f, sep='\t', index=False, header=False)
                write_txt_cost += time.perf_counter() - t_write

                total_features += len(region_features)

        if total_features == 0:
            t_empty = time.perf_counter()
            self._write_empty_chunk_output(output_file)
            write_empty_cost = time.perf_counter() - t_empty

            total_cost = time.perf_counter() - t_total

            return total_features

        t0 = time.perf_counter()
        df = pd.read_csv(output_file, sep='\t', header=0)
        read_csv_cost = time.perf_counter() - t0

        if df.empty:
            t_empty = time.perf_counter()
            self._write_empty_chunk_output(output_file)
            write_empty_cost = time.perf_counter() - t_empty

            total_cost = time.perf_counter() - t_total

            return total_features

        t0 = time.perf_counter()
        df = df.set_index(['#chrom', 'pos', 'ref', 'alt'])
        set_index_cost = time.perf_counter() - t0

        parquet_file = str(output_file).replace('.txt', '.parquet')

        t0 = time.perf_counter()
        df.to_parquet(parquet_file, index=True, engine='pyarrow', compression='snappy')
        write_parquet_cost = time.perf_counter() - t0

        total_cost = time.perf_counter() - t_total
        return total_features

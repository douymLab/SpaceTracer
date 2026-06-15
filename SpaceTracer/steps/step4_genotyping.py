#!/usr/bin/env python3
"""
Genotyping Step - DataFrame optimized version
"""

from typing import Dict
import os
from functools import partial
import gc

import pandas as pd
import psutil

from SpaceTracer.steps.base import BaseStep
from SpaceTracer.utils.logger import get_logger
from SpaceTracer.utils.parallel import parallel_map
from SpaceTracer.utils.utils import str2bool, load_manifest_tsv, save_manifest_tsv

from SpaceTracer.cores.genotyping_01_combine import (
    UMICombiner_from_spot,
    ClusterAlleleFilter,
    UMICombiner_from_cluster,
)
from SpaceTracer.cores.genotyping_02_ind_genotype import (
    ClusterVAFCalculator_from_df,
    IndGenoCalculator,
)
from SpaceTracer.cores.genotyping_03_spot_genotype import SpotGenoCalculator


def _fmt_bytes(n):
    if n is None:
        return "NA"
    units = ["B", "KB", "MB", "GB", "TB"]
    n = float(n)
    for u in units:
        if n < 1024 or u == units[-1]:
            return f"{n:.2f} {u}"
        n /= 1024


def profile_df(name, df):
    proc = psutil.Process(os.getpid())
    rss = proc.memory_info().rss

    if df is None:
        print(f"[PROFILE] {name}: df=None | RSS={_fmt_bytes(rss)}", flush=True)
        return

    if not isinstance(df, pd.DataFrame):
        print(f"[PROFILE] {name}: type={type(df)} | RSS={_fmt_bytes(rss)}", flush=True)
        return

    mem = df.memory_usage(deep=True).sum()
    print(
        f"[PROFILE] {name}: shape={df.shape}, df_mem={_fmt_bytes(mem)}, RSS={_fmt_bytes(rss)}",
        flush=True
    )


def profile_obj(name, obj):
    proc = psutil.Process(os.getpid())
    rss = proc.memory_info().rss
    try:
        n = len(obj)
    except Exception:
        n = "NA"

    print(
        f"[PROFILE] {name}: type={type(obj)}, len={n}, RSS={_fmt_bytes(rss)}",
        flush=True
    )


def profile_gc(tag=""):
    proc = psutil.Process(os.getpid())
    before = proc.memory_info().rss
    gc.collect()
    after = proc.memory_info().rss
    # print(
    #     f"[PROFILE] gc {tag}: RSS before={_fmt_bytes(before)}, after={_fmt_bytes(after)}",
    #     flush=True
    # )


model_name = __name__
logger = get_logger(model_name)


# [ADDED]
def _empty_prior_df():
    return pd.DataFrame(columns=["identifier", "fA", "fT", "fC", "fG"])


# [ADDED]
def _load_prior_region(prior_file, chrom, start_pos, end_pos, chunksize=200000):
    if prior_file == "" or prior_file is None:
        return _empty_prior_df()

    chrom = str(chrom)
    start_pos = int(start_pos)
    end_pos = int(end_pos)

    cols = ["chrom", "pos", "ref", "fA", "fT", "fC", "fG"]
    parts = []
    entered_target_chrom = False
    chunk_counter = 0
    hit_rows = 0

    reader = pd.read_csv(
        prior_file,
        sep="\t",
        header=None,
        comment="#",
        names=cols,
        chunksize=chunksize,
        dtype={
            "chrom": "string",
            "pos": "int64",
            "ref": "string",
            "fA": "float64",
            "fT": "float64",
            "fC": "float64",
            "fG": "float64",
        }
    )

    for chunk_df in reader:
        chunk_counter += 1

        chrom_series = chunk_df["chrom"]
        has_target = (chrom_series == chrom).any()

        if not entered_target_chrom:
            if not has_target:
                continue
            entered_target_chrom = True
        else:
            if not has_target:
                break

        sub = chunk_df[chrom_series == chrom]
        if sub.empty:
            continue

        chunk_min_pos = int(sub["pos"].iloc[0])
        chunk_max_pos = int(sub["pos"].iloc[-1])

        if chunk_max_pos < start_pos:
            continue

        if chunk_min_pos > end_pos:
            break

        sub = sub[(sub["pos"] >= start_pos) & (sub["pos"] <= end_pos)]
        if not sub.empty:
            parts.append(sub)
            hit_rows += len(sub)

        if chunk_max_pos > end_pos:
            break

    if not parts:
        print(
            f"[PROFILE] prior region chrom={chrom} start={start_pos} end={end_pos} "
            f"read_chunks={chunk_counter} hit_rows=0",
            flush=True
        )
        return _empty_prior_df()

    df = pd.concat(parts, ignore_index=True)
    df["identifier"] = df["chrom"].astype(str) + "_" + df["pos"].astype(str)
    df = df[["identifier", "fA", "fT", "fC", "fG"]]

    print(
        f"[PROFILE] prior region chrom={chrom} start={start_pos} end={end_pos} "
        f"read_chunks={chunk_counter} hit_rows={hit_rows}",
        flush=True
    )
    return df


class GenotypingStep(BaseStep):
    def get_inputs(self, context: Dict) -> Dict[str, str]:
        spot_count_file = context.get("spot_count_file", "")
        prior_file = context.get("prior_file")
        manifest_path = context.get("manifest_path")
        return {
            "manifest_path": manifest_path,
            "spot_count_file": spot_count_file,
            "prior_file": prior_file,
        }

    def optional_parameters(self, context: Dict) -> Dict[str, str]:
        related_files = {}
        cluster_file = context.get("cluster_file", "")
        seq_type = self.config.get("sequence_type")

        cluster_df = self._load_cluster_file(cluster_file, seq_type)
        related_files["cluster"] = cluster_df
        return related_files

    def _load_cluster_file(self, cluster_file: str, seq_type: str) -> pd.DataFrame:
        if not os.path.exists(cluster_file):
            return pd.DataFrame()

        def format_cluster(val):
            if pd.isnull(val):
                return "NA"
            if isinstance(val, float) and val.is_integer():
                return str(int(val))
            return str(val)

        if seq_type == "stereo":
            df = pd.read_csv(
                cluster_file,
                sep="\t",
                header=None,
                names=["x", "y", "cluster"],
                na_values=[]
            )
            df["barcode"] = df["x"].astype(str) + "_" + df["y"].astype(str)
            df = df.drop(columns=["x", "y"])

        elif seq_type in ["visium", "visium-HD"]:
            df = pd.read_csv(cluster_file, sep="\t", header=0, na_values=[])

            if len(df.columns) == 1 or "cluster" not in df.columns:
                df = pd.read_csv(
                    cluster_file,
                    sep="\t",
                    header=None,
                    names=["barcode", "cluster"],
                    na_values=[]
                )
            else:
                df = df.rename(columns={df.columns[0]: "barcode"})
        else:
            return pd.DataFrame()

        df["cluster"] = df["cluster"].apply(format_cluster)
        return df

    def _get_chunk_output_paths(self, context: Dict, chunk: str) -> Dict[str, str]:
        chunk_dir = os.path.join(self.step_dir, chunk)
        os.makedirs(chunk_dir, exist_ok=True)

        return {
            "cluster_count_file": os.path.join(chunk_dir, "cluster_count.out"),
            "cluster_count_filter_file": os.path.join(chunk_dir, "cluster_count_filter.out"),
            "ind_count_filter_file": os.path.join(chunk_dir, "ind_count_filter.out"),
            "ind_geno_file": os.path.join(chunk_dir, "ind_geno.out"),
            "ind_geno_filter_file": os.path.join(chunk_dir, "ind_geno_filter.out"),
            "ind_geno_filter_mutation_list": os.path.join(chunk_dir, "ind_geno_filter.out.mutation.list"),
            "germline_file": os.path.join(chunk_dir, "germline.out"),
            "cluster_vaf_file": os.path.join(chunk_dir, "cluster_vaf.out"),
            "spot_geno_file": os.path.join(chunk_dir, "spot_geno.parquet"),
        }

    def get_outputs(self, context: Dict) -> Dict[str, str]:
        out_dir = self.step_dir
        os.makedirs(out_dir, exist_ok=True)
        return {
            "genotype_results": os.path.join(out_dir, "genotyping_chunk_manifest.tsv")
        }

    def _run_one_chunk(
        self,
        row: Dict[str, str],
        cluster_df,
        prior_file,  # [CHANGED] 原 prior_df -> prior_file
        epsQ: int,
        alpha: float,
        epsAF: float,
        mu: float,
        thr_dp: int,
        pop_vaf: float,
        filter_oneallele: bool,
        cell_num,
        bins,
        ind_geno_workers: int,
        spot_geno_workers: int,
        context: Dict,
    ) -> Dict[str, str]:
        chunk = row["chunk"]
        spot_count = row["parquet_file"]

        if not os.path.exists(spot_count):
            raise FileNotFoundError(f"spot_count parquet not found: {spot_count}")

        outputs = self._get_chunk_output_paths(context, chunk)

        cluster_count_file = outputs["cluster_count_file"]
        cluster_count_filter_file = outputs["cluster_count_filter_file"]
        ind_count_filter_file = outputs["ind_count_filter_file"]
        ind_geno_file = outputs["ind_geno_file"]
        ind_geno_filter_file = outputs["ind_geno_filter_file"]
        germline_file = outputs["germline_file"]
        cluster_vaf_file = outputs["cluster_vaf_file"]
        spot_geno_file = outputs["spot_geno_file"]

        combiner_spot = UMICombiner_from_spot(epsQ)
        allele_filter = ClusterAlleleFilter(alpha, epsAF)
        combiner_cluster = UMICombiner_from_cluster(epsQ)
        ind_geno_calc = IndGenoCalculator()

        # 1. load + cluster-level combine
        spot_df = combiner_spot.load_df(spot_count)
        if self.seq_type=="stereo":
            cluster_count_df = combiner_spot.combine_bins_df(spot_df, cluster_df)
        else:
            cluster_count_df = combiner_spot.combine_cluster_df(spot_df, cluster_df)
        combiner_spot.save_df(cluster_count_df, cluster_count_file)

        # 2. allele filter
        cluster_count_filter_df = allele_filter.filter_df(cluster_count_df)

        cluster_count_filter_out = cluster_count_filter_df.rename(columns={"chrom": "#chrom"})
        cluster_count_filter_out.to_csv(cluster_count_filter_file, sep="\t", index=None, na_rep="NA")
        del cluster_count_filter_out

        del cluster_count_df

        # 3. individual-level combine
        ind_count_filter_df = combiner_cluster.combine_ind_df(cluster_count_filter_df)

        ind_count_filter_out = ind_count_filter_df.rename(columns={"chrom": "#chrom"})
        ind_count_filter_out.to_csv(ind_count_filter_file, sep="\t", index=None, na_rep="NA")
        del ind_count_filter_out

        # 4. genotype
        chunk_chrom = row["chrom"]
        chunk_start = int(row["start_pos"])
        chunk_end = int(row["end_pos"])

        prior_sub_df = _load_prior_region(
            prior_file=prior_file,
            chrom=chunk_chrom,
            start_pos=chunk_start,
            end_pos=chunk_end
        )

        geno_df, geno_filter_df, germ_df = ind_geno_calc.calculate_individual_genotype_df(
            ind_count_df=ind_count_filter_df,
            prior_df=prior_sub_df,  # [CHANGED]
            geno_file=ind_geno_file,
            geno_filter_file=ind_geno_filter_file,
            germline_file=germline_file,
            mu=mu,
            thr_dp=thr_dp,
            pop_vaf=pop_vaf,
            filter_oneallele=filter_oneallele,
            max_workers=ind_geno_workers
        )

        del ind_count_filter_df

        has_data = not geno_filter_df.empty

        if has_data:
            cluster_count_filter_for_vaf = cluster_count_filter_df.rename(columns={"chrom": "#chrom"})

            cluster_vaf_df = ClusterVAFCalculator_from_df(
                geno_filter_df,
                cluster_count_filter_for_vaf,
                outfile_path=cluster_vaf_file
            )

            del cluster_count_filter_for_vaf

            SpotGenoCalculator(
                bins,
                epsQ,
                thr_dp,
                pop_vaf,
                cell_num,
                max_workers=spot_geno_workers
            ).run_from_df(
                spot_count_df=spot_df,
                ind_geno_df=geno_filter_df,
                cluster_df=cluster_df,
                cluster_vaf_df=cluster_vaf_df,
                output_file=spot_geno_file
            )

            del cluster_vaf_df
            del cluster_count_filter_df
            del geno_df
            del geno_filter_df
            del germ_df
            del spot_df
            del prior_sub_df  # [ADDED]
            profile_gc(f"{chunk} end has_data")

            return {
                "chunk": chunk,
                "spot_count_file": spot_count,
                **outputs,
            }

        else:
            del cluster_count_filter_df
            del geno_df
            del geno_filter_df
            del germ_df
            del spot_df
            del prior_sub_df  # [ADDED]
            return {}

    def _run(self, context: Dict):
        related_files = self.optional_parameters(context)
        cluster_df = related_files["cluster"]

        parameters = self.get_step_config()
        alpha = float(parameters["alpha"])
        epsQ = int(parameters["epsQ"])
        epsAF = float(parameters["epsAF"])
        mu = float(parameters["mu"])
        thr_dp = int(parameters["thr_dp"])
        pop_vaf = float(parameters["pop_vaf"])
        filter_oneallele = str2bool(parameters["filter_oneallele"])

        inputs = self.get_inputs(context)
        spot_count_manifest = inputs["spot_count_file"]
        prior_file = inputs["prior_file"]

        manifest_df = pd.read_csv(inputs["manifest_path"], sep="\t")
        manifest_df = manifest_df.rename(columns={"chunk_id": "chunk"})  # [CHANGED]

        if not os.path.exists(spot_count_manifest):
            raise FileNotFoundError(f"spot_count manifest not found: {spot_count_manifest}")

        if isinstance(prior_file, int) and prior_file == 0:
            prior_file = ""
        elif isinstance(prior_file, str):
            if prior_file != "" and not os.path.exists(prior_file):
                raise FileNotFoundError(f"prior_file not found: {prior_file}")
        else:
            raise RuntimeError(f"Wrong prior input: {prior_file}")

        rows = load_manifest_tsv(spot_count_manifest)
        if not rows:
            raise ValueError(f"No chunk records found in manifest: {spot_count_manifest}")

        # [ADDED] merge manifest region info into rows
        rows_df = pd.DataFrame(rows)
        rows_df = rows_df.merge(
            manifest_df[["chunk", "chrom", "start_pos", "end_pos"]],
            on="chunk",
            how="left"
        )

        missing_region = rows_df[["chrom", "start_pos", "end_pos"]].isna().any(axis=1)
        if missing_region.any():
            bad_chunks = rows_df.loc[missing_region, "chunk"].tolist()
            raise ValueError(f"Missing region info for chunks: {bad_chunks[:10]}")

        rows = rows_df.to_dict("records")

        outputs = self.get_outputs(context)
        result_manifest = outputs["genotype_results"]

        cell_num = context.get("cell_num")
        seq_type = self.config.get("sequence_type")
        if seq_type == "stereo":
            bins = self.config.get("bin_size", 100)
        else:
            bins = None

        max_workers = self.config.get("runtime", {}).get("max_parallel", self.threads)
        parallel_backend = self.config.get("runtime", {}).get("parallel_backend", "process")

        ind_geno_workers = 1
        spot_geno_workers = 1

        worker = partial(
            self._run_one_chunk,
            cluster_df=cluster_df,
            prior_file=prior_file,  
            epsQ=epsQ,
            alpha=alpha,
            epsAF=epsAF,
            mu=mu,
            thr_dp=thr_dp,
            pop_vaf=pop_vaf,
            filter_oneallele=filter_oneallele,
            cell_num=cell_num,
            bins=bins,
            ind_geno_workers=ind_geno_workers,
            spot_geno_workers=spot_geno_workers,
            context=context,
        )

        chunk_results = parallel_map(
            rows,
            worker_fn=worker,
            max_workers=max_workers,
            desc="genotyping",
            raise_on_error=True,
            backend=parallel_backend,
        )

        save_manifest_tsv(chunk_results, result_manifest)

    def get_step_config(self) -> Dict:
        return self.config.get("steps", {}).get("genotyping", {})


# [KEPT] 保留原函数，方便回滚；当前版本不再使用它
def _load_prior_file(file):
    if file == "":
        df = pd.DataFrame()
    else:
        df = pd.read_csv(file, sep="\t", header=None, comment="#")
        columns = ["chrom", "pos", "ref", "fA", "fT", "fC", "fG"]
        df.columns = columns
        df["identifier"] = df["chrom"] + "_" + df["pos"].astype(str)
        df = df[["identifier", "fA", "fT", "fC", "fG"]]
    return df

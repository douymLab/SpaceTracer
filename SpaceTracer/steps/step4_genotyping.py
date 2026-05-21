#!/usr/bin/env python3
"""
Genotyping Step - DataFrame optimized version
"""

from typing import Dict
import os
from functools import partial

import pandas as pd

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


model_name = __name__
logger = get_logger(model_name)


class GenotypingStep(BaseStep):
    def get_inputs(self, context: Dict) -> Dict[str, str]:
        spot_count_file = context.get('spot_count_file', '')
        prior_file = context.get('prior_file')
        return {
            'spot_count_file': spot_count_file,
            'prior_file': prior_file
        }

    def optional_parameters(self, context: Dict) -> Dict[str, str]:
        related_files = {}

        cluster_file = context.get('cluster_file', '')
        if os.path.exists(cluster_file):
            cluster_df = pd.read_csv(cluster_file, sep="\t", header=0, na_values=[])

            if len(cluster_df.columns) == 1 or 'cluster' not in cluster_df.columns:
                cluster_df = pd.read_csv(
                    cluster_file,
                    sep="\t",
                    header=None,
                    names=['barcode', 'cluster'],
                    na_values=[]
                )
            else:
                cluster_df = cluster_df.rename(columns={cluster_df.columns[0]: 'barcode'})

            cluster_df['cluster'] = cluster_df['cluster'].apply(
                lambda x: str(int(x)) if isinstance(x, float) and x.is_integer()
                else str(x) if pd.notnull(x) else "NA"
            )
        else:
            cluster_df = pd.DataFrame()

        related_files["cluster"] = cluster_df
        return related_files

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
            "spot_geno_file": os.path.join(chunk_dir, "spot_geno.out"),
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
        prior_df,
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
        cluster_count_df = combiner_spot.combine_cluster_df(spot_df, cluster_df)
        combiner_spot.save_df(cluster_count_df, cluster_count_file)

        # 2. allele filter
        cluster_count_filter_df = allele_filter.filter_df(cluster_count_df)
        cluster_count_filter_out = cluster_count_filter_df.rename(columns={'chrom': '#chrom'})
        cluster_count_filter_out.to_csv(cluster_count_filter_file, sep="\t", index=None, na_rep='NA')

        # 3. individual-level combine
        ind_count_filter_df = combiner_cluster.combine_ind_df(cluster_count_filter_df)
        ind_count_filter_out = ind_count_filter_df.rename(columns={'chrom': '#chrom'})
        ind_count_filter_out.to_csv(ind_count_filter_file, sep="\t", index=None, na_rep='NA')

        # 4. genotype
        geno_df, geno_filter_df, germ_df = ind_geno_calc.calculate_individual_genotype_df(
            ind_count_df=ind_count_filter_df,
            prior_df=prior_df,
            geno_file=ind_geno_file,
            geno_filter_file=ind_geno_filter_file,
            germline_file=germline_file,
            mu=mu,
            thr_dp=thr_dp,
            pop_vaf=pop_vaf,
            filter_oneallele=filter_oneallele,
            max_workers=ind_geno_workers
        )

        has_data = not geno_filter_df.empty

        if has_data:
            cluster_count_filter_for_vaf = cluster_count_filter_df.rename(columns={'chrom': '#chrom'})
            cluster_vaf_df = ClusterVAFCalculator_from_df(
                geno_filter_df,
                cluster_count_filter_for_vaf,
                outfile_path=cluster_vaf_file
            )

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
            return {
                "chunk": chunk,
                "spot_count_file": spot_count,
                **outputs,
            }
        else:
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

        if not os.path.exists(spot_count_manifest):
            raise FileNotFoundError(f"spot_count manifest not found: {spot_count_manifest}")

        if isinstance(prior_file,int) and prior_file==0:
            prior_df=pd.DataFrame()
        elif isinstance(prior_file,str) and not os.path.exists(prior_file):
            raise FileNotFoundError(f"prior_file not found: {prior_file}")
        elif isinstance(prior_file,str) and os.path.exists(prior_file):
            prior_df = _load_prior_file(prior_file)
        else:
            raise RuntimeError(f"Wrong prior input: {prior_file}")

        rows = load_manifest_tsv(spot_count_manifest)
        if not rows:
            raise ValueError(f"No chunk records found in manifest: {spot_count_manifest}")

        outputs = self.get_outputs(context)
        result_manifest = outputs["genotype_results"]

        cell_num = self.config.get("cell_num")
        bins = self.config.get("bin_size")
        runtime_cfg = self.config.get("runtime", {})
        max_workers = runtime_cfg.get("max_parallel", self.threads)
        parallel_backend = runtime_cfg.get("parallel_backend", "thread")

        # Unified worker control for genotyping internals.
        # Backward-compatible behavior:
        # 1) explicit per-stage workers override
        # 2) shared geno_workers applies to both stages
        # 3) fallback to 1
        shared_geno_workers = int(runtime_cfg.get("geno_workers", 1))
        ind_geno_workers = int(runtime_cfg.get("ind_geno_workers", shared_geno_workers))
        if ind_geno_workers < 1:
            ind_geno_workers = 1
        spot_geno_workers = int(runtime_cfg.get("spot_geno_workers", shared_geno_workers))
        if spot_geno_workers < 1:
            spot_geno_workers = 1

        worker = partial(
            self._run_one_chunk,
            cluster_df=cluster_df,
            prior_df=prior_df,
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
        return self.config.get('steps', {}).get('genotyping', {})


def _load_prior_file(file):
    if file == "":
        df = pd.DataFrame()
    else:
        df = pd.read_csv(file, sep='\t', header=None, comment="#")
        columns = ["chrom", "pos", "ref", "fA", "fT", "fC", "fG"]
        df.columns = columns
        df['identifier'] = df['chrom'] + "_" + df['pos'].astype(str)
        df = df[["identifier", 'fA', 'fT', 'fC', 'fG']]
    return df

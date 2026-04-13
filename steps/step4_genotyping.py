#!/usr/bin/env python3
"""
Genotyping Step - This step will help you to get candidate sites by mpileup
"""

from typing import Dict
import os
import time
from functools import partial
import pandas as pd
from tqdm import tqdm
from SpaceTracer.cores.genotyping_02_ind_genotype import ClusterVAFCalculator, IndGenoCalculator
from SpaceTracer.cores.genotyping_03_spot_genotype import SpotGenoCalculator
from SpaceTracer.steps.base import BaseStep
from SpaceTracer.utils.logger import get_logger
from SpaceTracer.utils.parallel import parallel_map
from  SpaceTracer.utils.utils import str2bool,load_manifest_tsv,save_manifest_tsv
from SpaceTracer.cores.genotyping_01_combine import UMICombiner_from_spot,ClusterAlleleFilter,UMICombiner_from_cluster

model_name=__name__
logger = get_logger(model_name)

class GenotypingStep(BaseStep):

    def get_inputs(self, context: Dict) -> Dict[str, str]:
        """ input """
        spot_count_file=context.get('spot_count_file','')
        # spot_count_parquet=context.get('spot_count_parquet','')
        prior_file=context.get('prior_file')

        # if os.path.exists(spot_count_parquet):
        #     spot_count=spot_count_parquet
        # elif os.path.exists(spot_count_file):
        #         spot_count=spot_count_file
        # else: 
        #     raise FileNotFoundError(f'{spot_count_parquet}')
        inputs = {
            'spot_count_file': spot_count_file,
            'prior_file':prior_file
        }
        return inputs
    

    # def optional_parameters(self, context: Dict) -> Dict[str, str]:
    #     """ That's optional parameters """
    #     related_files={}

    #     cluster_file=context.get('cluster_file','')
    #     if os.path.exists(cluster_file):
    #         # cluster_df= pd.read_csv(cluster_file, sep="\t", header=None, names=['barcode', 'cluster'], na_values=[])
    #         cluster_df= pd.read_csv(cluster_file, sep="\t", header=0, na_values=[])
    #         cluster_df=cluster_df.rename(columns={'index':'barcode'})

    #         cluster_df['cluster'] = cluster_df['cluster'].apply(lambda x: str(int(x)) if isinstance(x, float) and x.is_integer() 
    #                                                         else str(x) if pd.notnull(x) else "NA")
    #     else:
    #         cluster_df=pd.DataFrame()
        
    #     related_files["cluster"]=cluster_df

    #     return related_files

    def optional_parameters(self, context: Dict) -> Dict[str, str]:
        related_files={}
        
        cluster_file=context.get('cluster_file','')
        if os.path.exists(cluster_file):
            # 尝试有表头读取
            cluster_df = pd.read_csv(cluster_file, sep="\t", header=0, na_values=[])
            
            # 如果读取后只有1列或列名不是预期的，尝试无表头读取
            if len(cluster_df.columns) == 1 or 'cluster' not in cluster_df.columns:
                cluster_df = pd.read_csv(cluster_file, sep="\t", header=None, 
                                        names=['barcode', 'cluster'], na_values=[])
            else:
                # 重命名第一列为barcode
                cluster_df = cluster_df.rename(columns={cluster_df.columns[0]: 'barcode'})
            
            # 统一处理cluster列
            cluster_df['cluster'] = cluster_df['cluster'].apply(lambda x: str(int(x)) if isinstance(x, float) and x.is_integer() 
                                                            else str(x) if pd.notnull(x) else "NA")
        else:
            cluster_df=pd.DataFrame()
        
        related_files["cluster"]=cluster_df
        
        return related_files

        
    def _load_prior_file(self,file):
        if file=="":
            df=pd.DataFrame()
        else:
            df = pd.read_csv(file,sep='\t',header=None,comment="#")
            columns=["chrom","pos","ref","fA","fT","fC","fG"]
            df.columns=columns
            df['identifier'] = df['chrom']+"_"+df['pos'].astype(str)
            df=df.set_index('identifier')
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
            "ind_geno_filter_mutation_list" : os.path.join(chunk_dir, "ind_geno_filter.out.mutation.list"),
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
        prior_file: str,
        epsQ: int,
        alpha: float,
        epsAF: float,
        mu: float,
        thr_dp: int,
        pop_vaf: float,
        filter_oneallele: bool,
        cell_num,
        bins,
        context: Dict,
    ) -> Dict[str, str]:
        t0 = time.perf_counter()
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

        # 1. cluster-level combine
        t1 = time.perf_counter()
        UMICombiner_from_spot(epsQ).combine_cluster(
            spot_count,
            cluster_count_file,
            cluster_df
        )
        t1 = time.perf_counter() - t1

        # 2. allele filter
        t2 = time.perf_counter()
        ClusterAlleleFilter(alpha, epsAF).filter(
            cluster_count_file,
            cluster_count_filter_file
        )
        t2 = time.perf_counter() - t2

        # 3. individual-level combine
        t3 = time.perf_counter()
        UMICombiner_from_cluster(epsQ).combine_ind(
            cluster_count_filter_file,
            ind_count_filter_file
        )
        t3 = time.perf_counter() - t3

        # 4. genotype
        t4 = time.perf_counter()
        IndGenoCalculator().calculate_individual_genotype(
            ind_count_filter_file,
            prior_file,
            ind_geno_file,
            ind_geno_filter_file,
            germline_file,
            mu,
            thr_dp,
            pop_vaf,
            filter_oneallele
        )
        t4 = time.perf_counter() - t4

        with open(ind_geno_filter_file, "r") as f:
            row_length = sum(1 for _ in f)

        has_data = row_length > 1
        if has_data:
            t5 = time.perf_counter()
            ClusterVAFCalculator(
                ind_geno_filter_file,
                cluster_count_filter_file,
                cluster_vaf_file
            )
            t5 = time.perf_counter() - t5

            # 6. spot genotype
            t6 = time.perf_counter()
            SpotGenoCalculator(
                bins,
                epsQ,
                thr_dp,
                pop_vaf,
                cell_num
            ).run(
                spot_count,
                ind_geno_filter_file,
                cluster_df,
                cluster_vaf_file,
                spot_geno_file
            )
            t6 = time.perf_counter() - t6
            total = time.perf_counter() - t0
            logger.info(
                "[genotyping chunk=%s] combine_cluster=%.2fs allele_filter=%.2fs "
                "combine_ind=%.2fs ind_genotype=%.2fs cluster_vaf=%.2fs "
                "spot_genotype=%.2fs total=%.2fs",
                chunk, t1, t2, t3, t4, t5, t6, total
            )
            return {
                "chunk": chunk,
                "spot_count_file": spot_count,
                **outputs,
            }
        else:
            total = time.perf_counter() - t0
            logger.info(
                "[genotyping chunk=%s] combine_cluster=%.2fs allele_filter=%.2fs "
                "combine_ind=%.2fs ind_genotype=%.2fs no-spot-output total=%.2fs",
                chunk, t1, t2, t3, t4, total
            )
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
        if not os.path.exists(prior_file):
            raise FileNotFoundError(f"prior_file not found: {prior_file}")

        rows = load_manifest_tsv(spot_count_manifest)
        if not rows:
            raise ValueError(f"No chunk records found in manifest: {spot_count_manifest}")

        outputs = self.get_outputs(context)
        result_manifest = outputs["genotype_results"]

        cell_num = self.config.get("cell_num")
        bins = self.config.get("bin_size")
        max_workers = self.config.get("runtime", {}).get("max_parallel", self.threads)
        parallel_backend = self.config.get("runtime", {}).get("parallel_backend", "thread")

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
            context=context,
        )

        chunk_results = parallel_map(
            rows,
            worker_fn=worker,
            max_workers=max_workers,
            desc=f"genotyping",
            raise_on_error=True,
            backend=parallel_backend,
        )
        
        save_manifest_tsv(chunk_results, result_manifest)


    
    def get_step_config(self) -> Dict:
        return self.config.get('steps', {}).get('genotyping', {})
    

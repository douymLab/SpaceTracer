#!/usr/bin/env python3
"""
Genotyping Step - This step will help you to get candidate sites by mpileup
"""

from typing import Dict
import os
import pandas as pd
from tqdm import tqdm
from SpaceTracer.cores.genotyping_02_ind_genotype import ClusterVAFCalculator, IndGenoCalculator
from SpaceTracer.cores.genotyping_03_spot_genotype import SpotGenoCalculator
from SpaceTracer.steps.base import BaseStep
from SpaceTracer.utils.logger import get_logger

from SpaceTracer.cores.genotyping_01_combine import UMICombiner_from_spot,ClusterAlleleFilter,UMICombiner_from_cluster

model_name=__name__
logger = get_logger(model_name)

class GenotypingStep(BaseStep):

    def get_inputs(self, context: Dict) -> Dict[str, str]:
        """ input """
        spot_count_file=context.get('spot_count_file','')
        spot_count_parquet=context.get('spot_count_parquet','')
        prior_file=context.get('prior_file')
        if os.path.exists(spot_count_parquet):
            spot_count=spot_count_parquet
        elif os.path.exists(spot_count_file):
                spot_count=spot_count_file
        inputs = {
            'spot_count': spot_count,
            'prior_file':prior_file
        }
        return inputs
    

    def optional_parameters(self, context: Dict) -> Dict[str, str]:
        """ That's optional parameters """
        related_files={}

        cluster_file=self.config.get('cluster','')
        if os.path.exists(cluster_file):
            cluster_df= pd.read_csv(cluster_file, sep="\t", header=None, names=['barcode', 'cluster'], na_values=[])
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

    def get_outputs(self,context: Dict) -> Dict[str, str]:
        """output"""
        if self.config.get('run').get('keep_intermediates'):
            cluster_count_file = os.path.join(self.step_dir, "cluster.count.out")
            cluster_count_filter_file = os.path.join(self.step_dir, "cluster_filter.count.out")
            ind_count_filter_file = os.path.join(self.step_dir, "ind_filter.count.out")
            ind_geno_file = os.path.join(self.step_dir, "ind_genotype.out")
            ind_geno_filter_file = os.path.join(self.step_dir, "ind_genotype.filter.out")
            ind_geno_filter_mutation_list = os.path.join(self.step_dir, "ind_genotype.filter.out.mutation.list")
            germline_file = os.path.join(self.step_dir, "germline.out")
            cluster_vaf_file = os.path.join(self.step_dir, "cluster_vaf.out")
            spot_geno_file = os.path.join(self.step_dir, "spot_genotype.out")
        
            return {
                'cluster_count_file': cluster_count_file,
                'cluster_count_filter_file': cluster_count_filter_file,
                'ind_count_filter_file': ind_count_filter_file,
                'ind_geno_file': ind_geno_file,
                'ind_geno_filter_file': ind_geno_filter_file,
                'ind_geno_filter_mutation_list': ind_geno_filter_mutation_list,
                'germline_file':germline_file,
                'cluster_vaf_file': cluster_vaf_file,
                'spot_geno_file': spot_geno_file
            }
        else:
            ind_geno_filter_file = os.path.join(self.step_dir, "ind_genotype.filter.out")
            ind_geno_filter_mutation_list = os.path.join(self.step_dir, "ind_genotype.filter.out.mutation.list")
            germline_file = os.path.join(self.step_dir, "germline.out")
            cluster_vaf_file = os.path.join(self.step_dir, "cluster_vaf.out")
            return {
                'ind_geno_filter_file': ind_geno_filter_file,
                'ind_geno_filter_mutation_list': ind_geno_filter_mutation_list,
                'germline_file':germline_file,
                'cluster_vaf_file': cluster_vaf_file,
                'spot_geno_file': spot_geno_file

            }


    
    def get_step_config(self) -> Dict:
        return self.config.get('steps', {}).get('genotyping', {})
    

    def _run(self, context: Dict):
        # parameters:
        related_files=self.optional_parameters(context)
        cluster_df=related_files['cluster']

        parameters=self.get_step_config()
        alpha=float(parameters['alpha'])
        epsQ=int(parameters['epsQ'])
        epsAF=float(parameters['epsAF'])
        mu=float(parameters['mu'])
        thr_dp=int(parameters['thr_dp'])
        pop_vaf=float(parameters['pop_vaf'])
        filter_oneallele=bool(parameters['filter_oneallele'])

        inputs=self.get_inputs(context)
        spot_count=inputs['spot_count']
        prior_file=inputs['prior_file']

        outputs=self.get_outputs(context)
        cluster_count_file = outputs["cluster_count_file"]
        cluster_count_filter_file =outputs["cluster_count_filter_file"]
        ind_count_filter_file=outputs["ind_count_filter_file"]
        ind_geno_file=outputs["ind_geno_file"]
        ind_geno_filter_file=outputs["ind_geno_filter_file"]
        germline_file=outputs["germline_file"]
        cluster_vaf_file=outputs["cluster_vaf_file"]
        spot_geno_file=outputs['spot_geno_file']

        # prior_df=self._load_prior_file(prior_file)

        # #genotype in ind level
        # cluster_count_df=UMICombiner_from_spot(epsQ).combine_cluster(spot_count, cluster_count_file, cluster_df)
        # cluster_allele_filter_df=ClusterAlleleFilter(alpha,epsAF).filter(cluster_count_df, cluster_count_filter_file)
        # ind_count_filter_df=UMICombiner_from_cluster(epsQ).combine_ind(cluster_count_filter_file,ind_count_filter_file)
        # ind_geno_filter_df=IndGenoCalculator().calculate_individual_genotype(ind_count_filter_df,prior_df, ind_geno_file,ind_geno_filter_file,germline_file, mu, thr_dp, pop_vaf, filter_oneallele)
        # print(spot_count)
        UMICombiner_from_spot(epsQ).combine_cluster(spot_count, cluster_count_file, cluster_df)
        ClusterAlleleFilter(alpha,epsAF).filter(cluster_count_file, cluster_count_filter_file)
        UMICombiner_from_cluster(epsQ).combine_ind(cluster_count_filter_file,ind_count_filter_file)
        IndGenoCalculator().calculate_individual_genotype(ind_count_filter_file,prior_file, ind_geno_file,ind_geno_filter_file,germline_file, mu, thr_dp, pop_vaf, filter_oneallele)

        # vaf in cluster
        ClusterVAFCalculator(ind_geno_filter_file,cluster_count_filter_file,cluster_vaf_file)
        
        # genotype in spot level
        cell_num=self.config.get('cell_num')
        bins=self.config.get('bin_size')
        SpotGenoCalculator(bins,epsQ,thr_dp,pop_vaf,cell_num).run(spot_count,ind_geno_filter_file,cluster_df,cluster_vaf_file,spot_geno_file)


from functools import partial
import multiprocessing
import os
from pathlib import Path
import pandas as pd
from tqdm import tqdm

from SpaceTracer.cores.spatial_features import handle_per_line
from SpaceTracer.steps.base import BaseStep
from SpaceTracer.utils.read_files import handle_barcode, load_spot_genotypes_data
from SpaceTracer.utils.utils import check_dir


class SpatialFeatureStep(BaseStep):
    def get_inputs(self, context):
        inputs={'ind_geno_filter_file': context.get('ind_geno_filter_file',''),
                'ind_geno_filter_mutation_list': context.get('ind_geno_filter_mutation_list',''),
                'spot_geno_file':context.get('spot_geno_file','')
                }
        return inputs
    
    def get_outputs(self, context):
        return {'spatial_feature': os.path.join(self.step_dir,"spatial_feature.txt")}
    
    def get_step_config(self):
        return self.config.get('steps', {}).get('spatial_feature', {})
    
    def _run(self,context):
        tissue_positions=self.config["tissue_position"]
        barcode_dir=Path(self.step_dir)/"barcode_dir"
        check_dir(barcode_dir)
        threads=self.threads

        inputs=self.get_inputs(context)
        spot_genotype_file=inputs["spot_geno_file"]
        ind_geno_filter_mutation_list=inputs["ind_geno_filter_mutation_list"]

        outputs=self.get_outputs(context)
        out_spatial_features=outputs["spatial_feature"]

        parameters=self.get_step_config()
        alpha=float(parameters['alpha'])
        thr_r2=float(parameters['thr_r2'])
        thr_prob=float(parameters['thr_prob'])
        thr_likelihood=float(parameters['thr_likelihood'])
        thr_vaf=float(parameters['thr_vaf'])

        plot_supp=bool(parameters['plot_supp'])
        fig_size=int(parameters['fig_size'])
        method=str(parameters['method'])
        num_directions=int(parameters['num_directions'])
        in_name=None
        
        print(alpha,thr_r2,thr_prob,thr_likelihood,thr_vaf, \
                            plot_supp, fig_size,method,num_directions)
        mutation_identifier_list=pd.read_csv(ind_geno_filter_mutation_list,header=None).iloc[:, 0].tolist()
        spot_geno_df=load_spot_genotypes_data(spot_genotype_file)
        barcode_dict=handle_barcode(tissue_positions)

        # uncommant below running with multiprocessing
        partial_func=partial(handle_per_line,barcode_dir,self.step_dir,in_name,spot_geno_df,barcode_dict, \
                            alpha,thr_r2,thr_prob,thr_likelihood,thr_vaf, \
                            plot_supp, fig_size,method,num_directions)
        with multiprocessing.Pool(threads) as pool:
            results = list(tqdm(pool.imap(partial_func, mutation_identifier_list, chunksize=10), total=len(mutation_identifier_list)))


        colnames = [
            '#chrom', 'pos', 'ref', 'alt', 'test_sig',
            'early_mutation', 'late_mutation', 'verylate_mutation',
            'ks_stat', 'ks_pvalue',
            'moranI_stat', 'moranI_pvalue',
            'mutant_rate', 'mutant_rate_prob', 'mutant_rate_likelihood', 'mutant_rate_vaf',
            'mean_vaf', 'max_vaf',
            'r_squared', 'wilcoxon_stat', 'wilcoxon_pvalue',
            'outlier_clusters', 'outlier_vaf', 'outlier_moranI_stat', 'outlier_moranI_pvalue'
        ]
        
        with open(out_spatial_features, 'w') as f:
            f.write('\t'.join(colnames) + '\n')
            for values in results:
                if values:
                    f.write('\t'.join(values) + '\n')
        
        print(f"✓ TSV written: {out_spatial_features}")
        
        df = pd.read_csv(
            out_spatial_features,
            sep='\t',header=0
        )
        
        df['identifier'] = (
            df['#chrom'] + '_' +
            df['pos'].astype(str) + '_' +
            df['ref'] + '_' +
            df['alt']
        )
        df = df.set_index('identifier')
        
        parquet_file = str(out_spatial_features).replace('.txt', '.parquet')
        df.to_parquet(parquet_file, index=True, compression='snappy')



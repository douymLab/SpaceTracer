
import os
from pathlib import Path
from typing import Dict

import numpy as np
import pandas as pd

from SpaceTracer.cores.get_cellNum import get_cb_ub
from SpaceTracer.steps.base import BaseStep

from SpaceTracer.utils.logger import get_logger
model_name=__name__
logger = get_logger(model_name)


class CellNumStep(BaseStep):
    def get_inputs(self, context):
        inputs = {
            'in_filter_bam': context.get('in_filter_bam')
        }
        return inputs

    
    def get_outputs(self, context):
        cell_num=context.get('cell_num',0)
        if isinstance(cell_num,int):
            if cell_num==0: # 0 means we'll run cell number function.
                return {'cell_num': os.path.join(self.work_dir, "refined_cell_num.txt")}
            else: # other int means the data was cell/sub-cell level.
                return {'cell_num': cell_num}
        elif Path(cell_num).exists(): # also file is allowed
            return {'cell_num': cell_num}
        else:
            raise ValueError(f'Wrong cell number input {cell_num}')

    def optional_parameters(self, context: Dict) -> Dict[str, str]:
        """ That's optional parameters """
        
        cluster_file=self.context.get('cluster_file')
        if os.path.exists(cluster_file):
            cluster_df= pd.read_csv(cluster_file, sep="\t", header=None, names=['barcode', 'cluster'], na_values=[])
            cluster_df['cluster'] = cluster_df['cluster'].apply(lambda x: str(int(x)) if isinstance(x, float) and x.is_integer() 
                                                            else str(x) if pd.notnull(x) else "NA")
        else:
            cluster_df=pd.DataFrame()
        
        return cluster_df
    
    def _run(self,context):
        cell_num=context.get('cell_num',0)
        if isinstance(cell_num,int) and cell_num==0:
            save_dir=Path(self.work_dir)

            # input:
            bam_file=self.get_inputs(context)['in_filter_bam']

            #output:
            outputs=self.get_outputs(context)
            count_file=save_dir/"raw_cell_num.txt"
            cell_num_file=outputs['cell_num']

            #parameter:
            rerun=True
            seq_type=self.config.get('sequence_type')
            bins=self.config.get('bins')

            if rerun or not os.path.exists(count_file):
                cb_ub_df=get_cb_ub(bam_file,count_file,seq_type,bins)
            else:
                cb_ub_df = pd.read_csv(count_file, sep='\t', header=0)
            
            cluster_df=self.optional_parameters(context)

            if cluster_df.empty:                        
                file_merged = pd.merge(cluster_df, cb_ub_df, on='barcode')
            else:
                file_merged=cb_ub_df.copy()
                file_merged['cluster']='bulk'

            cluster_sums = file_merged.groupby('cluster')['nUMI'].sum()
            max_cluster = cluster_sums.idxmax()
            max_cluster_data = file_merged[file_merged['cluster'] == max_cluster]
            median_nUMI = max_cluster_data['nUMI'].median()
            file_merged['calculated_cell_num'] = np.ceil(20 * file_merged['nUMI'] / median_nUMI)
            file_merged['refined_cell_num'] = np.where(file_merged['calculated_cell_num'] > 25, 25, file_merged['calculated_cell_num'])

            final_output_df = file_merged[['barcode', 'cluster', 'nUMI', 'refined_cell_num']]
            final_output_df=final_output_df.rename(columns={'nUMI': 'UMI_counts'})

            # final_output_df = file_merged[['barcode', 'cluster', 'nUMI', 'nREAD', 'refined_cell_num']]
            final_output_df.to_csv(cell_num_file, index=False,sep="\t")
        else:
            pass


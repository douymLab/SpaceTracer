import os
from pathlib import Path

import numpy as np
import pandas as pd

from SpaceTracer.cores.get_cellNum import get_cb_ub
from SpaceTracer.steps.base import BaseStep
from SpaceTracer.utils.logger import get_logger

model_name = __name__
logger = get_logger(model_name)


class CellNumStep(BaseStep):
    def get_inputs(self, context):
        return {
            'in_filter_bam': context.get('in_filter_bam'),
            'cluster_file': context.get('cluster_file')
        }

    def get_outputs(self, context):
        cell_num = context.get('cell_num', self.config.get('cell_num', 0))

        if isinstance(cell_num, int):
            if cell_num == 0:
                return {'cell_num': os.path.join(self.work_dir, "refined_cell_num.txt")}
            else:
                return {'cell_num': cell_num}

        elif isinstance(cell_num, str) and Path(cell_num).exists():
            return {'cell_num': cell_num}

        else:
            raise ValueError(f'Wrong cell number input: {cell_num}')

    def _load_cluster_df(self, context) -> pd.DataFrame:
        """
        Load cluster file from context.
        If cluster file does not exist or is empty, return an empty dataframe.
        Empty dataframe means cell number estimation will run in bulk mode.
        """
        cluster_file = context.get('cluster_file')

        if not cluster_file:
            logger.info("No cluster_file found in context. Cell number will be estimated in bulk mode.")
            return pd.DataFrame(columns=['barcode', 'cluster'])

        if not os.path.exists(cluster_file):
            logger.info(f"cluster_file does not exist: {cluster_file}. Cell number will be estimated in bulk mode.")
            return pd.DataFrame(columns=['barcode', 'cluster'])

        if os.path.getsize(cluster_file) == 0:
            logger.info(f"cluster_file is empty: {cluster_file}. Cell number will be estimated in bulk mode.")
            return pd.DataFrame(columns=['barcode', 'cluster'])

        try:
            cluster_df = pd.read_csv(
                cluster_file,
                sep="\t",
                header=None,
                names=['barcode', 'cluster'],
                na_values=[]
            )

            if cluster_df.empty:
                logger.info(f"cluster_file is empty after parsing: {cluster_file}. Using bulk mode.")
                return pd.DataFrame(columns=['barcode', 'cluster'])

            cluster_df['cluster'] = cluster_df['cluster'].apply(
                lambda x: str(int(x)) if isinstance(x, float) and x.is_integer()
                else str(x) if pd.notnull(x) else "NA"
            )

            return cluster_df

        except Exception as e:
            logger.exception(
                f"Failed to read cluster_file={cluster_file}. Cell number will be estimated in bulk mode. Error: {e}"
            )
            return pd.DataFrame(columns=['barcode', 'cluster'])

    def _run(self, context):
        cell_num = context.get('cell_num', self.config.get('cell_num', 0))

        # Case 1: constant integer
        if isinstance(cell_num, int) and cell_num != 0:
            logger.info(f"Constant cell_num is provided: {cell_num}. Skip automatic cell number estimation.")
            return

        # Case 2: user provided file
        if isinstance(cell_num, str):
            if Path(cell_num).exists():
                logger.info(f"Using user-provided cell_num file: {cell_num}")
                return
            else:
                raise FileNotFoundError(f"Provided cell_num file does not exist: {cell_num}")

        # Case 3: auto estimation
        if not (isinstance(cell_num, int) and cell_num == 0):
            raise ValueError(f"Unsupported cell_num setting: {cell_num}")

        save_dir = Path(self.work_dir)

        # input
        bam_file = self.get_inputs(context)['in_filter_bam']
        if not bam_file or not os.path.exists(bam_file):
            raise FileNotFoundError(f"in_filter_bam not found: {bam_file}")

        # output
        outputs = self.get_outputs(context)
        count_file = save_dir / "raw_cell_num.txt"
        cell_num_file = outputs['cell_num']

        # parameters
        rerun = True
        seq_type = self.config.get('sequence_type')
        bins = self.config.get('bins')

        # step 1. get barcode-level UMI counts
        if rerun or not os.path.exists(count_file):
            logger.info("Generating raw cell number counts from BAM.")
            cb_ub_df = get_cb_ub(bam_file, count_file, seq_type, bins)
        else:
            logger.info(f"Loading existing raw cell number counts: {count_file}")
            cb_ub_df = pd.read_csv(count_file, sep='\t', header=0)

        if cb_ub_df.empty:
            logger.warning("Raw cell number dataframe is empty. Writing an empty output file.")
            pd.DataFrame(columns=['barcode', 'cluster', 'UMI_counts', 'refined_cell_num']).to_csv(
                cell_num_file, index=False, sep="\t"
            )
            return

        # step 2. load cluster information
        cluster_df = self._load_cluster_df(context)

        # step 3. merge cluster info or use bulk mode
        if not cluster_df.empty:
            logger.info("Estimating cell number with cluster-aware mode.")
            file_merged = pd.merge(cluster_df, cb_ub_df, on='barcode', how='inner')

            if file_merged.empty:
                logger.warning(
                    "Merged dataframe between cluster_file and BAM-derived barcode counts is empty. "
                    "Falling back to bulk mode."
                )
                file_merged = cb_ub_df.copy()
                file_merged['cluster'] = 'bulk'
        else:
            logger.info("Estimating cell number with bulk mode.")
            file_merged = cb_ub_df.copy()
            file_merged['cluster'] = 'bulk'

        # step 4. calculate refined cell number
        cluster_sums = file_merged.groupby('cluster')['nUMI'].sum()
        max_cluster = cluster_sums.idxmax()
        max_cluster_data = file_merged[file_merged['cluster'] == max_cluster]
        median_nUMI = max_cluster_data['nUMI'].median()

        if pd.isna(median_nUMI) or median_nUMI == 0:
            logger.warning("median_nUMI is invalid (NaN or 0). Setting refined_cell_num to 1 for all barcodes.")
            file_merged['refined_cell_num'] = 1
        else:
            file_merged['calculated_cell_num'] = np.ceil(20 * file_merged['nUMI'] / median_nUMI)
            file_merged['refined_cell_num'] = np.where(
                file_merged['calculated_cell_num'] > 25,
                25,
                file_merged['calculated_cell_num']
            )

        final_output_df = file_merged[['barcode', 'cluster', 'nUMI', 'refined_cell_num']].copy()
        final_output_df = final_output_df.rename(columns={'nUMI': 'UMI_counts'})

        final_output_df.to_csv(cell_num_file, index=False, sep="\t")
        logger.info(f"Cell number result written to: {cell_num_file}")

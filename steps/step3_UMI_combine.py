#!/usr/bin/env python3
"""
UMI combine
"""

from functools import partial
import gc
import multiprocessing
import os
from typing import Dict

import pandas as pd
from SpaceTracer.steps.base import BaseStep
from SpaceTracer.utils.logger import get_logger

from SpaceTracer.cores.UMI_combine import combine_UMI_spot_for_both_mosaic_and_error

model_name=__name__
logger = get_logger(model_name)

# ── Column definitions ────────────────────────────────────────────────────────

# Matches the barcode_list structure built inside the worker function:
# [chrom, pos, ".", ref, alt, barcode, consensus_read_count_str,
#  qual_A, qual_T, qual_C, qual_G]
COLUMNS_MAIN = [
    "#chrom", "pos", "ID", "ref", "alt",
    "barcode", "consensus_read_count", "qA", "qT", "qC", "qG",
]

# error_allele file: one row per allele observation at each position
COLUMNS_ERROR_ALLELE = ["#chrom", "pos", "ref", "alt", "strand"]

# ── Step class ────────────────────────────────────────────────────────────────
class UMICombineStep(BaseStep):
    def get_inputs(self, context: Dict) -> Dict[str, str]:
        """input"""
        inputs = {
            'in_filter_bam': context.get('in_filter_bam'),
            'reference': self.config.get('genome_fasta'),
            'filter_mpileup_file':context.get('filter_mpileup_file')
        }

        return inputs
    
    def get_outputs(self,context: Dict) -> Dict[str, str]:
        """output"""
        return {
            'spot_count_file': os.path.join(self.step_dir, 'spot.count.txt'),
            'spot_count_parquet': os.path.join(self.step_dir, 'spot.count.parquet'),
            'error_count_file': os.path.join(self.step_dir, 'error.count.txt')
        }
    
    def _run(self,context: Dict):
        inputs=self.get_inputs(context)
        bam_file=inputs["in_filter_bam"]
        mpileup_file=inputs["filter_mpileup_file"]
        
        seq_type=self.config.get('sequence_type')
        thread=self.threads

        outputs=self.get_outputs(context)
        spot_count_file=outputs['spot_count_file']
        spot_count_parquet=outputs['spot_count_parquet']
        error_count_file=outputs['error_count_file']

        df = pd.read_csv(mpileup_file, sep="\t")
        
        worker_args = []
        for _, row in df.iterrows():
            types = {t.strip() for t in str(row["type"]).split(",")}

            # Pass the actual type string if present, else empty string
            check_mosaic = True if "candidate_somatic" in types else False
            check_error  = True if "candidate_error"  in types else False

            # Everything the worker needs to locate the position in the BAM.
            # Adjust fields to match your function's actual `identifier` usage.
            identifier = "_".join([str(row["#chrom"]), str(row["pos"]), row["ref"], row["alt1"]])
            worker_args.append((check_mosaic, check_error, identifier))

        # partial_func=partial(combine_UMI_spot_for_both_mosaic_and_error,bam_file,seq_type)
        partial_func = partial(_worker_wrapper, bam_file, seq_type)
        with (
                multiprocessing.Pool(thread) as pool,
                open(spot_count_file,  "w") as f_main,
                open(error_count_file, "w") as f_err,
            ):
                # Write headers once before any results arrive
                pd.DataFrame(columns=COLUMNS_MAIN).to_csv(
                    f_main, header=True, index=False, sep="\t"
                )
                pd.DataFrame(columns=COLUMNS_ERROR_ALLELE).to_csv(
                    f_err, header=True, index=False, sep="\t"
                )

                for result in pool.imap(partial_func, worker_args, chunksize=10):

                    # Worker returned None → unrecoverable error, skip position
                    if result is None:
                        continue

                    spot_count_list, error_list = result

                    # new_list == [] → skip this position entirely
                    if not spot_count_list:
                        continue

                    # Write UMI rows (main file)
                    pd.DataFrame(spot_count_list, columns=COLUMNS_MAIN).to_csv(
                        f_main, header=False, index=False, sep="\t"
                    )

                    # Write error allele row.
                    # error_allele and strand are single characters (one per position).
                    # chrom/pos/ref/alt are taken from identifier to avoid
                    # re-parsing new_list.
                    if error_list:
                        error_row = [error_list]
                        pd.DataFrame(error_row, columns=COLUMNS_ERROR_ALLELE).to_csv(
                            f_err, header=False, index=False, sep="\t"
                        )

                    del result
                    gc.collect()

        df_main = pd.read_csv(spot_count_file, sep='\t')
        df_main.to_parquet(spot_count_parquet, index=False)

        logger.info(
            f"[UmiCombine] Done → {spot_count_file}, {error_count_file}"
        )

        return {
            "umi_combined_file": spot_count_file,
            "error_file": error_count_file,
        }
    
    # ── Module-level worker (must NOT be inside the class) ───────────────────────

def _worker_wrapper(bam_file,seq_type, args):
    """
    Called by each pool worker for one position.

    Parameters (via partial + imap)
    ────────────────────────────────
    bam_file    str   fixed for the whole run
    seq_type    str   fixed for the whole run
    args        tuple (check_mosaic, check_error, run_type, identifier)

    Returns
    ───────
    (new_list, error_allele, strand)   on success
    None                               on error (position is skipped upstream)
    """
    check_mosaic, check_error, identifier = args
    try:
        mosaic_spot_list, error_list = combine_UMI_spot_for_both_mosaic_and_error(
            bam_file,
            check_mosaic,   # "" when this position is not candidate_mosaic
            check_error,    # "" when this position is not candidate_error
            seq_type,
            identifier,
        )
        # new_list  : list of rows (multiple barcodes per position)
        # error_allele : single char, e.g. "A"  — or "" if not applicable
        # strand       : single char, e.g. "+"  — or "" if not applicable
        return mosaic_spot_list, error_list
    
    except Exception as exc:
        logger.warning(f"[UmiCombine] worker failed for {identifier}: {exc}")
        return None
#!/usr/bin/env python3
"""
PhyloSOLID: Robust phylogeny reconstruction from single-cell data despite inherent error and sparsity

This pipeline constructs phylogenetic trees from spatial transcriptomics data with:
- Germline variant filtering (optional)
- Scaffold tree construction with CV threshold optimization
- Full-resolved tree building with mutation integration
- Artifact removal and discordance metric calculation

Author: Qing
Date: 2025/09/16
Update: 2025/10/13
Latest: 2026/07/29

Usage:
    python [this script] -s SAMPLE_ID -i /path/to/data -o /path/to/output

Output Structure:
    outputpath/
    ├── 01_germline_filter/
    ├── 02_scaffold_builder/
    │   ├── CV0.3/
    │   ├── CV0.4/
    │   └── ...
    ├── 03_mutation_integrator/
    │   ├── CV0.3/
    │   ├── CV0.4/
    │   └── ...
    ├── 04_final_results/
    └── cv_threshold_search_results.csv

"""


import time
start_time = time.perf_counter()

import os
os.environ['PYTHONHASHSEED'] = '42'

import warnings
warnings.filterwarnings("ignore", category=FutureWarning)


################################################################################################
########################################## PhyloSOLID ##########################################
################################################################################################
import logging
import copy
import random
import pandas as pd
import numpy as np
from tqdm import tqdm
from copy import deepcopy
import re
import json
import sys
import shutil

logger = logging.getLogger(__name__)

from src.data_loader import load_all
from src.scrna_classifier import real_time_classifier_predict
from src.germline_filter import identify_germline_variants
from src.germline_filter import *
from src.scaffold_builder import build_scaffold_tree
from src.scaffold_builder import *
from src.mutation_integrator import *
from src.full_tree_builder import build_fully_resolved_tree


# ------------------------------
# Configure logging
# ------------------------------
logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")

# ------------------------------
# Project parameters and file paths
# ------------------------------
import multiprocessing as mp
import argparse
from argparse import ArgumentParser
parser = argparse.ArgumentParser()

parser.add_argument("-s", "--sampleid", default="", type=str, help="The sampleid you can set and check.")
parser.add_argument("-i", "--inputpath", default="yourpath/data", type=str, help="The inputpath contains the preprocessing results from raw posterior-reads data.")
parser.add_argument("-o", "--outputpath", default="yourpath/results", type=str, help="The outputpath you want to save results.")
parser.add_argument("-c", "--celltype_file", default=None, type=str, help="The celltype_file you should provide. If you can't generate this file, please set 'None'.")
parser.add_argument("--is_predict_germ", default="no", choices=["yes", "no"], type=str, help="Select 'yes' or 'no' to determine whether to predict germline mutations.")
parser.add_argument("--is_detect_passtree_by_dp", default="no", choices=["yes", "no"], type=str, help="Select 'yes' or 'no' to determine whether to run Dynamic programing step.")
parser.add_argument("--is_filter_quality", default="yes", choices=["yes", "no"], type=str, help="Select 'yes' or 'no' to determine whether to filter mutations in scaffold steps by coverage quality.")
parser.add_argument("--cv_rank_thresh", default="0.3", type=str, 
                    help="""CV rank threshold for coverage-based filtration.
                    Options:
                      - Single value: '0.3' (run once with this value)
                      - Comma-separated: '0.3,0.5,0.7' (search over these values)
                      - Range: '0.3-0.7:0.1' (from 0.3 to 0.7 with step 0.1)
                      - 'auto': use default presets [0.3, 0.4, 0.5, 0.6, 0.7]""")
parser.add_argument("--remove_artifact_mutations", default="yes", choices=["yes", "no"], type=str, help="Select 'yes' or 'no' to determine whether to permanently remove artifact mutations.")
parser.add_argument("--seed", default=42, type=int, help="Random seed for reproducibility")

args = parser.parse_args()

# Set all random seeds
from src.reproducibility import set_seed, deterministic_choice
set_seed(args.seed)


# ============================================================
# Helper function: Parse cv_rank_thresh parameter
# ============================================================
def parse_cv_thresholds(cv_rank_thresh_str):
    """
    Parse cv_rank_thresh parameter, returning a list of thresholds to search.
    
    Examples:
        '0.3'           -> [0.3]
        '0.3,0.5,0.7'   -> [0.3, 0.5, 0.7]
        '0.3-0.7:0.1'   -> [0.3, 0.4, 0.5, 0.6, 0.7]
        '0.3-0.7'       -> [0.3, 0.4, 0.5, 0.6, 0.7]  (default step 0.1)
        'auto'          -> [0.3, 0.4, 0.5, 0.6, 0.7]
    """
    cv_rank_thresh_str = cv_rank_thresh_str.strip()
    
    # Case 1: 'auto' -> use default presets
    if cv_rank_thresh_str.lower() == 'auto':
        return [0.3, 0.4, 0.5, 0.6, 0.7]
    
    # Case 2: Range format: '0.3-0.7:0.1'
    range_pattern = r'^([\d.]+)\s*-\s*([\d.]+)(?:\s*:\s*([\d.]+))?$'
    match = re.match(range_pattern, cv_rank_thresh_str)
    if match:
        start = float(match.group(1))
        end = float(match.group(2))
        step = float(match.group(3)) if match.group(3) else 0.1
        
        # Generate sequence and round to 2 decimal places
        values = []
        current = start
        while current <= end + 1e-9:
            values.append(round(current, 2))
            current += step
        
        # Ensure end value is included
        if round(end, 2) not in values:
            values.append(round(end, 2))
            values = sorted(set(values))
        
        return values
    
    # Case 3: Comma-separated: '0.3,0.5,0.7'
    if ',' in cv_rank_thresh_str:
        values = [float(v.strip()) for v in cv_rank_thresh_str.split(',') if v.strip()]
        return values
    
    # Case 4: Single value
    try:
        return [float(cv_rank_thresh_str)]
    except ValueError:
        raise ValueError(f"Invalid cv_rank_thresh format: '{cv_rank_thresh_str}'. "
                         f"Supported formats: '0.3', '0.3,0.5,0.7', '0.3-0.7:0.1', 'auto'")


# get parameters
sampleid = args.sampleid
inputpath = args.inputpath
outputpath = args.outputpath
celltype_file = args.celltype_file
is_predict_germ = args.is_predict_germ
is_detect_passtree_by_dp = args.is_detect_passtree_by_dp
is_filter_quality = args.is_filter_quality
cv_rank_thresh_str = args.cv_rank_thresh
remove_artifact_mutations = args.remove_artifact_mutations

# Parse CV thresholds
cv_thresholds = parse_cv_thresholds(cv_rank_thresh_str)
is_search_mode = len(cv_thresholds) > 1

# Create main output directories
outputpath_01 = os.path.join(outputpath, "01_germline_filter")
outputpath_02 = os.path.join(outputpath, "02_scaffold_builder")
outputpath_03 = os.path.join(outputpath, "03_mutation_integrator")
outputpath_04 = os.path.join(outputpath, "04_final_results")

os.makedirs(outputpath_01, exist_ok=True)
os.makedirs(outputpath_02, exist_ok=True)
os.makedirs(outputpath_03, exist_ok=True)
os.makedirs(outputpath_04, exist_ok=True)


# Display parameters for verification
logger.info(f"sampleid: {sampleid}")
logger.info(f"inputpath: {inputpath}")
logger.info(f"outputpath: {outputpath}")
logger.info(f"celltype_file: {celltype_file}")
logger.info(f"is_predict_germ: {is_predict_germ}")
logger.info(f"is_detect_passtree_by_dp: {is_detect_passtree_by_dp}")
logger.info(f"is_filter_quality: {is_filter_quality}")
logger.info(f"cv_rank_thresh: {cv_rank_thresh_str}")
logger.info(f"  -> parsed as: {cv_thresholds}")
logger.info(f"  -> mode: {'SEARCH' if is_search_mode else 'SINGLE'}")
logger.info(f"remove_artifact_mutations: {remove_artifact_mutations}")
logger.info("")
logger.info("Directory structure:")
logger.info(f"  01_germline_filter: {outputpath_01}")
logger.info(f"  02_scaffold_builder: {outputpath_02}")
logger.info(f"  03_mutation_integrator: {outputpath_03}")
logger.info(f"  04_final_results: {outputpath_04}")


# ------------------------------
# Parameter settings
# ------------------------------

SETTING_PARAMS = {
    "models_path": "phylosolid/models/scdna",
    
    # 1 data loader
    "p_thresh": 0.5,
    
    # 2 germline filter
    "mcf_cutoff": 0.05,
    "mcn_cutoff": 5,
    
    # Pairwise correlation criteria (Section 2.1)
    "pair_N11_min": 0,
    "jaccard_thresh": 0.2,
    "jaccard_low": 0.1,
    "fraction_parent_child_thresh": 0.9,
    
    # 3.1 Initial filtration
    "posterior_threshold": 0.5,
    "maf_max_threshold": 0.3,
    "maf_mean_threshold": 0.1,
    
    # 3.2 Coverage-based filtration
    "na_prop_thresh_global": 0.95,
    "cv_thresh": 6.0,
    # "cv_rank_thresh" will be set dynamically
    
    # 3.3 Consensus correlation graph
    "consensus_runs": 100,
    "consensus_clone_freq_thresh": 0.1,
    "resolution_of_graph": 1,
        
    # 3.4 Penalty-based placement
    "general_weight_NA": 0.001,
    "fnfp_ratio": 0.1,
    "phi": 1.0,
    
    # 4.1 Dynamic programming
    "pass_tree_cutoff": 0.9,
    "unpass_tree_cutoff": 0.1,
    
    # 4.2 fp_ratio and fn_ratio
    "fp_ratio_cutoff_across_tree": 0.2,
    "fn_ratio_cutoff_across_tree": 0.9,
    "fp_ratio_cutoff_within_subclone": 0.1,
    "fp_ratio_persite_cutoff": 0.1,
    "fp_count_persite_cutoff": 0,
    
    "fp_ratio_per_mutation_cross_all_cells_cutoff": 0.2,
    "fp_count_per_mutation_cross_all_cells_cutoff": 5,
    "fp_ratio_per_cell_cross_all_muts_cutoff": 0.5,
    
    "intersection_vs_fn_flipping_ratio_cutoff": 0.2,
    "intersection_cell_count_on_mutation_cutoff": 5,
    "intersection_cell_ratio_on_mutation_cutoff": 0.2,
    "intersection_count_per_cells_cutoff": 1,
    "flipping_count_fn_per_cells_cutoff": 1,
    "flipping_to_1_count_per_cells_cutoff": 2
}

params = SETTING_PARAMS


# ------------------------------
# Step 1: Load data
# ------------------------------
logger.info("===== Step1: Loading data ...")
# Load raw sequencing data and build initial binary genotype matrix.
data = load_all(inputpath)
P_raw, V_raw, C_raw, A_raw = data["P"], data["V"], data["C"], data["A"]
df_features = data['features']
df_reads_raw = data['df_reads']
I_raw = build_binary_I(P_raw, V_raw, C_raw, params["p_thresh"])

logger.info(f"Loaded data: {len(P_raw)} cells, {len(I_raw.columns)} mutations")


##### Remove cells with no mutations (all zeros)
I_filtered = I_raw[I_raw.eq(1).any(axis=1)]
I = reorder_columns_by_mutant_stats(I_filtered, df_features)[0]
all_mutations = list(I.columns)

cells_in_I = I.index
cols_in_I = I.columns
P = P_raw.loc[cells_in_I, cols_in_I]
V = V_raw.loc[cells_in_I, cols_in_I]
C = C_raw.loc[cells_in_I, cols_in_I]
A = A_raw.loc[cells_in_I, cols_in_I]
bulk_row = df_reads_raw.loc[['bulk']]
df_reads_cells = df_reads_raw.loc[cells_in_I, cols_in_I]
df_reads = pd.concat([bulk_row, df_reads_cells])

df_features_new, empty_mutations = update_features_matrix(I, df_reads, df_features, params["mcf_cutoff"])


# ------------------------------
# Step 2: Classifier for identifying mosaic mutations (SKIPPED)
# ------------------------------
logger.info("===== Step2: Classifier ... SKIPPED (using all mutations as candidates)")

candidate_mutations = list(I_raw.columns)

P_candidate = P[candidate_mutations].copy()
V_candidate = V[candidate_mutations].copy()
A_candidate = A[candidate_mutations].copy()
C_candidate = C[candidate_mutations].copy()
I_candidate = I[candidate_mutations].copy()
df_reads_candidate = df_reads[candidate_mutations].copy()


# ------------------------------
# Step 3. Germline filtering
# ------------------------------
logger.info("===== Step3: Predict germline mutations ...")

if is_predict_germ == "yes":
    logging.info("Running germline filtering ...")
    stats_df, germline_mutations = identify_germline_variants(
        P=P, V=V, C=C, df_reads=df_reads, df_features_new=df_features_new, 
        p_thresh=params["p_thresh"],
        mcf_cutoff=params["mcf_cutoff"],
        mcn_cutoff=params["mcn_cutoff"],
        outputpath=outputpath_01,
        sampleid=sampleid
    )
else:
    germline_mutations = set()

predicted_germline_mutations = list(germline_mutations)
rescued_germline_mutations = [mut for mut in predicted_germline_mutations 
                             if mut in df_features_new.columns 
                             and df_features_new.loc['mutant_cell_fraction_detected', mut] < 0.5]

removed_germline_mutations = [i for i in predicted_germline_mutations if i not in rescued_germline_mutations]


logging.info(f"Identified {len(predicted_germline_mutations)} germline variants")


plot_heatmap_with_germline_mutations(I, predicted_germline_mutations, 
                                     os.path.join(outputpath_01, sampleid + ".heatmap_with_predicted_germline_mutations_and_histograms.pdf"))
removed_artifact_mutations = []
somatic_mutations_init = [i for i in all_mutations if i not in removed_germline_mutations and i not in removed_artifact_mutations]
somatic_mutations = list((reorder_columns_by_mutant_stats(I[somatic_mutations_init], df_features_new)[0]).columns)
P_somatic = P[somatic_mutations].copy()
V_somatic = V[somatic_mutations].copy()
A_somatic = A[somatic_mutations].copy()
C_somatic = C[somatic_mutations].copy()
I_somatic = I[somatic_mutations].copy()
df_reads_somatic = df_reads[somatic_mutations].copy()

I_somatic_withNA3 = I_somatic.replace({np.nan: 3}).astype(int)
I_somatic_withNA3.to_csv(os.path.join(outputpath_01, "I_somatic_withNA3.txt"), sep="\t")
df_features_new = add_mutation_proportions_to_features(df_features_new, I_somatic)


# ------------------------------
# Step 4: Scaffold builder (with dynamic cv_rank_thresh)
# ------------------------------
logger.info("===== Step4: Construct scaffold tree ...")

if celltype_file is None or celltype_file == "None":
    barcodes = df_reads_somatic.index.tolist()
    df_celltype = pd.DataFrame({
        "barcode": barcodes,
        "cell_type": ["default_type"] * len(barcodes)
    })
else:
    df_celltype = pd.read_csv(celltype_file, sep="\t")

df_celltype.to_csv(os.path.join(outputpath_02, "df_celltype.txt"), sep="\t")
logger.info(f"Celltype data loaded: {df_celltype.shape[0]} cells")

logging.info("Running scaffold building ...")
immune_mutations = []


# ============================================================
# Core function: Run full pipeline for a given cv_rank_thresh
# ============================================================
def run_pipeline_for_cv_threshold(cv_value, output_dir):
    """
    Run the full PhyloSOLID pipeline for a single CV threshold value.
    
    Parameters
    ----------
    cv_value : float
        The CV rank threshold to use
    output_dir : str
        Base output directory for this run
    
    Returns
    -------
    dict : {
        'cv_value': float,
        'omega_pre_qc': float,
        'omega_final': float,
        'omega_sum': float,
        'scaffold_count': int,
        'success': bool,
        'error': str or None,
        'scaffold_dir': str,
        'mutint_dir': str,
        'T_current': TreeNode,
        'M_current': pd.DataFrame,
        'root_mutations': List[str],
        'all_conflict_mutations': List[str],
        'I_attached': pd.DataFrame,
        'attached_mutations': List[str],
        'scaffold_mutations': List[str],
        'somatic_mutations': List[str],
        'no_group_mutations': List[str],
        'to_be_removed_cells': List[str],
        'identified_doublet_cells': List[str],
        'to_be_removed_mutations_by_fp_mutations_cross_all_cells': List[str],
        'final_remained_mutations': List[str],
        'final_conflict_mutations': List[str],
    }
    """
    # Create a copy of params with the current cv_rank_thresh
    params_local = copy.deepcopy(params)
    params_local['cv_rank_thresh'] = cv_value
    
    # Create run-specific output directories under 02_scaffold_builder and 03_mutation_integrator
    cv_label = str(cv_value).replace('.', '_')
    run_outputpath_scaffold = os.path.join(outputpath_02, f"CV{cv_label}")
    run_outputpath_full = os.path.join(outputpath_03, f"CV{cv_label}")
    os.makedirs(run_outputpath_scaffold, exist_ok=True)
    os.makedirs(run_outputpath_full, exist_ok=True)
    
    logger.info(f"  Running with cv_rank_thresh = {cv_value}")
    logger.info(f"  Scaffold output: {run_outputpath_scaffold}")
    logger.info(f"  Mutation integrator output: {run_outputpath_full}")
    
    try:
        # ---- Scaffold builder ----
        results_of_scaffold = build_scaffold_tree(
            P_somatic=P_somatic, 
            V_somatic=V_somatic, 
            A_somatic=A_somatic, 
            C_somatic=C_somatic, 
            I_somatic=I_somatic,
            df_reads_somatic=df_reads_somatic,
            df_features_new=df_features_new,
            params=params_local,
            is_filter_quality=is_filter_quality,
            outputpath_scaffold=run_outputpath_scaffold,
            sampleid=sampleid,
            df_celltype=df_celltype,
            immune_mutations=immune_mutations
        )
        
        # build_scaffold_tree now returns 14 values
        (T_scaffold, M_scaffold, df_flipping_spots, df_total_flipping_count, 
         final_cleaned_I_selected_withNA3, final_cleaned_M_scaffold, 
         backbone_mutations, mutation_group, spots_to_split, group_mutations, 
         no_group_mutations, remained_mutations, conflict_mutations, root_mutations) = results_of_scaffold
        
        scaffold_mutations = list(M_scaffold.columns)
        
        # ---- Step 5: DP pass tree & prepare data ----
        if is_detect_passtree_by_dp == "yes":
            pass_tree_cutoff = params_local['pass_tree_cutoff']
            unpass_tree_cutoff = params_local['unpass_tree_cutoff']
            p_thresh = params_local["p_thresh"]
            
            df_DP_results, passtree_mutations, onecell_mutations = run_dp_pass_tree(
                data=data, 
                df_features_new=df_features_new, 
                M_scaffold=M_scaffold, 
                outputpath_full=run_outputpath_full, 
                scaffold_mutations=scaffold_mutations,
                p_thresh=p_thresh,
                pass_tree_cutoff=pass_tree_cutoff,
                unpass_tree_cutoff=unpass_tree_cutoff,
                is_log_value_for_likelihoods=True
            )
        else:
            passtree_mutations = all_mutations
        
        attached_mutations = [i for i in passtree_mutations 
                             if i not in scaffold_mutations 
                             and i not in removed_germline_mutations 
                             and i not in removed_artifact_mutations]
        
        # Prepare data for full-resolved tree building
        I_attached_selected = I[scaffold_mutations + attached_mutations]
        I_attached_selected_sorted = I_attached_selected[
            I_attached_selected.apply(lambda col: (col == 1).sum(), axis=0).sort_values(ascending=False).index
        ]
        I_attached_sorted_non_empty = I_attached_selected_sorted[
            I_attached_selected_sorted.eq(1).any(axis=1)
        ]
        
        P_attached_sorted_non_empty = P.loc[
            I_attached_sorted_non_empty.index,
            I_attached_sorted_non_empty.columns
        ]
        
        I_attached_split, P_attached_split = split_spots_by_immune_mutations(
            spots_to_split, 
            [i for i in immune_mutations if i in I_attached_sorted_non_empty.columns], 
            I_attached_sorted_non_empty, 
            P_attached_sorted_non_empty
        )
        I_attached, sorting_stats_of_I_attached = reorder_columns_by_mutant_stats(
            I_attached_split, 
            df_features_new,
            min_cell_threshold=30,
            bin_size=5,
            descending=True
        )
        P_attached = P_attached_split[I_attached.columns]
        
        all_conflict_mutations = conflict_mutations.copy()
        
        # ---- Step 6-8: Build fully resolved tree ----
        (T_current, M_current, root_mutations, all_conflict_mutations,
         omega_before_qc, to_be_removed_cells, identified_doublet_cells,
         to_be_removed_mutations_by_fp_mutations_cross_all_cells,
         final_remained_mutations, final_conflict_mutations) = build_fully_resolved_tree(
            T_scaffold=T_scaffold,
            M_scaffold=M_scaffold,
            scaffold_mutations=scaffold_mutations,
            I_attached=I_attached,
            P_attached=P_attached,
            df_features_new=df_features_new,
            params=params_local,
            outputpath_full=run_outputpath_full,
            sampleid=sampleid,
            attached_mutations=attached_mutations,
            immune_mutations=immune_mutations,
            spots_to_split=spots_to_split,
            conflict_mutations=conflict_mutations,
            remove_artifact_mutations=remove_artifact_mutations,
            logger_obj=logger
        )
        
        # ---- Compute omega_final (post-QC) ----
        M_for_omega = M_current.drop(columns=['ROOT'], errors='ignore')
        mutations_on_tree_for_omega = M_for_omega.columns.to_series().apply(lambda x: x.split("|")).explode().unique().tolist()
        M_for_omega = split_merged_columns(M_for_omega, mutations_on_tree_for_omega)
        
        M_for_omega_clean = M_for_omega.loc[:, (M_for_omega != 0).any(axis=0)]
        M_for_omega_clean = M_for_omega_clean.loc[(M_for_omega_clean != 0).any(axis=1)]
        
        I_for_omega = I_attached.loc[M_for_omega_clean.index, M_for_omega_clean.columns].replace({np.nan: 3}).astype(int)
        
        N_deltaFP = ((I_for_omega == 1) & (M_for_omega_clean == 0)).sum().sum()
        N_deltaFN = ((I_for_omega == 0) & (M_for_omega_clean == 1)).sum().sum()
        
        omega_final = N_deltaFP + params_local['fnfp_ratio'] * N_deltaFN
        
        logger.info(f"  cv_rank_thresh={cv_value}:")
        logger.info(f"    Omega (pre-QC): {omega_before_qc:.4f}")
        logger.info(f"    Omega (final): {omega_final:.4f}")
        logger.info(f"    Omega (sum): {omega_before_qc + omega_final:.4f}")
        
        return {
            'cv_value': cv_value,
            'omega_pre_qc': omega_before_qc,
            'omega_final': omega_final,
            'omega_sum': omega_before_qc + omega_final,
            'scaffold_count': len(scaffold_mutations),
            'success': True,
            'error': None,
            'scaffold_dir': run_outputpath_scaffold,
            'mutint_dir': run_outputpath_full,
            'T_current': T_current,
            'M_current': M_current,
            'root_mutations': root_mutations,
            'all_conflict_mutations': all_conflict_mutations,
            'I_attached': I_attached,
            'attached_mutations': attached_mutations,
            'scaffold_mutations': scaffold_mutations,
            'somatic_mutations': somatic_mutations,
            'no_group_mutations': no_group_mutations,
            'to_be_removed_cells': to_be_removed_cells,
            'identified_doublet_cells': identified_doublet_cells,
            'to_be_removed_mutations_by_fp_mutations_cross_all_cells': to_be_removed_mutations_by_fp_mutations_cross_all_cells,
            'final_remained_mutations': final_remained_mutations,
            'final_conflict_mutations': final_conflict_mutations,
        }
        
    except Exception as e:
        logger.error(f"Error for cv_rank_thresh={cv_value}: {e}")
        import traceback
        traceback.print_exc()
        return {
            'cv_value': cv_value,
            'omega_pre_qc': float('inf'),
            'omega_final': float('inf'),
            'omega_sum': float('inf'),
            'scaffold_count': 0,
            'success': False,
            'error': str(e),
            'scaffold_dir': None,
            'mutint_dir': None,
            'T_current': None,
            'M_current': None,
            'root_mutations': [],
            'all_conflict_mutations': [],
            'I_attached': None,
            'attached_mutations': [],
            'scaffold_mutations': [],
            'somatic_mutations': [],
            'no_group_mutations': [],
            'to_be_removed_cells': [],
            'identified_doublet_cells': [],
            'to_be_removed_mutations_by_fp_mutations_cross_all_cells': [],
            'final_remained_mutations': [],
            'final_conflict_mutations': [],
        }


# ============================================================
# Main execution: Single or search mode
# ============================================================

if is_search_mode:
    # ============================================================
    # SEARCH MODE: Run pipeline for all CV thresholds
    # ============================================================
    logger.info("=" * 80)
    logger.info("=== CV THRESHOLD SEARCH MODE ===")
    logger.info(f"Searching over {len(cv_thresholds)} values: {cv_thresholds}")
    logger.info("=" * 80)
    
    search_results = []
    
    for cv_val in tqdm(cv_thresholds, desc="Searching CV thresholds"):
        logger.info(f"\n--- Testing cv_rank_thresh = {cv_val} ---")
        result = run_pipeline_for_cv_threshold(cv_val, outputpath)
        search_results.append(result)
    
    # ---- Summarize results ----
    df_results = pd.DataFrame(search_results)
    
    # Find optimal threshold (minimize omega_sum)
    df_valid = df_results[df_results['omega_sum'] != float('inf')]
    
    logger.info("=" * 80)
    logger.info("=== CV THRESHOLD SEARCH RESULTS ===")
    logger.info("=" * 80)
    logger.info("\n" + df_results.to_string())
    
    if not df_valid.empty:
        best_idx = df_valid['omega_sum'].idxmin()
        best_cv = df_valid.loc[best_idx, 'cv_value']
        best_omega_sum = df_valid.loc[best_idx, 'omega_sum']
        
        logger.info("=" * 80)
        logger.info(f"BEST CV RANK THRESHOLD: {best_cv}")
        logger.info(f"  Omega (pre-QC): {df_valid.loc[best_idx, 'omega_pre_qc']:.4f}")
        logger.info(f"  Omega (final): {df_valid.loc[best_idx, 'omega_final']:.4f}")
        logger.info(f"  Omega sum: {best_omega_sum:.4f}")
        logger.info(f"  Scaffold count: {df_valid.loc[best_idx, 'scaffold_count']}")
        logger.info("=" * 80)
        
        # Save search results
        df_results.to_csv(os.path.join(outputpath, "cv_threshold_search_results.csv"), index=False)
        
        # ---- Extract all variables from the best result ----
        best_result = df_valid.loc[best_idx]
        
        T_current = best_result.get('T_current')
        M_current = best_result.get('M_current')
        root_mutations = best_result.get('root_mutations', [])
        all_conflict_mutations = best_result.get('all_conflict_mutations', [])
        omega_before_qc = best_result.get('omega_pre_qc', 0.0)
        
        # Step 9 required variables
        I_attached = best_result.get('I_attached')
        attached_mutations = best_result.get('attached_mutations', [])
        scaffold_mutations = best_result.get('scaffold_mutations', [])
        somatic_mutations = best_result.get('somatic_mutations', [])
        no_group_mutations = best_result.get('no_group_mutations', [])
        to_be_removed_cells = best_result.get('to_be_removed_cells', [])
        identified_doublet_cells = best_result.get('identified_doublet_cells', [])
        to_be_removed_mutations_by_fp_mutations_cross_all_cells = best_result.get('to_be_removed_mutations_by_fp_mutations_cross_all_cells', [])
        final_remained_mutations = best_result.get('final_remained_mutations', [])
        final_conflict_mutations = best_result.get('final_conflict_mutations', [])
        best_mutint_dir = best_result.get('mutint_dir')
        
        # Store best_cv for summary
        optimal_cv = best_cv
        
        # ---- Copy best results to 04_final_results ----
        logger.info("=" * 80)
        logger.info(f"Copying best results (CV={best_cv}) to 04_final_results/")
        logger.info("=" * 80)
        
        if best_mutint_dir and os.path.exists(best_mutint_dir):
            for f in os.listdir(best_mutint_dir):
                src = os.path.join(best_mutint_dir, f)
                dst = os.path.join(outputpath_04, f)
                if os.path.isfile(src):
                    shutil.copy2(src, dst)
                    logger.info(f"  Copied: {f}")
            
            # Also copy the scaffold directory summary
            best_scaffold_dir = best_result.get('scaffold_dir')
            if best_scaffold_dir and os.path.exists(best_scaffold_dir):
                scaffold_dst = os.path.join(outputpath_04, "scaffold_summary")
                os.makedirs(scaffold_dst, exist_ok=True)
                for f in os.listdir(best_scaffold_dir):
                    src = os.path.join(best_scaffold_dir, f)
                    dst = os.path.join(scaffold_dst, f)
                    if os.path.isfile(src):
                        shutil.copy2(src, dst)
        
        # Create a README file in final_results
        with open(os.path.join(outputpath_04, "README.txt"), 'w') as f:
            f.write("=" * 80 + "\n")
            f.write("PhyloSOLID FINAL RESULTS\n")
            f.write("=" * 80 + "\n\n")
            f.write(f"Sample ID: {sampleid}\n")
            f.write(f"Optimal CV threshold: {best_cv}\n")
            f.write(f"Omega (pre-QC): {df_valid.loc[best_idx, 'omega_pre_qc']:.4f}\n")
            f.write(f"Omega (final): {df_valid.loc[best_idx, 'omega_final']:.4f}\n")
            f.write(f"Omega (combined): {best_omega_sum:.4f}\n")
            f.write(f"Scaffold count: {df_valid.loc[best_idx, 'scaffold_count']}\n\n")
            f.write("This directory contains the best results only.\n")
            f.write("For detailed outputs for all CV thresholds, see:\n")
            f.write(f"  - {outputpath_02}/\n")
            f.write(f"  - {outputpath_03}/\n")
            f.write("=" * 80 + "\n")
        
    else:
        logger.error("No valid results from threshold search!")
        sys.exit(1)
        
else:
    # ============================================================
    # SINGLE MODE: Run with the single CV threshold
    # ============================================================
    cv_value = cv_thresholds[0]
    logger.info("=" * 80)
    logger.info(f"SINGLE CV THRESHOLD MODE: cv_rank_thresh = {cv_value}")
    logger.info("=" * 80)
    
    # Run the full pipeline with the single CV threshold
    result = run_pipeline_for_cv_threshold(cv_value, outputpath)
    
    if not result['success']:
        logger.error(f"Pipeline failed for cv_rank_thresh={cv_value}")
        sys.exit(1)
    
    # ---- Extract all variables from result ----
    T_current = result.get('T_current')
    M_current = result.get('M_current')
    root_mutations = result.get('root_mutations', [])
    all_conflict_mutations = result.get('all_conflict_mutations', [])
    omega_before_qc = result.get('omega_pre_qc', 0.0)
    
    # Step 9 required variables
    I_attached = result.get('I_attached')
    attached_mutations = result.get('attached_mutations', [])
    scaffold_mutations = result.get('scaffold_mutations', [])
    somatic_mutations = result.get('somatic_mutations', [])
    no_group_mutations = result.get('no_group_mutations', [])
    to_be_removed_cells = result.get('to_be_removed_cells', [])
    identified_doublet_cells = result.get('identified_doublet_cells', [])
    to_be_removed_mutations_by_fp_mutations_cross_all_cells = result.get('to_be_removed_mutations_by_fp_mutations_cross_all_cells', [])
    final_remained_mutations = result.get('final_remained_mutations', [])
    final_conflict_mutations = result.get('final_conflict_mutations', [])
    
    # Store cv_value for summary
    optimal_cv = cv_value
    
    # ---- Copy results to 04_final_results ----
    best_mutint_dir = result.get('mutint_dir')
    if best_mutint_dir and os.path.exists(best_mutint_dir):
        for f in os.listdir(best_mutint_dir):
            src = os.path.join(best_mutint_dir, f)
            dst = os.path.join(outputpath_04, f)
            if os.path.isfile(src):
                shutil.copy2(src, dst)
    
    # Save cv threshold search results (single value)
    df_single = pd.DataFrame([{
        'cv_value': cv_value,
        'omega_pre_qc': result['omega_pre_qc'],
        'omega_final': result['omega_final'],
        'omega_sum': result['omega_sum'],
        'scaffold_count': result['scaffold_count'],
        'scaffold_dir': result.get('scaffold_dir'),
        'mutint_dir': result.get('mutint_dir'),
    }])
    df_single.to_csv(os.path.join(outputpath, "cv_threshold_search_results.csv"), index=False)
    
    # Create README
    with open(os.path.join(outputpath_04, "README.txt"), 'w') as f:
        f.write("=" * 80 + "\n")
        f.write("PhyloSOLID FINAL RESULTS\n")
        f.write("=" * 80 + "\n\n")
        f.write(f"Sample ID: {sampleid}\n")
        f.write(f"CV threshold: {cv_value}\n")
        f.write(f"Omega (pre-QC): {result['omega_pre_qc']:.4f}\n")
        f.write(f"Omega (final): {result['omega_final']:.4f}\n")
        f.write(f"Omega (combined): {result['omega_sum']:.4f}\n")
        f.write(f"Scaffold count: {result['scaffold_count']}\n\n")
        f.write("This directory contains the final results.\n")
        f.write("For detailed outputs, see:\n")
        f.write(f"  - {outputpath_02}/CV{str(cv_value).replace('.', '_')}/\n")
        f.write(f"  - {outputpath_03}/CV{str(cv_value).replace('.', '_')}/\n")
        f.write("=" * 80 + "\n")


# ============================================================
# Step 9: Post-processing & output (using the final tree)
# ============================================================

logger.info("=" * 80)
logger.info("Step 9: Post-processing & output")
logger.info("=" * 80)
logger.info("")

# ----------------------------------------------------------------------------
# Step 9.1: Prepare final data
# ----------------------------------------------------------------------------
logger.info("STEP 9.1: Prepare final data")
logger.info("-" * 80)

# ---- Drop ROOT column and add root mutations back ----
M_current_filtered = M_current.drop(columns=['ROOT'], errors='ignore')

for mut_on_root in root_mutations:
    M_current_filtered.insert(0, mut_on_root, 1)

# ---- Expand merged columns ----
mutations_on_T_current = M_current_filtered.columns.to_series().apply(lambda x: x.split("|")).explode().unique().tolist()

T_full = copy.deepcopy(T_current)
M_full = split_merged_columns(M_current_filtered, mutations_on_T_current)

logger.info("Final full-resolved tree:")
print_tree_logger(T_full)

logger.info(f"  Final tree cells: {M_full.shape[0]}")
logger.info(f"  Final tree mutations: {M_full.shape[1]}")
logger.info("")


# ----------------------------------------------------------------------------
# Step 9.2: Output results to final_results
# ----------------------------------------------------------------------------
logger.info("STEP 9.2: Output result files")
logger.info("-" * 80)

# ---- Create output directory inside final_results ----
phylo_dir = os.path.join(outputpath_04, "phylo")
os.makedirs(phylo_dir, exist_ok=True)

# ---- Export I matrix with NA=3 ----
I_full_withNA3 = I_attached.replace({np.nan: 3}).astype(int)
I_full_withNA3.to_csv(os.path.join(phylo_dir, "I_full_withNA3.txt"), sep="\t")

# ---- Export M matrix ----
WriteTfile(os.path.join(phylo_dir, "M_full_basedPivots.filtered_sites_inferred"), M_full, M_full.index.tolist(), M_full.columns.tolist(), judge="yes")

# ---- Clean and export final matrices ----
final_cleaned_M_full = M_full.loc[:, (M_full != 0).any(axis=0)]
final_cleaned_M_full = final_cleaned_M_full.loc[(final_cleaned_M_full != 0).any(axis=1)]

kept_rows = final_cleaned_M_full.index
kept_cols = final_cleaned_M_full.columns

final_cleaned_I_full_withNA3 = I_full_withNA3.loc[kept_rows, kept_cols]

WriteTfile(os.path.join(phylo_dir, "final_cleaned_M_full_basedPivots.filtered_sites_inferred"), 
           final_cleaned_M_full, final_cleaned_M_full.index.tolist(), final_cleaned_M_full.columns.tolist(), judge="yes")
final_cleaned_I_full_withNA3.to_csv(os.path.join(phylo_dir, "final_cleaned_I_full_withNA3_for_circosPlot.txt"), sep="\t")

logger.info(f"  Output directory: {phylo_dir}")
logger.info("")


# ----------------------------------------------------------------------------
# Step 9.3: Identify flipping spots
# ----------------------------------------------------------------------------
logger.info("STEP 9.3: Identify flipping spots")
logger.info("-" * 80)

df_bin_withNA3_for_flipping = final_cleaned_I_full_withNA3.copy()
df_phylogeny = final_cleaned_M_full.copy()

# ---- Compute flipping spots ----
false_negative_flipping_spots = df_bin_withNA3_for_flipping.apply(
    lambda col: find_flipping_spots(col, df_phylogeny[col.name], condition_in_bin=0, condition_phylogeny=1)
)
false_positive_flipping_spots = df_bin_withNA3_for_flipping.apply(
    lambda col: find_flipping_spots(col, df_phylogeny[col.name], condition_in_bin=1, condition_phylogeny=0)
)
NAto1_flipping_spots = df_bin_withNA3_for_flipping.apply(
    lambda col: find_flipping_spots(col, df_phylogeny[col.name], condition_in_bin=3, condition_phylogeny=1)
)
NAto0_flipping_spots = df_bin_withNA3_for_flipping.apply(
    lambda col: find_flipping_spots(col, df_phylogeny[col.name], condition_in_bin=3, condition_phylogeny=0)
)

# ---- Handle empty results ----
if false_negative_flipping_spots.empty:
    false_negative_flipping_spots = {col: [] for col in df_bin_withNA3_for_flipping.columns}

if false_positive_flipping_spots.empty:
    false_positive_flipping_spots = {col: [] for col in df_bin_withNA3_for_flipping.columns}

if NAto1_flipping_spots.empty:
    NAto1_flipping_spots = {col: [] for col in df_bin_withNA3_for_flipping.columns}

if NAto0_flipping_spots.empty:
    NAto0_flipping_spots = {col: [] for col in df_bin_withNA3_for_flipping.columns}

# ---- Build flipping spots dataframe ----
df_flipping_spots = pd.DataFrame({
    'Mutation': df_bin_withNA3_for_flipping.columns,
    'delta_FN_spots': [', '.join(false_negative_flipping_spots.get(col, [])) for col in df_bin_withNA3_for_flipping.columns],
    'delta_FP_spots': [', '.join(false_positive_flipping_spots.get(col, [])) for col in df_bin_withNA3_for_flipping.columns],
    'NA_to_1_spots': [', '.join(NAto1_flipping_spots.get(col, [])) for col in df_bin_withNA3_for_flipping.columns],
    'NA_to_0_spots': [', '.join(NAto0_flipping_spots.get(col, [])) for col in df_bin_withNA3_for_flipping.columns]
})
df_flipping_spots.to_csv(os.path.join(phylo_dir, "df_flipping_spots.txt"), sep="\t", index=False)

logger.info("")


# ----------------------------------------------------------------------------
# Step 9.4: Calculate total flipping counts
# ----------------------------------------------------------------------------
logger.info("STEP 9.4: Calculate total flipping counts")
logger.info("-" * 80)

# ---- Compute total discordance counts ----
total_FN_flipping = ((df_bin_withNA3_for_flipping == 0) & (df_phylogeny == 1)).sum().sum()
total_FP_flipping = ((df_bin_withNA3_for_flipping == 1) & (df_phylogeny == 0)).sum().sum()
total_NAto0 = ((df_bin_withNA3_for_flipping == 3) & (df_phylogeny == 0)).sum().sum()
total_NAto1 = ((df_bin_withNA3_for_flipping == 3) & (df_phylogeny == 1)).sum().sum()

omega_final = total_FP_flipping + params['fnfp_ratio'] * total_FN_flipping

logger.info("")
logger.info("  ┌─────────────────────────────────────────────────────────────────────┐")
logger.info("  │              WEIGHTED DISCORDANCE INDEX (FINAL)                    │")
logger.info("  ├─────────────────────────────────────────────────────────────────────┤")
logger.info(f"  │  Weighted Discordance Index (Omega)        : {omega_final:>10.4f}      │")
logger.info(f"  │    - delta_FP discordance                  : {total_FP_flipping:>10}          │")
logger.info(f"  │    - delta_FN discordance                  : {total_FN_flipping:>10}          │")
logger.info(f"  │    - NA->0 imputations                     : {total_NAto0:>10}       │")
logger.info(f"  │    - NA->1 imputations                     : {total_NAto1:>10}       │")
logger.info(f"  │    - FN/FP weight (lambda)                 : {params['fnfp_ratio']:>10.1f}      │")
logger.info("  └─────────────────────────────────────────────────────────────────────┘")

# ---- Save total flipping counts ----
df_total_flipping_count = pd.DataFrame({
    'total_delta_FP': [total_FP_flipping],
    'total_delta_FN': [total_FN_flipping],
    'total_NA_to_0': [total_NAto0],
    'total_NA_to_1': [total_NAto1],
    'weighted_discordance_index_Omega': [omega_final]
})
df_total_flipping_count.to_csv(os.path.join(phylo_dir, "df_total_flipping_count.txt"), sep="\t", index=False)

# ---- Save per-site flipping counts ----
df_flip_counts_tree = calculate_flip_counts_per_site(df_bin_withNA3_for_flipping, df_phylogeny)
df_flip_counts_tree.to_csv(os.path.join(phylo_dir, "df_flipping_count_for_each_mut.txt"), sep="\t", index=True)

logger.info(f"The shape of final_cleaned_M_full.shape: {final_cleaned_M_full.shape}")
logger.info("")


# ----------------------------------------------------------------------------
# Step 9.5: Tree format and clone information
# ----------------------------------------------------------------------------
logger.info("STEP 9.5: Tree format and clone information")
logger.info("-" * 80)

# ---- Export tree as JSON ----
tree_dict = tree_to_dict(T_full)

with open(os.path.join(phylo_dir, 'final_cleaned_tree_node.json'), 'w') as f:
    json.dump(tree_dict, f, indent=4)

# ---- Export tree as text ----
T_full.save_to_file(os.path.join(phylo_dir, 'final_cleaned_tree_node.txt'))

# ---- Assign clone labels to cells ----
mutation_clones = get_mutation_clone_and_backbone_mut_as_keys_by_first_level_with_frequency(T_full, I_attached)
df_barcode_clones = assign_clone_labels(M_full, mutation_clones)

df_barcode_clones.to_csv(os.path.join(phylo_dir, "df_barcode_clones_from_phylo_tree.csv"), sep=',', index=False)

# ---- Export tree in Newick format ----
try:
    newick_str = tree_to_newick(T_full)
    with open(os.path.join(phylo_dir, 'final_cleaned_tree.newick'), 'w') as f:
        f.write(newick_str + ';')
    logger.info("  Exported Newick format: final_cleaned_tree.newick")
except NameError:
    logger.warning("  tree_to_newick function not available, skipping Newick export")

logger.info("")


# ----------------------------------------------------------------------------
# Step 9.6: Save run summary to final_results
# ----------------------------------------------------------------------------
logger.info("STEP 9.6: Save run summary")
logger.info("-" * 80)

run_summary_file = os.path.join(outputpath_04, "run_summary.txt")

with open(run_summary_file, 'w') as f:
    f.write("=" * 80 + "\n")
    f.write("PhyloSOLID RUN SUMMARY\n")
    f.write("=" * 80 + "\n\n")
    f.write(f"Sample ID: {sampleid}\n")
    f.write(f"Run date: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
    f.write(f"Seed: {args.seed}\n")
    f.write(f"Output path: {outputpath_04}\n\n")
    
    f.write("-" * 80 + "\n")
    f.write("PARAMETERS\n")
    f.write("-" * 80 + "\n")
    f.write(f"  is_predict_germ: {is_predict_germ}\n")
    f.write(f"  is_detect_passtree_by_dp: {is_detect_passtree_by_dp}\n")
    f.write(f"  is_filter_quality: {is_filter_quality}\n")
    f.write(f"  cv_rank_thresh: {cv_rank_thresh_str}\n")
    if is_search_mode:
        f.write(f"  cv_thresholds_searched: {cv_thresholds}\n")
        f.write(f"  optimal_cv_threshold: {optimal_cv if 'optimal_cv' in locals() else 'N/A'}\n")
    else:
        f.write(f"  cv_threshold: {cv_value if 'cv_value' in locals() else 'N/A'}\n")
    f.write(f"  remove_artifact_mutations: {remove_artifact_mutations}\n")
    f.write(f"  fnfp_ratio: {params['fnfp_ratio']}\n")
    f.write(f"  fp_ratio_cutoff_within_subclone: {params['fp_ratio_cutoff_within_subclone']}\n")
    f.write(f"  fp_ratio_cutoff_across_tree: {params['fp_ratio_cutoff_across_tree']}\n")
    f.write(f"  fn_ratio_cutoff_across_tree: {params['fn_ratio_cutoff_across_tree']}\n")
    f.write(f"  fp_ratio_persite_cutoff: {params['fp_ratio_persite_cutoff']}\n\n")
    
    f.write("-" * 80 + "\n")
    f.write("INPUT STATISTICS\n")
    f.write("-" * 80 + "\n")
    f.write(f"  Total mutations: {len(all_mutations)}\n")
    f.write(f"  Germline removed: {len(predicted_germline_mutations)}\n")
    f.write(f"  Somatic mutations: {len(somatic_mutations) if somatic_mutations else 'N/A'}\n")
    f.write(f"  Scaffold mutations: {len(scaffold_mutations) if scaffold_mutations else 'N/A'}\n")
    f.write(f"  Accessory mutations: {len(attached_mutations) if attached_mutations else 'N/A'}\n")
    f.write(f"  No-group mutations: {len(no_group_mutations) if no_group_mutations else 'N/A'}\n")
    f.write(f"  Cells: {M_current.shape[0] if M_current is not None else 'N/A'}\n\n")
    
    f.write("-" * 80 + "\n")
    f.write("MUTATION CLASSIFICATION\n")
    f.write("-" * 80 + "\n")
    f.write(f"  M_scaffold: {len(scaffold_mutations) if scaffold_mutations else 'N/A'}\n")
    f.write(f"  M_accessory (integrated): {len(attached_mutations) if attached_mutations else 'N/A'}\n")
    f.write(f"  M_artifact (REMOVED): {len(to_be_removed_mutations_by_fp_mutations_cross_all_cells) if to_be_removed_mutations_by_fp_mutations_cross_all_cells else 'N/A'}\n")
    f.write(f"  M_root (root-assigned): {len(root_mutations) if root_mutations else 'N/A'}\n")
    f.write(f"  M_ambiguous (conflict): {len(all_conflict_mutations) if all_conflict_mutations else 'N/A'}\n\n")
    
    f.write("-" * 80 + "\n")
    f.write("CELL CLASSIFICATION\n")
    f.write("-" * 80 + "\n")
    f.write(f"  C_resolved (final): {final_cleaned_M_full.shape[0]}\n")
    f.write(f"  C_orphan (REMOVED): {len(to_be_removed_cells) if to_be_removed_cells else 'N/A'}\n")
    f.write(f"  C_chimeric (REMOVED): {len(identified_doublet_cells) if identified_doublet_cells else 'N/A'}\n\n")
    
    f.write("-" * 80 + "\n")
    f.write("DISCORDANCE METRICS\n")
    f.write("-" * 80 + "\n")
    f.write(f"  Ω (pre-QC): {omega_before_qc:.4f}\n")
    f.write(f"  Ω (final): {omega_final:.4f}\n")
    f.write(f"  Ω (reduction): {omega_before_qc - omega_final:.4f}\n")
    f.write(f"  Ω (combined): {omega_before_qc + omega_final:.4f}\n\n")
    
    f.write("-" * 80 + "\n")
    f.write("OUTPUT FILES\n")
    f.write("-" * 80 + "\n")
    f.write(f"  Main results: {phylo_dir}/\n")
    f.write("  Files exported:\n")
    f.write("    - I_full_withNA3.txt\n")
    f.write("    - M_full_basedPivots.filtered_sites_inferred\n")
    f.write("    - final_cleaned_M_full_basedPivots.filtered_sites_inferred\n")
    f.write("    - final_cleaned_I_full_withNA3_for_circosPlot.txt\n")
    f.write("    - df_flipping_spots.txt\n")
    f.write("    - df_total_flipping_count.txt\n")
    f.write("    - df_flipping_count_for_each_mut.txt\n")
    f.write("    - final_cleaned_tree_node.json\n")
    f.write("    - final_cleaned_tree_node.txt\n")
    f.write("    - df_barcode_clones_from_phylo_tree.csv\n")
    f.write("=" * 80 + "\n")

logger.info(f"  Run summary saved to: {run_summary_file}")
logger.info("")


# ----------------------------------------------------------------------------
# Step 9 Summary
# ----------------------------------------------------------------------------
logger.info("=" * 80)
logger.info("Step 9 COMPLETED: All results exported successfully")
logger.info("-" * 80)
logger.info(f"  Output directory: {phylo_dir}")
logger.info("  Files exported:")
logger.info("    - I_full_withNA3.txt")
logger.info("    - M_full_basedPivots.filtered_sites_inferred")
logger.info("    - final_cleaned_M_full_basedPivots.filtered_sites_inferred")
logger.info("    - final_cleaned_I_full_withNA3_for_circosPlot.txt")
logger.info("    - df_flipping_spots.txt")
logger.info("    - df_total_flipping_count.txt")
logger.info("    - df_flipping_count_for_each_mut.txt")
logger.info("    - final_cleaned_tree_node.json")
logger.info("    - final_cleaned_tree_node.txt")
logger.info("    - df_barcode_clones_from_phylo_tree.csv")
logger.info("=" * 80)


# ============================================================================
# FINAL SUMMARY
# ============================================================================
logger.info("")
logger.info("=" * 80)
logger.info("PHYLOSOLID: PHYLOGENETIC RECONSTRUCTION COMPLETED")
logger.info("=" * 80)
logger.info("")
logger.info("  ┌─────────────────────────────────────────────────────────────────────┐")
logger.info("  │  MUTATION CLASSIFICATION SUMMARY                                   │")
logger.info("  ├─────────────────────────────────────────────────────────────────────┤")
logger.info(f"  │  M_scaffold                     : {len(scaffold_mutations) if scaffold_mutations else 'N/A':>10}│")
logger.info(f"  │  M_accessory (integrated)       : {len(attached_mutations) if attached_mutations else 'N/A':>10}│")
logger.info(f"  │  M_artifact (REMOVED)           : {len(to_be_removed_mutations_by_fp_mutations_cross_all_cells) if to_be_removed_mutations_by_fp_mutations_cross_all_cells else 'N/A':>10}│")
logger.info(f"  │  M_root (root-assigned)         : {len(root_mutations) if root_mutations else 'N/A':>10}│")
logger.info(f"  │  M_ambiguous (conflict)         : {len(all_conflict_mutations) if all_conflict_mutations else 'N/A':>10}│")
logger.info("  ├─────────────────────────────────────────────────────────────────────┤")
logger.info("  │  CELL CLASSIFICATION SUMMARY                                       │")
logger.info("  ├─────────────────────────────────────────────────────────────────────┤")
logger.info(f"  │  C_resolved (final)             : {final_cleaned_M_full.shape[0]:>10}│")
logger.info(f"  │  C_orphan (REMOVED)             : {len(to_be_removed_cells) if to_be_removed_cells else 'N/A':>10}│")
logger.info(f"  │  C_chimeric (REMOVED)           : {len(identified_doublet_cells) if identified_doublet_cells else 'N/A':>10}│")
logger.info("  ├─────────────────────────────────────────────────────────────────────┤")
logger.info("  │  DISCORDANCE METRICS                                               │")
logger.info("  ├─────────────────────────────────────────────────────────────────────┤")
logger.info(f"  │  Omega (pre-QC)                 : {omega_before_qc:.4f}      │")
logger.info(f"  │  Omega (final)                  : {omega_final:.4f}      │")
logger.info(f"  │  Omega (reduction)              : {omega_before_qc - omega_final:.4f}      │")
logger.info(f"  │  Omega (combined)               : {omega_before_qc + omega_final:.4f}      │")
logger.info("  │                                   (lower is better)               │")
logger.info("  └─────────────────────────────────────────────────────────────────────┘")
logger.info("")
logger.info("=" * 80)
logger.info("PhyloSOLID completed successfully!")
logger.info("=" * 80)


# ------------------------------
# End of Process
# ------------------------------
finish_time = time.perf_counter()
print("Program finished in {:.4f} seconds".format(finish_time - start_time))

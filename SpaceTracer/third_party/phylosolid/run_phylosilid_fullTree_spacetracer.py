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
Update: 2026/08/21 - Fixed parallel CV search determinism with spawn context
Update: 2026/09/08 - Added pruning_confidence criterion for optimal CV selection (aligned with scRNA version)

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
    ├── logs/
    ├── cv_search_summary.csv
    └── cv_selection_report.txt

"""


import time
start_time = time.perf_counter()

import os
os.environ['PYTHONHASHSEED'] = '42'
# Force matplotlib to use non-interactive backend (required for parallel execution)
os.environ['MPLBACKEND'] = 'Agg'
os.environ['QT_QPA_PLATFORM'] = 'offscreen'

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
import pickle
import multiprocessing as mp
from functools import partial

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
import sys

# Configure root logger format (includes PID for tracking)
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] [PID:%(process)d] %(message)s"
)

# Get root logger
root_logger = logging.getLogger()


def setup_main_log_file(outputpath):
    """Set up main log file for the root logger."""
    log_dir = os.path.join(outputpath, "logs")
    os.makedirs(log_dir, exist_ok=True)
    main_log_file = os.path.join(log_dir, "run_results.log")
    
    # Remove existing file handlers to avoid duplication
    for handler in root_logger.handlers[:]:
        if isinstance(handler, logging.FileHandler):
            root_logger.removeHandler(handler)
    
    # File handler for main log
    file_handler = logging.FileHandler(main_log_file, mode='w')
    file_handler.setLevel(logging.INFO)
    file_handler.setFormatter(logging.Formatter(
        "%(asctime)s [%(levelname)s] [PID:%(process)d] %(message)s"
    ))
    root_logger.addHandler(file_handler)
    
    # Console handler (keep terminal output)
    # Remove existing console handlers
    for handler in root_logger.handlers[:]:
        if isinstance(handler, logging.StreamHandler) and handler.stream == sys.stdout:
            root_logger.removeHandler(handler)
    
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(logging.INFO)
    console_handler.setFormatter(logging.Formatter(
        "%(asctime)s [%(levelname)s] %(message)s"
    ))
    root_logger.addHandler(console_handler)
    
    root_logger.info(f"Main log file: {main_log_file}")
    return main_log_file


# ------------------------------
# Project parameters and file paths
# ------------------------------
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
parser.add_argument("--cv_rank_thresh", default="auto", type=str, 
                    help="""CV rank threshold for coverage-based filtration.
                    Options:
                      - Single value: '0.3' (run once with this value)
                      - Comma-separated: '0.3,0.5,0.7' (search over these values)
                      - Range: '0.3-0.7:0.1' (from 0.3 to 0.7 with step 0.1)
                      - 'auto': use default presets [0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]""")
parser.add_argument("--remove_artifact_mutations", default="yes", choices=["yes", "no"], type=str, help="Select 'yes' or 'no' to determine whether to permanently remove artifact mutations.")
parser.add_argument("--seed", default=42, type=int, help="Random seed for reproducibility")
parser.add_argument("--max_workers", default=None, type=int, help="Maximum number of parallel workers for CV search. If not specified, uses CPU count - 1.")

# Pruning confidence threshold (fixed)
# pruning_ratio = omega_reduced / omega_final
# Confident if pruning_ratio < 10.0 (the "1/3 rule")
PRUNING_CONFIDENCE_THRESHOLD = 10.0

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
        'auto'          -> [0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    """
    cv_rank_thresh_str = cv_rank_thresh_str.strip()
    
    # Case 1: 'auto' -> use default presets
    if cv_rank_thresh_str.lower() == 'auto':
        return [0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    
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
max_workers = args.max_workers

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

# Setup main log file
main_log_file = setup_main_log_file(outputpath)

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
logger.info(f"max_workers: {max_workers if max_workers else 'auto (CPU count - 1)'}")
logger.info(f"pruning_confidence_threshold: {PRUNING_CONFIDENCE_THRESHOLD} (pruning_ratio < {PRUNING_CONFIDENCE_THRESHOLD} = confident)")
logger.info("")
logger.info("Directory structure:")
logger.info(f"  01_germline_filter: {outputpath_01}")
logger.info(f"  02_scaffold_builder: {outputpath_02}")
logger.info(f"  03_mutation_integrator: {outputpath_03}")
logger.info(f"  04_final_results: {outputpath_04}")
logger.info(f"  Main log file: {main_log_file}")


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
    # "cv_rank_thresh" will be set dynamically
    
    # 3.3 Consensus correlation graph
    "consensus_runs": 100,
    "consensus_clone_freq_thresh": 0.1,
    "resolution_of_graph": None,  # "None" indicates automatic calculation.
    "min_resolution": 0.5,        # The minimum value when the data is excellent
    "max_resolution": 2.0,        # The maximum value when the data range is extremely large
    
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


##### Output binary matrix
I_raw.to_csv(os.path.join(outputpath, "I_raw.txt"), sep="\t")
# NA=0
I_raw_withNA0 = I_raw.replace({np.nan: 0}).astype(int)
I_raw_withNA0.to_csv(os.path.join(outputpath, "I_raw_withNA0.txt"), sep="\t")
I_raw_withNA0_T = I_raw_withNA0.T
I_raw_withNA0_T.to_csv(os.path.join(outputpath, "I_raw_withNA0_T.txt"), sep="\t")
# NA=3
I_raw_withNA3 = I_raw.replace({np.nan: 3}).astype(int)
I_raw_withNA3.to_csv(os.path.join(outputpath, "I_raw_withNA3.txt"), sep="\t")
I_raw_withNA3_T = I_raw_withNA3.T
I_raw_withNA3_T.to_csv(os.path.join(outputpath, "I_raw_withNA3_T.txt"), sep="\t")


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
# Step 2: Classifier for identifying mosaic mutations (SKIPPED for spacetracer)
# ------------------------------
logger.info("===== Step2: Classifier ... SKIPPED (using all mutations as candidates)")

# For spacetracer mode, we skip the classifier and use all mutations as candidates
candidate_mutations = list(I_raw.columns)

P_candidate = P[candidate_mutations].copy()
V_candidate = V[candidate_mutations].copy()
A_candidate = A[candidate_mutations].copy()
C_candidate = C[candidate_mutations].copy()
I_candidate = I[candidate_mutations].copy()
df_reads_candidate = df_reads[candidate_mutations].copy()


# ------------------------------
# Step 3: Germline filtering
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
    dict : Result dictionary containing all outputs
    """
    # Create run-specific output directories
    cv_label = str(cv_value).replace('.', '_')
    run_outputpath_scaffold = os.path.join(outputpath_02, f"CV{cv_label}")
    run_outputpath_full = os.path.join(outputpath_03, f"CV{cv_label}")
    os.makedirs(run_outputpath_scaffold, exist_ok=True)
    os.makedirs(run_outputpath_full, exist_ok=True)
    
    # ============================================================
    # Setup CV-specific logging (does NOT propagate to root logger)
    # ============================================================
    log_dir = os.path.join(output_dir, "logs")
    os.makedirs(log_dir, exist_ok=True)
    log_file = os.path.join(log_dir, f"cv_{cv_label}.log")
    
    # Create CV-specific logger (does not propagate to root)
    cv_logger = logging.getLogger(f"CV_{cv_label}")
    cv_logger.setLevel(logging.INFO)
    cv_logger.propagate = False  # CRITICAL: prevents log messages from going to root logger
    
    # Clear existing handlers
    cv_logger.handlers = []
    
    # File handler - each CV has its own log file
    file_handler = logging.FileHandler(log_file, mode='w')
    file_handler.setLevel(logging.INFO)
    file_handler.setFormatter(logging.Formatter(
        "%(asctime)s [%(levelname)s] [PID:%(process)d] %(message)s"
    ))
    cv_logger.addHandler(file_handler)
    
    # Console handler - with CV label for easy identification
    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.INFO)
    console_handler.setFormatter(logging.Formatter(
        f"[CV{cv_label}] %(asctime)s [%(levelname)s] %(message)s"
    ))
    cv_logger.addHandler(console_handler)
    
    # Log CV start
    cv_logger.info("=" * 80)
    cv_logger.info(f"=== STARTING CV THRESHOLD SEARCH: cv_rank_thresh = {cv_value} ===")
    cv_logger.info("=" * 80)
    cv_logger.info(f"  Scaffold output: {run_outputpath_scaffold}")
    cv_logger.info(f"  Mutation integrator output: {run_outputpath_full}")
    cv_logger.info(f"  Log file: {log_file}")
    cv_logger.info("=" * 80)
    
    # Create params copy and set current CV value
    params_local = copy.deepcopy(params)
    params_local['cv_rank_thresh'] = cv_value
    
    # Store initial data dimensions for statistics
    initial_cells = I_somatic.shape[0]
    initial_muts = I_somatic.shape[1]
    initial_total = initial_cells + initial_muts
    
    try:
        # ---- Scaffold builder ----
        cv_logger.info("")
        cv_logger.info("STEP 4: Building scaffold tree...")
        cv_logger.info("-" * 40)
        
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
            immune_mutations=immune_mutations,
            logger_obj=cv_logger
        )
        
        # build_scaffold_tree returns 14 values
        (T_scaffold, M_scaffold, df_flipping_spots, df_total_flipping_count, 
         final_cleaned_I_selected_withNA3, final_cleaned_M_scaffold, 
         backbone_mutations, mutation_group, spots_to_split, group_mutations, 
         no_group_mutations, remained_mutations, conflict_mutations, root_mutations) = results_of_scaffold
        
        scaffold_mutations = list(M_scaffold.columns)
        cv_logger.info(f"  Scaffold tree built: {len(scaffold_mutations)} mutations")
        cv_logger.info("")
        
        # ---- Step 5: DP pass tree & prepare data ----
        cv_logger.info("STEP 5: Dynamic programming pass tree...")
        cv_logger.info("-" * 40)
        
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
        
        cv_logger.info(f"  Pass tree mutations: {len(passtree_mutations)}")
        cv_logger.info(f"  Attached mutations: {len(attached_mutations)}")
        cv_logger.info("")
        
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
        cv_logger.info("STEP 6-8: Building fully resolved tree...")
        cv_logger.info("-" * 40)
        
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
            logger_obj=cv_logger,
            cv_value=cv_value,
            export_phylo=True
        )
        
        cv_logger.info(f"  Fully resolved tree built: {M_current.shape[0]} cells, {M_current.shape[1]} mutations")
        cv_logger.info("")
        
        # ---- Compute omega_final (post-QC) and flipping counts ----
        cv_logger.info("Computing final metrics...")
        cv_logger.info("-" * 40)
        
        M_for_omega = M_current.drop(columns=['ROOT'], errors='ignore')
        mutations_on_tree_for_omega = M_for_omega.columns.to_series().apply(lambda x: x.split("|")).explode().unique().tolist()
        M_for_omega = split_merged_columns(M_for_omega, mutations_on_tree_for_omega)
        
        M_for_omega_clean = M_for_omega.loc[:, (M_for_omega != 0).any(axis=0)]
        M_for_omega_clean = M_for_omega_clean.loc[(M_for_omega_clean != 0).any(axis=1)]
        
        I_for_omega = I_attached.loc[M_for_omega_clean.index, M_for_omega_clean.columns].replace({np.nan: 3}).astype(int)
        
        # Count flipping stats
        N_deltaFP = ((I_for_omega == 1) & (M_for_omega_clean == 0)).sum().sum()
        N_deltaFN = ((I_for_omega == 0) & (M_for_omega_clean == 1)).sum().sum()
        N_NAto0 = ((I_for_omega == 3) & (M_for_omega_clean == 0)).sum().sum()
        N_NAto1 = ((I_for_omega == 3) & (M_for_omega_clean == 1)).sum().sum()
        
        omega_final = N_deltaFP + params_local['fnfp_ratio'] * N_deltaFN
        
        # ---- Calculate pruning metrics ----
        final_cells = M_for_omega_clean.shape[0]
        final_muts = M_for_omega_clean.shape[1]
        final_total = final_cells + final_muts
        
        omega_reduced = omega_before_qc - omega_final
        pruning_ratio = omega_reduced / omega_final if omega_final > 0 else float('inf')
        is_confident = pruning_ratio < PRUNING_CONFIDENCE_THRESHOLD
        retention_rate = final_total / initial_total if initial_total > 0 else 0.0
        
        # ---- Summary for this CV ----
        cv_logger.info("")
        cv_logger.info("=" * 80)
        cv_logger.info(f"CV THRESHOLD {cv_value} COMPLETED")
        cv_logger.info("=" * 80)
        cv_logger.info(f"  Omega (pre-QC):        {omega_before_qc:.4f}")
        cv_logger.info(f"  Omega (final):         {omega_final:.4f}")
        cv_logger.info(f"  Omega (reduced):       {omega_reduced:.4f}")
        cv_logger.info(f"  Pruning ratio:         {pruning_ratio:.4f} {'✓ Confident' if is_confident else '✗ Over-pruned'}")
        cv_logger.info("")
        cv_logger.info(f"  Flipping stats (final):")
        cv_logger.info(f"    delta_FP:            {N_deltaFP}")
        cv_logger.info(f"    delta_FN:            {N_deltaFN}")
        cv_logger.info(f"    NA->0:               {N_NAto0}")
        cv_logger.info(f"    NA->1:               {N_NAto1}")
        cv_logger.info("")
        cv_logger.info(f"  Cells:                 {M_current.shape[0]}")
        cv_logger.info(f"  Mutations:             {M_current.shape[1]}")
        cv_logger.info(f"  Retention rate:        {retention_rate:.4f}")
        cv_logger.info("=" * 80)
        cv_logger.info(f"  Log file: {log_file}")
        cv_logger.info("=" * 80)
        
        return {
            'cv_value': cv_value,
            'success': True,
            'error': None,
            'scaffold_dir': run_outputpath_scaffold,
            'mutint_dir': run_outputpath_full,
            # Omega metrics
            'omega_pre_qc': omega_before_qc,
            'omega_final': omega_final,
            'omega_sum': omega_before_qc + omega_final,
            'omega_reduced': omega_reduced,
            'pruning_ratio': pruning_ratio,
            'is_confident': is_confident,
            # Flipping counts
            'N_deltaFP': N_deltaFP,
            'N_deltaFN': N_deltaFN,
            'N_NAto0': N_NAto0,
            'N_NAto1': N_NAto1,
            # Data dimensions
            'initial_cells': initial_cells,
            'initial_muts': initial_muts,
            'initial_total': initial_total,
            'final_cells': final_cells,
            'final_muts': final_muts,
            'final_total': final_total,
            'retention_rate': retention_rate,
            'scaffold_count': len(scaffold_mutations),
            # Tree objects
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
        cv_logger.error(f"Error for cv_rank_thresh={cv_value}: {e}")
        import traceback
        traceback.print_exc()
        cv_logger.error(traceback.format_exc())
        
        return {
            'cv_value': cv_value,
            'success': False,
            'error': str(e),
            'scaffold_dir': None,
            'mutint_dir': None,
            'omega_pre_qc': float('inf'),
            'omega_final': float('inf'),
            'omega_sum': float('inf'),
            'omega_reduced': float('inf'),
            'pruning_ratio': float('inf'),
            'is_confident': False,
            'N_deltaFP': 0,
            'N_deltaFN': 0,
            'N_NAto0': 0,
            'N_NAto1': 0,
            'initial_cells': 0,
            'initial_muts': 0,
            'initial_total': 0,
            'final_cells': 0,
            'final_muts': 0,
            'final_total': 0,
            'retention_rate': 0.0,
            'scaffold_count': 0,
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
# Wrapper function for parallel execution (MUST be at module level)
# ============================================================
def run_single_cv(cv_value):
    """
    Wrapper function for parallel execution of a single CV threshold.
    """
    import random
    import numpy as np
    
    base_seed = 42
    specific_seed = int(base_seed + cv_value * 100000)
    random.seed(specific_seed)
    np.random.seed(specific_seed)
    os.environ['PYTHONHASHSEED'] = str(specific_seed)
    
    from src.reproducibility import set_seed
    set_seed(specific_seed, verbose=False)
    
    try:
        result = run_pipeline_for_cv_threshold(cv_value, outputpath)
        if result['success']:
            return result
        else:
            return {
                'cv_value': cv_value,
                'success': False,
                'error': result.get('error', 'Unknown error'),
                'scaffold_dir': None,
                'mutint_dir': None,
            }
    except Exception as e:
        logger.error(f"Critical error for CV={cv_value}: {e}")
        import traceback
        traceback.print_exc()
        return {
            'cv_value': cv_value,
            'success': False,
            'error': str(e),
            'scaffold_dir': None,
            'mutint_dir': None,
        }


# ============================================================
# Parallel execution function for CV threshold search
# ============================================================
def run_parallel_cv_search(cv_thresholds, max_workers=None):
    if max_workers is None:
        max_workers = max(1, mp.cpu_count() - 1)
    
    logger.info("=" * 80)
    logger.info(f"Starting parallel CV search with {max_workers} workers...")
    logger.info(f"Testing {len(cv_thresholds)} thresholds: {cv_thresholds}")
    logger.info(f"Pruning confidence threshold: pruning_ratio < {PRUNING_CONFIDENCE_THRESHOLD}")
    logger.info("=" * 80)
    
    ctx = mp.get_context('spawn')
    
    with ctx.Pool(processes=max_workers) as pool:
        results = []
        for i, result in enumerate(pool.imap_unordered(run_single_cv, cv_thresholds)):
            results.append(result)
            if result['success']:
                confident = "✓" if result.get('is_confident', False) else "⚠️"
                logger.info(f"[{i+1}/{len(cv_thresholds)}] {confident} CV={result['cv_value']}: "
                           f"Omega_final={result.get('omega_final', 'N/A'):.4f}, "
                           f"Pruning_ratio={result.get('pruning_ratio', 'N/A'):.4f}")
            else:
                logger.info(f"[{i+1}/{len(cv_thresholds)}] ✗ CV={result['cv_value']}: Failed")
    
    return results


# ============================================================
# Selection function: Choose best CV using pruning_confidence
# ============================================================
def select_best_cv(df_valid):
    """
    Select the best CV threshold using the pruning_confidence criterion.
    
    Logic:
        1. Calculate pruning_ratio = omega_reduced / omega_final
        2. Confident if pruning_ratio < PRUNING_CONFIDENCE_THRESHOLD
        3. If there are confident trees: select the one with LOWEST omega_final
        4. If NO confident trees: select the one with LOWEST omega_pre_qc
    
    Parameters
    ----------
    df_valid : pd.DataFrame
        DataFrame with valid results (success=True, finite values)
    
    Returns
    -------
    tuple : (best_idx, best_cv, best_score, df_sorted, selection_details, summary_stats)
    """
    if df_valid.empty:
        return None, None, None, None, None, None
    
    df_valid = df_valid.copy()
    
    # Calculate pruning metrics
    df_valid['omega_sum'] = df_valid['omega_pre_qc'] + df_valid['omega_final']
    df_valid['omega_reduced'] = df_valid['omega_pre_qc'] - df_valid['omega_final']
    df_valid['pruning_ratio'] = df_valid['omega_reduced'] / df_valid['omega_final']
    df_valid['is_confident'] = df_valid['pruning_ratio'] < PRUNING_CONFIDENCE_THRESHOLD
    
    # Separate confident and over-pruned
    df_confident = df_valid[df_valid['is_confident'] == True]
    df_overpruned = df_valid[df_valid['is_confident'] == False]
    
    n_confident = len(df_confident)
    n_overpruned = len(df_overpruned)
    
    selection_details = []
    selection_details.append("=" * 80)
    selection_details.append("CV SELECTION: Pruning Confidence Criterion")
    selection_details.append("-" * 80)
    selection_details.append(f"  Confidence rule: pruning_ratio = omega_reduced / omega_final < {PRUNING_CONFIDENCE_THRESHOLD}")
    selection_details.append("  discordance must remain in the final tree.")
    selection_details.append("")
    selection_details.append(f"  Found {n_confident} confident trees, {n_overpruned} over-pruned trees")
    selection_details.append("")
    selection_details.append("  ┌────────────────────────────────────────────────────────────────────────────────────────────┐")
    selection_details.append("  │  CV    Omega_pre    Omega_final    Omega_reduced   Pruning_ratio   Status                │")
    selection_details.append("  ├────────────────────────────────────────────────────────────────────────────────────────────┤")
    
    for idx, row in df_valid.iterrows():
        cv = row['cv_value']
        omega_pre = row['omega_pre_qc']
        omega_final = row['omega_final']
        omega_reduced = row['omega_reduced']
        pruning_ratio = row['pruning_ratio']
        is_confident = row['is_confident']
        status = "✓ Confident" if is_confident else "✗ Over-pruned"
        selection_details.append(
            f"  │  {cv:5.1f}  {omega_pre:>11.1f}  {omega_final:>13.1f}  {omega_reduced:>15.1f}  {pruning_ratio:>13.4f}  {status:>19} │"
        )
    selection_details.append("  └────────────────────────────────────────────────────────────────────────────────────────────┘")
    selection_details.append("")
    
    # -------- Step 1: If there are confident trees, select by LOWEST omega_final --------
    if not df_confident.empty:
        df_sorted = df_confident.sort_values(
            by=['omega_final', 'omega_sum'], 
            ascending=[True, True]
        )
        best_idx = df_sorted.index[0]
        best_cv = df_sorted.loc[best_idx, 'cv_value']
        best_score = df_sorted.loc[best_idx, 'omega_final']
        
        best_row = df_sorted.loc[best_idx]
        selection_details.append(f"  ✓ Selected from confident trees by LOWEST omega_final: CV = {best_cv}")
        selection_details.append(f"      omega_final     = {best_row['omega_final']:.4f}")
        selection_details.append(f"      omega_pre_qc    = {best_row['omega_pre_qc']:.4f}")
        selection_details.append(f"      omega_reduced   = {best_row['omega_reduced']:.4f}")
        selection_details.append(f"      pruning_ratio   = {best_row['pruning_ratio']:.4f} (confident ✓)")
        selection_details.append("")
        selection_details.append("=" * 80)
        
        summary_stats = {
            'selection_method': 'confident_trees_lowest_omega_final',
            'n_confident': n_confident,
            'n_overpruned': n_overpruned,
            'selected_cv': best_cv,
            'selected_omega_final': best_row['omega_final'],
            'selected_omega_pre_qc': best_row['omega_pre_qc'],
            'selected_omega_reduced': best_row['omega_reduced'],
            'selected_pruning_ratio': best_row['pruning_ratio'],
            'selected_retention_rate': best_row['retention_rate'],
            'selected_scaffold_count': best_row['scaffold_count'],
        }
        
        return best_idx, best_cv, best_score, df_sorted, selection_details, summary_stats
    
    # -------- Step 2: If NO confident trees, select by LOWEST omega_pre_qc --------
    else:
        df_sorted = df_valid.sort_values(
            by=['omega_pre_qc', 'omega_final'], 
            ascending=[True, True]
        )
        best_idx = df_sorted.index[0]
        best_cv = df_sorted.loc[best_idx, 'cv_value']
        best_score = df_sorted.loc[best_idx, 'omega_pre_qc']
        
        best_row = df_sorted.loc[best_idx]
        selection_details.append("  ⚠️ WARNING: No confident trees found (all have pruning_ratio >= 10.0)")
        selection_details.append("     This indicates potential over-pruning across all CV thresholds.")
        selection_details.append(f"  ✓ Falling back to LOWEST omega_pre_qc: CV = {best_cv}")
        selection_details.append(f"      omega_pre_qc    = {best_row['omega_pre_qc']:.4f}")
        selection_details.append(f"      omega_final     = {best_row['omega_final']:.4f}")
        selection_details.append(f"      omega_reduced   = {best_row['omega_reduced']:.4f}")
        selection_details.append(f"      pruning_ratio   = {best_row['pruning_ratio']:.4f} (over-pruned ⚠️)")
        selection_details.append("")
        selection_details.append("  Note: When all trees are over-pruned, we select the one with")
        selection_details.append("  the lowest pre-QC discordance, as it represents the least")
        selection_details.append("  aggressive initial model fit.")
        selection_details.append("")
        selection_details.append("=" * 80)
        
        summary_stats = {
            'selection_method': 'fallback_lowest_omega_pre_qc',
            'n_confident': n_confident,
            'n_overpruned': n_overpruned,
            'selected_cv': best_cv,
            'selected_omega_final': best_row['omega_final'],
            'selected_omega_pre_qc': best_row['omega_pre_qc'],
            'selected_omega_reduced': best_row['omega_reduced'],
            'selected_pruning_ratio': best_row['pruning_ratio'],
            'selected_retention_rate': best_row['retention_rate'],
            'selected_scaffold_count': best_row['scaffold_count'],
        }
        
        return best_idx, best_cv, best_score, df_sorted, selection_details, summary_stats


# ============================================================
# Helper function: Save comprehensive results (SIMPLIFIED)
# ============================================================
def save_comprehensive_results(search_results, outputpath, sampleid,
                               best_cv=None, best_score=None, 
                               selection_details=None, best_row=None):
    """
    Save comprehensive results - SIMPLIFIED VERSION.
    Only creates:
        - cv_search_summary.csv (all CV metrics)
        - cv_selection_report.txt (selection details)
    
    Parameters
    ----------
    search_results : list
        List of result dictionaries for all CVs
    outputpath : str
        Base output directory
    sampleid : str
        Sample ID
    best_cv : float, optional
        Best CV threshold
    best_score : float, optional
        Best score
    selection_details : list, optional
        List of selection detail strings
    best_row : pd.Series, optional
        Best result row from df_sorted
    """
    
    # ---- 1. Master summary (one file, all CVs) ----
    master_rows = []
    for result in search_results:
        if result.get('success', False):
            row = {
                'cv': result['cv_value'],
                'selected': result['cv_value'] == best_cv if best_cv is not None else False,
                'omega_pre_qc': result.get('omega_pre_qc', float('inf')),
                'omega_final': result.get('omega_final', float('inf')),
                'omega_sum': result.get('omega_sum', float('inf')),
                'omega_reduced': result.get('omega_reduced', float('inf')),
                'pruning_ratio': result.get('pruning_ratio', float('inf')),
                'confident': result.get('is_confident', False),
                'scaffold_count': result.get('scaffold_count', 0),
                'retention_rate': result.get('retention_rate', 0.0),
                'N_deltaFP': result.get('N_deltaFP', 0),
                'N_deltaFN': result.get('N_deltaFN', 0),
                'N_NAto0': result.get('N_NAto0', 0),
                'N_NAto1': result.get('N_NAto1', 0),
                'final_cells': result.get('final_cells', 0),
                'final_muts': result.get('final_muts', 0),
            }
            master_rows.append(row)
    
    df_master = pd.DataFrame(master_rows)
    if not df_master.empty:
        df_master = df_master.sort_values('cv')
        df_master.to_csv(os.path.join(outputpath, "cv_search_summary.csv"), index=False)
        logger.info(f"  CV summary saved to: cv_search_summary.csv")
    
    # ---- 2. Selection report ----
    with open(os.path.join(outputpath, "cv_selection_report.txt"), 'w') as f:
        f.write("=" * 80 + "\n")
        f.write("CV SELECTION REPORT\n")
        f.write("=" * 80 + "\n\n")
        f.write(f"Date: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"Sample ID: {sampleid}\n")
        f.write(f"Total CVs tested: {len(master_rows)}\n")
        f.write(f"Confident pruning threshold: pruning_ratio < {PRUNING_CONFIDENCE_THRESHOLD}\n")
        f.write(f"Selection criterion: pruning_confidence\n\n")
        
        # Selection details
        if selection_details:
            f.write("-" * 80 + "\n")
            f.write("SELECTION DETAILS\n")
            f.write("-" * 80 + "\n")
            for detail in selection_details:
                f.write(detail + "\n")
            f.write("\n")
        
        # Selected CV
        if best_cv is not None and best_row is not None:
            f.write("-" * 80 + "\n")
            f.write("SELECTED OPTIMAL CV\n")
            f.write("-" * 80 + "\n")
            f.write(f"  CV:                 {best_cv}\n")
            f.write(f"  omega_pre_qc:       {best_row.get('omega_pre_qc', 0):.4f}\n")
            f.write(f"  omega_final:        {best_row.get('omega_final', 0):.4f}\n")
            f.write(f"  omega_reduced:      {best_row.get('omega_reduced', 0):.4f}\n")
            f.write(f"  pruning_ratio:      {best_row.get('pruning_ratio', 0):.4f}\n")
            f.write(f"  scaffold_count:     {best_row.get('scaffold_count', 0)}\n")
            f.write(f"  retention_rate:     {best_row.get('retention_rate', 0):.4f}\n")
            f.write("=" * 80 + "\n")
    
    logger.info(f"  Selection report saved to: cv_selection_report.txt")


# ============================================================
# Main execution: Single or search mode
# ============================================================

if __name__ == '__main__':
    from multiprocessing import freeze_support
    freeze_support()
    
    if is_search_mode:
        # ============================================================
        # SEARCH MODE: Parallel execution
        # ============================================================
        logger.info("=" * 80)
        logger.info("=== CV THRESHOLD SEARCH MODE (PARALLEL) ===")
        logger.info(f"Searching over {len(cv_thresholds)} values: {cv_thresholds}")
        logger.info(f"Pruning confidence threshold: pruning_ratio < {PRUNING_CONFIDENCE_THRESHOLD}")
        logger.info("=" * 80)
        
        # Run parallel search
        search_results = run_parallel_cv_search(
            cv_thresholds=cv_thresholds,
            max_workers=max_workers
        )
        
        # ---- Summarize results ----
        df_results = pd.DataFrame(search_results)
        
        # Filter valid results
        df_valid = df_results[df_results['success'] == True]
        df_valid = df_valid[df_valid['omega_final'].notna()]
        df_valid = df_valid[df_valid['omega_final'] != float('inf')]
        
        logger.info("=" * 80)
        logger.info("=== CV THRESHOLD SEARCH RESULTS ===")
        logger.info("=" * 80)
        logger.info("CV    Omega_pre   Omega_final   Omega_red   Pruning_ratio   Confident   Scaffold")
        logger.info("-" * 95)
        for idx, row in df_results.iterrows():
            if row['success']:
                is_confident = "✓" if row.get('is_confident', False) else "✗"
                logger.info(f"{row['cv_value']:5.1f}  {row['omega_pre_qc']:10.4f}  {row['omega_final']:10.4f}  "
                           f"{row.get('omega_reduced', 0):10.4f}  {row.get('pruning_ratio', 0):13.4f}  "
                           f"{is_confident:>9}  {row.get('scaffold_count', 0):>8}")
            else:
                logger.info(f"{row['cv_value']:5.1f}  {'FAILED':>10}  {'FAILED':>10}  {'FAILED':>10}  {'FAILED':>13}  {'FAILED':>9}  {'N/A':>8}")
        logger.info("=" * 80)
        
        # ---- Select the best CV ----
        if not df_valid.empty:
            best_idx, best_cv, best_score, df_sorted, selection_details, summary_stats = select_best_cv(df_valid)
            
            if best_idx is not None:
                best_row = df_sorted.loc[best_idx]
                
                logger.info("=" * 80)
                logger.info("=== CV SELECTION DETAILS ===")
                for detail in selection_details:
                    logger.info(f"  {detail}")
                logger.info("=" * 80)
                
                # ---- Extract variables from best result ----
                best_result = best_row.to_dict()
                T_current = best_result.get('T_current')
                M_current = best_result.get('M_current')
                root_mutations = best_result.get('root_mutations', [])
                all_conflict_mutations = best_result.get('all_conflict_mutations', [])
                omega_before_qc = best_result.get('omega_pre_qc', 0.0)
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
                best_scaffold_dir = best_result.get('scaffold_dir')
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
                
                if best_scaffold_dir and os.path.exists(best_scaffold_dir):
                    scaffold_dst = os.path.join(outputpath_04, "scaffold_summary")
                    os.makedirs(scaffold_dst, exist_ok=True)
                    for f in os.listdir(best_scaffold_dir):
                        src = os.path.join(best_scaffold_dir, f)
                        dst = os.path.join(scaffold_dst, f)
                        if os.path.isfile(src):
                            shutil.copy2(src, dst)
                
                # ---- Save comprehensive results (SIMPLIFIED) ----
                save_comprehensive_results(
                    search_results=search_results,
                    outputpath=outputpath,
                    sampleid=sampleid,
                    best_cv=best_cv,
                    best_score=best_score,
                    selection_details=selection_details,
                    best_row=best_row
                )
            
            else:
                logger.error("No valid results from threshold search!")
                sys.exit(1)
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
        
        result = run_pipeline_for_cv_threshold(cv_value, outputpath)
        
        if not result['success']:
            logger.error(f"Pipeline failed for cv_rank_thresh={cv_value}")
            sys.exit(1)
        
        # Extract variables
        T_current = result.get('T_current')
        M_current = result.get('M_current')
        root_mutations = result.get('root_mutations', [])
        all_conflict_mutations = result.get('all_conflict_mutations', [])
        omega_before_qc = result.get('omega_pre_qc', 0.0)
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
        optimal_cv = cv_value
        
        # Copy results to final_results
        best_mutint_dir = result.get('mutint_dir')
        if best_mutint_dir and os.path.exists(best_mutint_dir):
            for f in os.listdir(best_mutint_dir):
                src = os.path.join(best_mutint_dir, f)
                dst = os.path.join(outputpath_04, f)
                if os.path.isfile(src):
                    shutil.copy2(src, dst)
        
        # Save single CV results
        df_single = pd.DataFrame([{
            'cv_value': cv_value,
            'omega_pre_qc': result['omega_pre_qc'],
            'omega_final': result['omega_final'],
            'omega_sum': result['omega_sum'],
            'omega_reduced': result.get('omega_reduced', float('inf')),
            'pruning_ratio': result.get('pruning_ratio', float('inf')),
            'is_confident': result.get('is_confident', False),
            'N_deltaFP': result.get('N_deltaFP', 0),
            'N_deltaFN': result.get('N_deltaFN', 0),
            'N_NAto0': result.get('N_NAto0', 0),
            'N_NAto1': result.get('N_NAto1', 0),
            'scaffold_count': result.get('scaffold_count', 0),
            'retention_rate': result.get('retention_rate', 0.0),
            'scaffold_dir': result.get('scaffold_dir'),
            'mutint_dir': result.get('mutint_dir'),
        }])
        df_single.to_csv(os.path.join(outputpath, "cv_threshold_search_results.csv"), index=False)
        
        # Save selection report for single mode
        with open(os.path.join(outputpath, "cv_selection_report.txt"), 'w') as f:
            f.write("=" * 80 + "\n")
            f.write("CV SELECTION REPORT (SINGLE MODE)\n")
            f.write("=" * 80 + "\n\n")
            f.write(f"Sample ID: {sampleid}\n")
            f.write(f"CV threshold: {cv_value}\n")
            f.write(f"Pruning confidence threshold: pruning_ratio < {PRUNING_CONFIDENCE_THRESHOLD}\n")
            f.write(f"omega_pre_qc: {result['omega_pre_qc']:.4f}\n")
            f.write(f"omega_final: {result['omega_final']:.4f}\n")
            f.write(f"omega_reduced: {result.get('omega_reduced', float('inf')):.4f}\n")
            f.write(f"pruning_ratio: {result.get('pruning_ratio', float('inf')):.4f}\n")
            f.write(f"is_confident: {result.get('is_confident', False)}\n")
            f.write(f"Scaffold count: {result['scaffold_count']}\n")
            f.write("=" * 80 + "\n")
        
        with open(os.path.join(outputpath_04, "README.txt"), 'w') as f:
            f.write("=" * 80 + "\n")
            f.write("PhyloSOLID FINAL RESULTS\n")
            f.write("=" * 80 + "\n\n")
            f.write(f"Sample ID: {sampleid}\n")
            f.write(f"CV threshold: {cv_value}\n")
            f.write(f"Pruning confidence threshold: pruning_ratio < {PRUNING_CONFIDENCE_THRESHOLD}\n")
            f.write(f"omega_pre_qc: {result['omega_pre_qc']:.4f}\n")
            f.write(f"omega_final: {result['omega_final']:.4f}\n")
            f.write(f"pruning_ratio: {result.get('pruning_ratio', float('inf')):.4f}\n")
            f.write(f"Scaffold count: {result['scaffold_count']}\n")
            f.write("=" * 80 + "\n")
    
    
    # ============================================================
    # Step 9: Post-processing & output (using the final tree)
    # ============================================================
    
    logger.info("=" * 80)
    logger.info("Step 9: Post-processing & output")
    logger.info("=" * 80)
    logger.info("")
    
    # ---- Step 9.1: Prepare final data ----
    logger.info("STEP 9.1: Prepare final data")
    logger.info("-" * 80)
    
    M_current_filtered = M_current.drop(columns=['ROOT'], errors='ignore')
    
    for mut_on_root in root_mutations:
        M_current_filtered.insert(0, mut_on_root, 1)
    
    mutations_on_T_current = M_current_filtered.columns.to_series().apply(lambda x: x.split("|")).explode().unique().tolist()
    
    T_full = copy.deepcopy(T_current)
    M_full = split_merged_columns(M_current_filtered, mutations_on_T_current)
    
    logger.info("Final full-resolved tree:")
    print_tree_logger(T_full)
    
    logger.info(f"  Final tree cells: {M_full.shape[0]}")
    logger.info(f"  Final tree mutations: {M_full.shape[1]}")
    logger.info("")
    
    # ---- Step 9.2: Output results ----
    logger.info("STEP 9.2: Output result files")
    logger.info("-" * 80)
    
    phylo_dir = os.path.join(outputpath_04, "phylo")
    os.makedirs(phylo_dir, exist_ok=True)
    
    I_full_withNA3 = I_attached.replace({np.nan: 3}).astype(int)
    I_full_withNA3.to_csv(os.path.join(phylo_dir, "I_full_withNA3.txt"), sep="\t")
    
    WriteTfile(os.path.join(phylo_dir, "M_full_basedPivots.filtered_sites_inferred"), 
               M_full, M_full.index.tolist(), M_full.columns.tolist(), judge="yes")
    
    final_cleaned_M_full = M_full.loc[:, (M_full != 0).any(axis=0)]
    final_cleaned_M_full = final_cleaned_M_full.loc[(final_cleaned_M_full != 0).any(axis=1)]
    
    kept_rows = final_cleaned_M_full.index
    kept_cols = final_cleaned_M_full.columns
    
    final_cleaned_I_full_withNA3 = I_full_withNA3.loc[kept_rows, kept_cols]
    
    WriteTfile(os.path.join(phylo_dir, "final_cleaned_M_full_basedPivots.filtered_sites_inferred"), 
               final_cleaned_M_full, final_cleaned_M_full.index.tolist(), 
               final_cleaned_M_full.columns.tolist(), judge="yes")
    final_cleaned_I_full_withNA3.to_csv(
        os.path.join(phylo_dir, "final_cleaned_I_full_withNA3_for_circosPlot.txt"), sep="\t"
    )
    
    logger.info(f"  Output directory: {phylo_dir}")
    logger.info("")
    
    # ---- Step 9.3: Identify flipping spots ----
    logger.info("STEP 9.3: Identify flipping spots")
    logger.info("-" * 80)
    
    df_bin_withNA3_for_flipping = final_cleaned_I_full_withNA3.copy()
    df_phylogeny = final_cleaned_M_full.copy()
    
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
    
    if false_negative_flipping_spots.empty:
        false_negative_flipping_spots = {col: [] for col in df_bin_withNA3_for_flipping.columns}
    if false_positive_flipping_spots.empty:
        false_positive_flipping_spots = {col: [] for col in df_bin_withNA3_for_flipping.columns}
    if NAto1_flipping_spots.empty:
        NAto1_flipping_spots = {col: [] for col in df_bin_withNA3_for_flipping.columns}
    if NAto0_flipping_spots.empty:
        NAto0_flipping_spots = {col: [] for col in df_bin_withNA3_for_flipping.columns}
    
    df_flipping_spots = pd.DataFrame({
        'Mutation': df_bin_withNA3_for_flipping.columns,
        'delta_FN_spots': [', '.join(false_negative_flipping_spots.get(col, [])) for col in df_bin_withNA3_for_flipping.columns],
        'delta_FP_spots': [', '.join(false_positive_flipping_spots.get(col, [])) for col in df_bin_withNA3_for_flipping.columns],
        'NA_to_1_spots': [', '.join(NAto1_flipping_spots.get(col, [])) for col in df_bin_withNA3_for_flipping.columns],
        'NA_to_0_spots': [', '.join(NAto0_flipping_spots.get(col, [])) for col in df_bin_withNA3_for_flipping.columns]
    })
    df_flipping_spots.to_csv(os.path.join(phylo_dir, "df_flipping_spots.txt"), sep="\t", index=False)
    
    logger.info("")
    
    # ---- Step 9.4: Calculate total flipping counts ----
    logger.info("STEP 9.4: Calculate total flipping counts")
    logger.info("-" * 80)
    
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
    
    df_total_flipping_count = pd.DataFrame({
        'total_delta_FP': [total_FP_flipping],
        'total_delta_FN': [total_FN_flipping],
        'total_NA_to_0': [total_NAto0],
        'total_NA_to_1': [total_NAto1],
        'weighted_discordance_index_Omega': [omega_final]
    })
    df_total_flipping_count.to_csv(os.path.join(phylo_dir, "df_total_flipping_count.txt"), sep="\t", index=False)
    
    df_flip_counts_tree = calculate_flip_counts_per_site(df_bin_withNA3_for_flipping, df_phylogeny)
    df_flip_counts_tree.to_csv(os.path.join(phylo_dir, "df_flipping_count_for_each_mut.txt"), sep="\t", index=True)
    
    logger.info(f"The shape of final_cleaned_M_full.shape: {final_cleaned_M_full.shape}")
    logger.info("")
    
    # ---- Step 9.5: Tree format and clone information ----
    logger.info("STEP 9.5: Tree format and clone information")
    logger.info("-" * 80)
    
    tree_dict = tree_to_dict(T_full)
    with open(os.path.join(phylo_dir, 'final_cleaned_tree_node.json'), 'w') as f:
        json.dump(tree_dict, f, indent=4)
    
    T_full.save_to_file(os.path.join(phylo_dir, 'final_cleaned_tree_node.txt'))
    
    mutation_clones = get_mutation_clone_and_backbone_mut_as_keys_by_first_level_with_frequency(T_full, I_attached)
    df_barcode_clones = assign_clone_labels(M_full, mutation_clones)
    df_barcode_clones.to_csv(os.path.join(phylo_dir, "df_barcode_clones_from_phylo_tree.csv"), sep=',', index=False)
    
    logger.info("")
    
    
    # ---- Final summary ----
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
    logger.info(f"  │  Omega (reduced)                : {omega_before_qc - omega_final:.4f}      │")
    logger.info(f"  │  Pruning ratio                  : {(omega_before_qc - omega_final) / omega_final if omega_final > 0 else float('inf'):.4f}      │")
    logger.info(f"  │  Retention rate                 : {(final_cleaned_M_full.shape[0] + final_cleaned_M_full.shape[1]) / (I_somatic.shape[0] + I_somatic.shape[1]):.4f}      │")
    logger.info("  │                                   (pruning ratio < 10.0 = confident)│")
    logger.info("  └─────────────────────────────────────────────────────────────────────┘")
    logger.info("")
    if is_search_mode:
        logger.info(f"  ✓ CV threshold search completed: {len(cv_thresholds)} values tested")
        logger.info(f"  ✓ Optimal CV threshold: {optimal_cv}")
    logger.info("")
    logger.info("=" * 80)
    logger.info("PhyloSOLID completed successfully!")
    logger.info("=" * 80)
    
    finish_time = time.perf_counter()
    logger.info("Program finished in {:.4f} seconds".format(finish_time - start_time))
#!/usr/bin/env python3

"""
scaffold_builder.py

Build scaffold tree for PhyloSOLID (Methods Section 3).

Implements:
- 3.1 Initial filtration based on genotype posterior and mutant reads count
- 3.2 Coverage-based filtration for ubiquitous expression
- 3.3 Identification of "backbone clones" (largest clones independent to each other)
- 3.4 Penalty-based placement of non-backbone mutations onto the scaffold tree

Inputs:
- P: posterior probability matrix (cells x mutations)
- M: mutant allele frequency matrix (cells x mutations)  
- C: coverage matrix (cells x mutations)
- I: binary indicator matrix {0,1,NA} (cells x mutations)
- df_reads
- df_celltype
- params: dictionary of parameters with defaults

Outputs:
- T_s_root: root node of scaffold tree
- scaffold_mutations: list of scaffold mutations
- backbone: list of backbone mutations
- placements: dictionary of non-backbone placements
- G_consensus: consensus correlation graph
- clusters: list of mutation clusters
"""

import os
import shutil
import numpy as np
import pandas as pd
import networkx as nx
import logging
import copy
import math
import random
from tqdm import tqdm
from copy import deepcopy
from itertools import combinations
from collections import defaultdict, Counter
from typing import Set, List, Dict, Optional, Tuple, Any, Union
from anytree import Node, RenderTree
import scphylo as scp
from scphylo.pl._helper import (
    _add_barplot,
    _add_chromplot,
    _clonal_cell_mutation_list,
    _get_tree,
    _newick_info2_mutation_list,
)
from concurrent.futures import ThreadPoolExecutor, as_completed
from src.germline_filter import reorder_columns_by_mutant_stats
from src.reproducibility import set_seed, deterministic_permutation, get_seed

# Try to import optional dependencies
try:
    import igraph as ig
    import leidenalg as la
    _HAS_LEIDEN = True
except ImportError:
    _HAS_LEIDEN = False
    try:
        import community as community_louvain
    except ImportError:
        community_louvain = None

_VERBOSE_TREE_UPDATES = os.environ.get("PHYLOSOLID_VERBOSE_TREES", "").lower() in {"1", "true", "yes", "on"}




def _log_tree_snapshot(active_logger, label, tree):
    if not _VERBOSE_TREE_UPDATES:
        return
    active_logger.info(label)
    # print_tree(tree)


def _as_bool_mask(series, index):
    aligned = pd.Series(series, index=index) if not isinstance(series, pd.Series) else series.reindex(index)
    values = aligned.to_numpy(copy=False)
    
    if np.issubdtype(values.dtype, np.bool_):
        return pd.Series(values, index=index, copy=False)
    
    if np.issubdtype(values.dtype, np.integer):
        return pd.Series(values != 0, index=index)
    
    if np.issubdtype(values.dtype, np.floating):
        return pd.Series(np.isfinite(values) & (values != 0), index=index)
    
    numeric = pd.to_numeric(aligned, errors="coerce")
    numeric_values = numeric.to_numpy(dtype=np.float64, na_value=np.nan)
    return pd.Series(np.isfinite(numeric_values) & (numeric_values != 0), index=index)


def _build_conflict_mask(matrix, conflict_nodes, index):
    vec_conflicts = pd.Series(False, index=index)
    for conflict in conflict_nodes:
        vec_conflicts |= _as_bool_mask(matrix[conflict], index)
    return vec_conflicts


def _build_scaffold_lineage_fill_update(final_position, final_imputed_vec, parent_dict):
    if final_position is None or final_imputed_vec is None:
        return {"mutation_chain": (), "cells": ()}
    
    mutation_chain = tuple(get_full_mutnode_chain_with_anchor_scaffold(final_position["anchor"], parent_dict))
    cells = tuple(final_imputed_vec.index[final_imputed_vec == 1].tolist())
    return {
        "mutation_chain": mutation_chain,
        "cells": cells,
    }


def _apply_scaffold_lineage_fill_to_matrix(M_current, lineage_fill_update):
    mutation_chain = lineage_fill_update.get("mutation_chain", ())
    cells = lineage_fill_update.get("cells", ())
    
    if not mutation_chain or not cells:
        return M_current
    
    cells_list = list(cells)
    for mutation in mutation_chain:
        if mutation not in M_current.columns:
            continue
        column_values = M_current.loc[cells_list, mutation]
        cells_to_update = column_values.index[column_values == 0]
        if len(cells_to_update) > 0:
            M_current.loc[cells_to_update, mutation] = 1
    
    return M_current


def get_random_chromosome_position(snp_name, logger_obj=None):
    """
    Parse or generate random chromosome and position information based on snp_name
    """
    
    log = logger_obj if logger_obj is not None else logging.getLogger(__name__)
    
    parts = snp_name.split('_')
    
    # If the format is correct, contains parts separated by underscores
    if len(parts) >= 2:
        chromosome = parts[0]
        position = parts[1]
        return chromosome, position
    else:
        # Generate random chromosome and position
        chromosome = str(random.randint(1, 22))  # Chromosomes 1-22
        position = str(random.randint(100000, 999999))  # 6-digit position
        log.info(f"Generated random position for {snp_name}: chromosome {chromosome}, position {position}")
        return chromosome, position


# -------------------------
# 3.1 Initial filtration
# -------------------------

def initial_filter(P: pd.DataFrame, V: pd.DataFrame, A: pd.DataFrame, C: pd.DataFrame, I: pd.DataFrame, 
                  params: Optional[Dict] = None) -> Tuple[List[str], List[str], pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Implement Section 3.1: Initial filtration based on genotype posterior and mutant reads count.
    Modified to compute average mutant AF only across mutant cells (V > 0), matching R logic.
    
    Args:
        P: Posterior probability matrix
        V: Variant/Mutant allele frequency matrix
        A: Alternative allele count matrix
        C: Coverage matrix
        I: Binary indicator matrix
        params: Parameters dictionary
        
    Returns:
        Tuple of (kept_cells, kept_mutations, P_sub, V_sub, C_sub, I_sub)
    """
    
    p_thresh = params["posterior_threshold"]
    maf_max_cut = params["maf_max_threshold"]
    maf_mean_cut = params["maf_mean_threshold"]
    
    
    # Cell filtering: keep cells with at least one supported mutation (posterior > threshold)
    cells_with_V = (V > 0).any(axis=1)
    cells_with_P = (P > p_thresh).any(axis=1)
    kept_cells_mask = cells_with_V & cells_with_P
    kept_cells = P.index[kept_cells_mask].tolist()
    
    if len(kept_cells) == 0:
        raise ValueError("No cells passed initial filtration")
        
    
    # Subset matrices to kept cells
    P_sub = P.loc[kept_cells].copy()
    V_sub = V.loc[kept_cells].copy()
    A_sub = A.loc[kept_cells].copy()
    C_sub = C.loc[kept_cells].copy()
    I_sub = I.loc[kept_cells].copy()    
    
    # Mutation filtering: compute average only across mutant cells (V>0)
    Mj_bar = {}
    max_Mi = {}
    
    for j in V_sub.columns:
        # Mutant cells only
        mutant_cells_mask = V_sub[j] > 0
        mutant_cells = V_sub.loc[mutant_cells_mask, j]
        
        Mj_bar[j] = mutant_cells.mean() if len(mutant_cells) > 0 else 0.0
        max_Mi[j] = V_sub[j].max()
    
    # Convert to series for filtering
    Mj_bar_series = pd.Series(Mj_bar)
    max_Mi_series = pd.Series(max_Mi)
    
    # Filter mutations
    mutation_mask = (max_Mi_series >= maf_max_cut) & (Mj_bar_series >= maf_mean_cut)
    kept_mutations = V_sub.columns[mutation_mask].tolist()
    
    if len(kept_mutations) == 0:
        raise ValueError("No mutations passed initial filtration")
        
    
    # Subset matrices to kept mutations
    P_sub = P_sub[kept_mutations].copy()
    V_sub = V_sub[kept_mutations].copy()
    A_sub = A_sub[kept_mutations].copy()
    C_sub = C_sub[kept_mutations].copy()
    I_sub = I_sub[kept_mutations].copy()
    
    return kept_cells, kept_mutations, P_sub, V_sub, A_sub, C_sub, I_sub


# -------------------------
# 3.2 Coverage-based filtration
# -------------------------
def filter_scaffold_muts_by_na_proportion(filtered_sites, df_reads, df_celltype, na_prop_thresh=0.9, logger_obj=None):
    """
    Identify high-confidence scaffold mutations (shared variants in relatively ubiquitously expressed genes)
    based on cross-cell-type coverage.
    
    Parameters
    ----------
    filtered_sites : list
        Mutations that have passed per-cell filters (MAF, coverage).
    df_reads : pd.DataFrame
        Rows = cells (first row can be 'bulk'), columns = mutations,
        values = 'mut_count/total_count' (string) or NaN
    df_celltype : pd.DataFrame
        DataFrame containing 'barcode' and 'cell_type' columns
    na_prop_thresh : float or None
        Uniform threshold for NA proportion across cell types.
        If None, use Q3 quantile from the data.
    logger_obj : logging.Logger, optional
        Logger instance for logging output. If None, uses global logger.
    
    Returns
    -------
    scaffold_mutations : list
        High-confidence shared mutations for building scaffold phylogeny
    NA_prop : pd.DataFrame
        NA proportion (mutation × cell type)
    """
    # Use provided logger or fall back to global logger
    log = logger_obj if logger_obj is not None else logging.getLogger(__name__)
    
    # --- 1. Drop bulk row ---
    reads = df_reads.drop(index='bulk', errors='ignore')
    
    # --- 2. Get all cell types ---
    cell_types = df_celltype['cell_type'].unique()
    
    # --- 3. Build coverage matrix (coverage = 1, no coverage or NA = 0) ---
    def has_coverage(val):
        if pd.isna(val):
            return 0
        try:
            _, total = val.split('/')
            return 1 if int(total) > 0 else 0
        except:
            return 0
    
    coverage_matrix = reads.applymap(has_coverage)
    
    # --- 4. Calculate NA proportion per mutation per cell type ---
    df_NA_prop = pd.DataFrame(index=filtered_sites, columns=cell_types, dtype=float)
    for mut in filtered_sites:
        for t in cell_types:
            cells_in_type = df_celltype.loc[df_celltype['cell_type'] == t, 'barcode']
            valid_cells = [c for c in cells_in_type if c in coverage_matrix.index]
            if len(valid_cells) == 0:
                df_NA_prop.loc[mut, t] = 1.0
            else:
                cov_values = coverage_matrix.loc[valid_cells, mut]
                df_NA_prop.loc[mut, t] = 1.0 - cov_values.sum() / len(valid_cells)
    
    # --- 5. Calculate cutoff ---
    if na_prop_thresh is not None:
        theta = pd.Series(na_prop_thresh, index=cell_types)
    else:
        theta = df_NA_prop.quantile(0.75, axis=0)
    
    # --- 6. Determine informative mutations ---
    informative = df_NA_prop.lt(theta, axis=1)
    
    # --- 7. Filter high-confidence scaffold mutations ---
    scaffold_mutations = []
    cell_prop = df_celltype['cell_type'].value_counts(normalize=True)
    dominant_ctypes = cell_prop[cell_prop > 0.9].index.tolist()
    has_dominant = len(dominant_ctypes) > 0
    
    for mut in informative.index:
        n_informative = informative.loc[mut].sum()
        n_celltypes = informative.shape[1]
        if n_celltypes == 1:
            if n_informative >= 1:
                scaffold_mutations.append(mut)
        elif has_dominant:
            if informative.loc[mut, dominant_ctypes].any():
                scaffold_mutations.append(mut)
        else:
            if n_informative >= 2:
                scaffold_mutations.append(mut)
    
    log.info(f"Step1 (coverage-based) mutations: {len(scaffold_mutations)}")
    return scaffold_mutations, df_NA_prop


def get_total_reads_withoutNAcells(x):
    if pd.isna(x):
        return np.nan
    try:
        _, total = x.split('/')
        return int(total)
    except:
        return np.nan

def get_total_reads_withNAcells(x):
    if pd.isna(x):
        return 0  # Treat NA as 0 reads
    try:
        _, total = x.split('/')
        return int(total)
    except:
        return 0  # In case of any parsing issues, treat as 0


def coverage_filters(kept_mutations, df_reads, df_celltype, params, outputpath, logger_obj=None):
    """
    High-confidence scaffold mutation filtering based on coverage.
    
    Step1: cross-cell-type coverage (NA proportion)
    Step2: CV filter (based on ranking)
    
    Filtering logic:
    - total_muts <= 10: keep all
    - total_muts > 10: keep max(10, ceil(total_muts * cv_rank_thresh))
      i.e., if the calculated number is less than 10, keep at least 10;
      otherwise keep the calculated proportion.
    
    Parameters
    ----------
    kept_mutations : list
        Candidate mutations
    df_reads : pd.DataFrame
        Rows = cells (first row is 'bulk'), columns = mutations, values = 'mut/total' or NaN
    df_celltype : pd.DataFrame
        DataFrame containing 'barcode' and 'cell_type' columns
    params : dict
        Contains 'na_prop_thresh_global' and 'cv_rank_thresh' (default 0.5)
    outputpath : str, optional
        Path to save summary CSV and visualization
    logger_obj : logging.Logger, optional
        Logger instance for logging output. If None, uses global logger.
    
    Returns
    -------
    final_scaffold_mutations : list
        High-confidence shared scaffold mutations
    summary_df : pd.DataFrame
        Statistics per mutation: median / CV / mean / std / pass_filter
    df_NA_prop : pd.DataFrame
        NA proportion from Step1
    """
    import numpy as np
    import pandas as pd
    
    # Use provided logger or fall back to global logger
    log = logger_obj if logger_obj is not None else logging.getLogger(__name__)
    
    if params is None:
        params = DEFAULT_PARAMS
    
    na_prop_thresh = params["na_prop_thresh_global"]
    cv_rank_thresh = params.get("cv_rank_thresh", 0.5)
    
    # --- Step1: coverage-based filter ---
    step1_mutations, df_NA_prop = filter_scaffold_muts_by_na_proportion(
        kept_mutations, df_reads, df_celltype, na_prop_thresh, logger_obj=log
    )
    log.info(f"=====> Step1 (coverage-based) mutations: {len(step1_mutations)}")
    
    # --- Step2: CV filter based on ranking ---
    # Convert the reads to total read counts
    reads_matrix_withoutNAcells = df_reads.drop(index='bulk', errors='ignore').applymap(get_total_reads_withoutNAcells)
    reads_matrix_withNAcells = df_reads.drop(index='bulk', errors='ignore').applymap(get_total_reads_withNAcells)
    reads_matrix_withNAcells = reads_matrix_withNAcells.applymap(
        lambda v: 1 if (not pd.isna(v) and v >= 1) else (0 if not pd.isna(v) else np.nan)
    )
    
    # Calculate CV for all mutations
    cv_dict = {}
    median_dict = {}
    mean_dict = {}
    std_dict = {}
    
    for mut in kept_mutations:
        values_for_median = reads_matrix_withoutNAcells[mut].dropna() if mut in reads_matrix_withoutNAcells else []
        values_for_cv = reads_matrix_withNAcells[mut].dropna()
        
        # median
        if len(values_for_median) == 0:
            median_dict[mut] = np.nan
            continue
        median_val = np.median(values_for_median)
        median_dict[mut] = median_val
        
        # cv
        mean_val = np.mean(values_for_cv)
        std_val = np.std(values_for_cv)
        cv = std_val / mean_val if mean_val > 0 else np.inf
        cv_dict[mut] = cv
        mean_dict[mut] = mean_val
        std_dict[mut] = std_val
    
    # Sort by CV value (ascending)
    cv_sorted = sorted([(mut, cv) for mut, cv in cv_dict.items() if not np.isnan(cv) and cv != np.inf], key=lambda x: x[1])
    cv_sorted_muts = [mut for mut, _ in cv_sorted]
    
    total_muts = len(cv_sorted_muts)
    
    # --- Determine number of mutations to keep ---
    # Logic: min(max(10, ceil(total_muts * cv_rank_thresh)), total_muts)
    # Special case: total_muts <= 10, keep all
    if total_muts == 0:
        step2_mutations = []
        num_to_keep = 0
    elif total_muts <= 10:
        # Total <= 10, keep all
        num_to_keep = total_muts
        step2_mutations = cv_sorted_muts[:num_to_keep]
    else:
        # Total > 10
        # Calculate number to keep based on rank threshold
        num_from_rank = int(np.ceil(total_muts * cv_rank_thresh))
        # Keep at least 10, but cannot exceed total
        num_to_keep = min(max(10, num_from_rank), total_muts)
        step2_mutations = cv_sorted_muts[:num_to_keep]
    
    log.info(f"=====> Step2 (CV rank filter): kept {num_to_keep} out of {total_muts} mutations "
             f"(top {num_to_keep/total_muts*100:.1f}% by CV, threshold={cv_rank_thresh:.0%})")
    
    # --- Intersection ---
    final_scaffold_mutations = list(set(step1_mutations) & set(step2_mutations))
    final_scaffold_mutations_sorted = [i for i in kept_mutations if i in final_scaffold_mutations]
    log.info(f"=====> Final scaffold mutations (intersection): {len(final_scaffold_mutations_sorted)}")
    
    # --- Summary ---
    df_cv_stats = pd.DataFrame({
        "cov_median": pd.Series(median_dict),
        "cov_CV": pd.Series(cv_dict),
        "cov_mean": pd.Series(mean_dict),
        "cov_std": pd.Series(std_dict)
    })
    df_cv_stats["pass_CV"] = df_cv_stats.index.isin(step2_mutations)
    df_cv_stats["pass_NA"] = df_cv_stats.index.isin(step1_mutations)
    df_cv_stats["pass_cov"] = df_cv_stats.index.isin(final_scaffold_mutations_sorted)
    
    df_summary = pd.concat([df_cv_stats, df_NA_prop], axis=1)
    
    # --- Save ---
    if outputpath is not None:
        os.makedirs(outputpath, exist_ok=True)
        df_summary.to_csv(os.path.join(outputpath, "Summary_df_in_scaffold_filtration.csv"))
        
        # Generate CV distribution visualization
        try:
            import matplotlib.pyplot as plt
            
            fig, axes = plt.subplots(1, 2, figsize=(14, 5))
            
            # Left: CV distribution
            cv_vals = df_summary['cov_CV'].dropna()
            cv_vals = cv_vals[cv_vals != np.inf]
            
            if len(cv_vals) > 0:
                axes[0].hist(cv_vals, bins=30, alpha=0.7, edgecolor='black')
                # Mark threshold line
                if num_to_keep > 0 and num_to_keep < total_muts:
                    threshold_cv = cv_sorted[num_to_keep-1][1] if num_to_keep > 0 else None
                    if threshold_cv is not None:
                        axes[0].axvline(threshold_cv, color='r', linestyle='--', 
                                      label=f'Selected threshold (CV={threshold_cv:.2f})')
                axes[0].set_xlabel('Coefficient of Variation (CV)')
                axes[0].set_ylabel('Number of Mutations')
                axes[0].set_title('CV Distribution')
                axes[0].legend()
            
            # Right: Filtering status
            if 'pass_cov' in df_summary.columns:
                df_sorted = df_summary.sort_values('cov_CV').dropna(subset=['cov_CV'])
                df_sorted = df_sorted[df_sorted['cov_CV'] != np.inf]
                colors = ['green' if x else 'red' for x in df_sorted['pass_cov']]
                axes[1].bar(range(len(df_sorted)), [1]*len(df_sorted), color=colors, alpha=0.6)
                axes[1].set_xlabel('Mutations (sorted by CV)')
                axes[1].set_ylabel('Kept')
                axes[1].set_title(f'Kept Mutations (green: kept, red: filtered)\n'
                                f'Kept: {sum(df_sorted["pass_cov"])}/{len(df_sorted)}')
                axes[1].set_ylim(0, 1.5)
            
            plt.tight_layout()
            plt.savefig(os.path.join(outputpath, "cv_filter_visualization.pdf"), dpi=150)
            plt.close()
        except Exception as e:
            log.warning(f"Could not generate CV plot: {e}")
    
    return final_scaffold_mutations_sorted, df_summary


def calculate_cv_for_subgrouping(df_reads_resolved, mutations, cv_threshold=6, logger_obj=None):
    """
    Calculate CV values for mutations to filter low-CV mutations for subgrouping.
    
    Parameters:
    -----------
    df_reads_resolved : pd.DataFrame
        Reads dataframe containing read count information (rows = cells, columns = mutations)
    mutations : list
        List of mutations to calculate CV for
    cv_threshold : float
        CV threshold for filtering (mutations with CV <= threshold are kept)
    logger_obj : logging.Logger, optional
        Logger instance for logging messages. If None, uses global logger.
        
    Returns:
    --------
    Tuple containing:
        - low_cv_mutations: List of mutations with CV <= threshold
        - cv_stats: DataFrame with CV statistics for all mutations
    """
    # Use provided logger or fall back to global logger
    log = logger_obj if logger_obj is not None else logging.getLogger(__name__)
    
    # Convert the reads to total read counts, treating NA as 0
    reads_matrix_withoutNAcells = df_reads_resolved.drop(index='bulk', errors='ignore').applymap(get_total_reads_withoutNAcells)
    reads_matrix_withNAcells = df_reads_resolved.drop(index='bulk', errors='ignore').applymap(get_total_reads_withNAcells)
    reads_matrix_withNAcells = reads_matrix_withNAcells.applymap(
        lambda v: 1 if (not pd.isna(v) and v >= 1) else (0 if not pd.isna(v) else np.nan)
    )
    
    low_cv_mutations = []
    median_dict = {}
    cv_dict = {}
    mean_dict = {}
    std_dict = {}
    
    for mut in mutations:
        values_for_median = reads_matrix_withoutNAcells[mut].dropna() if mut in reads_matrix_withoutNAcells else []
        values_for_cv = reads_matrix_withNAcells[mut].dropna()
        
        # Calculate median
        if len(values_for_median) == 0:
            median_dict[mut] = np.nan
            continue
        median_val = np.median(values_for_median)
        median_dict[mut] = median_val
        
        # Calculate CV
        mean_val = np.mean(values_for_cv)
        std_val = np.std(values_for_cv)
        cv = std_val / mean_val if mean_val > 0 else np.inf
        cv_dict[mut] = cv
        mean_dict[mut] = mean_val
        std_dict[mut] = std_val
        
        if cv <= cv_threshold:
            low_cv_mutations.append(mut)
    
    log.info(f"In this subgroup, low_CV mutations (CV <= {cv_threshold}): {len(low_cv_mutations)} out of {len(mutations)}")
    
    # --- Summary ---
    df_cv_stats = pd.DataFrame({
        "cov_median": pd.Series(median_dict),
        "cov_CV": pd.Series(cv_dict),
        "cov_mean": pd.Series(mean_dict),
        "cov_std": pd.Series(std_dict)
    })
    df_cv_stats["pass_CV"] = df_cv_stats.index.isin(low_cv_mutations)
    
    return low_cv_mutations, df_cv_stats



# -------------------------
# 3.3 Consensus correlation graph and clustering
# -------------------------

import numpy as np
import pandas as pd
import itertools
import networkx as nx
from collections import defaultdict
from tqdm import tqdm
from joblib import Parallel, delayed
from typing import List, Tuple, Dict
from src.germline_filter import pairwise_counts, jaccard_index, are_mutations_correlated

def compute_clone_conter(muts, corr_cache, n_shuffle=100):
    clone_counter = defaultdict(int)
    for ref in muts:
        other_muts = [m for m in muts if m != ref]
        for _ in range(n_shuffle):
            shuffled = list(deterministic_permutation(other_muts))
            remaining = [ref] + shuffled.copy()
            while remaining:
                curr_ref = remaining[0]
                current_clone = [curr_ref]
                next_remaining = []
                for m in remaining[1:]:
                    if corr_cache[(curr_ref, m)]:
                        current_clone.append(m)
                    else:
                        next_remaining.append(m)
                clone_counter[tuple(sorted(current_clone))] += 1
                remaining = next_remaining
    return clone_counter

def compute_clone_and_pair_weights(muts, corr_cache, n_shuffle=100, logger_obj=None):
    """
    Compute clone weights and pair weights based on correlation patterns.
    
    This function builds clones by iteratively grouping correlated mutations,
    then computes weights for each clone and each mutation pair.
    
    Parameters
    ----------
    muts : list of str
        List of all mutation IDs
    corr_cache : dict
        {(mut1, mut2): True/False} - Boolean indicating whether two mutations are correlated
    n_shuffle : int
        Number of shuffle iterations per mutation for permutation generation
    logger_obj : logging.Logger, optional
        Logger instance for logging messages. If None, uses global logger.
    
    Returns
    -------
    clone_weights : dict
        {tuple(mut_ids): weight} - Global weight for each clone, including singleton clones
    pair_weights : dict
        {tuple(m1, m2): weight} - Weight for each mutation pair (only for clones with length >= 2)
    """
    # Use provided logger or fall back to global logger
    log = logger_obj if logger_obj is not None else logging.getLogger(__name__)
    
    clone_weights = defaultdict(float)  # Accumulate global clone weights
    
    if n_shuffle <= 0:
        log.warning("n_shuffle <= 0, returning empty weights")
        return clone_weights, defaultdict(float)
    
    log.debug(f"Computing clone and pair weights for {len(muts)} mutations with {n_shuffle} shuffles")
    
    for ref in muts:
        other_muts = [m for m in muts if m != ref]
        ref_clone_counter = defaultdict(int)  # Count each clone occurrence for current reference
        
        # Note: deterministic_permutation returns the same order for the same input,
        # so repeating n_shuffle times is equivalent to running once and normalizing proportionally.
        shuffled = list(deterministic_permutation(other_muts))
        remaining = [ref] + shuffled.copy()
        
        while remaining:
            curr_ref = remaining[0]
            current_clone = [curr_ref]
            next_remaining = []
            
            for m in remaining[1:]:
                key1 = (curr_ref, m)
                key2 = (m, curr_ref)
                is_corr = corr_cache.get(key1, corr_cache.get(key2, False))
                if is_corr:
                    current_clone.append(m)
                else:
                    next_remaining.append(m)
            
            ref_clone_counter[tuple(sorted(current_clone))] += 1
            remaining = next_remaining
        
        for clone, count in ref_clone_counter.items():
            clone_weights[clone] += count
    
    # Step: Compute pair weights, only considering clones with length >= 2
    pair_weights = defaultdict(float)
    for clone, weight in clone_weights.items():
        if len(clone) > 1:
            for m1, m2 in itertools.combinations(sorted(clone), 2):
                pair_weights[(m1, m2)] += weight
    
    log.debug(f"Computed {len(clone_weights)} clone weights and {len(pair_weights)} pair weights")
    
    return clone_weights, pair_weights


from typing import List, Tuple, Dict
from collections import defaultdict
import itertools
import numpy as np
import pandas as pd

def _build_pairwise_correlation_cache_fast(
    I_S: pd.DataFrame,
    muts: List[str],
) -> Tuple[Dict[Tuple[str, str], bool], Dict[str, float], Dict[str, int]]:
    I_numeric = I_S.loc[:, muts].apply(pd.to_numeric, errors="coerce")
    values = I_numeric.to_numpy(copy=False)
    
    ones = np.asarray(values == 1, dtype=np.int32, order="C")
    zeros = np.asarray(values == 0, dtype=np.int32, order="C")
    
    n11 = ones.T @ ones
    n10 = ones.T @ zeros
    n01 = zeros.T @ ones
    
    mutant_cell_number_arr = ones.sum(axis=0).astype(np.int32, copy=False)
    covered_cell_number_arr = (ones + zeros).sum(axis=0).astype(np.int32, copy=False)
    
    denominator = n11 + n10 + n01
    jaccard = np.divide(
        n11,
        denominator,
        out=np.zeros_like(n11, dtype=np.float64),
        where=denominator > 0,
    )
    f_forward = np.divide(
        n11,
        mutant_cell_number_arr[:, None],
        out=np.zeros_like(n11, dtype=np.float64),
        where=mutant_cell_number_arr[:, None] > 0,
    )
    f_reverse = np.divide(
        n11,
        mutant_cell_number_arr[None, :],
        out=np.zeros_like(n11, dtype=np.float64),
        where=mutant_cell_number_arr[None, :] > 0,
    )
    
    corr_matrix = (n11 >= 1) & (
        (jaccard >= 0.08) |
        ((jaccard > 0) & (jaccard < 0.08) & (np.maximum(f_forward, f_reverse) >= 0.5))
    )
    np.fill_diagonal(corr_matrix, True)
    
    mutant_cell_fraction_arr = np.divide(
        mutant_cell_number_arr,
        covered_cell_number_arr,
        out=np.zeros_like(mutant_cell_number_arr, dtype=np.float64),
        where=covered_cell_number_arr > 0,
    )
    
    corr_cache = {}
    for i, u in enumerate(muts):
        for j, v in enumerate(muts):
            corr_cache[(u, v)] = bool(corr_matrix[i, j])
    
    mutant_cell_fraction = dict(zip(muts, mutant_cell_fraction_arr))
    mutant_cell_number = dict(zip(muts, mutant_cell_number_arr))
    return corr_cache, mutant_cell_fraction, mutant_cell_number

def get_correlation_graph_elements(I_S: pd.DataFrame, n_shuffle: int = 100, seed: int = 42, 
                                   cutoff_mcf_for_graph: float = 0.05, cutoff_mcn_for_graph: int = 5,
                                   logger_obj=None) -> Tuple[Dict[Tuple[str], float], Dict[Tuple[str,str], float]]:
    """
    Compute clone_weights and pair_weights, filtering out low-support singleton mutations.
    
    Parameters
    ----------
    I_S : pd.DataFrame
        Binary matrix (cells x mutations), values 0/1/NA
    n_shuffle : int
        Number of random permutations per reference mutation
    seed : int
        Random seed
    cutoff_mcf_for_graph : float
        Minimum mutant cell fraction threshold for singleton mutations
    cutoff_mcn_for_graph : int
        Minimum mutant cell number threshold for singleton mutations
    logger_obj : logging.Logger, optional
        Logger instance for logging messages. If None, uses global logger.
    
    Returns
    -------
    clone_weights : dict
        Dictionary mapping clone tuples to weight
    pair_weights : dict
        Dictionary mapping mutation pairs to weight
    """
    # Use provided logger or fall back to global logger
    log = logger_obj if logger_obj is not None else logging.getLogger(__name__)
    
    muts = list(I_S.columns)
    n_mut = len(muts)
    
    log.debug(f"Computing correlation graph elements for {n_mut} mutations")
    
    # Step 1: Precompute pairwise correlation cache
    corr_cache, mutant_cell_fraction, mutant_cell_number = _build_pairwise_correlation_cache_fast(I_S, muts)
    
    # Step 2: Compute clone weights and pair weights
    clone_weights, pair_weights = compute_clone_and_pair_weights(muts, corr_cache, n_shuffle=n_shuffle, logger_obj=log)
    
    # Step 3: Calculate mutation fraction and mutant cell number for each mutation
    # Step 4: Remove low-support singleton mutations from clone_weights
    filtered_count = 0
    for mut in muts:
        if clone_weights.get((mut,), 0) == n_mut:  # Singleton clone
            frac = mutant_cell_fraction.get(mut, 0)
            num = mutant_cell_number.get(mut, 0)
            if frac <= cutoff_mcf_for_graph or num <= cutoff_mcn_for_graph:
                filtered_count += 1
                log.info(f"Filtered out singleton low-support mutation: {mut}, frac={frac:.3f}, num={num}")
                clone_weights.pop((mut,), None)  # Remove this clone
    
    log.info(f"Filtered out {filtered_count} singleton low-support mutations")
    
    return clone_weights, pair_weights

def plot_mutation_graph(G_ig, mutation_group, pdf_file, figsize=(8,8), edge_scale=0.2, seed=42,
                        ordered_mutations=None, logger_obj=None):
    """
    Visualize the mutation graph with node colors representing groups and edge widths representing weights.
    
    Parameters
    ----------
    G_ig : igraph Graph
        The constructed igraph graph
    mutation_group : dict
        {mutation_id: group_id} - Mapping of each mutation to its group
    pdf_file : str
        Output PDF file path
    figsize : tuple
        Figure size (width, height)
    edge_scale : float
        Edge weight scaling factor for visualization
    seed : int
        Random seed for layout reproducibility
    ordered_mutations : list, optional
        Predefined list of mutation order. If provided, nodes will be drawn in this order.
        If not provided, nodes are sorted alphabetically for deterministic ordering.
    logger_obj : logging.Logger, optional
        Logger instance for logging messages. If None, uses global logger.
    
    Returns
    -------
    None
    """
    # Use provided logger or fall back to global logger
    log = logger_obj if logger_obj is not None else logging.getLogger(__name__)
    
    # 1. Convert to NetworkX graph
    G_nx = nx.Graph()
    for v in G_ig.vs:
        G_nx.add_node(v['name'])
    for e in G_ig.es:
        m1 = G_ig.vs[e.source]['name']
        m2 = G_ig.vs[e.target]['name']
        G_nx.add_edge(m1, m2, weight=float(e['weight']))
    
    # Sort nodes according to ordered_mutations
    if ordered_mutations is not None:
        ordered_mutations_set = set(ordered_mutations)
        nodes_in_graph = set(G_nx.nodes())
        sorted_nodes = [mut for mut in ordered_mutations if mut in nodes_in_graph]
        extra_nodes = sorted([mut for mut in nodes_in_graph if mut not in ordered_mutations_set])
        sorted_nodes = sorted_nodes + extra_nodes
    else:
        sorted_nodes = sorted(G_nx.nodes())
    
    # Node colors
    groups = [mutation_group[n] for n in sorted_nodes]
    unique_groups = sorted(set(groups))
    
    # Use tab20 colormap (compatible with both old and new matplotlib)
    try:
        color_map = plt.colormaps.get_cmap('tab20')
    except AttributeError:
        color_map = plt.cm.get_cmap('tab20')
    
    node_colors = [color_map(g / max(1, len(unique_groups) - 1)) for g in groups]
    
    # Edge widths
    sorted_edges = sorted(G_nx.edges())
    edge_weights = [G_nx[u][v]['weight'] for u, v in sorted_edges]
    edge_widths = [w * edge_scale for w in edge_weights]
    
    # Layout
    pos = nx.spring_layout(G_nx, seed=seed, k=0.5)
    
    # Draw
    plt.figure(figsize=figsize)
    nx.draw_networkx_nodes(G_nx, pos, nodelist=sorted_nodes, node_color=node_colors, node_size=200)
    nx.draw_networkx_edges(G_nx, pos, edgelist=sorted_edges, width=edge_widths, alpha=0.7)
    nx.draw_networkx_labels(G_nx, pos, font_size=10, font_color='black')
    plt.title("Mutation Graph with Leiden Groups")
    plt.axis('off')
    plt.margins(x=0.2, y=0.2)
    plt.tight_layout(pad=2.0)
    plt.savefig(pdf_file, dpi=300)
    plt.close()
    
    log.info(f"Mutation graph plot saved to: {pdf_file}")


def compute_adaptive_resolution(I_matrix: pd.DataFrame, 
                                min_resolution: float = 0.5,
                                max_resolution: float = 2.0,
                                logger_obj=None) -> Tuple[float, Dict]:
    """
    Compute adaptive resolution based on data sparsity

    Core logic:
        - Higher NA ratio → higher resolution
        - More rare mutations → higher resolution
        - Better data quality → lower resolution
    """
    if logger_obj is None:
        logger_obj = logging.getLogger(__name__)
    
    n_cells, n_muts = I_matrix.shape
    
    # 1. NA ratio
    na_ratio = I_matrix.isna().sum().sum() / (n_cells * n_muts)
    
    # 2. Rare mutation ratio (mutations present in <= 3 cells)
    mut_counts = (I_matrix.fillna(0) == 1).sum(axis=0)
    rare_ratio = (mut_counts <= 3).sum() / n_muts if n_muts > 0 else 0
    
    # 3. Effective mutation ratio (mutations with coverage >= 5)
    coverage = I_matrix.notna().sum(axis=0)
    effective_ratio = (coverage >= 5).sum() / n_muts if n_muts > 0 else 0
    
    # 4. Data Quality Index (DQI)
    # Higher DQI = better quality → lower resolution
    # Lower DQI = poorer quality → higher resolution
    
    # NA score (higher score for fewer NAs)
    na_score = 1 - min(1.0, na_ratio * 1.5)
    
    # Rare mutation score (higher score for fewer rare mutations)
    rare_score = 1 - min(1.0, rare_ratio * 2.0)
    
    # Effective mutation score (higher score for more effective mutations)
    effective_score = effective_ratio
    
    # Combined DQI
    dqi = (na_score * 0.4 + rare_score * 0.35 + effective_score * 0.25)
    
    # 5. DQI → Resolution mapping
    if dqi >= 0.7:
        resolution = min_resolution + (max_resolution - min_resolution) * 0.2  # 0.8
    elif dqi >= 0.5:
        resolution = min_resolution + (max_resolution - min_resolution) * 0.5  # 1.25
    elif dqi >= 0.3:
        resolution = min_resolution + (max_resolution - min_resolution) * 0.75  # 1.625
    else:
        resolution = max_resolution  # 2.0
    
    # Additional fine-tuning based on rare mutation ratio
    if rare_ratio > 0.3:
        resolution += 0.2
    elif rare_ratio > 0.2:
        resolution += 0.1
    
    resolution = round(min(max_resolution, max(min_resolution, resolution)), 1)
    
    metrics = {
        'na_ratio': na_ratio,
        'rare_ratio': rare_ratio,
        'effective_ratio': effective_ratio,
        'na_score': na_score,
        'rare_score': rare_score,
        'effective_score': effective_score,
        'dqi': dqi,
        'resolution': resolution,
        'n_cells': n_cells,
        'n_muts': n_muts
    }
    
    return resolution, metrics

import igraph as ig
import leidenalg
import random
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
def leiden_mutation_groups(clone_weights, pair_weights, pdf_file, resolution=1, seed=42, 
                           ordered_mutations=None, logger_obj=None):
    """
    Build a weighted co-occurrence graph from clone_weights and pair_weights,
    and use the Leiden algorithm to partition mutations into groups.
    """
    # Use provided logger or fall back to global logger
    log = logger_obj if logger_obj is not None else logging.getLogger(__name__)
    
    log.info("Running Leiden algorithm for mutation group detection...")
    
    # ============================================================
    # CRITICAL FIX: Force deterministic RNG state for Leiden algorithm
    # ============================================================
    # The Leiden algorithm uses randomness for optimization. Even though we pass
    # a 'seed' parameter, the underlying igraph and numpy RNG states may have
    # been perturbed by parallel execution or previous operations. We explicitly
    # reset both random and numpy random states here to ensure deterministic
    # community detection results across parallel workers.
    # ============================================================
    import random
    random.seed(seed)
    np.random.seed(seed)
    # ============================================================
    
    # 1. Collect all mutations
    all_mutations = set()
    for clone in clone_weights.keys():
        all_mutations.update(clone)
    
    log.debug(f"Collected {len(all_mutations)} unique mutations from clone_weights")
    
    # ===== Sort according to ordered_mutations =====
    if ordered_mutations is not None:
        ordered_mutations_set = set(ordered_mutations)
        all_mutations = [mut for mut in ordered_mutations if mut in all_mutations]
        extra_mutations = sorted([mut for mut in all_mutations if mut not in ordered_mutations_set])
        all_mutations = all_mutations + extra_mutations
    else:
        # If no ordered_mutations provided, sort alphabetically for deterministic ordering
        all_mutations = sorted(list(all_mutations))
    
    log.debug(f"Final mutation order: {len(all_mutations)} mutations")
    
    # 2. Build igraph graph with deterministic node order
    G_ig = ig.Graph()
    G_ig.add_vertices(all_mutations)
    
    # 3. Add edges with deterministic order
    sorted_pairs = sorted(pair_weights.items())
    edge_count = 0
    for (m1, m2), w in sorted_pairs:
        G_ig.add_edge(m1, m2, weight=float(w))
        edge_count += 1
    
    log.debug(f"Added {edge_count} edges to the graph")
    
    # 4. Set random seeds (redundant but ensures safety)
    random.seed(seed)
    np.random.seed(seed)
    
    # 5. Run Leiden algorithm
    log.info(f"Running Leiden algorithm with resolution={resolution} and seed={seed}")
    partition = leidenalg.find_partition(
        G_ig,
        leidenalg.RBConfigurationVertexPartition,
        weights='weight',
        resolution_parameter=resolution,
        seed=seed
    )
    
    num_groups = len(partition)
    log.info(f"Leiden algorithm found {num_groups} mutation groups")
    
    # 6. Build mutation -> group dictionary
    mutation_group = {}
    for idx, community in enumerate(partition):
        for v in community:
            mutation_group[G_ig.vs[v]['name']] = idx
        log.debug(f"Group {idx}: {len(community)} mutations")
    
    # 7. Plot graph
    log.info(f"Saving mutation graph plot to: {pdf_file}")
    plot_mutation_graph(G_ig, mutation_group, pdf_file, seed=seed, 
                       ordered_mutations=all_mutations, logger_obj=log)
    
    return mutation_group, partition, G_ig


def build_cluster_graph(G_ig, mutation_group, logger_obj=None):
    """
    Build a weighted cluster-level graph from the mutation graph.
    
    Parameters
    ----------
    G_ig : igraph Graph
        The constructed igraph graph
    mutation_group : dict
        {mutation_id: group_id} - Mapping of each mutation to its group
    logger_obj : logging.Logger, optional
        Logger instance for logging messages. If None, uses global logger.
    
    Returns
    -------
    cluster_graph : igraph Graph
        Weighted graph where nodes are clusters and edges represent inter-cluster connections
    """
    # Use provided logger or fall back to global logger
    log = logger_obj if logger_obj is not None else logging.getLogger(__name__)
    
    clusters = set(mutation_group.values())
    cluster_graph = ig.Graph()
    cluster_graph.add_vertices(list(clusters))
    
    # Calculate inter-cluster connection weights
    inter_cluster_weights = {}
    for edge in G_ig.es:
        source_mut = G_ig.vs[edge.source]['name']
        target_mut = G_ig.vs[edge.target]['name']
        
        source_cluster = mutation_group[source_mut]
        target_cluster = mutation_group[target_mut]
        
        if source_cluster != target_cluster:
            pair = tuple(sorted([source_cluster, target_cluster]))
            inter_cluster_weights[pair] = inter_cluster_weights.get(pair, 0) + edge['weight']
    
    # Add edges
    for (cluster1, cluster2), weight in inter_cluster_weights.items():
        cluster_graph.add_edge(cluster1, cluster2, weight=weight)
    
    log.debug(f"Built cluster graph with {len(clusters)} clusters and {len(inter_cluster_weights)} inter-cluster connections")
    
    return cluster_graph

def detect_hub_clusters(G_ig, mutation_group, logger_obj=None):
    """
    Detect hub clusters based on weighted degree centrality.
    
    Parameters
    ----------
    G_ig : igraph Graph
        The constructed igraph graph
    mutation_group : dict
        {mutation_id: group_id} - Mapping of each mutation to its group
    logger_obj : logging.Logger, optional
        Logger instance for logging messages. If None, uses global logger.
    
    Returns
    -------
    hub_clusters : list
        List of cluster IDs identified as hub clusters
    cluster_degrees : dict
        {cluster_id: weighted_degree} - Weighted degree for each cluster
    """
    # Use provided logger or fall back to global logger
    log = logger_obj if logger_obj is not None else logging.getLogger(__name__)
    
    # 1. Build cluster-level graph
    cluster_graph = build_cluster_graph(G_ig, mutation_group)
    
    # 2. Calculate weighted degree for each cluster
    cluster_degrees = {}
    for cluster_id in set(mutation_group.values()):
        weighted_degree = 0
        for edge in cluster_graph.es:
            source_cluster = cluster_graph.vs[edge.source]['name']
            target_cluster = cluster_graph.vs[edge.target]['name']
            
            if source_cluster == cluster_id or target_cluster == cluster_id:
                weighted_degree += edge['weight']
        
        cluster_degrees[cluster_id] = weighted_degree
    
    # 3. Identify hub clusters (threshold can be adjusted)
    hub_threshold = np.percentile(list(cluster_degrees.values()), 75)
    hub_clusters = [cluster_id for cluster_id, degree in cluster_degrees.items() 
                   if degree > hub_threshold]
    
    log.debug(f"Detected {len(hub_clusters)} hub clusters out of {len(cluster_degrees)} total clusters")
    log.debug(f"Hub threshold: {hub_threshold:.4f}")
    
    return hub_clusters, cluster_degrees


def resolved_spots_by_immune_mutations(I_scaffold, immune_mutations, P_scaffold, V_scaffold, A_scaffold, C_scaffold, df_reads_scaffold, p_threshold=0.5, logger_obj=None):
    """
    Split spots in I_scaffold matrix based on immune mutations, and simultaneously update
    P_scaffold, C_scaffold, V_scaffold, and A_scaffold matrices.
    
    (pandas DataFrame, rows are spots, columns are mutations; 
    split I_scaffold and process other matrices accordingly)
    
    Parameters
    ----------
    I_scaffold : pd.DataFrame
        Binary indicator matrix (spots x mutations)
    immune_mutations : list
        List of immune mutation IDs
    P_scaffold : pd.DataFrame
        Posterior probability matrix
    V_scaffold : pd.DataFrame
        Variant/Mutant allele frequency matrix
    A_scaffold : pd.DataFrame
        Alternative allele count matrix
    C_scaffold : pd.DataFrame
        Coverage matrix
    df_reads_scaffold : pd.DataFrame
        Reads info (alt/total) matrix
    p_threshold : float
        Threshold for determining if a mutation is present, default 0.5
    logger_obj : logging.Logger, optional
        Logger instance for logging messages. If None, uses global logger.
    
    Returns
    -------
    I_resolved : pd.DataFrame
        Processed I_scaffold matrix (Binary indicator matrix)
    P_resolved : pd.DataFrame
        Processed P_scaffold matrix (Posterior probability matrix)
    V_resolved : pd.DataFrame
        Processed V_scaffold matrix (Variant/Mutant allele frequency matrix)
    A_resolved : pd.DataFrame
        Processed A_scaffold matrix (Alternative allele count matrix)
    C_resolved : pd.DataFrame
        Processed C_scaffold matrix (Coverage matrix)
    df_reads_resolved : pd.DataFrame
        Processed df_reads_scaffold matrix (Reads info (alt/total) matrix)
    spots_to_split : list
        List of spot IDs that were split
    """
    # Use provided logger or fall back to global logger
    log = logger_obj if logger_obj is not None else logging.getLogger(__name__)
    
    # Find immune mutations that actually exist in I_scaffold columns
    actual_immune_mutations = [mut for mut in immune_mutations if mut in I_scaffold.columns]
    # log.info(f"Found {len(actual_immune_mutations)} immune mutations in the data")
    
    # Find non-immune mutations (columns in I_scaffold but not in immune_mutations)
    non_immune_mutations = [col for col in I_scaffold.columns if col not in immune_mutations]
    # log.info(f"Found {len(non_immune_mutations)} non-immune mutations")
    
    # Lists to store results
    resolved_rows = []
    resolved_index = []
    
    P_resolved_rows = []
    C_resolved_rows = []
    V_resolved_rows = []
    A_resolved_rows = []
    df_reads_resolved_rows = []
    
    split_count = 0
    spots_to_split = []  # Track spot IDs that need to be split
    
    # Iterate through each row (spot)
    for spot_name, row in I_scaffold.iterrows():
        # Convert '<NA>' to NaN for processing
        I_row_clean = row.replace('<NA>', np.nan)
        
        # Apply p_threshold to each mutation site
        I_row_clean = I_row_clean.apply(lambda x: 1 if pd.notna(x) and x > p_threshold else np.nan)
        
        # Check if any immune mutations are 1
        immune_has_1 = False
        if actual_immune_mutations:
            immune_values = I_row_clean[actual_immune_mutations]
            immune_has_1 = any(immune_values == 1)
        
        # Check if any non-immune mutations are 1
        non_immune_has_1 = False
        if non_immune_mutations:
            non_immune_values = I_row_clean[non_immune_mutations]
            non_immune_has_1 = any(non_immune_values == 1)
        
        # If both immune and non-immune mutations are 1, need to split
        if immune_has_1 and non_immune_has_1:
            split_count += 1
            spots_to_split.append(spot_name)  # Record spot that needs splitting
            
            # Create immune version of the row
            immune_row = row.copy()
            # Set non-immune mutations to NaN
            immune_row[non_immune_mutations] = np.nan
            resolved_rows.append(immune_row)
            resolved_index.append(f"{spot_name}-immune")
            
            # Create non-immune version of the row
            non_immune_row = row.copy()
            # Set immune mutations to NaN
            non_immune_row[actual_immune_mutations] = np.nan
            resolved_rows.append(non_immune_row)
            resolved_index.append(f"{spot_name}-non")
            
            # Process P_scaffold
            immune_P_row = P_scaffold.loc[spot_name].copy()
            immune_P_row[non_immune_mutations] = np.nan
            P_resolved_rows.append(immune_P_row)
            non_immune_P_row = P_scaffold.loc[spot_name].copy()
            non_immune_P_row[actual_immune_mutations] = np.nan
            P_resolved_rows.append(non_immune_P_row)
            
            # Process V_scaffold
            immune_M_row = V_scaffold.loc[spot_name].copy()
            immune_M_row[non_immune_mutations] = np.nan
            V_resolved_rows.append(immune_M_row)
            non_immune_M_row = V_scaffold.loc[spot_name].copy()
            non_immune_M_row[actual_immune_mutations] = np.nan
            V_resolved_rows.append(non_immune_M_row)
            
            # Process A_scaffold
            immune_A_row = A_scaffold.loc[spot_name].copy()
            immune_A_row[non_immune_mutations] = np.nan
            A_resolved_rows.append(immune_A_row)
            non_immune_A_row = A_scaffold.loc[spot_name].copy()
            non_immune_A_row[actual_immune_mutations] = np.nan
            A_resolved_rows.append(non_immune_A_row)
            
            # Process C_scaffold
            immune_C_row = C_scaffold.loc[spot_name].copy()
            immune_C_row[non_immune_mutations] = np.nan
            C_resolved_rows.append(immune_C_row)
            non_immune_C_row = C_scaffold.loc[spot_name].copy()
            non_immune_C_row[actual_immune_mutations] = np.nan
            C_resolved_rows.append(non_immune_C_row)
            
            # Process df_reads_scaffold
            immune_df_reads_row = df_reads_scaffold.loc[spot_name].copy()
            immune_df_reads_row[non_immune_mutations] = np.nan
            df_reads_resolved_rows.append(immune_df_reads_row)
            non_immune_df_reads_row = df_reads_scaffold.loc[spot_name].copy()
            non_immune_df_reads_row[actual_immune_mutations] = np.nan
            df_reads_resolved_rows.append(non_immune_df_reads_row)
                
        else:
            # No split needed, keep original row
            resolved_rows.append(row)
            resolved_index.append(spot_name)
            
            # Process P_scaffold
            P_resolved_rows.append(P_scaffold.loc[spot_name])
            V_resolved_rows.append(V_scaffold.loc[spot_name])
            A_resolved_rows.append(A_scaffold.loc[spot_name])
            C_resolved_rows.append(C_scaffold.loc[spot_name])
            df_reads_resolved_rows.append(df_reads_scaffold.loc[spot_name])
    
    # Create new DataFrames
    I_resolved = pd.DataFrame(resolved_rows, index=resolved_index, columns=I_scaffold.columns)
    P_resolved = pd.DataFrame(P_resolved_rows, index=resolved_index, columns=P_scaffold.columns)
    V_resolved = pd.DataFrame(V_resolved_rows, index=resolved_index, columns=V_scaffold.columns)
    A_resolved = pd.DataFrame(A_resolved_rows, index=resolved_index, columns=A_scaffold.columns)
    C_resolved = pd.DataFrame(C_resolved_rows, index=resolved_index, columns=C_scaffold.columns)
    
    # Add bulk row back to df_reads
    df_reads_resolved = pd.DataFrame(df_reads_resolved_rows, index=resolved_index, columns=df_reads_scaffold.columns)
    bulk_row = df_reads_scaffold.loc['bulk']
    df_reads_resolved_with_bulk = pd.concat([pd.DataFrame(bulk_row).T, df_reads_resolved])
    df_reads_resolved_with_bulk.index = ['bulk'] + df_reads_resolved.index.tolist()
    
    # log.info(f"Processing complete: {len(I_scaffold)} original spots, {len(I_resolved)} spots after processing")
    # log.info(f"Total {split_count} spots were split")
    # log.debug(f"Split spot IDs: {spots_to_split}")
    
    return I_resolved, P_resolved, V_resolved, A_resolved, C_resolved, df_reads_resolved_with_bulk, spots_to_split


# I_resolved, P_resolved, V_resolved, A_resolved, C_resolved = resolved_spots_by_immune_mutations(I_scaffold, immune_mutations, P_scaffold, V_scaffold, A_scaffold, C_scaffold, p_threshold=0.5)

def split_spots_by_immune_mutations_scaffold(
    spots_to_split: list,  # List of spots that need to be split
    immune_mutations: list,  # List of immune mutations
    I_process: pd.DataFrame,  # I mutation matrix
    P_process: pd.DataFrame,   # P matrix (posterior probability)
    logger_obj=None  # Add logger_obj parameter
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Split the given spots_to_split based on immune mutations into two versions:
    one with '-immune' suffix and one with '-non' suffix.
    Rows that do not need splitting are retained unchanged.
    
    Parameters
    ----------
    spots_to_split : list
        List of spot IDs (rows in I_process) that need to be split
    
    immune_mutations : list
        List of immune mutation IDs
    
    I_process : pd.DataFrame
        Mutation matrix (rows are spots, columns are mutations)
    
    P_process : pd.DataFrame
        Posterior probability matrix
    
    logger_obj : logging.Logger, optional
        Logger instance for logging messages. If None, uses global logger.
    
    Returns
    -------
    I_process_resolved : pd.DataFrame
        Processed mutation matrix after splitting
    
    P_process_resolved : pd.DataFrame
        Processed posterior probability matrix after splitting
    """
    # Use provided logger or fall back to global logger
    log = logger_obj if logger_obj is not None else logging.getLogger(__name__)
    
    # Lists to store split rows
    resolved_I_rows = []
    resolved_P_rows = []
    resolved_index = []
    
    # Iterate through each spot
    for spot_name, row in I_process.iterrows():
        # If this row needs to be split
        if spot_name in spots_to_split:
            # Get original rows for this spot
            I_row = row
            P_row = P_process.loc[spot_name]
            
            # Create -immune row: keep immune mutation columns, others set to NA
            immune_row = I_row.copy()
            immune_row[immune_mutations] = I_row[immune_mutations]  # Keep immune mutations
            immune_row[~I_row.index.isin(immune_mutations)] = np.nan  # Set non-immune to NA
            
            P_immune_row = P_row.copy()
            P_immune_row[~P_row.index.isin(immune_mutations)] = np.nan  # Set non-immune to NA
            
            # Create -non row: set immune mutation columns to NA, others unchanged
            non_immune_row = I_row.copy()
            non_immune_row[immune_mutations] = np.nan  # Set immune to NA
            
            P_non_immune_row = P_row.copy()
            P_non_immune_row[immune_mutations] = np.nan  # Set immune to NA
            
            # Append split rows to results
            resolved_I_rows.append(immune_row)
            resolved_P_rows.append(P_immune_row)
            resolved_index.append(f"{spot_name}-immune")
            
            resolved_I_rows.append(non_immune_row)
            resolved_P_rows.append(P_non_immune_row)
            resolved_index.append(f"{spot_name}-non")
        else:
            # Rows that don't need splitting are retained as-is
            resolved_I_rows.append(row)
            resolved_P_rows.append(P_process.loc[spot_name])
            resolved_index.append(spot_name)
    
    # Combine split rows into new DataFrames
    I_process_resolved = pd.DataFrame(resolved_I_rows, index=resolved_index, columns=I_process.columns)
    P_process_resolved = pd.DataFrame(resolved_P_rows, index=resolved_index, columns=P_process.columns)
    
    log.debug(f"Split {len(spots_to_split)} spots, resulting in {len(resolved_index)} rows")
    
    return I_process_resolved, P_process_resolved


def sort_I_hierarchical_freeze_ones_fixed(I_S, mutation_group, logger_obj=None):
    """
    Sort the mutation matrix I_S hierarchically:
    - Column order: first by group ascending, then by MCF (Mutation Cell Fraction) descending within each group
    - Row order: iteratively "freeze" cells with value 1 column by column, pushing cells with 1 to the front,
      and continuing with remaining cells
    
    Parameters
    ----------
    I_S : pd.DataFrame
        Mutation matrix, rows are cells, columns are mutations, values are {0, 1, NA}
    mutation_group : dict
        {mutation_id: group_id} - Mapping of each mutation to its group
    logger_obj : logging.Logger, optional
        Logger object for logging messages
    
    Returns
    -------
    I_sorted : pd.DataFrame
        Sorted mutation matrix
    mut_df_sorted : pd.DataFrame
        Sorted mutation information table containing mutation, group, and mcf
    group_to_muts : dict
        Mapping of each group to its list of mutations
    final_order : list
        Final order of cells (rows)
    """
    # Step 1: Keep only mutations that exist in I_S and are in mutation_group
    selected_muts = [m for m in I_S.columns if m in mutation_group]
    if len(selected_muts) == 0:
        raise ValueError("No keys from mutation_group found in I_S.columns.")
    
    if logger_obj:
        logger_obj.info(f"Starting hierarchical sorting with {len(selected_muts)} mutations")
    
    # Step 2: Calculate group and MCF (n1 / (n1+n0), NA is excluded) for each mutation
    mut_info = []
    for mut in selected_muts:
        grp = mutation_group[mut]
        col = I_S[mut]
        n1 = (col == 1).fillna(False).sum()
        n0 = (col == 0).fillna(False).sum()
        mcf = n1 / (n1 + n0) if (n1 + n0) > 0 else -1.0
        mut_info.append({"mutation": mut, "group": grp, "mcf": mcf})
    mut_df = pd.DataFrame(mut_info)
    
    # Step 3: Sort mutations by group ascending, then mcf descending within each group
    mut_df_sorted = mut_df.sort_values(
        by=["group", "mcf"], ascending=[True, False]
    ).reset_index(drop=True)
    sorted_mutations = mut_df_sorted["mutation"].tolist()
    
    if logger_obj:
        logger_obj.info(f"Sorted {len(sorted_mutations)} mutations into {mut_df_sorted['group'].nunique()} groups")
    
    # Step 4: Row sorting (iteratively freeze ones)
    remaining = list(I_S.index)   # Initial remaining cells (keep original order)
    final_order = []
    
    for mut in sorted_mutations:
        ones, zeros, nas = [], [], []
        col = I_S[mut]
        for cell in remaining:
            v = col.loc[cell]
            if pd.isna(v):
                nas.append(cell)
            elif v == 1:
                ones.append(cell)
            else:
                zeros.append(cell)
        # Freeze ones (append them directly to the front)
        final_order.extend(ones)
        # Remaining cells go to the next iteration
        remaining = zeros + nas
    
    # Step 5: Append any remaining cells to the end
    final_order.extend(remaining)
    
    # Step 6: Generate sorted DataFrame
    I_sorted = I_S.loc[final_order, sorted_mutations].copy()
    
    # Step 7: Verify sorting correctness (including NA matching)
    orig_sub = I_S.loc[final_order, sorted_mutations]
    mask_eq = orig_sub.eq(I_sorted)
    mask_both_na = orig_sub.isna() & I_sorted.isna()
    mask_ok = mask_eq.fillna(False) | mask_both_na
    if not mask_ok.all().all():
        bad = (~mask_ok).sum().sum()
        raise RuntimeError(f"Sorting mismatch: {bad} cell-mutation values do not match (should not happen).")
    
    # Step 8: Build group to mutations mapping
    group_to_muts = mut_df_sorted.groupby("group")["mutation"].apply(list).to_dict()
    
    if logger_obj:
        logger_obj.info("Hierarchical sorting completed successfully")
    
    return I_sorted, mut_df_sorted, group_to_muts, final_order

# I_selected_and_sorted, mut_df_sorted, group_to_muts, final_order = sort_I_hierarchical_freeze_ones_fixed(I_S, mutation_group)


def map_group_to_backbone_mutations(mutation_group, group_to_muts, backbone_mutations, logger_obj=None):
    """
    Replace group IDs in group_to_muts with corresponding backbone mutation IDs based on mutation_group mapping.
    
    Parameters
    ----------
    mutation_group : dict
        Mapping of each mutation to its group ID, format: {mutation: group}
    group_to_muts : dict
        Mapping of each group ID to its list of mutations, format: {group: [mutations]}
    backbone_mutations : list
        List of backbone mutation identifiers
    logger_obj : logging.Logger, optional
        Logger instance for logging messages. If None, uses global logger.
    
    Returns
    -------
    dict
        Dictionary with backbone mutations as keys and their corresponding mutation lists as values,
        format: {backbone_mutation: [mutations]}
    """
    # Use provided logger or fall back to global logger
    log = logger_obj if logger_obj is not None else logging.getLogger(__name__)
    
    # Step 1: Create a mapping from mutation to its group ID
    mutation_to_group = {mutation: group for mutation, group in mutation_group.items()}
    
    # Step 2: Use this mapping to replace keys in group_to_muts with corresponding backbone mutations
    group_to_muts_with_backbone = {}
    
    # Iterate through each group ID and rename with backbone mutation value
    for group, mutations in group_to_muts.items():
        # Find the backbone mutation that corresponds to this group ID
        new_key = [mutation for mutation in backbone_mutations if mutation_to_group.get(mutation) == group]
        if new_key:  # Ensure the corresponding backbone mutation was found
            group_to_muts_with_backbone[new_key[0]] = mutations
    
    log.debug(f"Mapped {len(group_to_muts)} groups to {len(group_to_muts_with_backbone)} backbone mutations")
    
    return group_to_muts_with_backbone

# # 使用示例
# group_to_muts_with_backbone = map_group_to_backbone_mutations(mutation_group, group_to_muts, backbone_mutations)

def compute_corr_cache_with_new_mut_scaffold(I_attached, existing_muts, new_mut, logger_obj=None):
    """
    Compute correlation cache including the new mutation.
    
    Parameters
    ----------
    I_attached : pd.DataFrame
        Mutation matrix
    existing_muts : list
        List of existing mutation IDs
    new_mut : str
        New mutation ID to be added
    logger_obj : logging.Logger, optional
        Logger instance for logging messages. If None, uses global logger.
    
    Returns
    -------
    corr_cache : dict
        Correlation cache for all mutation pairs
    """
    # Use provided logger or fall back to global logger
    log = logger_obj if logger_obj is not None else logging.getLogger(__name__)
    
    # All mutations (existing + new)
    all_muts = existing_muts + [new_mut]
    
    corr_cache = {}
    
    # Calculate correlations for all mutation pairs (including new mutation with existing ones)
    for u, v in itertools.combinations(all_muts, 2):
        corr = are_mutations_correlated(I_attached, u, v)
        corr_cache[(u, v)] = corr
        corr_cache[(v, u)] = corr
    
    # Set self-correlation to True
    for m in all_muts:
        corr_cache[(m, m)] = True
    
    return corr_cache


def compute_new_mut_clone_affinity_correct_scaffold(new_mut, mutation_clones_rescue, I_attached, n_shuffle=100, logger_obj=None):
    """
    Correct version: uses correlation cache that includes the new mutation.
    
    Parameters
    ----------
    new_mut : str
        New mutation ID
    mutation_clones_rescue : dict
        Dictionary of mutation clones
    I_attached : pd.DataFrame
        Mutation matrix
    n_shuffle : int
        Number of shuffles for permutation
    logger_obj : logging.Logger, optional
        Logger instance for logging messages. If None, uses global logger.
    
    Returns
    -------
    clone_affinity : dict
        Affinity scores for each clone
    detailed_scores : dict
        Detailed scores for each clone
    """
    # Use provided logger or fall back to global logger
    log = logger_obj if logger_obj is not None else logging.getLogger(__name__)
    
    # Get all existing mutations
    all_existing_muts = []
    for clone in mutation_clones_rescue.values():
        all_existing_muts.extend(clone)
    all_existing_muts = list(set(all_existing_muts))
    
    # Recalculate correlation cache including the new mutation
    # log.debug("Recalculating correlation cache with new mutation...")
    corr_cache = compute_corr_cache_with_new_mut_scaffold(I_attached, all_existing_muts, new_mut)
    
    clone_affinity = {}
    detailed_scores = {}
    
    log.info(f"Calculating correlation with tree mutations for: {new_mut}")
    
    for clone_rep, clone_muts in mutation_clones_rescue.items():
        clone_key = tuple(sorted(clone_muts))
        
        # Check correlation between new mutation and clone members
        direct_corr_count = 0
        correlation_details = []
        
        for existing_mut in clone_muts:
            key1 = (new_mut, existing_mut)
            key2 = (existing_mut, new_mut)
            
            is_corr = corr_cache.get(key1, False) or corr_cache.get(key2, False)
            
            if is_corr:
                direct_corr_count += 1
                # Get detailed count data for debugging
                counts = pairwise_counts(I_attached, new_mut, existing_mut)
                correlation_details.append(f"{existing_mut}: N11={counts['N11']}, J={jaccard_index(I_attached, new_mut, existing_mut):.3f}")
            else:
                counts = pairwise_counts(I_attached, new_mut, existing_mut)
                correlation_details.append(f"{existing_mut}: N11={counts['N11']} (not corr)")
        
        direct_ratio = direct_corr_count / len(clone_muts) if clone_muts else 0
        
        # log.debug(f"\nClone {clone_rep}: {len(clone_muts)} mutations")
        # log.debug(f"  Direct correlation: {direct_corr_count}/{len(clone_muts)} = {direct_ratio:.3f}")
        # 
        # # Show first few correlation details
        # for detail in correlation_details[:3]:
        #     log.debug(f"    {detail}")
        # if len(correlation_details) > 3:
        #     log.debug(f"    ... and {len(correlation_details) - 3} more mutations")
        
        # Calculate correlation with clone representative mutation
        rep_correlation = 0
        if clone_rep in clone_muts:
            key1 = (new_mut, clone_rep)
            key2 = (clone_rep, new_mut)
            rep_correlation = 1 if (corr_cache.get(key1, False) or corr_cache.get(key2, False)) else 0
            counts = pairwise_counts(I_attached, new_mut, clone_rep)
            # log.debug(f"  Representative correlation: {rep_correlation} (N11={counts['N11']}, J={jaccard_index(I_attached, new_mut, clone_rep):.3f})")
        
        # Calculate weighted score (considering clone size)
        base_score = direct_ratio
        size_factor = min(1.0, len(clone_muts) / 10)  # Avoid overly large clones from dominating
        final_score = base_score * (0.7 + 0.3 * size_factor)
        
        clone_affinity[clone_key] = final_score
        detailed_scores[clone_key] = {
            'direct_ratio': direct_ratio,
            'direct_correlations': direct_corr_count,
            'rep_correlation': rep_correlation,
            'clone_size': len(clone_muts)
        }
        
        # log.debug(f"  Final score: {final_score:.3f}")
    
    return clone_affinity, detailed_scores


def select_best_clone_scaffold(detailed_scores, logger_obj=None):
    """
    Select the best clone scaffold.
    
    Parameters
    ----------
    detailed_scores : dict
        Detailed scores dictionary for clones, keys are tuples containing mutations
    logger_obj : logging.Logger, optional
        Logger instance for logging messages. If None, uses global logger.
    
    Returns
    -------
    list
        Best clone(s) based on scoring criteria
    """
    # Use provided logger or fall back to global logger
    log = logger_obj if logger_obj is not None else logging.getLogger(__name__)
    
    # Step 1: Find the maximum rep_correlation value
    max_rep_correlation = max([score['rep_correlation'] for score in detailed_scores.values()])
    
    if max_rep_correlation == 0:
        return []
    
    # Step 2: Find all clones with the maximum rep_correlation
    max_rep_correlation_clones = [
        clone for clone, score in detailed_scores.items() if score['rep_correlation'] == max_rep_correlation
    ]
    
    if len(max_rep_correlation_clones) == 1:
        return max_rep_correlation_clones
    
    # Step 4: Compare direct_ratio
    max_direct_ratio = max([detailed_scores[clone]['direct_ratio'] for clone in max_rep_correlation_clones])
    
    if max_direct_ratio == 0:
        return []
    
    # Step 5: Find clones with maximum direct_ratio
    best_clones = [
        clone for clone in max_rep_correlation_clones 
        if detailed_scores[clone]['direct_ratio'] == max_direct_ratio
    ]
    
    if len(best_clones) == 1:
        return best_clones
    
    # Step 6: If still multiple, select the one with most mutations
    max_mutation_count = max([len(clone) for clone in best_clones])
    best_clones_with_most_mutations = [
        clone for clone in best_clones 
        if len(clone) == max_mutation_count
    ]
    
    if len(best_clones_with_most_mutations) == 1:
        return best_clones_with_most_mutations
    
    # Step 7: If mutation count is also the same, sort by first mutation name and take the first
    return [sorted(best_clones_with_most_mutations, key=lambda x: x[0])[0]]


def add_new_mutation_to_clone(mutation_group_with_non_group_mutations, assigned_clone, new_mut, logger_obj=None):
    """
    Add a new mutation to the appropriate group based on the assigned_clone.
    
    Parameters
    ----------
    mutation_group_with_non_group_mutations : dict
        Dictionary containing mutation to group mapping
    assigned_clone : list of tuple
        Each element is a tuple containing mutations that belong to the same group
    new_mut : str
        New mutation to be added to the group in assigned_clone
    logger_obj : logging.Logger, optional
        Logger instance for logging messages. If None, uses global logger.
    
    Returns
    -------
    dict
        Updated mutation_group_with_non_group_mutations dictionary
    """
    # Use provided logger or fall back to global logger
    log = logger_obj if logger_obj is not None else logging.getLogger(__name__)
    
    # Step 1: Find the group ID of new_mut
    group_of_new_mut = None
    for mutation in assigned_clone[0]:  # Assuming each tuple contains one group of mutations
        if mutation in mutation_group_with_non_group_mutations:
            group_of_new_mut = mutation_group_with_non_group_mutations[mutation]
            break  # Found group ID, exit loop
    
    if group_of_new_mut is not None:
        # Step 2: Add new_mut to that group
        mutation_group_with_non_group_mutations[new_mut] = group_of_new_mut
    else:
        log.error(f"Could not find group for mutation {new_mut}")
    
    return mutation_group_with_non_group_mutations




def plot_heatmap_with_celltype_by_your_sorting(I_raw, df_celltype, mutation_group, your_sorting_muts, pdf_file, logger_obj=None):
    """
    Plot a heatmap of the mutation matrix with cell type annotation bar,
    row/column statistics bar plots, and legends.
    
    Parameters
    ----------
    I_raw : pd.DataFrame
        Raw mutation matrix (cell x mutation), elements are {0, 1, NA}
    df_celltype : pd.DataFrame
        Cell type information for each cell (must contain "cell_type" column)
    mutation_group : dict
        {mutation_id: group_id} - Mapping of each mutation to its group
    your_sorting_muts : list
        Predefined mutation order (column order)
    pdf_file : str
        Output PDF file path for saving the figure
    logger_obj : logging.Logger, optional
        Logger instance for logging messages. If None, uses global logger.
    """
    # Use provided logger or fall back to global logger
    log = logger_obj if logger_obj is not None else logging.getLogger(__name__)
    
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap, BoundaryNorm
    from matplotlib.patches import Patch
    
    # -------------------
    # Step 1: Determine cell order (keep original order, filter needed cells)
    # -------------------
    ordered_leaf_nodes = list(I_raw.index)
    cells_to_use = [c for c in ordered_leaf_nodes if c in I_raw.index]
    I_sorted_cells = I_raw.loc[cells_to_use]
    
    # -------------------
    # Step 2: Determine mutation order (use externally provided ordering)
    # -------------------
    mutation_order = your_sorting_muts
    I_sorted = I_sorted_cells[mutation_order]
    
    # -------------------
    # Step 3: Convert matrix values (NA -> 2, for separate coloring)
    # -------------------
    I_numeric = I_sorted.fillna(np.nan).apply(pd.to_numeric, errors="coerce")
    I_plot = I_numeric.copy()
    I_plot = I_plot.where(~I_plot.isna(), 2)  # Fill NA with 2
    
    # -------------------
    # Step 4: Calculate row/column mutation counts (for bar plots)
    # -------------------
    row_sums = I_numeric.sum(axis=1, skipna=True)  # Mutations per cell
    col_sums = I_numeric.sum(axis=0, skipna=True)  # Cells supporting each mutation
    
    # -------------------
    # Step 5: Layout using GridSpec
    #   - Left: row bar plot
    #   - Center: heatmap
    #   - Top: column bar plot
    #   - Right: cell type annotation bar
    #   - Bottom: legend
    # -------------------
    fig = plt.figure(figsize=(12, 10))
    gs = fig.add_gridspec(5, 6,
                          width_ratios=[0.3, 0.3, 3, 0.05, 0.05, 0.05],
                          height_ratios=[0.5, 0.05, 3, 0.3, 0.3],
                          wspace=0.05, hspace=0.05)
    
    ax_row_bar = fig.add_subplot(gs[2, 0])   # Left row bar plot
    ax_heatmap = fig.add_subplot(gs[2, 2])   # Center heatmap
    ax_celltype_bar = fig.add_subplot(gs[2, 3])  # Right cell type annotation bar
    ax_col_bar = fig.add_subplot(gs[0, 2])   # Top column bar plot
    ax_dummy = fig.add_subplot(gs[0, 0])
    ax_dummy.axis("off")  # Placeholder
    
    # -------------------
    # Step 6: Draw heatmap
    #   - Color mapping: 0=light blue, 1=dark red, NA=white
    # -------------------
    cmap = ListedColormap(["#D4E8F0", "#7D2224", "white"])
    bounds = [0, 0.5, 1.5, 2.5]
    norm = BoundaryNorm(bounds, cmap.N)
    
    im = ax_heatmap.imshow(I_plot, aspect="auto", cmap=cmap,
                           interpolation="nearest", norm=norm)
    ax_heatmap.set_xlim(-0.5, I_plot.shape[1]-0.5)
    ax_heatmap.set_ylim(I_plot.shape[0]-0.5, -0.5)
    ax_heatmap.set_yticks([])
    
    # Set x-axis mutation labels
    ax_heatmap.set_xticks(range(len(I_plot.columns)))
    ax_heatmap.set_xticklabels(I_plot.columns, rotation=90, fontsize=6, ha='center')
    
    # -------------------
    # Step 7: Color x-axis mutation labels by group
    # -------------------
    unique_groups = sorted(set(mutation_group.values()))
    if len(unique_groups) > 0:
        cmap_groups = plt.cm.get_cmap("tab20", len(unique_groups))
        group_colors = {g: cmap_groups(i) for i, g in enumerate(unique_groups)}
        
        for label in ax_heatmap.get_xticklabels():
            mut_name = label.get_text()
            if mut_name in mutation_group:
                group_id = mutation_group[mut_name]
                label.set_color(group_colors.get(group_id, 'black'))
            else:
                label.set_color('black')
    
    # -------------------
    # Step 8: Column bar plot (cells supporting each mutation)
    # -------------------
    ax_col_bar.bar(range(len(col_sums)), col_sums.values,
                   color="#7D2224", alpha=0.7, align="center")
    ax_col_bar.set_xlim(ax_heatmap.get_xlim())
    ax_col_bar.set_xticks([])
    ax_col_bar.tick_params(axis="y", labelsize=8)
    ax_col_bar.set_ylabel("Cell Number\nper Mutation", fontsize=10)
    
    # -------------------
    # Step 9: Row bar plot (mutations per cell)
    # -------------------
    ax_row_bar.barh(range(len(row_sums)), row_sums.values,
                    color="#7D2224", alpha=0.7, align="center")
    ax_row_bar.set_ylim(ax_heatmap.get_ylim())
    ax_row_bar.set_yticks([])
    ax_row_bar.set_xlabel("Mutation\nBurden\nper Cell", fontsize=10)
    ax_row_bar.invert_xaxis()
    
    # -------------------
    # Step 10: Draw cell type annotation bar (right side)
    # -------------------
    celltypes = df_celltype["cell_type"].astype("category")
    type_codes = celltypes.cat.codes
    unique_types = celltypes.cat.categories
    
    if len(unique_types) > 0:
        cmap_types = plt.cm.get_cmap("tab20", len(unique_types))
        ax_celltype_bar.imshow(np.array(type_codes)[:, None], aspect="auto",
                               cmap=cmap_types, origin="upper")
    else:
        # If no cell types, fill with gray
        log.warning("No unique cell types found, using gray fill")
        ax_celltype_bar.imshow(np.zeros((len(type_codes), 1)), aspect="auto",
                               cmap="gray", origin="upper")
    ax_celltype_bar.set_xticks([])
    ax_celltype_bar.set_yticks([])
    
    # -------------------
    # Step 11: Add legends
    #   - Mutation values (0, 1, NA)
    #   - Cell types
    #   - Mutation groups
    # -------------------
    # Mutation values legend
    heatmap_handles = [Patch(facecolor=c, label=l) 
                       for c, l in zip(["#D4E8F0", "#7D2224", "white"],
                                       ["0 (No Mutation)", "1 (Mutation)", "NA (Missing)"])]
    fig.legend(handles=heatmap_handles, loc="lower center", ncol=3,
               bbox_to_anchor=(0.5, -0.03), frameon=False, fontsize=9,
               title="Mutation Values", title_fontsize=10)
    
    # Cell type legend (only if data exists)
    if len(unique_types) > 0:
        celltype_handles = [Patch(facecolor=cmap_types(i), label=label) 
                            for i, label in enumerate(unique_types)]
        fig.legend(handles=celltype_handles, loc="lower center",
                   ncol=min(len(unique_types), 5),
                   bbox_to_anchor=(0.5, -0.12), frameon=False, fontsize=9,
                   title="Cell Types", title_fontsize=10)
    else:
        log.warning("No cell types found, skipping cell type legend")
    
    # Mutation group legend (only if groups exist)
    if len(unique_groups) > 0:
        group_handles = [Patch(facecolor=color, label=f'Group {group_id}')
                         for group_id, color in group_colors.items()]
        fig.legend(handles=group_handles, loc="lower center",
                   ncol=min(len(group_handles), 5),
                   bbox_to_anchor=(0.5, -0.20), frameon=False, fontsize=9,
                   title="Mutation Groups", title_fontsize=10)
    else:
        log.warning("No mutation groups found, skipping group legend")
    
    # -------------------
    # Step 12: Save figure
    # -------------------
    plt.suptitle("Heatmap of Mutations with Cell Type Bar", fontsize=14, y=0.95)
    plt.tight_layout()
    # Dynamically adjust bottom space to avoid excessive blank space from empty legends
    plt.subplots_adjust(bottom=0.1 if len(unique_types) > 0 else 0.05)
    plt.savefig(pdf_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    log.info(f"Heatmap saved to: {pdf_file}")


def select_founder_mutations(I_selected_and_sorted: pd.DataFrame, mutation_group: dict, logger_obj=None):
    """
    Select the founder mutation for each group from I_selected_and_sorted based on mutation_group.
    The founder mutation is defined as the mutation with the highest Mutant Cell Fraction (MCF) within each group.
    
    Parameters
    ----------
    I_selected_and_sorted : pd.DataFrame
        Mutation matrix (cells x mutations), elements are {1, 0, NA}
    mutation_group : dict
        {mutation: group_id} - Mapping of each mutation to its group
    logger_obj : logging.Logger, optional
        Logger object for logging messages
    
    Returns
    -------
    dict
        {group_id: founder_mutation}
    """
    if logger_obj:
        logger_obj.info("Starting founder mutation selection...")
    
    # Result dictionary
    founder_dict = {}
    
    # Convert to Series for easier manipulation
    mut2group = pd.Series(mutation_group)
    
    if logger_obj:
        logger_obj.info(f"Processing {len(mut2group.groupby(mut2group))} mutation groups")
    
    for group_id, muts in mut2group.groupby(mut2group):
        muts = list(muts.index)  # Mutations in current group
        mcf_values = {}
        mcn_values = {}
        
        for mut in muts:
            col = I_selected_and_sorted[mut]
            num_mutant = (col == 1).sum()     # Number of mutant cells
            num_covered = col.notna().sum()   # Number of cells with sequencing coverage
            mcf = num_mutant / num_covered if num_covered > 0 else 0
            mcf_values[mut] = mcf
            mcn_values[mut] = num_mutant
        
        # Find the mutation with the highest MCF
        founder_mut = max(mcn_values, key=mcn_values.get)
        founder_dict[group_id] = founder_mut
        
        if logger_obj:
            logger_obj.debug(f"Group {group_id}: selected founder '{founder_mut}' with MCF={mcf_values[founder_mut]:.4f} (mutant cells={mcn_values[founder_mut]})")
    
    if logger_obj:
        logger_obj.info(f"Successfully selected {len(founder_dict)} founder mutations")
    
    return founder_dict


def calculate_intersection_counts_under_backbone_nodes_scaffold(mutation_list_under_backbone_nodes, M_current, I_attached, new_mut):
    """
    Calculate co-occurrence counts between new_mut and mutations in each backbone node's mutation list.
    
    Parameters
    ----------
    mutation_list_under_backbone_nodes : dict
        Mapping of backbone node to its mutation list
    M_current : pd.DataFrame
        Cell states on backbone nodes (0/1)
    I_attached : pd.DataFrame
        Cell states on other mutations (0/1/NaN)
    new_mut : str
        New mutation to check
    
    Returns
    -------
    dict
        Co-occurrence count for each backbone node
    """
    # Check if new_mut exists in I_attached
    if new_mut not in I_attached.columns:
        raise ValueError(f"Mutation {new_mut} is not in I_attached.columns")
    
    intersection_counts = {}
    
    for backbone_node, mutation_list in mutation_list_under_backbone_nodes.items():
        # 1. Find cells where current backbone node is 1
        backbone_cells = M_current[M_current[backbone_node] == 1].index
        
        # 2. Filter to cells that also exist in I_attached
        common_cells = backbone_cells.intersection(I_attached.index)
        
        if len(common_cells) == 0:
            intersection_counts[backbone_node] = 0
            continue
        
        # 3. Get the status of new_mut in these cells
        new_mut_status = I_attached.loc[common_cells, new_mut]
        
        # 4. Filter to cells where new_mut is 1
        new_mut_positive_cells = new_mut_status[new_mut_status == 1].index
        
        if len(new_mut_positive_cells) == 0:
            intersection_counts[backbone_node] = 0
            continue
        
        # 5. For each mutation in mutation_list, calculate co-occurrence with new_mut
        total_intersection = 0
        
        for mutation in mutation_list:
            if mutation in I_attached.columns:
                # Get current mutation status in new_mut-positive cells
                mut_status = I_attached.loc[new_mut_positive_cells, mutation]
                # Count cells where both are 1 (ignoring NaN values)
                intersection_count = (mut_status == 1).sum()
                total_intersection += intersection_count
        
        intersection_counts[backbone_node] = total_intersection
    
    return intersection_counts


def find_best_backbone_node_scaffold(mutation_list_under_backbone_nodes, M_current, intersection_counts_under_backbone_node):
    """
    Simplified version: returns only the best backbone node.
    
    Parameters
    ----------
    mutation_list_under_backbone_nodes : dict
        Mapping of backbone node to its mutation list
    M_current : pd.DataFrame
        Cell states on backbone nodes (0/1)
    intersection_counts_under_backbone_node : dict
        Co-occurrence counts for each backbone node
    
    Returns
    -------
    str
        Best backbone node
    """
    max_count = max(intersection_counts_under_backbone_node.values())
    max_nodes = [node for node, count in intersection_counts_under_backbone_node.items() 
                if count == max_count]
    
    if len(max_nodes) == 1:
        return max_nodes[0]
    
    # Tie-breaking: choose node with highest normalized count
    best_node = None
    best_normalized = -1
    
    for node in max_nodes:
        backbone_cells = M_current[M_current[node] == 1].index
        backbone_cell_count = len(backbone_cells)
        
        if backbone_cell_count > 0:
            normalized = intersection_counts_under_backbone_node[node] / backbone_cell_count
            if normalized > best_normalized:
                best_normalized = normalized
                best_node = node
        else:
            # If no cells, normalized value is 0
            if best_normalized < 0:  # No valid node found yet
                best_node = node
                best_normalized = 0
    
    return best_node


def find_best_backbone_for_new_mutation_scaffold(mutation_list_under_backbone_nodes, M_current, I_attached, new_mut, logger_obj=None):
    """
    Complete function: calculate co-occurrence counts and find the best backbone node.
    
    Parameters
    ----------
    mutation_list_under_backbone_nodes : dict
        Mapping of backbone node to its mutation list
    M_current : pd.DataFrame
        Cell states on backbone nodes (0/1)
    I_attached : pd.DataFrame
        Cell states on other mutations (0/1/NaN)
    new_mut : str
        New mutation to check
    logger_obj : logging.Logger, optional
        Logger object for logging messages
    
    Returns
    -------
    tuple
        (best_backbone, intersection_counts)
    """
    if logger_obj:
        logger_obj.debug(f"Finding best backbone for mutation: {new_mut}")
    
    # Step 1: Calculate co-occurrence counts
    intersection_counts = calculate_intersection_counts_under_backbone_nodes_scaffold(
        mutation_list_under_backbone_nodes, M_current, I_attached, new_mut
    )
    
    if logger_obj:
        logger_obj.debug(f"Intersection counts calculated for {len(intersection_counts)} backbone nodes")
    
    # Step 2: Find the best backbone node
    best_backbone = find_best_backbone_node_scaffold(
        mutation_list_under_backbone_nodes, M_current, intersection_counts
    )
    
    if logger_obj:
        logger_obj.debug(f"Best backbone node selected: {best_backbone}")
    
    return best_backbone, intersection_counts




# -------------------------
# 3.4 Penalty-based placement
# -------------------------

import copy
from typing import List, Dict, Tuple, Optional, Set
import numpy as np
import pandas as pd
import math
from collections import deque
from src.germline_filter import pairwise_counts

# -------------------------
# TreeNode: Minimal tree object
# -------------------------
class TreeNode:
    def __init__(self, name: str):
        self.name = name            # Usually mutation id or "ROOT"
        self.parent: Optional['TreeNode'] = None
        self.children: List['TreeNode'] = []
    
    def add_child(self, child: 'TreeNode'):
        child.parent = self
        self.children.append(child)
    
    def remove_child(self, child: 'TreeNode'):
        # Check if child exists
        if child in self.children:
            self.children.remove(child)
            child.parent = None
    
    def is_root(self):
        return self.parent is None
    
    def traverse(self):
        # Yield nodes in pre-order
        yield self
        for c in self.children:
            yield from c.traverse()
    
    def find(self, name: str) -> Optional['TreeNode']:
        for n in self.traverse():
            if n.name == name:
                return n
        return None
    
    def path_to_root(self) -> List['TreeNode']:
        node = self
        p = []
        while node is not None:
            p.append(node)
            node = node.parent
        return p[::-1]  # root->...->self
    
    def insert_on_edge(self, parent: 'TreeNode', child: 'TreeNode', new_node_name: str) -> 'TreeNode':
        """
        Insert a new node on the edge between parent and child.
        Assumes child is a child of parent.
        """
        if child not in parent.children:
            raise ValueError("child is not a child of parent")
        new_node = TreeNode(new_node_name)
        # Replace child with new_node under parent
        parent.remove_child(child)
        parent.add_child(new_node)
        # Attach child under new_node
        new_node.add_child(child)
        return new_node
    
    def add_new_parent_for_children(self, children_list: List['TreeNode'], new_node_name: str) -> 'TreeNode':
        """
        Move children_list from their original parent into a new parent node.
        The new parent is attached to the original parent of the first child.
        """
        parents = set(c.parent for c in children_list)
        if len(parents) != 1:
            # If children have different parents, define strategy here
            pass
        old_parent = children_list[0].parent
        new_node = TreeNode(new_node_name)
        if old_parent is None:
            # Attach to root (if allowed)
            pass
        else:
            # Replace these children in old_parent with new_node
            for c in children_list:
                old_parent.remove_child(c)
            old_parent.add_child(new_node)
        # Attach children to new_node
        for c in children_list:
            new_node.add_child(c)
        return new_node
    
    def add_leaf(self, parent: 'TreeNode', new_leaf_name: str) -> 'TreeNode':
        new_leaf = TreeNode(new_leaf_name)
        parent.add_child(new_leaf)
        return new_leaf
    
    def copy(self) -> 'TreeNode':
        """Deep copy the entire tree (preserving names)"""
        mapping = {}
        def _copy(node):
            n2 = TreeNode(node.name)
            mapping[node] = n2
            for c in node.children:
                n2.add_child(_copy(c))
            return n2
        return _copy(self)
    
    def all_nodes(self) -> List['TreeNode']:
        """Return all node objects in the tree"""
        return list(self.traverse())
    
    def all_names(self) -> List[str]:
        """Return all node names in the tree"""
        return [n.name for n in self.traverse()]
    
    def all_names_no_root(self) -> List[str]:
        """Return all node names excluding 'ROOT'"""
        return [n.name for n in self.traverse() if n.name != "ROOT"]
    
    def all_edges(self) -> List[tuple]:
        """Return all edges in the tree as (parent_name, child_name)"""
        return [(n.parent.name, n.name) for n in self.traverse() if n.parent]
    
    def to_string(self, level=0):
        """Convert tree structure to string"""
        indent = "  " * level
        result = f"{indent}└─ {self.name}\n"
        for child in self.children:
            result += child.to_string(level + 1)
        return result
    
    def save_to_file(self, filename):
        """Save tree structure to file"""
        with open(filename, 'w', encoding='utf-8') as f:
            f.write(self.to_string())
    
    def __str__(self):
        return self.to_string()


def print_tree_simple(node, level=0):
    """Simple version, using only print and not going through logging"""
    indent = "  " * level
    print(f"{indent}└─ {node.name}")
    for child in node.children:
        print_tree_simple(child, level + 1)


def print_tree(node, level=0, logger_obj=None):
    """
    Recursively print tree structure using logger.
    
    Parameters
    ----------
    node : TreeNode
        The root node of the tree to print
    level : int
        Current indentation level
    logger_obj : logging.Logger, optional
        Logger instance for logging messages. If None, uses global logger.
    """
    # Use provided logger or fall back to global logger
    log = logger_obj if logger_obj is not None else logging.getLogger(__name__)
    
    indent = "  " * level
    log.info(f"{indent}└─ {node.name}")
    for child in node.children:
        print_tree(child, level + 1, logger_obj)


def print_tree_logger(node, indent="", is_root=True, logger_obj=None):
    """
    Print tree structure with tree connecting lines.
    
    Parameters
    ----------
    node : TreeNode
        The node to print
    indent : str
        Indentation string
    is_root : bool
        Whether this is the root node
    logger_obj : logging.Logger, optional
        Logger to use. If None, uses global logger.
    """
    # Use provided logger or fall back to global logger
    log = logger_obj if logger_obj is not None else logging.getLogger(__name__)
    
    if is_root:
        log.info(node.name)
        for child in node.children:
            print_tree_logger(child, "  ", False, logger_obj)
    else:
        log.info(f"{indent}└─ {node.name}")
        new_indent = indent + "    "
        for child in node.children:
            print_tree_logger(child, new_indent, False, logger_obj)


##### Tree TreeNode Format Saving and Reading
import json

# Assume that T_full is your TreeNode object
def tree_to_dict(node):
    # Suppose you already have a TreeNode class or a similar structure
    tree_dict = {}
    tree_dict['name'] = node.name
    if node.children:
        tree_dict['children'] = [tree_to_dict(child) for child in node.children]
    return tree_dict


def print_tree_dict(tree, prefix=""):
    print(prefix + "└─ " + tree["name"])
    if "children" in tree:
        for child in tree["children"]:
            print_tree_dict(child, prefix + "  ")

# tree_json_file = "/storage/douyanmeiLab/yangqing/tools/PhyloMosaicGenie/pmg/src/phylosolid/samples/phylo_wrapped_up/sample_P6_merged/mutation_integrator.germline_no/phylo_v11.okVersion_tripleFPs_but_backbone.send_to_YD/final_cleaned_tree_node.json"

# with open(tree_json_file, 'r') as f:
#     tree_json = json.load(f)

# print_tree_dict(tree_json)




def find_paths_to_leaves_tree_node(start_node, logger_obj=None):
    """
    Find all paths from the given node to all leaf nodes (for TreeNode objects).
    
    Parameters
    ----------
    start_node : TreeNode
        The starting node
    logger_obj : logging.Logger, optional
        Logger instance for logging messages. If None, uses global logger.
    
    Returns
    -------
    list
        List containing all paths from the starting node to leaf nodes
    """
    # Use provided logger or fall back to global logger
    log = logger_obj if logger_obj is not None else logging.getLogger(__name__)
    
    def dfs(current_node, current_path, all_paths):
        """Depth-first search to find all paths to leaf nodes"""
        current_path = current_path + [current_node.name]
        
        if not current_node.children:  # If leaf node
            all_paths.append(current_path)
        else:
            for child in current_node.children:
                dfs(child, current_path, all_paths)
    
    all_paths = []
    dfs(start_node, [], all_paths)
    
    log.debug(f"Found {len(all_paths)} paths from '{start_node.name}' to leaf nodes")
    
    return all_paths


def get_all_descendants_tree_node(start_node, logger_obj=None):
    """
    Get all descendant nodes of the given node (for TreeNode objects).
    
    Parameters
    ----------
    start_node : TreeNode
        The starting node
    logger_obj : logging.Logger, optional
        Logger instance for logging messages. If None, uses global logger.
    
    Returns
    -------
    set
        Set of all descendant node names
    """
    # Use provided logger or fall back to global logger
    log = logger_obj if logger_obj is not None else logging.getLogger(__name__)
    
    def get_descendants(node):
        """Recursively get all descendant nodes"""
        descendants = set()
        for child in node.children:
            descendants.add(child.name)
            descendants |= get_descendants(child)
        return descendants
    
    result = get_descendants(start_node)
    log.debug(f"Found {len(result)} descendants of '{start_node.name}'")
    
    return result


def get_subtree_nodes_tree_node(start_node, logger_obj=None):
    """
    Get all node names in the subtree starting from the given node (for TreeNode objects).
    
    Parameters
    ----------
    start_node : TreeNode
        The starting node
    logger_obj : logging.Logger, optional
        Logger instance for logging messages. If None, uses global logger.
    
    Returns
    -------
    set
        Set of all node names in the subtree
    """
    # Use provided logger or fall back to global logger
    log = logger_obj if logger_obj is not None else logging.getLogger(__name__)
    
    def collect_nodes(node):
        nodes = {node.name}
        for child in node.children:
            nodes |= collect_nodes(child)
        return nodes
    
    result = collect_nodes(start_node)
    log.debug(f"Collected {len(result)} nodes in subtree rooted at '{start_node.name}'")
    
    return result




# -------------------------
# MutationMatrix: Packaging pandas DataFrame
# -------------------------
class MutationMatrix:
    def __init__(self, df_posteriors: pd.DataFrame):
        # df_posteriors: cells x mutations, values in [0,1] or NA/pd.NA
        self.P = df_posteriors.copy()
    
    def copy(self):
        return MutationMatrix(self.P.copy())
    
    def cells(self) -> List[str]:
        return list(self.P.index)
    
    def mutations(self) -> List[str]:
        return list(self.P.columns)
    
    def posterior(self, cell: str, mut: str):
        return self.P.at[cell, mut]
    
    def binary_call(self, threshold: float = 0.5) -> pd.DataFrame:
        return self.P.applymap(lambda x: 1 if (pd.notna(x) and x > threshold) else (0 if pd.notna(x) else pd.NA))
    
    def n_cells(self):
        return self.P.shape[0]


# -------------------------
# Build an initial backbone tree (T_B, M_B)
# -------------------------
def build_backbone_tree(backbone_mutations: List[str]) -> TreeNode:
    """
    Build initial backbone tree with root and backbone mutations as children.
    
    Args:
        backbone_mutations: List of backbone mutation IDs
        
    Returns:
        Root node of backbone tree
    """
    root = TreeNode("ROOT")
    for mut in backbone_mutations:
        mut_node = TreeNode(mut)
        root.add_child(mut_node)
    
    return root

def impute_backbone_clones(selected_df, backbone_mutations, mutation_group, logger_obj=None):
    """
    Calculate backbone clones based on group information.
    
    Parameters
    ----------
    selected_df : pd.DataFrame
        DataFrame containing mutation data (cells × mutations)
    backbone_mutations : list
        List of backbone mutations
    mutation_group : dict
        Mapping of mutation to group ID
    
    Returns
    -------
    exclusive_backbone_matrix : pd.DataFrame
        Binary matrix for each backbone clone (cells × backbone_clones)
    final_assignment : pd.Series
        Final backbone clone assignment for each cell
    """
    # Use provided logger or fall back to global logger
    log = logger_obj if logger_obj is not None else logging.getLogger(__name__)
    
    # Replace NA values with 0
    selected_df = selected_df.fillna(0)
    
    # Find the backbone mutation for each group
    group_to_backbone = {}
    for mut in backbone_mutations:
        group = mutation_group[mut]
        group_to_backbone[group] = mut
    
    # Calculate imputed clone vector for each backbone mutation
    backbone_vectors = {}
    sum_vectors = {}
    for group, backbone_mutation in group_to_backbone.items():
        # Find all mutations in this group
        group_muts = [mut for mut, g in mutation_group.items() if g == group]
        
        # Sum the vectors of these mutations
        group_df = selected_df[group_muts]
        sum_vector = group_df.sum(axis=1)
        sum_vectors[backbone_mutation] = sum_vector
        # Set values > 0 to 1 to get binary vector
        imputed_vector = (sum_vector > 0).astype(int)
        backbone_vectors[backbone_mutation] = imputed_vector
    
    # Create backbone matrix (cells × backbone_clones)
    backbone_matrix = pd.DataFrame(backbone_vectors)
    sum_matrix = pd.DataFrame(sum_vectors)
    
    # Handle cell assignment: determine which backbone clone each cell belongs to
    cell_assignments = {}
    
    for cell in selected_df.index:
        # Get the values for this cell across all backbone clones
        cell_backbone_values = backbone_matrix.loc[cell]
        
        # Find which backbone clones this cell belongs to (value = 1)
        belonging_clones = cell_backbone_values[cell_backbone_values == 1].index.tolist()
        
        if len(belonging_clones) == 0:
            # Does not belong to any backbone clone
            cell_assignments[cell] = None
        
        elif len(belonging_clones) == 1:
            # Belongs to only one backbone clone, assign directly
            cell_assignments[cell] = belonging_clones[0]
        
        else:
            # Belongs to multiple backbone clones, need to resolve conflict
            # Calculate mutation count for this cell in each conflicting clone
            clone_scores = {}
            for clone in belonging_clones:
                mutation_count = sum_vectors[clone].loc[cell]
                clone_scores[clone] = mutation_count
            
            # Find the highest score
            max_score = max(clone_scores.values())
            
            # Find all clones with the highest score
            top_clones = [clone for clone, score in clone_scores.items() if score == max_score]
            
            if len(top_clones) == 1:
                # Only one highest score, assign normally
                assigned_clone = top_clones[0]
            else:
                # Multiple clones have the same score, choose clone with strongest backbone expression
                backbone_expr = {clone: selected_df.loc[cell, clone] for clone in top_clones}
                assigned_clone = max(backbone_expr.items(), key=lambda x: x[1])[0]
            
            cell_assignments[cell] = assigned_clone
    
    # Create final backbone assignment vector
    final_assignment = pd.Series(cell_assignments, name='backbone_assignment')
    
    # Create exclusive backbone matrix (each cell belongs to only one clone)
    exclusive_backbone_matrix = pd.DataFrame(0, 
                                           index=selected_df.index, 
                                           columns=backbone_matrix.columns)
    
    for cell, assignment in cell_assignments.items():
        if assignment is not None:
            exclusive_backbone_matrix.loc[cell, assignment] = 1
    
    return exclusive_backbone_matrix, final_assignment


def impute_backbone_clones_weighted(selected_df, backbone_mutations, mutation_group, backbone_weight=2.0, logger_obj=None):
    """
    Calculate backbone clones based on group information with weighted scoring,
    allowing cells to be assigned across backbone clones.
    
    Parameters
    ----------
    selected_df : pd.DataFrame
        DataFrame containing mutation data (cells × mutations)
    backbone_mutations : list
        List of backbone mutations
    mutation_group : dict
        Mapping of mutation to group ID
    backbone_weight : float, default=2.0
        Weight coefficient for backbone mutations
        weight=1: all mutations equal
        weight>1: backbone mutations weighted higher
        weight<1: other mutations weighted higher
    
    Returns
    -------
    exclusive_backbone_matrix : pd.DataFrame
        Binary matrix with each cell assigned to only one backbone clone
    final_assignment : pd.Series
        Assigned backbone clone for each cell
    """
    # Use provided logger or fall back to global logger
    log = logger_obj if logger_obj is not None else logging.getLogger(__name__)
    
    # Replace NA values with 0
    selected_df = selected_df.fillna(0)
    
    # Find the backbone mutation for each group
    group_to_backbone = {}
    for mut in backbone_mutations:
        group = mutation_group[mut]
        group_to_backbone[group] = mut
    
    # Calculate weighted score for each backbone mutation
    weighted_scores = {}
    sum_vectors = {}  # Store original sum vectors for statistics
    
    for group, backbone_mutation in group_to_backbone.items():
        # Find all mutations in this group
        group_muts = [mut for mut, g in mutation_group.items() if g == group]
        group_df = selected_df[group_muts]
        
        # Calculate original sum vector (for statistics)
        sum_vector = group_df.sum(axis=1)
        sum_vectors[backbone_mutation] = sum_vector
        
        # Calculate weighted score: backbone mutation has higher weight
        weighted_vector = group_df.copy()
        if backbone_mutation in weighted_vector.columns:
            weighted_vector[backbone_mutation] = weighted_vector[backbone_mutation] * backbone_weight
        
        weighted_score = weighted_vector.sum(axis=1)
        weighted_scores[backbone_mutation] = weighted_score
    
    # Create score matrix
    weighted_matrix = pd.DataFrame(weighted_scores)
    
    # Assign cells to backbone clones
    cell_assignments = {}
    
    for cell in selected_df.index:
        # Get weighted scores for this cell across all backbone clones
        cell_scores = weighted_matrix.loc[cell]
        max_score = cell_scores.max()
        
        # If max score is 0, don't assign any clone
        if max_score == 0:
            cell_assignments[cell] = None
            continue
        
        # Find all clones with the highest score
        top_clones = cell_scores[cell_scores == max_score].index.tolist()
        
        if len(top_clones) == 1:
            # Only one highest score, assign directly
            assigned_clone = top_clones[0]
        else:
            # Multiple clones have the same score, choose clone with strongest backbone expression
            backbone_expr = {clone: selected_df.loc[cell, clone] for clone in top_clones}
            assigned_clone = max(backbone_expr.items(), key=lambda x: x[1])[0]
        
        cell_assignments[cell] = assigned_clone
    
    # Create final backbone assignment vector
    final_assignment = pd.Series(cell_assignments, name='backbone_assignment')
    
    # Create exclusive backbone matrix (each cell belongs to only one clone)
    exclusive_backbone_matrix = pd.DataFrame(0, 
                                           index=selected_df.index, 
                                           columns=backbone_mutations)
    
    for cell, assignment in cell_assignments.items():
        if assignment is not None:
            exclusive_backbone_matrix.loc[cell, assignment] = 1
    
    return exclusive_backbone_matrix, final_assignment




import pandas as pd
import numpy as np
from collections import defaultdict, Counter
from typing import Dict, List, Tuple

def find_intersecting_mutations(I: pd.DataFrame, backbone_mutation: str, min_N11: int = 1) -> List[str]:
    """Find all mutations that co-occur with the backbone mutation"""
    intersecting_muts = []
    
    for other_mut in I.columns:
        if other_mut == backbone_mutation:
            continue
            
        counts = pairwise_counts(I, backbone_mutation, other_mut)
        if counts['N11'] >= min_N11:
            intersecting_muts.append(other_mut)
    
    return intersecting_muts

def impute_backbone_clones_comprehensive(I_somatic: pd.DataFrame, 
                                       backbone_mutations: List[str], 
                                       min_N11: int = 1,
                                       remove_zeros: bool = False,
                                       logger_obj=None) -> Tuple[pd.DataFrame, pd.Series, dict]:
    """
    Perform comprehensive imputation based on the entire mutation matrix and backbone mutations.
    
    Parameters
    ----------
    I_somatic : pd.DataFrame
        DataFrame containing all mutation data (cells × mutations)
    backbone_mutations : list
        List of backbone mutations
    min_N11 : int, default=1
        Minimum co-occurrence threshold
    remove_zeros : bool, default=False
        Whether to automatically remove all-zero rows and columns
    logger_obj : logging.Logger, optional
        Logger instance for logging messages. If None, uses global logger.
    
    Returns
    -------
    exclusive_backbone_matrix : pd.DataFrame
        Assignment of each cell to backbone clones
    final_assignment : pd.Series
        Final backbone assignment for each cell
    intersection_sets : dict
        Intersection mutation sets for each backbone mutation
    """
    # Use provided logger or fall back to global logger
    log = logger_obj if logger_obj is not None else logging.getLogger(__name__)
    
    # Create a copy to avoid modifying original data
    I = I_somatic.copy()
    
    # Replace NA values with 0
    I = I.fillna(0)
    
    log.info(f"Original data shape: {I.shape}")
    log.info(f"Original cell count: {len(I.index)}")
    log.info(f"Original mutation count: {len(I.columns)}")
    
    # Optional: automatically remove all-zero rows and columns
    if remove_zeros:
        # Remove all-zero columns
        I = I.loc[:, (I != 0).any(axis=0)]
        # Remove all-zero rows
        I = I.loc[(I != 0).any(axis=1)]
        log.info(f"Cleaned data shape: {I.shape}")
        log.info(f"Cleaned cell count: {len(I.index)}")
        log.info(f"Cleaned mutation count: {len(I.columns)}")
    
    # Identify cells with mutations (at least one mutation value = 1)
    has_mutation = (I == 1).any(axis=1)
    mutated_cells = I.index[has_mutation]
    non_mutated_cells = I.index[~has_mutation]
    
    log.info(f"\nCell statistics:")
    log.info(f"  Cells with mutations: {len(mutated_cells)}")
    log.info(f"  Cells without mutations: {len(non_mutated_cells)}")
    
    # Only process cells with mutations for backbone assignment
    I_mutated = I.loc[mutated_cells]
    
    # Find intersecting mutations for each backbone mutation
    intersection_sets = {}
    backbone_vectors = {}
    sum_vectors = {}
    
    log.info("\nFinding intersecting mutations for each backbone mutation...")
    for backbone_mut in backbone_mutations:
        if backbone_mut not in I.columns:
            log.warning(f"Backbone mutation {backbone_mut} not in data columns, skipping")
            continue
            
        # Find all mutations co-occurring with this backbone mutation
        intersecting_muts = find_intersecting_mutations(I_mutated, backbone_mut, min_N11)
        intersection_sets[backbone_mut] = intersecting_muts
        
        log.info(f"Backbone {backbone_mut}: found {len(intersecting_muts)} intersecting mutations")
        
        # Always include backbone mutation itself + intersecting mutations
        all_relevant_muts = [backbone_mut] + intersecting_muts
        
        # Ensure all mutations are in the data
        available_muts = [mut for mut in all_relevant_muts if mut in I_mutated.columns]
        
        # Calculate sum vector of these mutations
        group_df = I_mutated[available_muts]
        sum_vector = group_df.sum(axis=1)
        sum_vectors[backbone_mut] = sum_vector
        
        # Binarize: if any relevant mutation is expressed, the backbone clone is active
        imputed_vector = (sum_vector > 0).astype(int)
        backbone_vectors[backbone_mut] = imputed_vector
    
    # Create backbone matrix (only on mutated cells)
    backbone_matrix = pd.DataFrame(backbone_vectors, index=mutated_cells)
    
    # Calculate mutant cell number for each mutation (for tie-breaking)
    mutation_prevalence = I_mutated.sum(axis=0)
    
    # Handle cell assignment
    cell_assignments = {}
    conflict_resolution_stats = defaultdict(int)
    unassigned_cells_mutations = {}  # Record mutation info for unassigned cells
    
    log.info(f"\nProcessing {len(mutated_cells)} mutated cells for assignment...")
    for i, cell in enumerate(mutated_cells):
        if i > 0 and i % 1000 == 0:
            log.info(f"  Processed {i}/{len(mutated_cells)} cells...")
        
        # Get this cell's values across all backbone clones
        cell_backbone_values = backbone_matrix.loc[cell]
        
        # Find which backbone clones this cell belongs to (value = 1)
        belonging_clones = cell_backbone_values[cell_backbone_values == 1].index.tolist()
        
        if len(belonging_clones) == 0:
            # Cell has mutations but doesn't belong to any backbone clone
            cell_assignments[cell] = None
            conflict_resolution_stats['mutated_but_unassigned'] += 1
            
            # Record this cell's mutations for debugging
            cell_mutations = I_mutated.columns[I_mutated.loc[cell] == 1].tolist()
            unassigned_cells_mutations[cell] = cell_mutations
        
        elif len(belonging_clones) == 1:
            # Belongs to only one backbone clone, assign directly
            cell_assignments[cell] = belonging_clones[0]
            conflict_resolution_stats['unique_assignment'] += 1
        
        else:
            # Belongs to multiple backbone clones, need to resolve conflict
            conflict_resolution_stats['conflicts'] += 1
            
            # Level 1: Count supporting mutations in each conflicting clone
            clone_scores = {}
            for clone in belonging_clones:
                support_count = sum_vectors[clone].loc[cell]
                clone_scores[clone] = support_count
            
            # Find the highest score
            max_score = max(clone_scores.values())
            
            # Find all clones with the highest score
            top_clones = [clone for clone, score in clone_scores.items() if score == max_score]
            
            if len(top_clones) == 1:
                # Only one highest score, assign normally
                assigned_clone = top_clones[0]
                conflict_resolution_stats['resolved_unique_max'] += 1
            else:
                # Multiple clones have the same score, need secondary judgment
                conflict_resolution_stats['ties'] += 1
                
                # Level 2: Compare prevalence of supporting mutations
                prevalence_scores = {}
                for clone in top_clones:
                    # Get all relevant mutations for this backbone clone
                    supporting_muts = [clone] + intersection_sets[clone]
                    available_muts = [mut for mut in supporting_muts if mut in mutation_prevalence.index]
                    # Calculate average prevalence of supporting mutations
                    avg_prevalence = mutation_prevalence[available_muts].mean()
                    prevalence_scores[clone] = avg_prevalence
                
                # Find clone with highest prevalence
                max_prevalence = max(prevalence_scores.values())
                top_prevalence_clones = [clone for clone, score in prevalence_scores.items() if score == max_prevalence]
                
                if len(top_prevalence_clones) == 1:
                    assigned_clone = top_prevalence_clones[0]
                    conflict_resolution_stats['resolved_by_prevalence'] += 1
                else:
                    # Level 3: Choose clone with largest imputed clone size
                    conflict_resolution_stats['final_ties'] += 1
                    clone_sizes = {}
                    for clone in top_prevalence_clones:
                        clone_size = backbone_vectors[clone].sum()
                        clone_sizes[clone] = clone_size
                    
                    max_size = max(clone_sizes.values())
                    max_size_clones = [clone for clone, size in clone_sizes.items() if size == max_size]
                    assigned_clone = max_size_clones[0]
                    conflict_resolution_stats['resolved_by_size'] += 1
            
            cell_assignments[cell] = assigned_clone
            conflict_resolution_stats[f'resolved_to_{assigned_clone}'] += 1
    
    # Cells without mutations remain unassigned
    for cell in non_mutated_cells:
        cell_assignments[cell] = None
        conflict_resolution_stats['no_mutation'] += 1
    
    # Critical fix: ensure all original cells are in the assignment dictionary
    all_original_cells = set(I_somatic.index)
    assigned_cells = set(cell_assignments.keys())
    missing_cells = all_original_cells - assigned_cells
    
    if missing_cells:
        log.warning(f"Found {len(missing_cells)} cells not assigned, setting them to unassigned")
        for cell in missing_cells:
            cell_assignments[cell] = None
            conflict_resolution_stats['missing_cells_set_to_unassigned'] += 1
    
    # Create final backbone assignment vector (all original cells)
    final_assignment = pd.Series(cell_assignments, name='backbone_assignment')
    
    # Create exclusive backbone matrix (all original cells)
    exclusive_backbone_matrix = pd.DataFrame(0, 
                                           index=I_somatic.index,  # Use original index
                                           columns=backbone_mutations)
    
    # Only set values for mutated cells that have assignments
    for cell, assignment in cell_assignments.items():
        if assignment is not None:
            exclusive_backbone_matrix.loc[cell, assignment] = 1
    
    # Analyze mutations with no intersection with backbone mutations
    log.info(f"\n=== Mutations with no intersection with backbone mutations ===")
    
    # Find all mutations that have no intersection with any backbone mutation
    all_mutations = set(I_mutated.columns)
    mutations_with_intersection = set()
    
    for backbone_mut, intersecting_muts in intersection_sets.items():
        mutations_with_intersection.add(backbone_mut)
        mutations_with_intersection.update(intersecting_muts)
    
    mutations_without_intersection = all_mutations - mutations_with_intersection
    
    log.info(f"Total mutations: {len(all_mutations)}")
    log.info(f"Mutations with intersection with backbone: {len(mutations_with_intersection)}")
    log.info(f"Mutations without intersection with backbone: {len(mutations_without_intersection)}")
    
    if mutations_without_intersection:
        log.info("\nMutations without intersection with backbone:")
        total_cells_in_isolated_muts = 0
        for mut in sorted(mutations_without_intersection):
            mut_cell_count = (I_mutated[mut] == 1).sum()
            total_cells_in_isolated_muts += mut_cell_count
            log.info(f"  {mut}: {mut_cell_count} cells")
        log.info(f"These isolated mutations affect {total_cells_in_isolated_muts} cells total")
    
    # Analyze mutation composition of unassigned cells
    if unassigned_cells_mutations:
        log.info(f"\n=== Detailed analysis of unassigned cells ===")
        log.info(f"Number of unassigned cells: {len(unassigned_cells_mutations)}")
        
        # Count mutation distribution in these cells
        all_unassigned_mutations = []
        for cell, muts in unassigned_cells_mutations.items():
            all_unassigned_mutations.extend(muts)
        
        mutation_counts = Counter(all_unassigned_mutations)
        
        log.info(f"Unassigned cells involve {len(mutation_counts)} different mutations")
        
        # Categorize by mutation type
        backbone_muts_in_unassigned = []
        intersecting_muts_in_unassigned = []
        isolated_muts_in_unassigned = []
        
        for mut, count in mutation_counts.items():
            if mut in backbone_mutations:
                backbone_muts_in_unassigned.append((mut, count))
            elif mut in mutations_with_intersection:
                intersecting_muts_in_unassigned.append((mut, count))
            else:
                isolated_muts_in_unassigned.append((mut, count))
        
        log.info(f"\nMutation types in unassigned cells:")
        log.info(f"  - Backbone mutations: {len(backbone_muts_in_unassigned)}")
        if backbone_muts_in_unassigned:
            for mut, count in sorted(backbone_muts_in_unassigned, key=lambda x: x[1], reverse=True):
                log.info(f"    {mut}: {count} cells")
        
        log.info(f"  - Mutations with backbone intersection: {len(intersecting_muts_in_unassigned)}")
        if intersecting_muts_in_unassigned:
            for mut, count in sorted(intersecting_muts_in_unassigned, key=lambda x: x[1], reverse=True)[:10]:
                log.info(f"    {mut}: {count} cells")
        
        log.info(f"  - Isolated mutations (no backbone intersection): {len(isolated_muts_in_unassigned)}")
        if isolated_muts_in_unassigned:
            log.info(f"    (showing first 20)")
            for mut, count in sorted(isolated_muts_in_unassigned, key=lambda x: x[1], reverse=True)[:20]:
                log.info(f"    {mut}: {count} cells")
        
        # Validation: check how many unassigned cells have only isolated mutations
        cells_with_only_isolated_muts = 0
        cells_with_backbone_or_intersecting = 0
        
        for cell, muts in unassigned_cells_mutations.items():
            has_backbone_or_intersecting = any(mut in mutations_with_intersection for mut in muts)
            if has_backbone_or_intersecting:
                cells_with_backbone_or_intersecting += 1
            else:
                cells_with_only_isolated_muts += 1
        
        log.info(f"\nUnassigned cells validation:")
        log.info(f"  - Cells with only isolated mutations: {cells_with_only_isolated_muts}")
        log.info(f"  - Cells with backbone or intersecting mutations: {cells_with_backbone_or_intersecting}")
        
        # If there are cells with backbone or intersecting mutations but still unassigned
        if cells_with_backbone_or_intersecting > 0:
            log.warning(f"\n{ cells_with_backbone_or_intersecting} cells have backbone or intersecting mutations but remain unassigned")
            log.warning("This may indicate imputation logic issues, please check:")
            
            problematic_cells = []
            for cell, muts in unassigned_cells_mutations.items():
                if any(mut in mutations_with_intersection for mut in muts):
                    problematic_cells.append((cell, muts))
            
            # Show first few problematic cells as examples
            for cell, muts in problematic_cells[:3]:
                backbone_muts_in_cell = [mut for mut in muts if mut in backbone_mutations]
                intersecting_muts_in_cell = [mut for mut in muts if mut in mutations_with_intersection and mut not in backbone_mutations]
                log.info(f"  Cell {cell}:")
                log.info(f"    Backbone mutations: {backbone_muts_in_cell}")
                log.info(f"    Intersecting mutations: {intersecting_muts_in_cell[:5]}...")  # Only show first 5
    
    # Final validation
    total_assigned_cells = exclusive_backbone_matrix.sum().sum()
    total_input_cells = len(I_somatic.index)
    total_mutated_cells = len(mutated_cells)
    total_final_assignment_cells = len(final_assignment)
    
    log.info(f"\n=== Final Validation ===")
    log.info(f"Original input cells: {total_input_cells}")
    log.info(f"Final assignment vector cells: {total_final_assignment_cells}")
    log.info(f"Cells with mutations: {total_mutated_cells}")
    log.info(f"Cells without mutations: {len(non_mutated_cells)}")
    log.info(f"Cells assigned to backbone: {total_assigned_cells}")
    log.info(f"Mutated but unassigned cells: {conflict_resolution_stats['mutated_but_unassigned']}")
    log.info(f"Assignment rate: {total_assigned_cells/total_mutated_cells*100:.2f}%")
    
    # Use assertions with friendly error messages
    if total_input_cells != total_final_assignment_cells:
        raise ValueError(f"Cell count mismatch! Input {total_input_cells} != Final assignment {total_final_assignment_cells}")
    
    if total_mutated_cells != total_assigned_cells + conflict_resolution_stats['mutated_but_unassigned']:
        raise ValueError(f"Mutated cell count mismatch! {total_mutated_cells} != {total_assigned_cells} + {conflict_resolution_stats['mutated_but_unassigned']}")
    
    log.info("\nConflict resolution statistics:")
    for stat, count in sorted(conflict_resolution_stats.items()):
        if count > 0:
            log.info(f"  {stat}: {count} cells")
    
    # Output final backbone clone sizes
    log.info("\nFinal backbone clone sizes:")
    total_backbone_cells = 0
    for backbone_mut in backbone_mutations:
        if backbone_mut in exclusive_backbone_matrix.columns:
            clone_size = exclusive_backbone_matrix[backbone_mut].sum()
            intersecting_count = len(intersection_sets.get(backbone_mut, []))
            log.info(f"  {backbone_mut}: {clone_size} cells, {intersecting_count} intersecting mutations")
            total_backbone_cells += clone_size
    
    log.info(f"Total cells across all backbone clones: {total_backbone_cells}")
    log.info(f"Mutated cells not assigned to any backbone: {conflict_resolution_stats['mutated_but_unassigned']}")
    
    return exclusive_backbone_matrix, final_assignment, intersection_sets




# -------------------------
# merge_identical_columns(matrix)
# Input: original mutation-by-cell matrix (each column is a mutation)
# Output: matrix with identical columns merged, column names become "mut1|mut2|mut3"
#         Also returns a mapping dictionary for reverse splitting
# split_merged_columns(matrix, mapping)
# Input: merged matrix + mapping dictionary
# Output: restored original mutation-by-cell matrix (each column is a mutation)
# -------------------------
def merge_duplicate_columns(df: pd.DataFrame) -> Tuple[pd.DataFrame, Dict[str, str]]:
    """
    Merge identical columns in the original df (cells x mutations) into 'a|b' format.
    
    Returns (merged_df, col_map):
      - merged_df: cells x merged_cols
      - col_map: original_col -> merged_col_name
    
    Note: pd.NA/np.nan values are treated as None in the signature.
    """
    sig_map = {}
    for col in df.columns:
        # signature: tuple of values where NaN -> None (hashable)
        sig = tuple(None if pd.isna(v) else v for v in df[col].to_list())
        sig_map.setdefault(sig, []).append(col)
    
    merged_cols = {}
    col_map = {}
    for cols in sig_map.values():
        merged_name = "|".join(sorted(cols))
        # Use the first column's data as representative
        merged_cols[merged_name] = df[cols[0]].copy()
        for c in cols:
            col_map[c] = merged_name
    
    merged_df = pd.DataFrame(merged_cols, index=df.index)
    
    return merged_df, col_map


def split_merged_columns(merged_matrix: pd.DataFrame, mut_list: list):
    """
    Split merged columns based on the mutation list.
    
    Parameters
    ----------
    merged_matrix : pd.DataFrame
        Dataframe with merged columns
    mut_list : list
        List of original mutation names to expand to
    
    Returns
    -------
    pd.DataFrame
        Dataframe with expanded columns (mut_list)
    """
    if merged_matrix.empty or not mut_list:
        return pd.DataFrame(index=merged_matrix.index, columns=mut_list)
    
    mutation_to_merged_col = {}
    for merged_col in merged_matrix.columns:
        for mutation in merged_col.split("|"):
            if mutation not in mutation_to_merged_col:
                mutation_to_merged_col[mutation] = merged_col
    
    result_data = {}
    na_series = pd.Series(pd.NA, index=merged_matrix.index)
    
    for mut in mut_list:
        if mut in mutation_to_merged_col:
            result_data[mut] = merged_matrix[mutation_to_merged_col[mut]]
        else:
            result_data[mut] = na_series
    
    return pd.DataFrame(result_data, index=merged_matrix.index)[mut_list]

# # Test
# df = pd.DataFrame({
#     'A': [1, 2, np.nan],
#     'B': [1, 2, np.nan],  # the same as A
#     'C': [1, 2, None],    # the same as A
#     'D': [1, 2, 3]        # different
# })
# mut_list = ['A', 'B', 'C', 'D']
# 
# merged = merge_duplicate_columns(df)
# recovered = split_merged_columns(merged, mut_list)




# -------------------------
# find potential positions that is intersected with new_mut based on backbone tree for scaffold tree
#   - intersection_nodes
#   - parent_dict
# -------------------------

def find_intersection_positions_within_group_directly_scaffold(
    T_current: TreeNode,
    new_mut: str,
    matrix,
    mutation_group: dict,
    min_overlap: int = 1,
    intersection_nodes: Optional[Set] = None,
    logger_obj=None
):
    """
    Optimized version based on intersection analysis, directly finding relevant positions,
    only considering nodes under clones that belong to the same group as the current mutation.
    """
    # Use provided logger or fall back to global logger
    log = logger_obj if logger_obj is not None else logging.getLogger(__name__)
    
    # 1. Find all nodes that intersect with the target mutation
    if intersection_nodes is None:
        intersection_nodes = find_all_intersect_muts_from_tree_by_matrix_scaffold(
            T_current, matrix, new_mut, min_overlap, logger_obj=log
        )
    
    if len(intersection_nodes) == 0:
        log.debug(f"No intersection nodes found for {new_mut}")
        return []
    
    # 2. Build tree parent-child dictionary
    tree_parent_dict = build_tree_parent_dict_scaffold(T_current)
    
    # 3. Get the group of the new mutation
    target_group = mutation_group[new_mut]
    
    # 4. Find all path nodes, filtered to only clones in the target group (returns sorted list)
    all_path_nodes = get_all_path_nodes_with_group_filter(intersection_nodes, tree_parent_dict, mutation_group, target_group)
    
    # 5. Pre-create a deep copy of the base tree
    base_tree_copy = deepcopy(T_current)
    
    # 6. Generate candidate positions only on relevant nodes
    candidate_positions = []
    
    # ============================================================
    # CRITICAL FIX: Iterate over sorted all_path_nodes to ensure deterministic order
    # ============================================================
    for node_name in all_path_nodes:
        if node_name == "ROOT":
            continue
            
        node = base_tree_copy.find(node_name)
        if node is None:
            log.warning(f"Node {node_name} not found in tree")
            continue
        
        candidate_positions.append(_create_on_node_candidate_fast_scaffold(base_tree_copy, node, new_mut))
        candidate_positions.append(_create_new_leaf_candidate_fast_scaffold(base_tree_copy, node, new_mut))
        
        for child in node.children:
            if child.name in all_path_nodes:
                candidate_positions.append(_create_on_edge_candidate_fast_scaffold(base_tree_copy, node, child, new_mut))
        
        if len(node.children) >= 2:
            path_children = [child for child in node.children if child.name in all_path_nodes]
            if len(path_children) >= 2:
                for r in range(2, min(4, len(path_children) + 1)):
                    for combo in combinations(path_children, r):
                        candidate_positions.append(_create_merge_candidate_fast_scaffold(base_tree_copy, node, combo, new_mut))
    
    return candidate_positions


def find_new_leaf_positions_for_target_node(T_current: TreeNode, new_mut: str, matrix, target_node, min_overlap=1, logger_obj=None):
    """
    Optimized version based on intersection analysis, but only returns new_leaf positions
    related to the target node.
    
    Parameters
    ----------
    T_current : TreeNode
        Current tree
    new_mut : str
        New mutation to place
    matrix : pd.DataFrame
        Mutation matrix
    target_node : TreeNode
        Target node to filter positions
    min_overlap : int, default=1
        Minimum overlap threshold
    logger_obj : logging.Logger, optional
        Logger instance for logging messages. If None, uses global logger.
    
    Returns
    -------
    list
        List of new_leaf candidate positions anchored at the target node
    """
    # Use provided logger or fall back to global logger
    log = logger_obj if logger_obj is not None else logging.getLogger(__name__)
    
    # 1. Find all nodes that intersect with the target mutation
    intersection_nodes = find_all_intersect_muts_from_tree_by_matrix_scaffold(
        T_current, matrix, new_mut, min_overlap, logger_obj=log
    )
    
    if len(intersection_nodes) == 0:
        log.debug(f"No intersection nodes found for {new_mut}")
        return []  # No intersection, return empty list
    
    # 2. Build tree parent-child dictionary
    tree_parent_dict = build_tree_parent_dict_scaffold(T_current)
    
    # 3. Find all path nodes
    all_path_nodes = find_all_path_nodes_scaffold(intersection_nodes, tree_parent_dict)
    
    # 4. Pre-create a deep copy of the base tree
    base_tree_copy = deepcopy(T_current)
    
    # 5. Generate candidate positions only on relevant nodes
    candidate_positions = []
    
    for node_name in all_path_nodes:
        if node_name == "ROOT":
            # Keep on_node type positions on ROOT
            node = base_tree_copy.find(node_name)
            if node is None:
                log.warning(f"Node {node_name} not found in tree")
                continue
            
            # Generate on_node candidate on ROOT
            candidate_positions.append(_create_on_node_candidate_fast_scaffold(base_tree_copy, node, new_mut))
            
            continue  # Skip other position types on ROOT
        
        node = base_tree_copy.find(node_name)
        if node is None:
            log.warning(f"Node {node_name} not found in tree")
            continue
        
        # --- 1) Place on node ---
        candidate_positions.append(_create_on_node_candidate_fast_scaffold(base_tree_copy, node, new_mut))
        
        # --- 2) New leaf ---
        candidate_positions.append(_create_new_leaf_candidate_fast_scaffold(base_tree_copy, node, new_mut))
        
        # --- 3) Place on each edge ---
        for child in node.children:
            if child.name in all_path_nodes:  # Only consider path children
                candidate_positions.append(_create_on_edge_candidate_fast_scaffold(base_tree_copy, node, child, new_mut))
        
        # --- 4) New parent merging multiple children ---
        if len(node.children) >= 2:
            # Only consider path children combinations
            path_children = [child for child in node.children if child.name in all_path_nodes]
            if len(path_children) >= 2:
                # Limit combinations to avoid explosion
                for r in range(2, min(4, len(path_children) + 1)):
                    for combo in combinations(path_children, r):
                        candidate_positions.append(_create_merge_candidate_fast_scaffold(base_tree_copy, node, combo, new_mut))
    
    # 6. Filter positions where anchor is target_node and placement_type is new_leaf
    target_positions = [
        pos for pos in candidate_positions 
        if pos['placement_type'] == 'new_leaf' and pos['anchor'] == target_node.name
    ]
    
    return target_positions


def find_intersection_positions_within_tree_directly_scaffold(
    T_current: TreeNode,
    new_mut: str,
    matrix,
    min_overlap: int = 1,
    intersection_nodes: Optional[Set] = None,
    logger_obj=None
):
    """
    Optimized version based on intersection analysis, directly finding relevant positions.
    """
    # Use provided logger or fall back to global logger
    log = logger_obj if logger_obj is not None else logging.getLogger(__name__)
    
    # 1. Find all nodes that intersect with the target mutation
    if intersection_nodes is None:
        intersection_nodes = find_all_intersect_muts_from_tree_by_matrix_scaffold(
            T_current, matrix, new_mut, min_overlap, logger_obj=log
        )
    
    if len(intersection_nodes) == 0:
        log.debug(f"No intersection nodes found for {new_mut}")
        return []
    
    # 2. Build tree parent-child dictionary
    tree_parent_dict = build_tree_parent_dict_scaffold(T_current)
    
    # 3. Find all path nodes (returns sorted list for deterministic iteration)
    all_path_nodes = find_all_path_nodes_scaffold(intersection_nodes, tree_parent_dict)
    
    # 4. Pre-create a deep copy of the base tree
    base_tree_copy = deepcopy(T_current)
    
    # 5. Generate candidate positions only on relevant nodes
    candidate_positions = []
    
    # ============================================================
    # CRITICAL FIX: Iterate over sorted all_path_nodes to ensure deterministic order
    # ============================================================
    for node_name in all_path_nodes:
        if node_name == "ROOT":
            node = base_tree_copy.find(node_name)
            if node is None:
                log.warning(f"Node {node_name} not found in tree")
                continue
            candidate_positions.append(_create_on_node_candidate_fast_scaffold(base_tree_copy, node, new_mut))
            continue
        
        node = base_tree_copy.find(node_name)
        if node is None:
            log.warning(f"Node {node_name} not found in tree")
            continue
        
        # --- 1) Place on node ---
        candidate_positions.append(_create_on_node_candidate_fast_scaffold(base_tree_copy, node, new_mut))
        
        # --- 2) New leaf ---
        candidate_positions.append(_create_new_leaf_candidate_fast_scaffold(base_tree_copy, node, new_mut))
        
        # --- 3) Place on each edge ---
        for child in node.children:
            if child.name in all_path_nodes:
                candidate_positions.append(_create_on_edge_candidate_fast_scaffold(base_tree_copy, node, child, new_mut))
        
        # --- 4) New parent merging multiple children ---
        if len(node.children) >= 2:
            path_children = [child for child in node.children if child.name in all_path_nodes]
            if len(path_children) >= 2:
                for r in range(2, min(4, len(path_children) + 1)):
                    for combo in combinations(path_children, r):
                        candidate_positions.append(_create_merge_candidate_fast_scaffold(base_tree_copy, node, combo, new_mut))
    
    return candidate_positions


def filter_candidate_positions_from_target_node(candidate_positions, subtree_nodes, logger_obj=None):
    """
    Filter candidate positions to only those anchored in the specified subtree nodes.
    
    Parameters
    ----------
    candidate_positions : list
        List of candidate positions
    subtree_nodes : set
        Set of subtree node names
    logger_obj : logging.Logger, optional
        Logger instance for logging messages. If None, uses global logger.
    
    Returns
    -------
    list
        Filtered list of candidate positions
    """
    # Use provided logger or fall back to global logger
    log = logger_obj if logger_obj is not None else logging.getLogger(__name__)
    
    filtered_positions = []
    
    for pos in candidate_positions:
        # Check if this position's anchor is in the subtree node set
        if pos.get('anchor') in subtree_nodes:
            filtered_positions.append(pos)
    
    log.info(f"Filtered {len(filtered_positions)} out of {len(candidate_positions)} positions within target subtree")
    return filtered_positions


def get_all_path_nodes_with_group_filter(intersection_nodes, tree_parent_dict, mutation_group, target_group):
    """
    Get path nodes related to intersection nodes, but only those belonging to the target_group
    in mutation_group.
    
    Parameters
    ----------
    intersection_nodes : set or list
        Set/list of intersection node names
    tree_parent_dict : dict
        Parent-child dictionary
    mutation_group : dict
        Mapping of mutation to group ID
    target_group : int or str
        Target group ID to filter nodes
    
    Returns
    -------
    list
        Sorted list of path nodes belonging to the target group
    """
    all_path_nodes = set()
    all_path_nodes.add('ROOT')  # Always include ROOT
    
    for node in intersection_nodes:
        path = get_path_to_root_scaffold(node, tree_parent_dict)
        all_path_nodes.update(path)
    
    # Paths between intersection nodes
    intersection_list = list(intersection_nodes)
    for i in range(len(intersection_list)):
        for j in range(i + 1, len(intersection_list)):
            path_between = get_path_between_nodes_scaffold(intersection_list[i], intersection_list[j], tree_parent_dict)
            all_path_nodes.update(path_between)
    
    # Only keep nodes belonging to target_group
    all_path_nodes = {node for node in all_path_nodes if valid_node_with_group(node, mutation_group, target_group)}
    
    # ============================================================
    # CRITICAL FIX: Return sorted list to ensure deterministic iteration order
    # ============================================================
    return sorted(all_path_nodes)


def valid_node_with_group(node, mutation_group, target_group):
    """
    Check if a node belongs to the target group in mutation_group.
    """
    node_muts = node.split("|")
    for mut in node_muts:
        if mut in mutation_group and mutation_group[mut] == target_group:
            return True
    return False


def build_tree_parent_dict_scaffold(tree):
    """Build parent-child dictionary directly from the tree."""
    parent_dict = {}
    for node in tree.traverse():
        for child in node.children:
            parent_dict[child.name] = node.name
    return parent_dict


def _create_on_node_candidate_fast_scaffold(base_tree_copy, node, new_mut):
    """Quickly create a candidate for placing on a node."""
    new_tree = deepcopy(base_tree_copy)
    anchor_node = new_tree.find(node.name)
    
    # Replace node logic
    new_node = TreeNode(new_mut)
    for child in anchor_node.children:
        new_node.add_child(deepcopy(child))
    
    if anchor_node.parent:
        parent = anchor_node.parent
        for i, child in enumerate(parent.children):
            if child.name == node.name:
                parent.children[i] = new_node
                new_node.parent = parent
                break
    else:
        new_tree = new_node
    
    return {
        "placement_type": "on_node",
        "anchor": node.name,
        "meta": {},
        "nodes": _extract_nodes_info_sacffold(new_tree),
        "edges": _extract_edges_info_scaffold(new_tree)
    }


def _create_new_leaf_candidate_fast_scaffold(base_tree_copy, node, new_mut):
    """Quickly create a candidate for adding a new leaf."""
    new_tree = deepcopy(base_tree_copy)
    anchor_node = new_tree.find(node.name)
    new_leaf = TreeNode(new_mut)
    anchor_node.add_child(new_leaf)
    
    return {
        "placement_type": "new_leaf",
        "anchor": node.name,
        "meta": {},
        "nodes": _extract_nodes_info_sacffold(new_tree),
        "edges": _extract_edges_info_scaffold(new_tree)
    }


def _create_on_edge_candidate_fast_scaffold(base_tree_copy, parent_node, child_node, new_mut):
    """Quickly create a candidate for placing on an edge."""
    new_tree = deepcopy(base_tree_copy)
    parent = new_tree.find(parent_node.name)
    child = new_tree.find(child_node.name)
    new_node = TreeNode(new_mut)
    parent.remove_child(child)
    parent.add_child(new_node)
    new_node.add_child(child)
    
    return {
        "placement_type": "on_edge",
        "anchor": parent_node.name,
        "meta": {"child": child_node.name},
        "nodes": _extract_nodes_info_sacffold(new_tree),
        "edges": _extract_edges_info_scaffold(new_tree)
    }


def _create_merge_candidate_fast_scaffold(base_tree_copy, parent_node, children_combo, new_mut):
    """Quickly create a candidate for merging child nodes."""
    new_tree = deepcopy(base_tree_copy)
    parent = new_tree.find(parent_node.name)
    combo_nodes = [new_tree.find(child.name) for child in children_combo]
    new_parent = TreeNode(new_mut)
    for cn in combo_nodes:
        parent.remove_child(cn)
        new_parent.add_child(cn)
    parent.add_child(new_parent)
    
    return {
        "placement_type": "new_parent_merge",
        "anchor": parent_node.name,
        "meta": {"merge_children": [child.name for child in children_combo]},
        "nodes": _extract_nodes_info_sacffold(new_tree),
        "edges": _extract_edges_info_scaffold(new_tree)
    }


def _extract_nodes_info_sacffold(tree):
    """Extract node information."""
    return [{"name": n.name,
             "parent": n.parent.name if n.parent else None,
             "children": [c.name for c in n.children]} 
            for n in tree.traverse()]


def _extract_edges_info_scaffold(tree):
    """Extract edge information."""
    return [(n.parent.name, n.name) for n in tree.traverse() if n.parent]


def build_parent_dict_from_candidates_scaffold(candidate_positions):
    """Build parent-child dictionary from candidate positions."""
    parent_dict = {}
    
    for pos in candidate_positions:
        ptype = pos["placement_type"]
        anchor = pos["anchor"]
        
        if ptype == "on_edge":
            child = pos["meta"]["child"]
            parent_dict[child] = anchor
        
        elif ptype == "new_parent_merge":
            merge_children = pos["meta"]["merge_children"]
            for child in merge_children:
                parent_dict[child] = anchor
    
    return parent_dict


def build_lineage_parent_dict_from_tree(tree, anchor):
    """
    Build lineage parent-child dictionary from the current tree from anchor to ROOT.
    
    Returns format:
    {'chr17_7389869_G_A': 'chr19_18929325_C_T', 'chr7_142870944_G_A': 'chr15_74920412_G_A', ...}
    Excludes ROOT.
    """
    lineage_parent_dict = {}
    
    def find_path(node, target, current_path):
        if node.name == target:
            # Build parent-child relationships on the path
            for i in range(len(current_path)-1, 0, -1):
                lineage_parent_dict[current_path[i]] = current_path[i-1]
            if current_path:  # Ensure anchor itself is in the relationship
                lineage_parent_dict[target] = current_path[-1] if current_path else None
            return True
        
        for child in node.children:
            if find_path(child, target, current_path + [node.name]):
                return True
        return False
    
    # Start searching from ROOT
    find_path(tree, anchor, [])
    
    # Filter out entries containing ROOT
    filtered_dict = {child: parent for child, parent in lineage_parent_dict.items() 
                    if child != "ROOT" and parent != "ROOT"}
    
    return filtered_dict


def find_all_intersect_muts_from_tree_by_matrix_scaffold(tree, matrix, target_mut, min_overlap=1, logger_obj=None):
    """
    Return all tree nodes that intersect with target_mut in the matrix
    with at least min_overlap cells co-occurring.
    """
    # Use provided logger or fall back to global logger
    log = logger_obj if logger_obj is not None else logging.getLogger(__name__)
    
    intersect_nodes = set()
    
    # Traverse all nodes (except ROOT)
    for node in tree.all_nodes():
        if node.name == "ROOT":
            continue
        
        has_intersection = False
        
        # Check if any mutation in the node intersects with the target mutation
        for mut in node.name.split("|"):
            if mut == target_mut:
                continue  # Skip self
            
            if mut not in matrix.columns or target_mut not in matrix.columns:
                continue  # Ensure mutations exist in the matrix
            
            # Calculate co-occurrence count
            mut_vec = matrix[mut].fillna(0)
            target_vec = matrix[target_mut].fillna(0)
            N11 = ((mut_vec == 1) & (target_vec == 1)).sum()
            
            # If any mutation has sufficient co-occurrence, mark this node
            if N11 >= min_overlap:
                has_intersection = True
                break  # One intersecting mutation is enough
        
        # If any mutation in the node intersects with the target, add the entire node
        if has_intersection:
            intersect_nodes.add(node.name)
    
    return intersect_nodes


def find_all_path_nodes_scaffold(intersection_nodes, tree_parent_dict):
    """
    Find all nodes on the paths connecting intersection nodes.
    
    Parameters
    ----------
    intersection_nodes : set or list
        Set/list of intersection node names
    tree_parent_dict : dict
        Parent-child dictionary (tree_parent_dict)
    
    Returns
    -------
    list
        Sorted list of all node names on paths between intersection nodes
    """
    all_path_nodes = set()
    all_path_nodes.add('ROOT')  # Always include ROOT
    
    # For each intersection node, find path from ROOT to that node
    for node in intersection_nodes:
        path_to_root = get_path_to_root_scaffold(node, tree_parent_dict)
        all_path_nodes.update(path_to_root)
    
    # Find paths connecting different intersection nodes
    intersection_list = list(intersection_nodes)
    for i in range(len(intersection_list)):
        for j in range(i + 1, len(intersection_list)):
            path_between = get_path_between_nodes_scaffold(intersection_list[i], intersection_list[j], tree_parent_dict)
            all_path_nodes.update(path_between)
    
    # ============================================================
    # CRITICAL FIX: Return sorted list to ensure deterministic iteration order
    # ============================================================
    # In parallel execution, set iteration order can vary across processes
    # due to different memory layouts or hash randomization. Sorting ensures
    # that all workers iterate over nodes in the exact same order.
    # ============================================================
    return sorted(all_path_nodes)


def get_path_to_root_scaffold(node, tree_parent_dict):
    """Find path from node to ROOT."""
    path = []
    current = node
    while current in tree_parent_dict:
        path.append(current)
        current = tree_parent_dict[current]
        if current == 'ROOT':
            path.append('ROOT')
            break
    return path


def get_path_between_nodes_scaffold(node1, node2, tree_parent_dict):
    """Find path between two nodes."""
    # Find path from node1 to ROOT
    path1 = get_path_to_root_scaffold(node1, tree_parent_dict)
    # Find path from node2 to ROOT
    path2 = get_path_to_root_scaffold(node2, tree_parent_dict)
    
    # Reverse paths to start from ROOT
    path1_from_root = list(reversed(path1))
    path2_from_root = list(reversed(path2))
    
    # Find Lowest Common Ancestor (LCA)
    lca = None
    for i in range(min(len(path1_from_root), len(path2_from_root))):
        if path1_from_root[i] == path2_from_root[i]:
            lca = path1_from_root[i]
        else:
            break
    
    if lca is None:
        return []
    
    # Build full path: node1 -> LCA -> node2
    lca_index1 = path1_from_root.index(lca)
    lca_index2 = path2_from_root.index(lca)
    
    path_node1_to_lca = path1_from_root[lca_index1:]
    path_lca_to_node2 = path2_from_root[lca_index2:]
    
    # Merge paths, removing duplicate LCA
    full_path = path_node1_to_lca + path_lca_to_node2[1:]
    return full_path




# -------------------------
# penalty function
# Expected genotype G generation: returns G (cells -> 0/1) for placement type
# placement_type: 'on_node' | 'as_child' | 'on_edge' | 'new_parent_merge' | 'new_leaf'
# Simplified strategy (can be replaced with stricter Bayesian inference):
#   - on_node: expected to be cells where any mutation in node subtree exists (using posterior > 0.5)
#   - as_child: expected to be subset of node subtree (use node subtree cells as superset, adjusted by posterior weights)
#   - on_edge: same as on_node (edge placement is treated as splitting child subset)
#   - new_parent_merge: union of children's cell sets
#   - new_leaf: subset of parent subtree (approx: empty set initially)
# -------------------------

def select_overlapping_with_minimum(results, logger_obj=None):
    """
    Select all positions that overlap with the minimum penalty candidates.
    
    Parameters
    ----------
    results : list
        List of result dictionaries containing penalty information
    logger_obj : logging.Logger, optional
        Logger instance for logging messages. If None, uses global logger.
    
    Returns
    -------
    list
        Selected positions
    """
    # Use provided logger or fall back to global logger
    log = logger_obj if logger_obj is not None else logging.getLogger(__name__)
    
    if not results:
        return []
    
    # Find the minimum penalty_min
    min_penalty = min(r['penalty_min'] for r in results)
    
    # Find all candidates with penalty_min equal to the minimum
    min_candidates = [r for r in results if r['penalty_min'] == min_penalty]
    
    # Calculate the maximum penalty_max among these minimum candidates
    min_group_max = max(r['penalty_max'] for r in min_candidates)
    
    # Select all positions that overlap with the minimum candidate group
    selected = []
    for result in results:
        # Select if candidate's penalty_min <= max penalty_max of the minimum group
        if result['penalty_min'] <= min_group_max:
            selected.append(result['position'])
    
    log.debug(f"Selected {len(selected)} positions out of {len(results)} candidates")
    
    return selected


def compute_binary_penalty_for_positions(new_mut, refined_positions, M_current, I_selected, logger_obj=None):
    """
    Compute binary penalty bounds for new_mut across all candidate positions.
    
    Parameters
    ----------
    new_mut : str
        Mutation to be added to the tree
    refined_positions : list of dict
        List of candidate positions, each containing placement_type, anchor, nodes, edges, etc.
    M_current : pd.DataFrame
        Backbone clone matrix (0/1)
    I_selected : pd.DataFrame
        Observed genotype matrix (contains NA)
    logger_obj : logging.Logger, optional
        Logger instance for logging messages. If None, uses global logger.
    
    Returns
    -------
    selected_positions : list
        Selected optimal positions
    df_penalty : pd.DataFrame
        Penalty min/max for each position
    """
    # Use provided logger or fall back to global logger
    log = logger_obj if logger_obj is not None else logging.getLogger(__name__)
    
    def compute_N_vectors(imputed_vec, observed_vec):
        """Calculate N_10 and N_01 counts."""
        imputed_vec = imputed_vec.astype(int)
        obs_vec = observed_vec.copy()
        obs_vec = obs_vec.fillna(-1).astype(int)
        
        N_10 = np.sum((obs_vec == 1) & (imputed_vec == 0))
        N_01 = np.sum((obs_vec == 0) & (imputed_vec == 1))
        return N_10, N_01
    
    log.debug(f"Computing binary penalty for mutation: {new_mut}")
    log.debug(f"Number of candidate positions: {len(refined_positions)}")
    
    results = []
    new_mut_bin_vector = I_selected[new_mut].replace({pd.NA: np.nan}).fillna(0).astype(int)
    
    for idx, pos in enumerate(refined_positions):
        placement_type = pos['placement_type']
        anchor = pos['anchor']
        
        # Generate imputed vector
        imputed_vec = pd.Series(0, index=M_current.index)  # Default all zeros
        
        if placement_type == 'on_node':
            # new_mut matches the cells of the anchor node
            if anchor not in M_current.columns:
                raise ValueError(f"Anchor {anchor} not in M_current")
            imputed_vec = M_current[anchor]
        
        elif placement_type == 'new_leaf':
            # Get parent and child information
            parent = anchor
            if parent == 'ROOT':
                vec_parent = pd.Series(1, index=M_current.index)  # Root is 1
            else:
                if parent not in M_current.columns:
                    raise ValueError(f"Parent {parent} not in M_current")
                vec_parent = M_current[parent]
            
            # New mutation intersection: ensure new_mut is only 1 under the parent
            imputed_vec = new_mut_bin_vector & vec_parent
            
            # Get child nodes
            child_nodes = [node['name'] for node in pos['nodes'] if node['parent'] == parent]
            child_nodes = [i for i in child_nodes if i != new_mut]
            
            # Ensure new_mut is mutually exclusive with children
            for child in child_nodes:
                if child not in M_current.columns:
                    raise ValueError(f"Child {child} not in M_current")
                vec_child = M_current[child]
                imputed_vec &= ~vec_child  # Mutually exclusive with children
            
            # Ensure that only where new_mut is 1, the parent chain is 1
            imputed_vec &= new_mut_bin_vector
            
        elif placement_type == 'on_edge':
            parent = anchor
            child = pos['meta']['child']
            
            # Process parent vector
            if parent == 'ROOT':
                vec_parent = pd.Series(1, index=M_current.index)
            else:
                if parent not in M_current.columns:
                    raise ValueError(f"Parent {anchor} not in M_current")
                vec_parent = M_current[anchor]
            
            # Process target child vector
            if child not in M_current.columns:
                raise ValueError(f"Child {child} not in M_current")
            
            vec_child = M_current[child]
            
            # Get all other children (siblings) under parent
            sibling_nodes = [node['name'] for node in pos['nodes'] 
                            if node['parent'] == parent and node['name'] != child]
            sibling_nodes = [i for i in sibling_nodes if i != new_mut]
            
            # Initialize imputed_vec: include target child, within parent range
            imputed_vec = vec_child.copy()
            
            # Add new mutation within parent range
            imputed_vec = (imputed_vec | (new_mut_bin_vector & vec_parent)).astype(int)
            
            # Key: exclude all sibling node cells (mutually exclusive)
            for sibling in sibling_nodes:
                if sibling not in M_current.columns:
                    raise ValueError(f"Sibling {sibling} not in M_current")
                vec_sibling = M_current[sibling]
                imputed_vec = (imputed_vec & ~vec_sibling).astype(int)
            
            # Final ensure within parent range
            imputed_vec = (imputed_vec & vec_parent).astype(int)
        
        elif placement_type == 'new_parent_merge':
            parent = anchor
            merge_children = pos['meta']['merge_children']
            
            # Process parent vector
            if parent == 'ROOT':
                vec_parent = pd.Series(1, index=M_current.index)
            else:
                if parent not in M_current.columns:
                    raise ValueError(f"Parent {parent} not in M_current for edge placement")
                vec_parent = M_current[parent]
            
            # Process merged children vectors
            vec_children = pd.Series(0, index=M_current.index)
            for c in merge_children:
                if c not in M_current.columns:
                    raise ValueError(f"Child {c} not in M_current for edge placement")
                vec_children = ((vec_children == 1) | (M_current[c] == 1)).astype(int)
            
            # Initialize imputed_vec: include merged children, within parent range
            imputed_vec = vec_children.copy()
            
            # Add new mutation within parent range
            imputed_vec = (imputed_vec | (new_mut_bin_vector & vec_parent)).astype(int)
            
            # Key: exclude all other sibling node cells (mutually exclusive)
            # Get all other children (siblings) under parent, excluding merged children
            sibling_nodes = [node['name'] for node in pos['nodes'] 
                            if node['parent'] == parent and node['name'] not in merge_children]
            sibling_nodes = [i for i in sibling_nodes if i != new_mut]
            
            for sibling in sibling_nodes:
                if sibling not in M_current.columns:
                    raise ValueError(f"Sibling {sibling} not in M_current")
                vec_sibling = M_current[sibling]
                imputed_vec = (imputed_vec & ~vec_sibling).astype(int)
            
            # Final ensure within parent range
            imputed_vec = (imputed_vec & vec_parent).astype(int)
        
        else:
            raise ValueError(f"Unknown placement_type: {placement_type}")
        
        # Corresponding observed vector
        if new_mut not in I_selected.columns:
            observed_vec = pd.Series(np.nan, index=M_current.index)  # No observed data, all NA
        else:
            observed_vec = I_selected[new_mut]
        
        N_10, N_01 = compute_N_vectors(imputed_vec, observed_vec)
        
        # Calculate penalty
        penalty_min = 0.5 * N_10 + 0.05 * N_01
        penalty_max = 1.0 * N_10 + 0.1 * N_01
        
        results.append({
            'position_index': idx,
            'placement_type': placement_type,
            'anchor': anchor,
            'N_10': N_10,
            'N_01': N_01,
            'penalty_min': penalty_min,
            'penalty_max': penalty_max,
            'position': pos
        })
    
    df_penalty = pd.DataFrame(results)
    
    # Select all potentially optimal positions
    # Filtering logic: select candidates with low lower bound and upper bound not exceeding other candidates' lower bounds
    selected_positions = select_overlapping_with_minimum(results, logger_obj=log)
    
    log.debug(f"Selected {len(selected_positions)} optimal positions out of {len(results)} candidates")
    
    return selected_positions, df_penalty

# selected_positions, df_penalty = compute_binary_penalty_for_positions('chr17_7389869_G_A', refined_positions, M_current, I_selected)
# df_penalty[['position_index', 'placement_type', 'anchor', 'N_10', 'N_01', 'penalty_min', 'penalty_max']]
# #    position_index placement_type             anchor  N_10  N_01  penalty_min  penalty_max
# # 0               0        on_node  chr1_39034563_T_A     0    25         1.25          2.5
# # 1               1       new_leaf  chr1_39034563_T_A     0     0         0.00          0.0


##### Bayesian penalty
def compute_bayesian_penalty_each_pos(input_binary_vec, posterior_vec, imputed_vec, fnfp_ratio=0.1, omega_NA=0.001, logger_obj=None):
    """
    Bayesian penalty calculation based on data prior distribution.
    
    Parameters
    ----------
    input_binary_vec : array-like
        Observed binary vector (0/1/NA)
    posterior_vec : array-like
        Posterior probability vector
    imputed_vec : array-like
        Imputed genotype vector (0/1)
    fnfp_ratio : float, default=0.1
        False negative to false positive ratio
    omega_NA : float, default=0.001
        NA penalty weight
    logger_obj : logging.Logger, optional
        Logger instance for logging messages. If None, uses global logger.
    
    Returns
    -------
    tuple
        (total_penalty, weight_na_to_1, weight_na_to_0)
    """
    # Use provided logger or fall back to global logger
    log = logger_obj if logger_obj is not None else logging.getLogger(__name__)
    
    # 1. Calculate prior probabilities
    non_nan_values = input_binary_vec[~np.isnan(input_binary_vec)]
    count_1 = np.sum(non_nan_values == 1)
    count_0 = np.sum(non_nan_values == 0)
    total_count = count_1 + count_0
    
    if total_count > 0:
        p_mutation = count_1 / total_count      # Mutation prior probability
        p_wildtype = count_0 / total_count      # Wildtype prior probability
    else:
        p_mutation = p_wildtype = 0.5
    
    # 2. Set weights based on prior probabilities
    # NA→1 penalty: proportional to wildtype probability (more wildtype = more unreasonable to impute 1)
    weight_na_to_1 = omega_NA * p_wildtype  
    # NA→0 penalty: proportional to mutation probability (more mutations = more unreasonable to impute 0)
    weight_na_to_0 = omega_NA * p_mutation
    
    # 3. Count various cases
    fp_count = np.sum((input_binary_vec == 1) & (imputed_vec == 0))
    fn_count = np.sum((input_binary_vec == 0) & (imputed_vec == 1))
    na_to_1_count = np.sum(np.isnan(input_binary_vec) & (imputed_vec == 1))
    na_to_0_count = np.sum(np.isnan(input_binary_vec) & (imputed_vec == 0))
    
    # 4. Calculate penalties
    fp_penalty = fp_count
    fn_penalty = fn_count * fnfp_ratio
    na_to_1_penalty = na_to_1_count * weight_na_to_1
    na_to_0_penalty = na_to_0_count * weight_na_to_0
    
    # 5. Enhanced FP penalty calculation
    total_penalty = fp_penalty + fn_penalty + na_to_1_penalty + na_to_0_penalty
    
    return total_penalty, weight_na_to_1, weight_na_to_0


def compute_bayesian_penalty_each_chain_mut_by_pos(input_binary_vec, posterior_vec, imputed_vec, fnfp_ratio, weight_na_to_1, weight_na_to_0, logger_obj=None):
    """
    Bayesian penalty calculation for chain mutations by position.
    
    Parameters
    ----------
    input_binary_vec : array-like
        Observed binary vector (0/1/NA)
    posterior_vec : array-like
        Posterior probability vector
    imputed_vec : array-like
        Imputed genotype vector (0/1)
    fnfp_ratio : float
        False negative to false positive ratio
    weight_na_to_1 : float
        NA to 1 penalty weight
    weight_na_to_0 : float
        NA to 0 penalty weight
    logger_obj : logging.Logger, optional
        Logger instance for logging messages. If None, uses global logger.
    
    Returns
    -------
    float
        Total penalty
    """
    # Use provided logger or fall back to global logger
    log = logger_obj if logger_obj is not None else logging.getLogger(__name__)
    
    # 3. Count various cases
    fp_count = np.sum((input_binary_vec == 1) & (imputed_vec == 0))
    fn_count = np.sum((input_binary_vec == 0) & (imputed_vec == 1))
    na_to_1_count = np.sum(np.isnan(input_binary_vec) & (imputed_vec == 1))
    na_to_0_count = np.sum(np.isnan(input_binary_vec) & (imputed_vec == 0))
    
    # 4. Calculate penalties
    fp_penalty = fp_count
    fn_penalty = fn_count * fnfp_ratio
    na_to_1_penalty = na_to_1_count * weight_na_to_1
    na_to_0_penalty = na_to_0_count * weight_na_to_0
    
    total_penalty = fp_penalty + fn_penalty + na_to_1_penalty + na_to_0_penalty
    
    return total_penalty


# ============================================================
# Compute penalties for all candidate positions
# ============================================================
def compute_bayesian_penalty_for_all_positions_scaffold(
    new_mut, selected_positions, T_current, M_current, I_selected, P_selected, parent_dict, intersection_nodes, 
    ω_NA=0.001, fnfp_ratio=0.1, φ=1, logger_obj=None
):
    """
    Compute penalties for all candidate positions without modifying M_current.
    
    Optimization strategies:
    1. Cache all reusable computation results (conflict masks, mutation chains, etc.)
    2. Use NumPy vectorization instead of loops
    3. Automatically decide whether to use parallel execution based on number of positions:
       - <= 20: serial (avoid parallel overhead)
       - > 20: parallel (utilize multi-core)
    
    Parameters
    ----------
    new_mut : str
        New mutation to be processed
    selected_positions : list
        List of candidate positions, each containing placement_type, anchor, etc.
    T_current : dict
        Current phylogenetic tree structure
    M_current : pandas.DataFrame
        Current mutation matrix, index: cells, columns: mutations
    I_selected : pandas.DataFrame
        Mutation presence matrix, index: cells, columns: mutations
    P_selected : pandas.DataFrame
        Posterior probability matrix, index: cells, columns: mutations
    parent_dict : dict
        Parent node dictionary for building mutation chains
    intersection_nodes : list
        List of intersection nodes
    ω_NA : float, default=0.001
        NA weight parameter
    fnfp_ratio : float, default=0.1
        False negative to false positive ratio
    φ : float, default=1
        Bayesian penalty parameter
    logger_obj : logging.Logger, optional
        Logger instance for logging messages. If None, uses global logger.
    
    Returns
    -------
    df_penalty : pandas.DataFrame
        DataFrame containing penalties for all candidate positions
        Each row corresponds to one candidate position
        Includes detailed penalty components and imputed vector
    """
    import pandas as pd
    import numpy as np
    import logging
    from joblib import Parallel, delayed
    
    # Use provided logger or fall back to global logger
    log = logger_obj if logger_obj is not None else logging.getLogger(__name__)
    
    # If no candidate positions, return empty DataFrame
    if len(selected_positions) == 0:
        log.warning(f"No selected positions for mutation {new_mut}")
        return pd.DataFrame()
    
    # ============================================================
    # 1. Pre-compute common data
    # ============================================================
    
    matrix_index = M_current.index
    matrix_columns = M_current.columns.tolist()
    
    # Get new_mut presence vector
    new_mut_bin_vector = I_selected[new_mut].replace({pd.NA: np.nan}).fillna(0).astype(int)
    new_mut_mask = new_mut_bin_vector.reindex(matrix_index, fill_value=0).to_numpy(dtype=bool)
    
    # Calculate mutation statistics
    input_binary_vec_full = I_selected[new_mut].replace({pd.NA: np.nan})
    na_ratio = input_binary_vec_full.isna().mean()
    mut_ratio = input_binary_vec_full.fillna(0).mean()
    N_nodes_beforeT = len(T_current.all_nodes())
    
    # ============================================================
    # 2. Build caching system
    # ============================================================
    
    # Cache M_current columns
    m_current_cache = {}
    for col in matrix_columns:
        if col != 'ROOT':
            m_current_cache[col] = M_current[col].reindex(matrix_index, fill_value=0).to_numpy(dtype=np.int8)
    
    # Cache I_selected columns
    i_selected_cache = {}
    for col in I_selected.columns:
        i_selected_cache[col] = I_selected[col].reindex(matrix_index, fill_value=np.nan).to_numpy(dtype=np.float64)
    
    # Cache conflict masks
    conflict_cache = {}
    
    def get_conflict_mask(anchor, sibling_nodes, exclude_nodes=None):
        """Get conflict mask with caching"""
        if exclude_nodes is None:
            exclude_nodes = []
        
        cache_key = (anchor, tuple(sorted(sibling_nodes)), tuple(sorted(exclude_nodes)))
        
        if cache_key in conflict_cache:
            return conflict_cache[cache_key]
        
        lineage_parent = build_lineage_parent_dict_from_tree(T_current, anchor)
        lineage_conflict_nodes = get_all_conflict_nodes_outside_lineage_scaffold(
            anchor, lineage_parent, matrix_columns, exclude_nodes=exclude_nodes
        )
        
        all_conflict_nodes = list(set(sibling_nodes + lineage_conflict_nodes))
        
        if not all_conflict_nodes:
            conflict_mask = np.zeros(len(matrix_index), dtype=bool)
        else:
            conflict_mask = np.zeros(len(matrix_index), dtype=bool)
            for node in all_conflict_nodes:
                if node in m_current_cache:
                    conflict_mask = conflict_mask | (m_current_cache[node] == 1)
                elif node in matrix_columns:
                    node_array = M_current[node].reindex(matrix_index, fill_value=0).to_numpy(dtype=np.int8)
                    conflict_mask = conflict_mask | (node_array == 1)
        
        conflict_cache[cache_key] = conflict_mask
        return conflict_mask
    
    # Cache mutation chains
    chain_cache = {}
    def get_full_chain(anchor):
        if anchor not in chain_cache:
            chain_cache[anchor] = get_full_mutnode_chain_with_anchor_scaffold(anchor, parent_dict)
        return chain_cache[anchor]
    
    # Cache node mutations
    node_mutations_cache = {}
    def get_node_mutations(node_name):
        if node_name not in node_mutations_cache:
            node_mutations_cache[node_name] = node_name.split("|")
        return node_mutations_cache[node_name]
    
    # Cache children
    children_cache = {}
    def get_children(node):
        if node not in children_cache:
            children_cache[node] = find_children_of_node_scaffold(node, matrix_columns, parent_dict)
        return children_cache[node]
    
    # ============================================================
    # 3. Single position computation function (core logic)
    # ============================================================
    
    def compute_single_position(idx, pos):
        """Compute penalty for a single candidate position"""
        placement_type = pos['placement_type']
        anchor = pos['anchor']
        
        imputed_array = np.zeros(len(matrix_index), dtype=np.int8)
        merge_penalty = 0
        N_nodes = N_nodes_beforeT
        
        # Calculate imputed vector based on placement type
        if placement_type == 'on_node':
            sibling_nodes = [n['name'] for n in pos['nodes'] if n['parent'] == anchor and n['name'] != new_mut]
            conflict_mask = get_conflict_mask(anchor, sibling_nodes)
            new_mut_cleaned = new_mut_mask & ~conflict_mask
            
            if anchor in m_current_cache:
                anchor_array = m_current_cache[anchor]
            else:
                anchor_array = M_current[anchor].reindex(matrix_index, fill_value=0).to_numpy(dtype=np.int8)
            
            imputed_array = ((anchor_array == 1) | new_mut_cleaned).astype(np.int8)
            N_nodes = N_nodes_beforeT + 1
            
        elif placement_type == 'new_leaf':
            sibling_nodes = [n['name'] for n in pos['nodes'] if n['parent'] == anchor and n['name'] != new_mut]
            conflict_mask = get_conflict_mask(anchor, sibling_nodes)
            new_mut_cleaned = new_mut_mask & ~conflict_mask
            imputed_array = new_mut_cleaned.astype(np.int8)
            N_nodes = N_nodes_beforeT + 2
            
        elif placement_type == 'on_edge':
            child = pos['meta']['child']
            sibling_nodes = [n['name'] for n in pos['nodes'] if n['parent'] == anchor and n['name'] not in [child, new_mut]]
            conflict_mask = get_conflict_mask(anchor, sibling_nodes)
            new_mut_cleaned = new_mut_mask & ~conflict_mask
            
            if child in m_current_cache:
                child_array = m_current_cache[child]
            else:
                child_array = M_current[child].reindex(matrix_index, fill_value=0).to_numpy(dtype=np.int8)
            
            imputed_array = ((child_array == 1) | new_mut_cleaned).astype(np.int8)
            N_nodes = N_nodes_beforeT + 2
            
        elif placement_type == 'new_parent_merge':
            merge_children = pos['meta']['merge_children']
            sibling_nodes = [n['name'] for n in pos['nodes'] if n['parent'] == anchor and n['name'] not in merge_children + [new_mut]]
            conflict_mask = get_conflict_mask(anchor, sibling_nodes, merge_children)
            new_mut_cleaned = new_mut_mask & ~conflict_mask
            
            children_union = np.zeros(len(matrix_index), dtype=bool)
            for c in merge_children:
                if c in m_current_cache:
                    children_union = children_union | (m_current_cache[c] == 1)
                else:
                    child_array = M_current[c].reindex(matrix_index, fill_value=0).to_numpy(dtype=np.int8)
                    children_union = children_union | (child_array == 1)
            
            imputed_array = (children_union | new_mut_cleaned).astype(np.int8)
            merge_penalty = np.log(len(merge_children)) * 0.5
            N_nodes = N_nodes_beforeT + 2
            
        else:
            raise ValueError(f"Unknown placement_type: {placement_type}")
        
        imputed_vec = pd.Series(imputed_array, index=matrix_index)
        imputed_bool = imputed_array.astype(bool)
        
        # Calculate penalty
        full_mutnode_chain = get_full_chain(anchor)
        
        posterior_vec = P_selected[new_mut]
        input_binary_vec = I_selected[new_mut]
        
        new_mut_penalty, actual_na_flip_ratio, refined_ω_NA, φ_adjusted, weight_na_to_1, weight_na_to_0 = compute_dynamic_penalty_scaffold(
            input_binary_vec, posterior_vec, imputed_vec, fnfp_ratio, ω_NA, φ,
            na_ratio, mut_ratio, placement_type, N_nodes, logger_obj=log
        )
        
        # Calculate chain penalty
        chain_penalty = 0
        chain_mutations_count = 0
        
        for node in full_mutnode_chain:
            if node == 'ROOT':
                continue
            
            mutations_on_node = get_node_mutations(node)
            scored_mutations = [m for m in mutations_on_node if m != new_mut]
            chain_mutations_count += len(scored_mutations)
            
            if not scored_mutations:
                continue
            
            if node in m_current_cache:
                node_array = m_current_cache[node]
            else:
                node_array = M_current[node].reindex(matrix_index, fill_value=0).to_numpy(dtype=np.float64)
            
            cells_to_flip_mask = imputed_bool & ((node_array == 0) | np.isnan(node_array))
            
            if not cells_to_flip_mask.any():
                continue
            
            flip_indices = np.where(cells_to_flip_mask)[0]
            imputed_ones = np.ones(len(flip_indices), dtype=np.int8)
            
            for mutation in scored_mutations:
                if mutation in i_selected_cache:
                    mut_input_array = i_selected_cache[mutation]
                else:
                    mut_input_array = I_selected[mutation].reindex(matrix_index, fill_value=np.nan).to_numpy(dtype=np.float64)
                
                mut_input_subset = mut_input_array[flip_indices]
                
                mut_penalty = compute_bayesian_penalty_each_chain_mut_by_pos(
                    mut_input_subset,
                    None,
                    imputed_ones,
                    fnfp_ratio,
                    weight_na_to_1,
                    weight_na_to_0,
                    logger_obj=log
                )
                chain_penalty += mut_penalty
        
        total_chain_penalty = new_mut_penalty + chain_penalty
        
        log_N_nodes_penalty = np.log(N_nodes)
        BIC_penalty = φ_adjusted * np.log(N_nodes)
        
        root_penalty = 0
        if anchor == 'ROOT':
            root_penalty = np.log(N_nodes) * 0.5
        
        base_total_penalty = total_chain_penalty + log_N_nodes_penalty + BIC_penalty + merge_penalty + root_penalty
        
        intersection_penalty = compute_intersection_based_penalty_scaffold(
            new_mut, pos, intersection_nodes, M_current, I_selected, na_ratio, mut_ratio, actual_na_flip_ratio, logger_obj=log
        )
        
        hierarchy_penalty = compute_hierarchy_penalty_scaffold(
            new_mut, pos, M_current, I_selected, parent_dict, na_ratio, mut_ratio, actual_na_flip_ratio, logger_obj=log
        )
        
        total_penalty = base_total_penalty + intersection_penalty + hierarchy_penalty
        
        return {
            'position_index': idx,
            'placement_type': placement_type,
            'anchor': anchor,
            'new_mut_penalty': new_mut_penalty,
            'chain_penalty': chain_penalty,
            'total_chain_penalty': total_chain_penalty,
            'N_nodes': N_nodes,
            'BIC_penalty': BIC_penalty,
            'log_N_nodes_penalty': log_N_nodes_penalty,
            'merge_penalty': merge_penalty,
            'root_penalty': root_penalty,
            'base_total_penalty': base_total_penalty,
            'intersection_penalty': intersection_penalty,
            'hierarchy_penalty': hierarchy_penalty,
            'total_penalty': total_penalty,
            'position': pos,
            'imputed_vec': imputed_vec,
            'na_ratio': na_ratio,
            'mut_ratio': mut_ratio,
            'actual_na_flip_ratio': actual_na_flip_ratio,
            'chain_mutations_count': chain_mutations_count,
            'weight_na_to_1': weight_na_to_1,
            'weight_na_to_0': weight_na_to_0,
            'refined_ω_NA': refined_ω_NA,
            'φ_adjusted': φ_adjusted,
            'base_ω_NA': ω_NA,
            'base_φ': φ,
            'full_mutnode_chain': full_mutnode_chain
        }
    
    # ============================================================
    # 4. Automatic execution mode selection
    # ============================================================
    
    n_positions = len(selected_positions)
    
    # Auto-threshold: based on empirically optimized values
    # - <= 20: serial (parallel overhead > benefit)
    # - > 20: parallel (utilize multi-core)
    if n_positions <= 20:
        # Serial execution (with cache optimization)
        results = [compute_single_position(idx, pos) for idx, pos in enumerate(selected_positions)]
    else:
        # Parallel execution (with cache optimization)
        n_jobs = -1  # Use all available CPU cores
        log.info(f"Automatically using parallel mode for {n_positions} positions")
        results = Parallel(n_jobs=n_jobs, verbose=10 if n_positions > 100 else 0)(
            delayed(compute_single_position)(idx, pos)
            for idx, pos in enumerate(selected_positions)
        )
    
    # ============================================================
    # 5. Return results
    # ============================================================
    
    df_penalty = pd.DataFrame(results)
    return df_penalty


# ============================================================
# Apply specific position to tree and matrix
# ============================================================

def apply_position_to_tree_scaffold(
    new_mut, position, imputed_vec, T_current, M_current, I_selected, parent_dict, logger_obj=None
):
    """
    Apply the specific position to trees and matrices.
    
    Parameters
    ----------
    new_mut : str
        New mutation to add
    position : dict
        Position dictionary containing placement_type, anchor, etc.
    imputed_vec : pd.Series
        Imputed genotype vector
    T_current : TreeNode
        Current tree
    M_current : pd.DataFrame
        Current mutation matrix
    I_selected : pd.DataFrame
        Selected mutation matrix
    parent_dict : dict
        Parent node dictionary
    logger_obj : logging.Logger, optional
        Logger instance for logging messages. If None, uses global logger.
    
    Returns
    -------
    tuple
        (T_updated, M_updated)
    """
    # Use provided logger or fall back to global logger
    log = logger_obj if logger_obj is not None else logging.getLogger(__name__)
    
    M_updated = M_current.copy()
    T_updated = T_current.copy()
    
    anchor = position['anchor']
    placement_type = position['placement_type']
    
    # Get full mutation chain from anchor to ROOT
    full_mutnode_chain = get_full_mutnode_chain_with_anchor_scaffold(anchor, parent_dict)
    
    # Ensure imputed_vec has the same index as M_updated
    if not imputed_vec.index.equals(M_updated.index):
        imputed_vec = imputed_vec.reindex(M_updated.index, fill_value=0)
    
    # Get cells where imputed_vec = 1
    cells_with_final_one = imputed_vec[imputed_vec == 1].index.tolist()
    if len(cells_with_final_one) > 0:
        # Propagate mutations along the chain
        for cell in cells_with_final_one:
            for mutation in full_mutnode_chain:
                if M_updated.loc[cell, mutation] == 0:
                    M_updated.loc[cell, mutation] = 1
    
    # Add new mutation column to M_updated (consider 'on_node')
    if placement_type == 'on_node':
        new_colname = anchor + "|" + new_mut
        M_updated.drop(anchor, axis=1, inplace=True)
        M_updated[new_colname] = imputed_vec
    else:
        M_updated[new_mut] = imputed_vec
    
    # Update tree structure
    T_updated = add_new_mutation_to_tree_independent(new_mut, T_updated, position)
    
    log.debug(f"Applied position: {placement_type} on anchor {anchor} for mutation {new_mut}")
    
    return T_updated, M_updated


def get_all_conflict_nodes_outside_lineage_scaffold(anchor, parent_dict, all_columns, exclude_nodes=None, logger_obj=None):
    """
    Get all possible conflict nodes outside the lineage.
    
    Parameters
    ----------
    anchor : str
        Current anchor node
    parent_dict : dict
        Parent node dictionary
    all_columns : list
        All mutation nodes
    exclude_nodes : list, optional
        Nodes to exclude (e.g., merge_children)
    logger_obj : logging.Logger, optional
        Logger instance for logging messages. If None, uses global logger.
    
    Returns
    -------
    list
        Conflict nodes outside the lineage
    """
    # Use provided logger or fall back to global logger
    log = logger_obj if logger_obj is not None else logging.getLogger(__name__)
    
    if exclude_nodes is None:
        exclude_nodes = []
    
    # Get all nodes in the current lineage (path from anchor to ROOT)
    lineage_nodes = set()
    current = anchor
    while current is not None:
        lineage_nodes.add(current)
        parent = parent_dict.get(current, None)
        if parent is None or parent == 'ROOT':
            break
        current = parent
    lineage_nodes.add('ROOT')
    
    # Get all nodes not in the current lineage
    conflict_nodes = []
    for node in all_columns:
        if node == 'ROOT':
            continue
        if node not in lineage_nodes and node not in exclude_nodes:
            conflict_nodes.append(node)
    
    log.debug(f"Found {len(conflict_nodes)} conflict nodes outside lineage for anchor {anchor}")
    
    return conflict_nodes


def get_full_mutnode_chain_with_anchor_scaffold(anchor, parent_dict, logger_obj=None):
    """
    Get the full mutation chain from anchor to ROOT (including anchor itself).
    
    Parameters
    ----------
    anchor : str
        Anchor node
    parent_dict : dict
        Parent node dictionary
    logger_obj : logging.Logger, optional
        Logger instance for logging messages. If None, uses global logger.
    
    Returns
    -------
    list
        Mutation chain from anchor to ROOT
    """
    # Use provided logger or fall back to global logger
    log = logger_obj if logger_obj is not None else logging.getLogger(__name__)
    
    mutation_chain = [anchor]
    current = anchor
    
    while True:
        parent = parent_dict.get(current, None)
        if parent is None or parent == 'ROOT':
            break
        mutation_chain.append(parent)
        current = parent
    
    log.debug(f"Mutation chain for anchor {anchor}: {mutation_chain}")
    
    return mutation_chain


def compute_dynamic_penalty_scaffold(input_binary_vec, posterior_vec, imputed_vec, fnfp_ratio, base_ω_NA, base_φ, na_ratio, mut_ratio, placement_type, N_nodes, logger_obj=None):
    """
    Dynamically adjust all penalty weights based on actual NA→1 flip ratio.
    
    Parameters
    ----------
    input_binary_vec : pd.Series
        Observed binary vector (0/1/NA)
    posterior_vec : pd.Series
        Posterior probability vector
    imputed_vec : pd.Series
        Imputed genotype vector
    fnfp_ratio : float
        False negative to false positive ratio
    base_ω_NA : float
        Base NA penalty weight
    base_φ : float
        Base BIC penalty parameter
    na_ratio : float
        NA ratio in data
    mut_ratio : float
        Mutation ratio in data
    placement_type : str
        Type of placement
    N_nodes : int
        Number of nodes in the tree
    logger_obj : logging.Logger, optional
        Logger instance for logging messages. If None, uses global logger.
    
    Returns
    -------
    tuple
        (penalty, actual_na_flip_ratio, refined_ω_NA, φ_adjusted, weight_na_to_1, weight_na_to_0)
    """
    # Use provided logger or fall back to global logger
    log = logger_obj if logger_obj is not None else logging.getLogger(__name__)
    
    # Calculate actual NA→1 flip ratio
    na_mask = input_binary_vec.isna()
    na_to_one_count = ((imputed_vec == 1) & na_mask).sum()
    total_na_count = na_mask.sum()
    
    if total_na_count > 0:
        actual_na_flip_ratio = na_to_one_count / total_na_count
    else:
        actual_na_flip_ratio = 0
    
    # Calculate w_fn value
    w_fn = fnfp_ratio * 1  # Default w_fn = 0.1
    
    # Fine-tuned dynamic adjustment: w_na approaches w_fn when flip ratio is high
    if actual_na_flip_ratio > 0.8:  # Over 80% of NAs flipped to 1
        refined_ω_NA = w_fn * 0.9  # Almost equal to w_fn
    elif actual_na_flip_ratio > 0.6:  # Over 60% of NAs flipped to 1
        refined_ω_NA = w_fn * 0.7  # 70% of w_fn
    elif actual_na_flip_ratio > 0.4:  # Over 40% of NAs flipped to 1
        refined_ω_NA = w_fn * 0.5  # 50% of w_fn
    elif actual_na_flip_ratio > 0.2:  # Over 20% of NAs flipped to 1
        refined_ω_NA = base_ω_NA * 4.0  # Significantly increase penalty
    elif actual_na_flip_ratio > 0.1:  # Over 10% of NAs flipped to 1
        refined_ω_NA = base_ω_NA * 2.0
    elif actual_na_flip_ratio > 0.05:  # Over 5% of NAs flipped to 1
        refined_ω_NA = base_ω_NA * 1.5
    else:
        refined_ω_NA = base_ω_NA  # Few flips, keep original penalty
    
    # Ensure it doesn't exceed w_fn (NA flips are inherently less reliable than FN)
    refined_ω_NA = min(refined_ω_NA, w_fn * 0.95)
    
    # Dynamically adjust BIC penalty based on NA characteristics and actual flip ratio
    φ_adjusted = compute_dynamic_bic_penalty_scaffold(
        base_φ, na_ratio, mut_ratio, actual_na_flip_ratio, placement_type, N_nodes, logger_obj=log
    )
    
    # Calculate penalty using adjusted weights
    penalty, weight_na_to_1, weight_na_to_0 = compute_bayesian_penalty_each_pos(
        input_binary_vec, posterior_vec, imputed_vec, fnfp_ratio, refined_ω_NA, logger_obj=log
    )
    
    log.debug(f"Dynamic penalty: NA flip ratio={actual_na_flip_ratio:.3f}, refined_ω_NA={refined_ω_NA:.4f}, φ_adjusted={φ_adjusted:.3f}")
    
    return penalty, actual_na_flip_ratio, refined_ω_NA, φ_adjusted, weight_na_to_1, weight_na_to_0


def compute_dynamic_bic_penalty_scaffold(base_φ, na_ratio, mut_ratio, actual_na_flip_ratio, placement_type, N_nodes, logger_obj=None):
    """
    Dynamically adjust BIC penalty based on NA characteristics and actual flip ratio.
    
    Parameters
    ----------
    base_φ : float
        Base BIC penalty parameter
    na_ratio : float
        NA ratio in data
    mut_ratio : float
        Mutation ratio in data
    actual_na_flip_ratio : float
        Actual NA→1 flip ratio
    placement_type : str
        Type of placement
    N_nodes : int
        Number of nodes in the tree
    logger_obj : logging.Logger, optional
        Logger instance for logging messages. If None, uses global logger.
    
    Returns
    -------
    float
        Adjusted BIC penalty
    """
    # Use provided logger or fall back to global logger
    log = logger_obj if logger_obj is not None else logging.getLogger(__name__)
    
    # Base adjustment: high NA, low mutation reduces BIC penalty
    if na_ratio > 0.7 and mut_ratio < 0.1:
        base_adjustment = 0.3
    elif na_ratio > 0.5 and mut_ratio < 0.2:
        base_adjustment = 0.6
    else:
        base_adjustment = 1.0
    
    # Further adjustment based on actual flips
    if actual_na_flip_ratio > 0.4:
        # If actual flip ratio is high, this placement may be unreliable
        # For new node creation, be more lenient (lower penalty)
        # For existing node placement, increase penalty (may cause many unreliable flips)
        if placement_type in ['new_leaf', 'on_edge']:
            flip_adjustment = 0.8  # New node creation, further reduce penalty
        else:
            flip_adjustment = 1.2  # Existing node placement, increase penalty
    elif actual_na_flip_ratio > 0.2:
        if placement_type in ['new_leaf', 'on_edge']:
            flip_adjustment = 0.9
        else:
            flip_adjustment = 1.1
    else:
        flip_adjustment = 1.0
    
    # Node count consideration: more nodes means more lenient for new node creation
    node_adjustment = 1.0
    if N_nodes > 20 and placement_type in ['new_leaf', 'on_edge']:
        node_adjustment = 0.9
    elif N_nodes > 50 and placement_type in ['new_leaf', 'on_edge']:
        node_adjustment = 0.8
    
    φ_adjusted = base_φ * base_adjustment * flip_adjustment * node_adjustment
    
    log.debug(f"Dynamic BIC: base={base_φ}, adjusted={φ_adjusted:.3f}")
    
    return max(φ_adjusted, 0.1)  # Ensure it doesn't go too low


def compute_intersection_based_penalty_scaffold(new_mut, position, intersection_nodes, M_current, I_selected, na_ratio, mut_ratio, actual_na_flip_ratio, logger_obj=None):
    """
    Adjust penalty based on intersection patterns and actual NA flips.
    
    Parameters
    ----------
    new_mut : str
        New mutation
    position : dict
        Position dictionary
    intersection_nodes : set
        Set of intersection nodes
    M_current : pd.DataFrame
        Current mutation matrix
    I_selected : pd.DataFrame
        Selected mutation matrix
    na_ratio : float
        NA ratio in data
    mut_ratio : float
        Mutation ratio in data
    actual_na_flip_ratio : float
        Actual NA→1 flip ratio
    logger_obj : logging.Logger, optional
        Logger instance for logging messages. If None, uses global logger.
    
    Returns
    -------
    float
        Extra penalty
    """
    # Use provided logger or fall back to global logger
    log = logger_obj if logger_obj is not None else logging.getLogger(__name__)
    
    extra_penalty = 0
    placement_type = position['placement_type']
    anchor = position['anchor']
    
    # Get new_mut's mutant cells
    mut_cells = set(I_selected.index[I_selected[new_mut].fillna(0) == 1])
    
    if placement_type == 'on_node':
        anchor_cells = set(M_current.index[M_current[anchor] == 1])
        
        # Case 1: new_mut only strongly intersects with one node, but placed on a different early node
        if len(intersection_nodes) == 1:
            sole_intersection = list(intersection_nodes)[0]
            if anchor != sole_intersection and anchor in [n for n in M_current.columns if n != 'ROOT']:
                # Adjust based on actual flip ratio
                flip_based_multiplier = 1.0 + actual_na_flip_ratio  # More flips = heavier penalty
                extra_penalty += np.log(len(M_current.columns)) * 0.8 * flip_based_multiplier
        
        # Case 2: anchor covers far more cells than new_mut's mutant cells
        if len(anchor_cells) > len(mut_cells) * 3:
            # Adjust FN penalty based on actual flip ratio
            fn_penalty = np.log(len(anchor_cells) - len(mut_cells)) * 0.3
            extra_penalty += fn_penalty * (1.0 + actual_na_flip_ratio)
            
        # Case 3: High NA mutation placed on non-intersecting node with high flip ratio
        if na_ratio > 0.7 and anchor not in intersection_nodes and actual_na_flip_ratio > 0.3:
            extra_penalty += np.log(len(M_current.columns)) * 0.5
    
    elif placement_type in ['new_leaf', 'on_edge']:
        # For high NA, low mutation with low flip ratio, reduce penalty for opening new nodes
        if na_ratio > 0.7 and mut_ratio < 0.1 and actual_na_flip_ratio < 0.2:
            bonus = -np.log(len(M_current.columns)) * 0.4 * (1.0 - actual_na_flip_ratio)
            extra_penalty += bonus
            
        # For mutations intersecting with multiple nodes, reduce penalty for opening new nodes
        if len(intersection_nodes) >= 2:
            bonus = -np.log(len(M_current.columns)) * 0.3
            extra_penalty += bonus
    
    log.debug(f"Intersection penalty for {new_mut}: {extra_penalty:.3f}")
    
    return extra_penalty


def compute_hierarchy_penalty_scaffold(new_mut, position, M_current, I_selected, parent_dict, 
                            na_ratio, mut_ratio, actual_na_flip_ratio, logger_obj=None):
    """
    Calculate hierarchy rationality penalty, considering actual NA flips.
    
    Parameters
    ----------
    new_mut : str
        New mutation
    position : dict
        Position dictionary
    M_current : pd.DataFrame
        Current mutation matrix
    I_selected : pd.DataFrame
        Selected mutation matrix
    parent_dict : dict
        Parent node dictionary
    na_ratio : float
        NA ratio in data
    mut_ratio : float
        Mutation ratio in data
    actual_na_flip_ratio : float
        Actual NA→1 flip ratio
    logger_obj : logging.Logger, optional
        Logger instance for logging messages. If None, uses global logger.
    
    Returns
    -------
    float
        Hierarchy penalty
    """
    # Use provided logger or fall back to global logger
    log = logger_obj if logger_obj is not None else logging.getLogger(__name__)
    
    hierarchy_penalty = 0
    placement_type = position['placement_type']
    anchor = position['anchor']
    
    mut_cells = set(I_selected.index[I_selected[new_mut].fillna(0) == 1])
    
    if placement_type == 'on_node':
        anchor_cells = set(M_current.index[M_current[anchor] == 1])
        
        # Check if this anchor has children
        anchor_children = find_children_of_node_scaffold(anchor, M_current.columns, parent_dict)
        
        if anchor_children:
            # If anchor has children and new_mut has a specific pattern with children
            children_intersection = False
            for child in anchor_children:
                child_cells = set(M_current.index[M_current[child] == 1])
                if mut_cells.issubset(child_cells) and len(mut_cells) < len(child_cells) * 0.8:
                    children_intersection = True
                    break
            
            if children_intersection:
                # Adjust hierarchy penalty based on actual flip ratio
                flip_multiplier = 1.0 + actual_na_flip_ratio
                hierarchy_penalty += np.log(len(anchor_children) + 1) * 0.4 * flip_multiplier
    
    elif placement_type == 'new_leaf':
        parent = anchor
        parent_cells = set(M_current.index[M_current[parent] == 1])
        
        # If new_mut's cells are a proper subset of parent cells, this is a reasonable subclone
        if mut_cells.issubset(parent_cells) and len(mut_cells) < len(parent_cells) * 0.8:
            # Give reward, reward magnitude depends on actual flip ratio (fewer flips = more reward)
            bonus_multiplier = 1.0 - actual_na_flip_ratio  # Fewer flips = more reward
            hierarchy_penalty -= np.log(len(parent_cells) - len(mut_cells)) * 0.2 * bonus_multiplier
            
        # For high NA mutations with low flip ratio, further reward reasonable subclone placement
        if na_ratio > 0.7 and mut_cells.issubset(parent_cells) and actual_na_flip_ratio < 0.2:
            bonus_multiplier = 1.0 - actual_na_flip_ratio
            hierarchy_penalty -= np.log(len(parent_cells)) * 0.3 * bonus_multiplier
    
    log.debug(f"Hierarchy penalty for {new_mut}: {hierarchy_penalty:.3f}")
    
    return hierarchy_penalty


def find_children_of_node_scaffold(node, all_columns, parent_dict, logger_obj=None):
    """
    Find direct children of a node.
    
    Parameters
    ----------
    node : str
        Node name
    all_columns : list
        All column names
    parent_dict : dict
        Parent node dictionary
    logger_obj : logging.Logger, optional
        Logger instance for logging messages. If None, uses global logger.
    
    Returns
    -------
    list
        List of direct children
    """
    # Use provided logger or fall back to global logger
    log = logger_obj if logger_obj is not None else logging.getLogger(__name__)
    
    children = []
    for col in all_columns:
        if col == node or col == 'ROOT':
            continue
        if parent_dict.get(col) == node:
            children.append(col)
    
    log.debug(f"Children of {node}: {children}")
    
    return children




def add_new_mutation_to_tree_independent(new_mut, T_current, final_position, logger_obj=None):
    """
    Add a new mutation to the tree independently based on the placement type.
    
    Parameters
    ----------
    new_mut : str
        New mutation to add
    T_current : TreeNode
        Current tree
    final_position : dict
        Final position dictionary containing placement_type, anchor, meta
    logger_obj : logging.Logger, optional
        Logger instance for logging messages. If None, uses global logger.
    
    Returns
    -------
    TreeNode
        Updated tree
    """
    # Use provided logger or fall back to global logger
    log = logger_obj if logger_obj is not None else logging.getLogger(__name__)
    
    placement_type = final_position['placement_type']
    anchor = final_position['anchor']
    meta = final_position['meta']
    
    # Get anchor node
    anchor_node = T_current.find(anchor)
    if not anchor_node:
        raise ValueError(f"Anchor node {anchor} not found in the tree.")
    
    # Select insertion strategy based on placement_type
    if placement_type == 'on_node':
        # For on_node, just append |new_mut to the anchor node name
        new_mut_name = new_mut
        anchor_node.name = f"{anchor_node.name}|{new_mut_name}"
        log.debug(f"Added {new_mut} to node {anchor}")
        return T_current
    
    elif placement_type == 'new_leaf':
        # For new_leaf, add a new leaf node under the anchor
        new_leaf_name = new_mut
        new_leaf = TreeNode(new_leaf_name)
        anchor_node.add_child(new_leaf)
        log.debug(f"Added new leaf {new_mut} under {anchor}")
        return T_current
    
    elif placement_type == 'on_edge':
        # For on_edge, find the child in meta and insert a new node
        child_name = meta.get('child')
        child_node = T_current.find(child_name)
        if not child_node:
            raise ValueError(f"Child node {child_name} not found in the tree.")
        
        # Create new node and insert it between anchor and child
        new_edge_name = new_mut
        T_current.insert_on_edge(anchor_node, child_node, new_edge_name)
        log.debug(f"Inserted {new_mut} on edge between {anchor} and {child_name}")
        return T_current
    
    elif placement_type == 'new_parent_merge':
        # For new_parent_merge, merge multiple children under a new parent
        children_to_merge = meta.get('merge_children', [])
        children_nodes = [T_current.find(child) for child in children_to_merge]
        if None in children_nodes:
            raise ValueError(f"One or more merge children not found in the tree.")
        
        new_parent_name = new_mut
        new_parent = anchor_node.add_new_parent_for_children(children_nodes, new_parent_name)
        log.debug(f"Merged children {children_to_merge} under new parent {new_mut} attached to {anchor}")
        return T_current
    
    else:
        raise ValueError(f"Unknown placement type: {placement_type}")


def add_new_mutation_to_tree_conflict_free(new_mut, T_current, pos, M_current, parent_dict, I_selected, logger_obj=None):
    """
    Add a new mutation to the tree, ensuring no conflicts with the existing tree structure.
    Based on the logic from compute_bayesian_penalty_for_positions_scaffold.
    
    Parameters
    ----------
    new_mut : str
        New mutation to add
    T_current : TreeNode
        Current tree
    pos : dict
        Position dictionary
    M_current : pd.DataFrame
        Current mutation matrix
    parent_dict : dict
        Parent node dictionary
    I_selected : pd.DataFrame
        Selected mutation matrix
    logger_obj : logging.Logger, optional
        Logger instance for logging messages. If None, uses global logger.
    
    Returns
    -------
    tuple
        (T_current, M_current, imputed_vec)
    """
    # Use provided logger or fall back to global logger
    log = logger_obj if logger_obj is not None else logging.getLogger(__name__)
    
    placement_type = pos['placement_type']
    anchor = pos['anchor']
    meta = pos.get('meta', {})
    
    # Get anchor node
    anchor_node = T_current.find(anchor)
    if not anchor_node:
        raise ValueError(f"Anchor node {anchor} not found in the tree.")
    
    # Get binary vector of new mutation
    new_mut_bin_vector = I_selected[new_mut].replace({pd.NA: np.nan}).fillna(0).astype(int)
    
    # Default imputed vector
    imputed_vec = pd.Series(0, index=M_current.index)
    
    # Calculate conflict-free imputed_vec based on placement_type
    if placement_type == 'on_node':
        # Place on node, use anchor node's vector directly
        imputed_vec = M_current[anchor].astype(int)
        
        # Update node name
        anchor_node.name = f"{anchor_node.name}|{new_mut}"
        log.debug(f"Added {new_mut} to node {anchor} (on_node)")
        
    elif placement_type == 'new_leaf':
        parent = anchor
        vec_parent = M_current[parent] if parent != 'ROOT' else pd.Series(1, index=M_current.index)
        
        # Get direct sibling conflict nodes
        sibling_nodes = [n['name'] for n in pos['nodes'] if n['parent'] == parent and n['name'] != new_mut]
        
        # Get all conflict nodes outside lineage
        lineage_conflict_nodes = get_all_conflict_nodes_outside_lineage_scaffold(
            parent, build_lineage_parent_dict_from_tree(T_current, anchor), M_current.columns, logger_obj=log
        )
        
        # Merge all conflict nodes (deduplicate)
        all_conflict_nodes = list(set(sibling_nodes + lineage_conflict_nodes))
        
        # Build conflict vector
        vec_conflicts = pd.Series(0, index=M_current.index)
        for conflict in all_conflict_nodes:
            vec_conflicts |= M_current[conflict]
        
        # Correct logic: first exclude all conflicts, then intersect with parent
        new_mut_cleaned = new_mut_bin_vector & ~vec_conflicts
        imputed_vec = new_mut_cleaned.astype(int)
        
        # Add new leaf node
        new_leaf = TreeNode(new_mut)
        anchor_node.add_child(new_leaf)
        log.debug(f"Added new leaf {new_mut} under {anchor} (new_leaf)")
        
    elif placement_type == 'on_edge':
        parent = anchor
        child = pos['meta']['child']
        vec_parent = M_current[parent] if parent != 'ROOT' else pd.Series(1, index=M_current.index)
        vec_child = M_current[child]
        
        # Get direct sibling conflict nodes
        sibling_nodes = [n['name'] for n in pos['nodes'] if n['parent'] == parent and n['name'] not in [child, new_mut]]
        
        # Get all conflict nodes outside lineage
        lineage_conflict_nodes = get_all_conflict_nodes_outside_lineage_scaffold(
            parent, build_lineage_parent_dict_from_tree(T_current, anchor), M_current.columns, logger_obj=log
        )
        
        # Merge all conflict nodes (deduplicate)
        all_conflict_nodes = list(set(sibling_nodes + lineage_conflict_nodes))
        
        # Build conflict vector
        vec_conflicts = pd.Series(0, index=M_current.index)
        for conflict in all_conflict_nodes:
            vec_conflicts |= M_current[conflict]
        
        # Correct logic: child ∪ (cleaned new_mut ∩ parent)
        new_mut_cleaned = new_mut_bin_vector & ~vec_conflicts
        imputed_vec = (vec_child | new_mut_cleaned).astype(int)
        
        # Get child node for insertion
        child_node = T_current.find(child)
        if not child_node:
            raise ValueError(f"Child node {child} not found in the tree.")
        
        # Insert new node on edge
        T_current.insert_on_edge(anchor_node, child_node, new_mut)
        log.debug(f"Inserted {new_mut} on edge between {anchor} and {child} (on_edge)")
        
    elif placement_type == 'new_parent_merge':
        parent = anchor
        merge_children = pos['meta']['merge_children']
        vec_parent = M_current[parent] if parent != 'ROOT' else pd.Series(1, index=M_current.index)
        
        # Build union vector of children
        vec_children = pd.Series(0, index=M_current.index)
        for c in merge_children:
            vec_children |= M_current[c]
        
        # Get direct sibling conflict nodes
        sibling_nodes = [n['name'] for n in pos['nodes'] if n['parent'] == parent and n['name'] not in merge_children + [new_mut]]
        
        # Get all conflict nodes outside lineage (excluding merge_children)
        lineage_conflict_nodes = get_all_conflict_nodes_outside_lineage_scaffold(
            parent, build_lineage_parent_dict_from_tree(T_current, anchor), M_current.columns, 
            exclude_nodes=merge_children, logger_obj=log
        )
        
        # Merge all conflict nodes (deduplicate)
        all_conflict_nodes = list(set(sibling_nodes + lineage_conflict_nodes))
        
        # Build conflict vector
        vec_conflicts = pd.Series(0, index=M_current.index)
        for conflict in all_conflict_nodes:
            vec_conflicts |= M_current[conflict]
        
        # Correct logic: children ∪ (cleaned new_mut ∩ parent)
        new_mut_cleaned = new_mut_bin_vector & ~vec_conflicts
        imputed_vec = (vec_children | new_mut_cleaned).astype(int)
        
        # Get child nodes for merging
        children_nodes = [T_current.find(child) for child in merge_children]
        if None in children_nodes:
            raise ValueError(f"One or more merge children not found in the tree.")
        
        # Create new parent and merge children
        new_parent = anchor_node.add_new_parent_for_children(children_nodes, new_mut)
        log.debug(f"Merged children {merge_children} under new parent {new_mut} attached to {anchor} (new_parent_merge)")
        
    else:
        raise ValueError(f"Unknown placement type: {placement_type}")
    
    # Update mutation matrix M_current
    # First add the new column
    M_current[new_mut] = imputed_vec
    
    # For cells where imputed_vec is 1, ensure their parent chain is also 1
    if placement_type != 'on_node':  # on_node doesn't need this as it already inherits parent state
        full_mutnode_chain = get_full_mutnode_chain_with_anchor_scaffold(anchor, parent_dict, logger_obj=log)
        cells_with_final_one = imputed_vec[imputed_vec == 1].index.tolist()
        
        if len(cells_with_final_one) > 0:
            for cell in cells_with_final_one:
                for mutation in full_mutnode_chain:
                    if M_current.loc[cell, mutation] == 0:
                        M_current.loc[cell, mutation] = 1
    
    log.info(f"Successfully added mutation {new_mut} to tree with placement type: {placement_type}")
    
    return T_current, M_current, imputed_vec


def WriteTfile(out_prefix, matrix, rownames, colnames, judge, logger_obj=None):
    """
    Write matrix output as an integer matrix in the format specified in the documentation.
    
    Parameters
    ----------
    out_prefix : str
        Output file prefix
    matrix : pd.DataFrame
        Matrix to write
    rownames : list
        Row names
    colnames : list
        Column names
    judge : str
        Whether to check for conflict-free status ('yes' or 'no')
    logger_obj : logging.Logger, optional
        Logger instance for logging messages. If None, uses global logger.
    """
    # Use provided logger or fall back to global logger
    log = logger_obj if logger_obj is not None else logging.getLogger(__name__)
    
    matrix_output = matrix.astype(int)
    df_output = pd.DataFrame(matrix_output)
    df_output.index = rownames
    df_output.columns = colnames
    df_output.index.name = "cellIDxmutID"
    is_cf = scp.ul.is_conflict_free_gusfield(df_output)
    
    if judge == "yes":
        if is_cf:
            log.info("Current tree is conflict-free, outputting binary matrix and plotting phylogenetic tree")
            df_output.to_csv(out_prefix + ".CFMatrix", sep="\t")
            tree = scp.ul.to_tree(df_output)
            scp.pl.clonal_tree(tree, output_file=out_prefix + ".tree_scphylo.pdf")
        else:
            log.warning("Current tree is NOT conflict-free!")
    else:
        log.info("Only outputting binary matrix (no conflict-free check)")
        df_output.to_csv(out_prefix + ".CFMatrix", sep="\t")


def find_flipping_spots(series_in_bin, series_phylogeny, condition_in_bin, condition_phylogeny, logger_obj=None):
    """
    Find a list of eligible line names (Spots) based on conditions.
    
    Parameters
    ----------
    series_in_bin : pd.Series
        Binary series
    series_phylogeny : pd.Series
        Phylogeny series
    condition_in_bin : int
        Condition for binary series
    condition_phylogeny : int
        Condition for phylogeny series
    logger_obj : logging.Logger, optional
        Logger instance for logging messages. If None, uses global logger.
    
    Returns
    -------
    list
        List of eligible line names (spots)
    """
    # Use provided logger or fall back to global logger
    log = logger_obj if logger_obj is not None else logging.getLogger(__name__)
    
    result = series_in_bin[(series_in_bin == condition_in_bin) & (series_phylogeny == condition_phylogeny)].index.tolist()
    
    log.debug(f"Found {len(result)} flipping spots for conditions: bin={condition_in_bin}, phylo={condition_phylogeny}")
    
    return result


def integrate_mutations_to_scaffold_within_group(
    sorted_attached_mutations, T_current, M_current, I_attached, P_attached, 
    mutation_group, ω_NA, fnfp_ratio, φ, logger_obj=None, max_retries=None,
    show_progress=False
):
    """
    Process external mutations and integrate them into the scaffold phylogenetic tree (with rollback and retry support).
    
    Parameters
    ----------
    sorted_attached_mutations : list
        Sorted list of mutations to process
    T_current : TreeNode
        Current phylogenetic tree structure
    M_current : pd.DataFrame
        Current mutation matrix
    I_attached : pd.DataFrame
        Attached mutation information
    P_attached : pd.DataFrame
        Attached probability information
    mutation_group : dict
        Mutation group information
    ω_NA : float
        NA weight parameter
    fnfp_ratio : float
        False negative to false positive ratio
    φ : float
        Bayesian penalty parameter
    logger_obj : logging.Logger, optional
        Logger instance for logging messages. If None, uses global logger.
    max_retries : int, optional
        Maximum number of candidate positions to try. None means try all.
    show_progress : bool, default=False
        Whether to show progress bars (disabled in parallel environments)
    
    Returns
    -------
    tuple : (external_mutations, conflict_mutations, T_current, M_current)
    """
    # Use provided logger or fall back to global logger
    log = logger_obj if logger_obj is not None else logging.getLogger(__name__)
    
    external_mutations = []
    conflict_mutations = []
    
    # Determine whether to show progress (only in interactive environments)
    use_progress = show_progress and sys.stdout.isatty()
    
    mut_iterator = sorted_attached_mutations
    if use_progress:
        mut_iterator = tqdm(sorted_attached_mutations, desc="Processing mutations", unit="mutation")
    
    for new_mut in mut_iterator:
        
        # 1. Find intersection nodes
        intersection_nodes = find_all_intersect_muts_from_tree_by_matrix_scaffold(
            T_current, I_attached, new_mut, logger_obj=log
        )
        if len(intersection_nodes) == 0:
            external_mutations.append(new_mut)
            continue
        
        # 2. Get candidate positions
        refined_positions = find_intersection_positions_within_group_directly_scaffold(
            T_current, new_mut, I_attached, mutation_group, min_overlap=1, intersection_nodes=intersection_nodes, logger_obj=log
        )
        parent_dict = build_parent_dict_from_candidates_scaffold(refined_positions)
        
        if len(refined_positions) == 0:
            external_mutations.append(new_mut)
            log.warning(f"Mutation {new_mut} added to external_mutations (no candidate positions found)")
            continue
        
        # 3. Backup current state
        M_backup = M_current.copy()
        T_backup = T_current.copy()
        
        # 4. Calculate penalties for all candidate positions
        df_penalty = compute_bayesian_penalty_for_all_positions_scaffold(
            new_mut, refined_positions, T_current, M_current, I_attached, P_attached, 
            parent_dict, intersection_nodes, ω_NA=ω_NA, fnfp_ratio=fnfp_ratio, φ=φ,
            logger_obj=log
        )
        
        if df_penalty.empty:
            external_mutations.append(new_mut)
            log.warning(f"Mutation {new_mut} added to external_mutations (no valid penalty scores)")
            continue
        
        # 5. Sort by penalty, filter out positions where imputed_vec is all zeros
        df_valid = df_penalty[df_penalty['imputed_vec'].apply(lambda x: x.sum() > 0)]
        if df_valid.empty:
            external_mutations.append(new_mut)
            log.warning(f"Mutation {new_mut} added to external_mutations (all imputed_vec are zero)")
            continue
        
        df_sorted = df_valid.sort_values('total_penalty')
        
        # 6. Determine number of candidates to try
        if max_retries is None:
            candidates_to_try = df_sorted
        else:
            candidates_to_try = df_sorted.head(max_retries)
        
        # 7. Iterate through candidate positions
        placed = False
        for attempt, (idx, row) in enumerate(candidates_to_try.iterrows()):
            
            # Restore backup
            M_current = M_backup.copy()
            T_current = T_backup.copy()
            
            try:
                # Apply this position to tree and matrix
                T_current, M_current = apply_position_to_tree_scaffold(
                    new_mut, row['position'], row['imputed_vec'], 
                    T_current, M_current, I_attached, parent_dict, logger_obj=log
                )
                
                # Check for conflicts
                if scp.ul.is_conflict_free_gusfield(M_current):
                    placed = True
                    log.debug(f"Mutation {new_mut} successfully placed at position {row['position_index']}")
                    break
                else:
                    log.warning(f"✗ Position {row['position_index']} caused conflict, trying next candidate")
                    
            except Exception as e:
                log.error(f"Error placing mutation at position {row['position_index']}: {e}")
                continue
        
        # 8. Handle result
        if not placed:
            M_current = M_backup.copy()
            T_current = T_backup.copy()
            conflict_mutations.append(new_mut)
            log.warning(f"Mutation {new_mut} added to conflict_mutations (all {len(candidates_to_try)} candidates failed)")
    
    return external_mutations, conflict_mutations, T_current, M_current


def process_subtree_mutations_to_specific_node(
    subtree_groups, target_node_names, T_current, M_current, I_selected, P_selected, 
    ω_NA, fnfp_ratio, φ, logger_obj=None, root_mutations=None, mutation_group=None,
    show_progress=False
):
    """
    Process external mutations through subtree groups, supporting both multi-mutation groups
    and single-mutation groups, attaching them to the specified node.
    
    Parameters
    ----------
    subtree_groups : list
        List of subtree groups
    target_node_names : str or list
        Target node name(s) for attachment
    T_current : TreeNode
        Current phylogenetic tree
    M_current : pd.DataFrame
        Current mutation matrix
    I_selected : pd.DataFrame
        Selected mutation information
    P_selected : pd.DataFrame
        Selected probability information
    ω_NA : float
        NA weight parameter
    fnfp_ratio : float
        False negative to false positive ratio
    φ : float
        Bayesian penalty parameter
    logger_obj : logging.Logger, optional
        Logger instance for logging messages. If None, uses global logger.
    root_mutations : list, optional
        List of root mutations
    mutation_group : dict, optional
        Mutation group information
    show_progress : bool, default=False
        Whether to show progress bars (disabled in parallel environments)
    
    Returns
    -------
    tuple : (external_mutations, conflict_mutations, T_current, M_current, root_mutations)
    """
    # Use provided logger or fall back to global logger
    log = logger_obj if logger_obj is not None else logging.getLogger(__name__)
    
    if root_mutations is None:
        root_mutations = []
    
    if mutation_group is None:
        mutation_group = {}
    
    external_mutations = []
    conflict_mutations = []
    
    # Find target node
    if target_node_names is None:
        raise ValueError(f"Target node '{target_node_names}' not found in the tree")
    
    # Determine whether to show progress (only in interactive environments)
    use_progress = show_progress and sys.stdout.isatty()
    
    # Separate subtree groups into multi-mutation and singleton groups
    multi_mut_subtree_groups = [g for g in subtree_groups if len(g) > 1]
    singleton_subtree_groups = [g for g in subtree_groups if len(g) == 1]
    
    # ============================================================
    # Process multi-mutation groups
    # ============================================================
    if use_progress:
        group_iterator = tqdm(multi_mut_subtree_groups, desc="Processing multiple subtrees")
    else:
        group_iterator = multi_mut_subtree_groups
    
    for group_idx, group in enumerate(group_iterator):
        
        sorted_group = sorted(group, key=lambda mut: I_selected[mut].sum(), reverse=True)
        viable_mutations = [mut for mut in sorted_group 
                           if len(find_all_intersect_muts_from_tree_by_matrix_scaffold(T_current, I_selected, mut, logger_obj=log)) > 0]
        non_viable_mutations = [mut for mut in sorted_group if mut not in viable_mutations]
        
        if not viable_mutations:
            external_mutations.extend(sorted_group)
            continue
        
        sorted_group = viable_mutations + non_viable_mutations
        current_group_start_node = None
        reattached_mutations = []
        
        # Process mutations within group
        if use_progress:
            mut_iterator = tqdm(sorted_group, desc="Processing mutations in group")
        else:
            mut_iterator = sorted_group
        
        for idx, subtree_mut in enumerate(mut_iterator):
            
            if idx == 0:
                # First mutation: use the integrate function
                ext_temp, conf_temp, T_current, M_current = integrate_mutations_to_scaffold_within_group(
                    sorted_attached_mutations=[subtree_mut],
                    T_current=T_current,
                    M_current=M_current,
                    I_attached=I_selected,
                    P_attached=P_selected,
                    mutation_group=mutation_group,
                    ω_NA=ω_NA,
                    fnfp_ratio=fnfp_ratio,
                    φ=φ,
                    logger_obj=log,
                    max_retries=None,
                    show_progress=False
                )
                external_mutations.extend(ext_temp)
                conflict_mutations.extend(conf_temp)
                
                if subtree_mut in ext_temp or subtree_mut in conf_temp:
                    log.warning(f"First mutation {subtree_mut} failed to place, skipping group")
                    reattached_mutations.extend(sorted_group[idx+1:])
                    continue
                
                current_group_start_node = T_current.find(subtree_mut)
            
            else:
                if current_group_start_node is None:
                    reattached_mutations.append(subtree_mut)
                    continue
                
                # Subsequent mutations: use the integrate function
                ext_temp, conf_temp, T_current, M_current = integrate_mutations_to_scaffold_within_group(
                    sorted_attached_mutations=[subtree_mut],
                    T_current=T_current,
                    M_current=M_current,
                    I_attached=I_selected,
                    P_attached=P_selected,
                    mutation_group=mutation_group,
                    ω_NA=ω_NA,
                    fnfp_ratio=fnfp_ratio,
                    φ=φ,
                    logger_obj=log,
                    max_retries=None,
                    show_progress=False
                )
                external_mutations.extend(ext_temp)
                conflict_mutations.extend(conf_temp)
                
                if subtree_mut in ext_temp or subtree_mut in conf_temp:
                    reattached_mutations.append(subtree_mut)
                    continue
        
        # Process reattached mutations
        if reattached_mutations:
            for subtree_mut in reattached_mutations:
                # Use the integrate function
                ext_temp, conf_temp, T_current, M_current = integrate_mutations_to_scaffold_within_group(
                    sorted_attached_mutations=[subtree_mut],
                    T_current=T_current,
                    M_current=M_current,
                    I_attached=I_selected,
                    P_attached=P_selected,
                    mutation_group=mutation_group,
                    ω_NA=ω_NA,
                    fnfp_ratio=fnfp_ratio,
                    φ=φ,
                    logger_obj=log,
                    max_retries=None,
                    show_progress=False
                )
                external_mutations.extend(ext_temp)
                conflict_mutations.extend(conf_temp)
    
    # ============================================================
    # Process singleton mutation groups
    # ============================================================
    if use_progress:
        singleton_iterator = tqdm(singleton_subtree_groups, desc="Processing singleton subtrees")
    else:
        singleton_iterator = singleton_subtree_groups
    
    for group in singleton_iterator:
        subtree_mut = group[0]
        
        # Use the integrate function
        ext_temp, conf_temp, T_current, M_current = integrate_mutations_to_scaffold_within_group(
            sorted_attached_mutations=[subtree_mut],
            T_current=T_current,
            M_current=M_current,
            I_attached=I_selected,
            P_attached=P_selected,
            mutation_group=mutation_group,
            ω_NA=ω_NA,
            fnfp_ratio=fnfp_ratio,
            φ=φ,
            logger_obj=log,
            max_retries=None,
            show_progress=False
        )
        external_mutations.extend(ext_temp)
        conflict_mutations.extend(conf_temp)
    
    log.info(f"Processed {len(subtree_groups)} subtree groups: {len(multi_mut_subtree_groups)} multi-mutation, {len(singleton_subtree_groups)} singleton")
    log.info(f"External mutations: {len(external_mutations)}, Conflict mutations: {len(conflict_mutations)}")
    
    return external_mutations, conflict_mutations, T_current, M_current, root_mutations


# -------------------------
# Sub-grouping within Backbone Groups
# -------------------------

def find_latest_hub_node(root, hub_mutations, logger_obj=None):
    """
    Find all latest nodes containing hub_mutations, and return these nodes along with the hub_mutations they contain.
    
    Parameters
    ----------
    root : TreeNode
        Root node of the tree
    hub_mutations : List[str]
        List of hub mutations to search for
    logger_obj : logging.Logger, optional
        Logger instance for logging messages. If None, uses global logger.
    
    Returns
    -------
    Tuple[List[TreeNode], List[List[str]]]
        (list of latest nodes, list of hub_mutations in each corresponding node)
    """
    # Use provided logger or fall back to global logger
    log = logger_obj if logger_obj is not None else logging.getLogger(__name__)
    
    hub_nodes = []
    
    # Traverse all non-ROOT nodes
    for node in root.traverse():
        if node.name == "ROOT":
            continue
            
        # Parse mutations in the node
        if "|" in node.name:
            node_mutations = node.name.split("|")
        else:
            node_mutations = [node.name]
        
        # Find hub_mutations contained in this node
        hub_in_node = [mut for mut in node_mutations if mut in hub_mutations]
        
        if hub_in_node:
            hub_nodes.append((node, hub_in_node))
    
    if not hub_nodes:
        log.debug("No hub nodes found")
        return [], []
    
    # Calculate depth of all nodes
    def calculate_depth(node):
        """Calculate node depth (ROOT depth is 0)"""
        depth = 0
        current = node
        while current.parent:
            depth += 1
            current = current.parent
        return depth
    
    # Calculate depth for all hub nodes
    hub_nodes_with_depth = [(node, hub_muts, calculate_depth(node)) 
                           for node, hub_muts in hub_nodes]
    
    # Find the maximum depth
    max_depth = max(depth for _, _, depth in hub_nodes_with_depth)
    
    # Find all nodes with depth equal to the maximum depth
    latest_nodes = [(node, hub_muts) for node, hub_muts, depth in hub_nodes_with_depth 
                   if depth == max_depth]
    
    # Separate nodes list and mutations list
    if latest_nodes:
        nodes, mutations = zip(*latest_nodes)
        log.debug(f"Found {len(nodes)} latest hub nodes at depth {max_depth}")
        return list(nodes), list(mutations)
    else:
        return [], []


# # Usage example:
# latest_nodes, found_mutations_list = find_latest_hub_node(T_current, hub_mutations)

# if latest_nodes:
#     log.info(f"Found {len(latest_nodes)} latest nodes containing hub mutations:")
#     for i, (node, mutations) in enumerate(zip(latest_nodes, found_mutations_list)):
#         log.info(f"Node {i+1}: {node.name}")
#         log.info(f"  Hub mutations contained: {mutations}")
# else:
#     log.info("No nodes found containing hub mutations")


# More concise version (using vectorized operations)
def calculate_group_cooccurrence_fraction(sorted_parallel_groups, backbone_mut, I_selected, logger_obj=None):
    """
    Simplified version using vectorized operations for improved efficiency.
    
    Parameters
    ----------
    sorted_parallel_groups : dict
        Dictionary of sorted parallel groups
    backbone_mut : str
        Backbone mutation name
    I_selected : pd.DataFrame
        Selected mutation matrix
    logger_obj : logging.Logger, optional
        Logger instance for logging messages. If None, uses global logger.
    
    Returns
    -------
    dict
        {group_id: fraction} for each group
    """
    # Use provided logger or fall back to global logger
    log = logger_obj if logger_obj is not None else logging.getLogger(__name__)
    
    # Find cells where backbone mutation is 1
    backbone_pos_mask = I_selected[backbone_mut] == 1
    backbone_cells = I_selected[backbone_pos_mask].index
    
    if len(backbone_cells) == 0:
        log.debug(f"No cells with backbone mutation {backbone_mut} found")
        return {group_id: 0.0 for group_id in sorted_parallel_groups.keys()}
    
    group_fractions = {}
    
    for group_id, mutation_list in sorted_parallel_groups.items():
        # Check if mutations exist
        valid_mutations = [mut for mut in mutation_list if mut in I_selected.columns]
        if not valid_mutations:
            group_fractions[group_id] = 0.0
            continue
        
        # Extract data and calculate
        group_data = I_selected.loc[backbone_cells, valid_mutations]
        
        # Use max(axis=1) to find if each cell has at least one mutation = 1 (ignoring NA)
        # Fill NA with 0 first, then find max
        group_data_filled = group_data.fillna(0)
        has_mutation = (group_data_filled == 1).any(axis=1)
        
        fraction = has_mutation.mean() if len(has_mutation) > 0 else 0.0
        group_fractions[group_id] = fraction
    
    log.debug(f"Calculated co-occurrence fractions for {len(group_fractions)} groups")
    
    return group_fractions


# # Usage example:
# group_fractions = calculate_group_cooccurrence_fraction(
#     sorted_parallel_groups, backbone_mut, I_selected
# )
# # Example output:
# # {0: 0.85, 1: 0.72, 2: 0.63}


def find_max_fraction_group(group_fractions, logger_obj=None):
    """
    Find the group key with the maximum fraction.
    
    Parameters
    ----------
    group_fractions : dict
        {group_id: fraction}
    logger_obj : logging.Logger, optional
        Logger instance for logging messages. If None, uses global logger.
    
    Returns
    -------
    tuple
        (max_key, max_value)
    """
    # Use provided logger or fall back to global logger
    log = logger_obj if logger_obj is not None else logging.getLogger(__name__)
    
    if not group_fractions:
        log.debug("No group fractions provided")
        return None, 0.0
    
    max_key = max(group_fractions, key=group_fractions.get)
    max_value = group_fractions[max_key]
    
    log.debug(f"Maximum fraction group: {max_key}, fraction: {max_value:.3f}")
    
    return max_key, max_value


# # Usage example
# max_key, max_value = find_max_fraction_group(group_fractions)
# log.info(f"Maximum fraction group: {max_key}, fraction: {max_value:.3f}")
# # Output: Maximum fraction group: 1, fraction: 0.069


def calculate_cv_for_single_mutation(df_reads_group, mutation, cv_thresh, logger_obj=None):
    """
    Calculate CV statistics for a single mutation.
    
    Parameters
    ----------
    df_reads_group : pd.DataFrame
        DataFrame containing coverage data
    mutation : str
        Mutation name
    cv_thresh : float
        CV threshold
    logger_obj : logging.Logger, optional
        Logger instance for logging messages. If None, uses global logger.
    
    Returns
    -------
    dict
        Dictionary containing CV statistics
    """
    # Use provided logger or fall back to global logger
    log = logger_obj if logger_obj is not None else logging.getLogger(__name__)
    
    try:
        cov_data = df_reads_group[mutation].dropna()
        if len(cov_data) == 0:
            log.debug(f"No coverage data for mutation {mutation}")
            return {'cov_median': 0, 'cov_CV': float('inf'), 'cov_mean': 0, 'cov_std': 0, 'pass_CV': False}
        
        cov_median = np.median(cov_data)
        cov_mean = np.mean(cov_data)
        cov_std = np.std(cov_data)
        cov_CV = cov_std / cov_mean if cov_mean > 0 else float('inf')
        
        result = {
            'cov_median': cov_median,
            'cov_CV': cov_CV,
            'cov_mean': cov_mean,
            'cov_std': cov_std,
            'pass_CV': cov_CV <= cv_thresh
        }
        
        log.debug(f"CV for {mutation}: median={cov_median:.2f}, CV={cov_CV:.2f}, pass={result['pass_CV']}")
        
        return result
        
    except Exception as e:
        log.warning(f"Error calculating CV for mutation {mutation}: {str(e)}")
        return {'cov_median': 0, 'cov_CV': float('inf'), 'cov_mean': 0, 'cov_std': 0, 'pass_CV': False}


def perform_subgrouping_within_backbone_groups_and_build_initial_scaffold_tree(
    sorted_I_resolved, sorted_P_resolved, T_current, M_current, mutation_group, 
    backbones_of_group, df_reads_resolved, df_features_new, outputpath, sampleid, 
    logger_obj, params, M_B, root_mutations, cutoff_mcf_for_graph, cutoff_mcn_for_graph,
    show_progress=False
):
    """
    Perform sub-grouping within each backbone group to identify sub-clones.
    
    Parameters
    ----------
    sorted_I_resolved : pd.DataFrame
        Resolved binary mutation matrix
    sorted_P_resolved : pd.DataFrame
        Resolved posterior mutation matrix
    T_current : TreeNode
        Current tree structure
    M_current : pd.DataFrame
        Current mutation matrix
    mutation_group : dict
        Mapping of mutations to main groups
    backbones_of_group : dict
        Mapping of group IDs to backbone mutations
    df_reads_resolved : pd.DataFrame
        Reads dataframe for coverage information
    outputpath : str
        Output directory path
    sampleid : str
        Sample ID
    logger_obj : logging.Logger
        Logger instance
    params : dict
        Parameters including general settings
    M_B : pd.DataFrame
        Mutation backbone matrix
    root_mutations : list
        List of root mutations
    cutoff_mcf_for_graph : float
        MCF cutoff for graph construction
    cutoff_mcn_for_graph : int
        MCN cutoff for graph construction
    show_progress : bool, default=False
        Whether to show progress bars
    
    Returns
    -------
    Tuple containing:
        - complete_mutation_hierarchy: Complete hierarchy information for all mutations
        - subgroup_backbone_mutations: Backbone mutations for each subgroup
        - subgroup_details: Detailed information about each subgroup
        - T_current: Updated tree structure
        - M_current: Updated mutation matrix
        - root_mutations: Updated root mutations
        - external_mutations: Mutations that failed to be added to the tree
        - conflict_mutations: Mutations that caused conflicts during placement
    """
    # Use provided logger or fall back to global logger
    log = logger_obj if logger_obj is not None else logging.getLogger(__name__)
    
    # Store final sub-grouping results
    mutation_subgroups = {}
    subgroup_backbone_mutations = {}
    subgroup_details = {}
    external_mutations = []
    conflict_mutations = []
    
    # Define global parameters
    ω_NA = params.get('general_weight_NA', 0.001)
    fnfp_ratio = params.get('fnfp_ratio', 0.1)
    φ = params.get('phi', 1.0)
    
    log.info(f"Starting sub-grouping within backbone groups for {len(backbones_of_group)} groups")
    
    # Determine whether to show progress (only in interactive environments)
    use_progress = show_progress and sys.stdout.isatty()
    
    group_iterator = backbones_of_group.items()
    if use_progress:
        group_iterator = tqdm(backbones_of_group.items(), desc="Processing backbone groups")
    
    for group_id, backbone_mut in group_iterator:
        
        dir_subgroup = os.path.join(outputpath, f"backbone_clone_{group_id}.{backbone_mut}")
        os.makedirs(dir_subgroup, exist_ok=True)
        
        # Get all mutations in this group (excluding the backbone mutation)
        group_muts = [mut for mut, gid in mutation_group.items()
                     if gid == group_id and mut != backbone_mut]
        
        # Generate subgroup data matrix
        I_group = sorted_I_resolved.loc[list(M_B[M_B[backbone_mut]==1].index), group_muts]
        df_reads_group = df_reads_resolved.loc[['bulk'] + list(M_B[M_B[backbone_mut]==1].index), group_muts]
        
        if len(group_muts) < 1:
            # condition_0: This backbone group has only the backbone mutation, no other mutations to attach
            log.debug(f"Group {group_id}: No additional mutations to process")
            continue
        
        elif len(group_muts) == 1:
            # condition_1: Only one mutation, attach it directly
            log.debug(f"Group {group_id}: Single mutation {group_muts[0]}, attaching directly")
            
            mutation_subgroups[group_muts[0]] = f"{group_id}_0"
            subgroup_backbone_mutations[f"{group_id}_0"] = group_muts[0]
            subgroup_details[f"{group_id}_0"] = {
                'main_group': group_id,
                'backbone_mutation': backbone_mut,
                'mutations': group_muts,
                'is_trivial': True
            }
            
            # Attach single mutation directly to tree
            subtree_groups = [group_muts]
            latest_nodes, found_mutations_list = find_latest_hub_node(T_current, [backbone_mut], logger_obj=log)
            target_node_names = [node.name for node in latest_nodes]
            external_mutations_temp, conflict_mutations_temp, T_current, M_current, root_mutations = process_subtree_mutations_to_specific_node(
                subtree_groups, target_node_names, T_current, M_current, sorted_I_resolved, sorted_P_resolved, 
                ω_NA, fnfp_ratio, φ, log, root_mutations, mutation_group=mutation_group, show_progress=False
            )            
            external_mutations.extend(external_mutations_temp)
            conflict_mutations.extend(conflict_mutations_temp)
            continue
        
        # Use all mutations directly, no CV filtering
        I_group_all = I_group[group_muts].copy()
        
        # Further filter out mutations that are all 0 or all NA
        valid_mutations = []
        for mut in group_muts:
            mut_data = I_group_all[mut]
            if (mut_data == 1).any() and (mut_data == 0).any():
                valid_mutations.append(mut)
        
        # Add invalid mutations to external_mutations
        invalid_mutations = [mut for mut in group_muts if mut not in valid_mutations]
        external_mutations.extend(invalid_mutations)
        
        I_group_final = I_group_all[valid_mutations]
        
        if I_group_final.shape[1] < 2:
            # condition_2: Not enough mutations for graph, attach by MAF order
            log.debug(f"Group {group_id}: Not enough mutations ({I_group_final.shape[1]}) for graph, attaching individually")
            
            for i, mut in enumerate(valid_mutations):
                subgroup_id = f"{group_id}_{i}"
                mutation_subgroups[mut] = subgroup_id
                subgroup_backbone_mutations[subgroup_id] = mut
                subgroup_details[subgroup_id] = {
                    'main_group': group_id,
                    'backbone_mutation': backbone_mut,
                    'mutations': [mut],
                    'is_trivial': True
                }
            
            if valid_mutations:
                subtree_groups = [[mut] for mut in valid_mutations]
                latest_nodes, found_mutations_list = find_latest_hub_node(T_current, [backbone_mut], logger_obj=log)
                target_node_names = [node.name for node in latest_nodes]
                external_mutations_temp, conflict_mutations_temp, T_current, M_current, root_mutations = process_subtree_mutations_to_specific_node(
                    subtree_groups, target_node_names, T_current, M_current, sorted_I_resolved, sorted_P_resolved, 
                    ω_NA, fnfp_ratio, φ, log, root_mutations, mutation_group=mutation_group, show_progress=False
                )
                external_mutations.extend(external_mutations_temp)
                conflict_mutations.extend(conflict_mutations_temp)
            continue
        
        try:
            # Step 1: Calculate correlation weights between mutations
            log.debug(f"Group {group_id}: Computing correlation graph for {len(valid_mutations)} mutations")
            clone_weights_sub, pair_weights_sub = get_correlation_graph_elements(
                I_group_final, 100, 42, cutoff_mcf_for_graph, cutoff_mcn_for_graph, logger_obj=log
            )
            
            if not pair_weights_sub or len(pair_weights_sub) < 2:
                # condition_2: Insufficient graph structure for community detection
                if not pair_weights_sub:
                    log.warning(f"Group {group_id}: pair_weights is empty, all clones are singleton mutations")
                else:
                    log.warning(f"Group {group_id}: insufficient edges ({len(pair_weights_sub)}) for community detection")
                
                for i, mut in enumerate(valid_mutations):
                    subgroup_id = f"{group_id}_{i}"
                    mutation_subgroups[mut] = subgroup_id
                    subgroup_backbone_mutations[subgroup_id] = mut
                    subgroup_details[subgroup_id] = {
                        'main_group': group_id,
                        'backbone_mutation': backbone_mut,
                        'mutations': [mut],
                        'is_trivial': True,
                        'community_detection_skipped': True,
                        'skip_reason': 'insufficient_graph_structure'
                    }
                
                subtree_groups = [[mut] for mut in valid_mutations]
                latest_nodes, found_mutations_list = find_latest_hub_node(T_current, [backbone_mut], logger_obj=log)
                target_node_names = [node.name for node in latest_nodes]
                external_mutations_temp, conflict_mutations_temp, T_current, M_current, root_mutations = process_subtree_mutations_to_specific_node(
                    subtree_groups, target_node_names, T_current, M_current, sorted_I_resolved, sorted_P_resolved, 
                    ω_NA, fnfp_ratio, φ, log, root_mutations, mutation_group=mutation_group, show_progress=False
                )
                external_mutations.extend(external_mutations_temp)
                conflict_mutations.extend(conflict_mutations_temp)
                continue
            
            else:
                # Step 2: Use Leiden algorithm for sub-grouping
                log.debug(f"Group {group_id}: Running Leiden algorithm for community detection")
                mutation_subgroup, partition_sub, G_ig_sub = leiden_mutation_groups(
                    clone_weights_sub, pair_weights_sub, 
                    os.path.join(dir_subgroup, f"{sampleid}.group_{group_id}_subgraph.pdf"), 
                    params.get('resolution_of_graph', 1.0),
                    seed=get_seed(),
                    ordered_mutations=I_group_final.columns.tolist(),
                    logger_obj=log
                )
                external_mutations_outgroup = [i for i in valid_mutations if i not in list(mutation_subgroup.keys())]
                external_mutations.extend(external_mutations_outgroup)
                
                # Step 3: Detect hub clusters in the graph
                hub_clusters, cluster_degrees = detect_hub_clusters(G_ig_sub, mutation_subgroup, logger_obj=log)
                
                if len(hub_clusters) == 0:
                    # condition_3: Can group but no hub detected, attach by MAF order
                    log.debug(f"Group {group_id}: No hub clusters detected")
                    if len(set(mutation_subgroup.values())) <= 2:
                        sorted_subgroup_mutations_but_backbone = [i for i in sorted_I_resolved.columns if i in list(mutation_subgroup.keys())]
                        external_mutations_temp, conflict_mutations_temp, T_current, M_current = integrate_mutations_to_scaffold_within_group(
                            sorted_attached_mutations=sorted_subgroup_mutations_but_backbone,
                            T_current=T_current,
                            M_current=M_current,
                            I_attached=sorted_I_resolved,
                            P_attached=sorted_P_resolved,
                            mutation_group=mutation_group,
                            ω_NA=ω_NA,
                            fnfp_ratio=fnfp_ratio,
                            φ=φ,
                            logger_obj=log,
                            show_progress=False
                        )
                        external_mutations.extend(external_mutations_temp)
                        conflict_mutations.extend(conflict_mutations_temp)
                    else:
                        sorted_parallel_groups = {v: sorted([k for k in mutation_subgroup if mutation_subgroup[k] == v],
                                                           key=lambda x: list(sorted_I_resolved.columns).index(x))
                                                for v in set(mutation_subgroup.values())}
                        subtree_groups = list(sorted_parallel_groups.values())
                        latest_nodes, found_mutations_list = find_latest_hub_node(T_current, [backbone_mut], logger_obj=log)
                        target_node_names = [node.name for node in latest_nodes]
                        external_mutations_temp, conflict_mutations_temp, T_current, M_current, root_mutations = process_subtree_mutations_to_specific_node(
                            subtree_groups, target_node_names, T_current, M_current, sorted_I_resolved, sorted_P_resolved, 
                            ω_NA, fnfp_ratio, φ, log, root_mutations, mutation_group=mutation_group, show_progress=False
                        )
                        external_mutations.extend(external_mutations_temp)
                        conflict_mutations.extend(conflict_mutations_temp)
                
                elif len(hub_clusters) == 1:
                    # condition_4: One hub cluster found, attach hub first, then process remaining
                    log.debug(f"Group {group_id}: Single hub cluster detected")
                    hub_mutations = [i for i, m in mutation_subgroup.items() if m == hub_clusters[0]]
                    sorted_subgroup_mutations_which_hub = [i for i in sorted_I_resolved.columns if i in hub_mutations]
                    external_mutations_temp, conflict_mutations_temp, T_current, M_current = integrate_mutations_to_scaffold_within_group(
                        sorted_attached_mutations=sorted_subgroup_mutations_which_hub,
                        T_current=T_current,
                        M_current=M_current,
                        I_attached=sorted_I_resolved,
                        P_attached=sorted_P_resolved,
                        mutation_group=mutation_group,
                        ω_NA=ω_NA,
                        fnfp_ratio=fnfp_ratio,
                        φ=φ,
                        logger_obj=log,
                        show_progress=False
                    )
                    external_mutations.extend(external_mutations_temp)
                    conflict_mutations.extend(conflict_mutations_temp)
                    
                    I_group_dehub = I_group_final.drop(columns=hub_mutations, errors='ignore').copy()
                    
                    if I_group_dehub.shape[1] >= 2:
                        clone_weights_dehub, pair_weights_dehub = get_correlation_graph_elements(
                            I_group_dehub, 100, 42, cutoff_mcf_for_graph, cutoff_mcn_for_graph, logger_obj=log
                        )
                        
                        if pair_weights_dehub and len(pair_weights_dehub) >= 2:
                            mutation_dehubgroup, partition_dehub, G_ig_dehub = leiden_mutation_groups(
                                clone_weights_dehub, pair_weights_dehub, 
                                os.path.join(dir_subgroup, f"{sampleid}.group_{group_id}_for_group_exclude_hub_cluster.pdf"), 
                                params.get('resolution_of_graph', 1.0),
                                seed=get_seed(),
                                ordered_mutations=I_group_dehub.columns.tolist(),
                                logger_obj=log
                            )
                            hub_clusters_dehub, cluster_degrees_dehub = detect_hub_clusters(G_ig_dehub, mutation_dehubgroup, logger_obj=log)
                            
                            if all(np.array(list(cluster_degrees_dehub.values())) == 0):
                                sorted_parallel_groups = {v: sorted([k for k in mutation_dehubgroup if mutation_dehubgroup[k] == v],
                                                                   key=lambda x: list(sorted_I_resolved.columns).index(x))
                                                        for v in set(mutation_dehubgroup.values())}
                                subtree_groups = list(sorted_parallel_groups.values())
                                latest_nodes, found_mutations_list = find_latest_hub_node(T_current, hub_mutations, logger_obj=log)
                                target_node_names = [node.name for node in latest_nodes]
                                external_mutations_temp, conflict_mutations_temp, T_current, M_current, root_mutations = process_subtree_mutations_to_specific_node(
                                    subtree_groups, target_node_names, T_current, M_current, sorted_I_resolved, sorted_P_resolved, 
                                    ω_NA, fnfp_ratio, φ, log, root_mutations, mutation_group=mutation_group, show_progress=False
                                )
                                external_mutations.extend(external_mutations_temp)
                                conflict_mutations.extend(conflict_mutations_temp)
                            else:
                                external_mutations.extend(list(I_group_dehub.columns))
                
                else:
                    # condition_5: Multiple independent clusters, attach in parallel
                    log.debug(f"Group {group_id}: Multiple hub clusters detected ({len(hub_clusters)})")
                    if all(np.array(list(cluster_degrees.values())) == 0):
                        sorted_parallel_groups = {v: sorted([k for k in mutation_subgroup if mutation_subgroup[k] == v],
                                                           key=lambda x: list(sorted_I_resolved.columns).index(x))
                                                for v in set(mutation_subgroup.values())}
                        
                        group_fractions = calculate_group_cooccurrence_fraction(
                            sorted_parallel_groups, backbone_mut, I_group_final, logger_obj=log
                        )
                        max_key, max_value = find_max_fraction_group(group_fractions, logger_obj=log)
                        
                        if max_value > 0.9:
                            log.debug(f"Group {group_id}: Group {max_key} has high co-occurrence ({max_value:.3f})")
                            I_group_maxsub = I_group_final[sorted_parallel_groups[max_key]].copy()
                            
                            if I_group_maxsub.shape[1] >= 2:
                                clone_weights_maxsub, pair_weights_maxsub = get_correlation_graph_elements(
                                    I_group_maxsub, 100, 42, cutoff_mcf_for_graph, cutoff_mcn_for_graph, logger_obj=log
                                )
                                
                                if pair_weights_maxsub and len(pair_weights_maxsub) >= 2:
                                    mutation_maxsubgroup, partition_maxsub, G_ig_maxsub = leiden_mutation_groups(
                                        clone_weights_maxsub, pair_weights_maxsub, 
                                        os.path.join(dir_subgroup, f"{sampleid}.group_{group_id}_for_more90percent_group.pdf"), 
                                        params.get('resolution_of_graph', 1.0),
                                        seed=get_seed(),
                                        ordered_mutations=I_group_maxsub.columns.tolist(),
                                        logger_obj=log
                                    )
                                    hub_clusters_maxsub, cluster_degrees_maxsub = detect_hub_clusters(G_ig_maxsub, mutation_maxsubgroup, logger_obj=log)
                                    
                                    if all(np.array(list(cluster_degrees_maxsub.values())) == 0) or len(cluster_degrees_maxsub) == 1:
                                        sorted_parallel_groups_inner = {v: sorted([k for k in mutation_maxsubgroup if mutation_maxsubgroup[k] == v],
                                                                                  key=lambda x: list(sorted_I_resolved.columns).index(x))
                                                                      for v in set(mutation_maxsubgroup.values())}
                                        subtree_groups = list(sorted_parallel_groups_inner.values())
                                        latest_nodes, found_mutations_list = find_latest_hub_node(T_current, [backbone_mut], logger_obj=log)
                                        target_node_names = [node.name for node in latest_nodes]
                                        external_mutations_temp, conflict_mutations_temp, T_current, M_current, root_mutations = process_subtree_mutations_to_specific_node(
                                            subtree_groups, target_node_names, T_current, M_current, sorted_I_resolved, sorted_P_resolved, 
                                            ω_NA, fnfp_ratio, φ, log, root_mutations, mutation_group=mutation_group, show_progress=False
                                        )
                                        external_mutations.extend(external_mutations_temp)
                                        conflict_mutations.extend(conflict_mutations_temp)
                        
                        else:
                            log.debug(f"Group {group_id}: No group with high co-occurrence, attaching in parallel")
                            subtree_groups = list(sorted_parallel_groups.values())
                            latest_nodes, found_mutations_list = find_latest_hub_node(T_current, [backbone_mut], logger_obj=log)
                            target_node_names = [node.name for node in latest_nodes]
                            external_mutations_temp, conflict_mutations_temp, T_current, M_current, root_mutations = process_subtree_mutations_to_specific_node(
                                subtree_groups, target_node_names, T_current, M_current, sorted_I_resolved, sorted_P_resolved, 
                                ω_NA, fnfp_ratio, φ, log, root_mutations, mutation_group=mutation_group, show_progress=False
                            )
                            external_mutations.extend(external_mutations_temp)
                            conflict_mutations.extend(conflict_mutations_temp)
            
            # Step 4: Sort subgroup matrix
            log.debug(f"Group {group_id}: Sorting subgroup matrix")
            I_group_sorted, mut_df_sorted_sub, subgroup_to_muts, final_order_sub = sort_I_hierarchical_freeze_ones_fixed(
                I_group_final, mutation_subgroup, logger_obj=log
            )
            
            # Step 5: Select founder mutation for each subgroup
            subgroup_backbones = select_founder_mutations(I_group_sorted, mutation_subgroup, logger_obj=log)
            
            # Store subgroup results
            for mut, subgroup_id in mutation_subgroup.items():
                full_subgroup_id = f"{group_id}_{subgroup_id}"
                mutation_subgroups[mut] = full_subgroup_id
            
            # Store subgroup backbone mutations
            for subgroup_id, sub_backbone in subgroup_backbones.items():
                full_subgroup_id = f"{group_id}_{subgroup_id}"
                subgroup_backbone_mutations[full_subgroup_id] = sub_backbone
            
            # Store subgroup details
            for subgroup_id in set(mutation_subgroup.values()):
                full_subgroup_id = f"{group_id}_{subgroup_id}"
                subgroup_muts = [mut for mut, sid in mutation_subgroup.items() if sid == subgroup_id]
                subgroup_details[full_subgroup_id] = {
                    'main_group': group_id,
                    'backbone_mutation': backbone_mut,
                    'subgroup_backbone': subgroup_backbones.get(subgroup_id),
                    'mutations': subgroup_muts,
                    'is_trivial': False
                }
        
        except Exception as e:
            log.error(f"Error in sub-grouping for group {group_id}: {str(e)}")
            # If sub-grouping fails, create independent subgroups for each mutation
            for i, mut in enumerate(valid_mutations):
                subgroup_id = f"{group_id}_{i}"
                mutation_subgroups[mut] = subgroup_id
                subgroup_backbone_mutations[subgroup_id] = mut
                subgroup_details[subgroup_id] = {
                    'main_group': group_id,
                    'backbone_mutation': backbone_mut,
                    'mutations': [mut],
                    'is_trivial': True
                }
            
            if valid_mutations:
                subtree_groups = [[mut] for mut in valid_mutations]
                latest_nodes, found_mutations_list = find_latest_hub_node(T_current, [backbone_mut], logger_obj=log)
                target_node_names = [node.name for node in latest_nodes]
                external_mutations_temp, conflict_mutations_temp, T_current, M_current, root_mutations = process_subtree_mutations_to_specific_node(
                    subtree_groups, target_node_names, T_current, M_current, sorted_I_resolved, sorted_P_resolved, 
                    ω_NA, fnfp_ratio, φ, log, root_mutations, mutation_group=mutation_group, show_progress=False
                )
                external_mutations.extend(external_mutations_temp)
                conflict_mutations.extend(conflict_mutations_temp)
    
    # Create complete mutation hierarchy
    complete_mutation_hierarchy = {}
    
    # Add backbone mutations to hierarchy
    for group_id, backbone_mut in backbones_of_group.items():
        complete_mutation_hierarchy[backbone_mut] = {
            'main_group': group_id,
            'sub_group': f"{group_id}_backbone",
            'is_backbone': True,
            'is_sub_backbone': False
        }
    
    # Add subgroup backbone mutations
    for full_subgroup_id, sub_backbone in subgroup_backbone_mutations.items():
        main_group_id = int(full_subgroup_id.split('_')[0])
        complete_mutation_hierarchy[sub_backbone] = {
            'main_group': main_group_id,
            'sub_group': full_subgroup_id,
            'is_backbone': False,
            'is_sub_backbone': True
        }
    
    # Add other mutations
    for mut, full_subgroup_id in mutation_subgroups.items():
        if mut not in complete_mutation_hierarchy:
            main_group_id = int(full_subgroup_id.split('_')[0])
            complete_mutation_hierarchy[mut] = {
                'main_group': main_group_id,
                'sub_group': full_subgroup_id,
                'is_backbone': False,
                'is_sub_backbone': False
            }
    
    # Save grouping results
    df_mutation_hierarchy = pd.DataFrame.from_dict(complete_mutation_hierarchy, orient='index')
    df_mutation_hierarchy.reset_index(inplace=True)
    df_mutation_hierarchy.rename(columns={'index': 'mutation'}, inplace=True)
    df_mutation_hierarchy.to_csv(os.path.join(outputpath, f"{sampleid}.mutation_hierarchy.csv"), index=False)
    
    log.info(f"Sub-grouping completed: {len(complete_mutation_hierarchy)} mutations processed")
    log.info(f"  External mutations: {len(external_mutations)}")
    log.info(f"  Conflict mutations: {len(conflict_mutations)}")
    
    return complete_mutation_hierarchy, subgroup_backbone_mutations, subgroup_details, T_current, M_current, root_mutations, external_mutations, conflict_mutations




def generate_new_leaf_on_root_scaffold(T_current: TreeNode, new_mut: str, logger_obj=None):
    """
    Generate a new leaf position attached to the ROOT node.
    
    Parameters
    ----------
    T_current : TreeNode
        Current tree structure
    new_mut : str
        New mutation to add
    logger_obj : logging.Logger, optional
        Logger instance for logging messages. If None, uses global logger.
    
    Returns
    -------
    dict
        Position dictionary for the new leaf
    """
    # Use provided logger or fall back to global logger
    log = logger_obj if logger_obj is not None else logging.getLogger(__name__)
    
    # Deep copy the tree to avoid modifying the original
    new_tree = deepcopy(T_current)
    
    # Find the root node (assumed to be named "ROOT")
    root_node = new_tree.find("ROOT")
    if not root_node:
        raise ValueError("ROOT node not found in the tree")
    
    # Create a new leaf node
    new_leaf = TreeNode(new_mut)
    
    # Add the new leaf node to the root's children
    root_node.add_child(new_leaf)
    
    # Generate the corresponding candidate position
    new_leaf_position = {
        "placement_type": "new_leaf",
        "anchor": "ROOT",
        "meta": {},
        "nodes": [{"name": n.name,
                   "parent": n.parent.name if n.parent else None,
                   "children": [c.name for c in n.children]} for n in new_tree.traverse()],
        "edges": [(n.parent.name, n.name) for n in new_tree.traverse() if n.parent]
    }
    
    log.debug(f"Generated new leaf position for {new_mut} under ROOT")
    
    return new_leaf_position


def check_mutation_has_uncovered_cells(M_current, subtree_mut, I_attached, logger_obj=None):
    """
    Check if a mutation still has cells with value 1 that are not covered by the tree.
    
    Parameters
    ----------
    M_current : pd.DataFrame
        Current mutation matrix
    subtree_mut : str
        Mutation to check
    I_attached : pd.DataFrame
        Attached mutation information
    logger_obj : logging.Logger, optional
        Logger instance for logging messages. If None, uses global logger.
    
    Returns
    -------
    tuple
        (has_uncovered, uncovered_cells_with_mut)
    """
    # Use provided logger or fall back to global logger
    log = logger_obj if logger_obj is not None else logging.getLogger(__name__)
    
    # Get all columns in M_current except ROOT (tree mutation nodes)
    tree_mut_columns = [col for col in M_current.columns if col != 'ROOT']
    
    # Find cells that are all 0 on tree nodes (uncovered cells)
    if len(tree_mut_columns) == 0:
        # If there are no other mutation nodes on the tree, all cells are uncovered
        uncovered_cells = M_current.index.tolist()
    else:
        uncovered_cells = M_current[
            M_current[tree_mut_columns].sum(axis=1) == 0
        ].index.tolist()
    
    if len(uncovered_cells) == 0:
        log.debug(f"No uncovered cells found for mutation {subtree_mut}")
        return False, []
    
    # Check if the current mutation has 1s in these uncovered cells
    mut_values = I_attached.loc[uncovered_cells, subtree_mut]
    cells_with_mut = mut_values[mut_values == 1.0].index.tolist()
    
    has_uncovered = len(cells_with_mut) > 0
    log.debug(f"Mutation {subtree_mut}: {len(cells_with_mut)} uncovered cells with mutation")
    
    return has_uncovered, cells_with_mut


def process_misassigned_mutations_direct_to_root(
    subtree_groups, T_current, M_current, I_attached, P_attached, 
    ω_NA, fnfp_ratio, φ, logger_obj=None, root_mutations=None, max_retries=None,
    show_progress=False
):
    """
    Process external mutations through subtree groups, supporting both multi-mutation groups
    and single-mutation groups separately (with rollback and retry support).
    
    Parameters
    ----------
    subtree_groups : list
        List of subtree groups
    T_current : TreeNode
        Current tree structure
    M_current : pd.DataFrame
        Current mutation matrix
    I_attached : pd.DataFrame
        Attached mutation information
    P_attached : pd.DataFrame
        Attached probability information
    ω_NA : float
        NA weight parameter
    fnfp_ratio : float
        False negative to false positive ratio
    φ : float
        Bayesian penalty parameter
    logger_obj : logging.Logger, optional
        Logger instance for logging messages. If None, uses global logger.
    root_mutations : list, optional
        List of root mutations
    max_retries : int, optional
        Maximum number of candidates to try
    show_progress : bool, default=False
        Whether to show progress bars
    
    Returns
    -------
    tuple
        (remained_mutations, conflict_mutations, T_current, M_current, root_mutations)
    """
    # Use provided logger or fall back to global logger
    log = logger_obj if logger_obj is not None else logging.getLogger(__name__)
    
    if root_mutations is None:
        root_mutations = []
    
    remained_mutations = []
    conflict_mutations = []
    
    # Separate subtree groups into multi-mutation and singleton groups
    multi_mut_subtree_groups = [g for g in subtree_groups if len(g) > 1]
    singleton_subtree_groups = [g for g in subtree_groups if len(g) == 1]
    
    # Determine whether to show progress (only in interactive environments)
    use_progress = show_progress and sys.stdout.isatty()
    
    # ============================================================
    # Process multi-mutation groups (length > 1)
    # ============================================================
    if use_progress:
        group_iterator = tqdm(multi_mut_subtree_groups, desc="Processing multiple subtrees")
    else:
        group_iterator = multi_mut_subtree_groups
    
    for group_idx, group in enumerate(group_iterator):
        
        # Sort by number of 1s in I_attached (descending)
        sorted_group = sorted(group, key=lambda subtree_mut: I_attached[subtree_mut].sum(), reverse=True)
        
        # Add mutations one by one, first attached to ROOT
        reattached_mutations = []
        
        if use_progress:
            mut_iterator = tqdm(sorted_group, desc="Processing mutations in group")
        else:
            mut_iterator = sorted_group
        
        for idx, subtree_mut in enumerate(mut_iterator):
            T_rollback = copy.deepcopy(T_current)
            M_rollback = M_current.copy()
            
            # ---- Check: does this mutation still have 1s in uncovered cells? ----
            has_uncovered, uncovered_cells_with_mut = check_mutation_has_uncovered_cells(
                M_current, subtree_mut, I_attached, logger_obj=log
            )
            
            if not has_uncovered:
                # All 1 cells are already covered by tree mutations
                if subtree_mut not in remained_mutations:
                    remained_mutations.append(subtree_mut)
                log.info(f"Mutation {subtree_mut} added to remained_mutations (no uncovered cells with mutation)")
                continue  # Skip this mutation
            
            log.debug(f"Mutation {subtree_mut} has {len(uncovered_cells_with_mut)} uncovered cells with mutation: {uncovered_cells_with_mut[:5]}...")
            
            if idx == 0:
                # First mutation: find intersection nodes and place it
                intersection_nodes = find_all_intersect_muts_from_tree_by_matrix_scaffold(
                    T_current, I_attached, subtree_mut, logger_obj=log
                )
                
                if len(intersection_nodes) == 0:
                    # No intersection, attach directly to ROOT
                    root_new_leaf_position = generate_new_leaf_on_root_scaffold(T_current, subtree_mut, logger_obj=log)
                    selected_positions = [root_new_leaf_position]
                else:
                    # Get candidate positions
                    potential_positions = find_intersection_positions_within_tree_directly_scaffold(
                        T_current, subtree_mut, I_attached, min_overlap=1, intersection_nodes=intersection_nodes, logger_obj=log
                    )
                    root_new_leaf_position = generate_new_leaf_on_root_scaffold(T_current, subtree_mut, logger_obj=log)
                    selected_positions = potential_positions + [root_new_leaf_position]
                
                parent_dict = build_parent_dict_from_candidates_scaffold(selected_positions)
                
                # Backup current state
                M_backup = M_current.copy()
                T_backup = T_current.copy()
                
                # Calculate penalties for all candidate positions
                df_penalty = compute_bayesian_penalty_for_all_positions_scaffold(
                    subtree_mut, selected_positions, T_current, M_current, I_attached, P_attached, 
                    parent_dict, intersection_nodes, ω_NA=ω_NA, fnfp_ratio=fnfp_ratio, φ=φ,
                    logger_obj=log
                )
                
                if df_penalty.empty:
                    reattached_mutations.append(subtree_mut)
                    log.warning(f"Mutation {subtree_mut} added to reattached_mutations (no valid penalty scores)")
                    continue
                
                df_valid = df_penalty[df_penalty['imputed_vec'].apply(lambda x: x.sum() > 0)]
                if df_valid.empty:
                    reattached_mutations.append(subtree_mut)
                    log.warning(f"Mutation {subtree_mut} added to reattached_mutations (all imputed_vec are zero)")
                    continue
                
                df_sorted = df_valid.sort_values('total_penalty')
                
                # Check if best candidate is on ROOT
                best_row = df_sorted.iloc[0]
                best_position = best_row['position']
                
                if (best_position.get('anchor') == 'ROOT' and 
                    best_position.get('placement_type') == 'on_node'):
                    # This mutation should be added directly to root_mutations
                    if subtree_mut not in root_mutations:
                        root_mutations.append(subtree_mut)
                    log.info(f"Mutation {subtree_mut} added to root_mutations (best position is on ROOT)")
                    continue
                
                if max_retries is None:
                    candidates_to_try = df_sorted
                else:
                    candidates_to_try = df_sorted.head(max_retries)
                
                # Try candidate positions
                placed = False
                for attempt, (idx_row, row) in enumerate(candidates_to_try.iterrows()):
                    
                    M_current = M_backup.copy()
                    T_current = T_backup.copy()
                    
                    try:
                        T_current, M_current = apply_position_to_tree_scaffold(
                            subtree_mut, row['position'], row['imputed_vec'], 
                            T_current, M_current, I_attached, parent_dict, logger_obj=log
                        )
                        
                        if scp.ul.is_conflict_free_gusfield(M_current):
                            placed = True
                            log.debug(f"Mutation {subtree_mut} successfully placed at position {row['position_index']}")
                            break
                        else:
                            log.warning(f"✗ Position {row['position_index']} caused conflict, trying next candidate")
                            
                    except Exception as e:
                        log.error(f"Error placing mutation at position {row['position_index']}: {e}")
                        continue
                
                if not placed:
                    M_current = M_backup.copy()
                    T_current = T_backup.copy()
                    reattached_mutations.append(subtree_mut)
                    log.warning(f"Mutation {subtree_mut} added to reattached_mutations (all candidates failed)")
            
            else:
                # Subsequent mutations: attach to the growing group
                # Find intersection nodes
                intersection_nodes = find_all_intersect_muts_from_tree_by_matrix_scaffold(
                    T_current, I_attached, subtree_mut, logger_obj=log
                )
                
                if len(intersection_nodes) == 0:
                    root_new_leaf_position = generate_new_leaf_on_root_scaffold(T_current, subtree_mut, logger_obj=log)
                    selected_positions = [root_new_leaf_position]
                else:
                    potential_positions = find_intersection_positions_within_tree_directly_scaffold(
                        T_current, subtree_mut, I_attached, min_overlap=1, intersection_nodes=intersection_nodes, logger_obj=log
                    )
                    root_new_leaf_position = generate_new_leaf_on_root_scaffold(T_current, subtree_mut, logger_obj=log)
                    selected_positions = potential_positions + [root_new_leaf_position]
                
                parent_dict = build_parent_dict_from_candidates_scaffold(selected_positions)
                
                M_backup = M_current.copy()
                T_backup = T_current.copy()
                
                df_penalty = compute_bayesian_penalty_for_all_positions_scaffold(
                    subtree_mut, selected_positions, T_current, M_current, I_attached, P_attached, 
                    parent_dict, intersection_nodes, ω_NA=ω_NA, fnfp_ratio=fnfp_ratio, φ=φ,
                    logger_obj=log
                )
                
                if df_penalty.empty:
                    reattached_mutations.append(subtree_mut)
                    log.warning(f"Mutation {subtree_mut} added to reattached_mutations (no valid penalty scores)")
                    continue
                
                df_valid = df_penalty[df_penalty['imputed_vec'].apply(lambda x: x.sum() > 0)]
                if df_valid.empty:
                    reattached_mutations.append(subtree_mut)
                    log.warning(f"Mutation {subtree_mut} added to reattached_mutations (all imputed_vec are zero)")
                    continue
                
                df_sorted = df_valid.sort_values('total_penalty')
                
                best_row = df_sorted.iloc[0]
                best_position = best_row['position']
                
                if (best_position.get('anchor') == 'ROOT' and 
                    best_position.get('placement_type') == 'on_node'):
                    if subtree_mut not in root_mutations:
                        root_mutations.append(subtree_mut)
                    log.info(f"Mutation {subtree_mut} added to root_mutations (best position is on ROOT)")
                    continue
                
                if max_retries is None:
                    candidates_to_try = df_sorted
                else:
                    candidates_to_try = df_sorted.head(max_retries)
                
                placed = False
                for attempt, (idx_row, row) in enumerate(candidates_to_try.iterrows()):
                    
                    M_current = M_backup.copy()
                    T_current = T_backup.copy()
                    
                    try:
                        T_current, M_current = apply_position_to_tree_scaffold(
                            subtree_mut, row['position'], row['imputed_vec'], 
                            T_current, M_current, I_attached, parent_dict, logger_obj=log
                        )
                        
                        if scp.ul.is_conflict_free_gusfield(M_current):
                            placed = True
                            log.debug(f"Mutation {subtree_mut} successfully placed at position {row['position_index']}")
                            break
                        else:
                            log.warning(f"✗ Position {row['position_index']} caused conflict, trying next candidate")
                            
                    except Exception as e:
                        log.error(f"Error placing mutation at position {row['position_index']}: {e}")
                        continue
                
                if not placed:
                    M_current = M_backup.copy()
                    T_current = T_backup.copy()
                    reattached_mutations.append(subtree_mut)
                    log.warning(f"Mutation {subtree_mut} added to reattached_mutations (all candidates failed)")
        
        # Handle reattached mutations
        if reattached_mutations:
            for subtree_mut in reattached_mutations:
                # Check if mutation has uncovered cells
                has_uncovered, uncovered_cells_with_mut = check_mutation_has_uncovered_cells(
                    M_current, subtree_mut, I_attached, logger_obj=log
                )
                
                if not has_uncovered:
                    if subtree_mut not in remained_mutations:
                        remained_mutations.append(subtree_mut)
                    log.info(f"Mutation {subtree_mut} added to remained_mutations (no uncovered cells with mutation)")
                    continue
                
                # Try to place as new leaf under ROOT
                root_new_leaf_position = generate_new_leaf_on_root_scaffold(T_current, subtree_mut, logger_obj=log)
                selected_positions = [root_new_leaf_position]
                
                parent_dict = build_parent_dict_from_candidates_scaffold(selected_positions)
                intersection_nodes = {'ROOT'}
                
                M_backup = M_current.copy()
                T_backup = T_current.copy()
                
                df_penalty = compute_bayesian_penalty_for_all_positions_scaffold(
                    subtree_mut, selected_positions, T_current, M_current, I_attached, P_attached, 
                    parent_dict, intersection_nodes, ω_NA=ω_NA, fnfp_ratio=fnfp_ratio, φ=φ,
                    logger_obj=log
                )
                
                if df_penalty.empty or df_penalty[df_penalty['imputed_vec'].apply(lambda x: x.sum() > 0)].empty:
                    reattached_mutations.append(subtree_mut)
                    log.warning(f"Mutation {subtree_mut} could not be placed, adding to reattached_mutations")
                    continue
                
                df_valid = df_penalty[df_penalty['imputed_vec'].apply(lambda x: x.sum() > 0)]
                df_sorted = df_valid.sort_values('total_penalty')
                best_row = df_sorted.iloc[0]
                
                if (best_row['position'].get('anchor') == 'ROOT' and 
                    best_row['position'].get('placement_type') == 'on_node'):
                    if subtree_mut not in root_mutations:
                        root_mutations.append(subtree_mut)
                    log.info(f"Mutation {subtree_mut} added to root_mutations")
                    continue
                
                try:
                    T_current, M_current = apply_position_to_tree_scaffold(
                        subtree_mut, best_row['position'], best_row['imputed_vec'], 
                        T_current, M_current, I_attached, parent_dict, logger_obj=log
                    )
                    log.info(f"Mutation {subtree_mut} successfully placed as reattached mutation")
                except Exception as e:
                    log.error(f"Failed to place mutation {subtree_mut} as reattached: {e}")
                    M_current = M_backup.copy()
                    T_current = T_backup.copy()
                    conflict_mutations.append(subtree_mut)
                    log.warning(f"Mutation {subtree_mut} added to conflict_mutations")
    
    # ============================================================
    # Process singleton mutation groups (length = 1)
    # ============================================================
    if use_progress:
        singleton_iterator = tqdm(singleton_subtree_groups, desc="Processing singleton subtrees")
    else:
        singleton_iterator = singleton_subtree_groups
    
    for group in singleton_iterator:
        T_rollback = copy.deepcopy(T_current)
        M_rollback = M_current.copy()
        
        subtree_mut = group[0]
        
        # ---- Check: does this mutation still have 1s in uncovered cells? ----
        has_uncovered, uncovered_cells_with_mut = check_mutation_has_uncovered_cells(
            M_current, subtree_mut, I_attached, logger_obj=log
        )
        
        if not has_uncovered:
            # All 1 cells are already covered by tree mutations
            if subtree_mut not in remained_mutations:
                remained_mutations.append(subtree_mut)
            log.info(f"Mutation {subtree_mut} added to remained_mutations (no uncovered cells with mutation)")
            continue  # Skip this mutation
        
        log.debug(f"Mutation {subtree_mut} has {len(uncovered_cells_with_mut)} uncovered cells with mutation: {uncovered_cells_with_mut[:5]}...")
        
        # Find intersection nodes
        intersection_nodes = find_all_intersect_muts_from_tree_by_matrix_scaffold(
            T_current, I_attached, subtree_mut, logger_obj=log
        )
        
        if len(intersection_nodes) == 0:
            intersection_nodes = {'ROOT'}
        
        # Get candidate positions using optimized method
        potential_positions = find_intersection_positions_within_tree_directly_scaffold(
            T_current, subtree_mut, I_attached, min_overlap=1, intersection_nodes=intersection_nodes, logger_obj=log
        )
        root_new_leaf_position = generate_new_leaf_on_root_scaffold(T_current, subtree_mut, logger_obj=log)
        parent_dict = build_parent_dict_from_candidates_scaffold(potential_positions)
        selected_positions = [p for p in potential_positions if p['anchor'] == 'ROOT'] + [root_new_leaf_position]
        
        # ---- Backup current state ----
        M_backup = M_current.copy()
        T_backup = T_current.copy()
        
        # ---- Calculate penalties for all candidate positions ----
        df_penalty = compute_bayesian_penalty_for_all_positions_scaffold(
            subtree_mut, selected_positions, T_current, M_current, I_attached, P_attached, 
            parent_dict, intersection_nodes, ω_NA=ω_NA, fnfp_ratio=fnfp_ratio, φ=φ,
            logger_obj=log
        )
        
        if df_penalty.empty:
            conflict_mutations.append(subtree_mut)
            log.warning(f"Mutation {subtree_mut} added to conflict_mutations (no valid penalty scores)")
            continue
        
        df_valid = df_penalty[df_penalty['imputed_vec'].apply(lambda x: x.sum() > 0)]
        if df_valid.empty:
            conflict_mutations.append(subtree_mut)
            log.warning(f"Mutation {subtree_mut} added to conflict_mutations (all imputed_vec are zero)")
            continue
        
        df_sorted = df_valid.sort_values('total_penalty')
        
        # ---- Check if the best candidate position is on ROOT ----
        best_row = df_sorted.iloc[0]
        best_position = best_row['position']
        
        # Check if it's an on_node placement on ROOT
        if (best_position.get('anchor') == 'ROOT' and 
            best_position.get('placement_type') == 'on_node'):
            # This mutation should be added directly to root_mutations
            if subtree_mut not in root_mutations:
                root_mutations.append(subtree_mut)
            log.info(f"Mutation {subtree_mut} added to root_mutations (best position is on ROOT)")
            continue  # Skip this mutation, don't add to tree
        
        if max_retries is None:
            candidates_to_try = df_sorted
        else:
            candidates_to_try = df_sorted.head(max_retries)
        
        # ---- Try candidate positions ----
        placed = False
        for attempt, (idx_row, row) in enumerate(candidates_to_try.iterrows()):
            
            M_current = M_backup.copy()
            T_current = T_backup.copy()
            
            try:
                T_current, M_current = apply_position_to_tree_scaffold(
                    subtree_mut, row['position'], row['imputed_vec'], 
                    T_current, M_current, I_attached, parent_dict, logger_obj=log
                )
                
                if scp.ul.is_conflict_free_gusfield(M_current):
                    placed = True
                    log.debug(f"Mutation {subtree_mut} successfully placed at position {row['position_index']}")
                    break
                else:
                    log.warning(f"✗ Position {row['position_index']} caused conflict, trying next candidate")
                    
            except Exception as e:
                log.error(f"Error placing mutation at position {row['position_index']}: {e}")
                continue
        
        if not placed:
            M_current = M_backup.copy()
            T_current = T_backup.copy()
            conflict_mutations.append(subtree_mut)
            log.warning(f"Mutation {subtree_mut} added to conflict_mutations (all candidates failed)")
    
    log.info(f"Processed {len(subtree_groups)} subtree groups: {len(multi_mut_subtree_groups)} multi-mutation, {len(singleton_subtree_groups)} singleton")
    log.info(f"Remained mutations: {len(remained_mutations)}, Conflict mutations: {len(conflict_mutations)}")
    
    return remained_mutations, conflict_mutations, T_current, M_current, root_mutations


def cluster_external_mutations_by_intersection_scaffold(I_selected, external_mutations, min_shared=1, logger_obj=None):
    """
    Cluster external mutations into subtree groups based on intersection,
    and provide a reasonable addition order within each group.
    
    Parameters
    ----------
    I_selected : pd.DataFrame
        Mutation × cells matrix (rows are mutations, columns are cells)
    external_mutations : list[str]
        List of mutations to process
    min_shared : int, default=1
        Minimum number of cells where two mutations co-occur (both 1) to be considered intersecting
    logger_obj : logging.Logger, optional
        Logger instance for logging messages. If None, uses global logger.
    
    Returns
    -------
    list[list[str]]
        Each subtree (mutation group) with internal order that can be used directly
        for add_new_mutation_to_tree_independent()
    """
    # Use provided logger or fall back to global logger
    log = logger_obj if logger_obj is not None else logging.getLogger(__name__)
    
    log.debug(f"Clustering {len(external_mutations)} external mutations with min_shared={min_shared}")
    
    # 1. Build intersection graph
    G = nx.Graph()
    G.add_nodes_from(external_mutations)
    
    edge_count = 0
    for i, m1 in enumerate(external_mutations):
        for m2 in external_mutations[i+1:]:
            if m1 not in I_selected.columns or m2 not in I_selected.columns:
                continue
            v1 = I_selected[m1].fillna(0).astype(int)
            v2 = I_selected[m2].fillna(0).astype(int)
            inter = (v1 & v2).sum()
            if inter >= min_shared:
                G.add_edge(m1, m2)
                edge_count += 1
    
    log.debug(f"Built intersection graph with {len(external_mutations)} nodes and {edge_count} edges")
    
    # 2. Find connected components (each component is a subtree)
    subtree_groups = []
    components = list(nx.connected_components(G))
    log.debug(f"Found {len(components)} connected components")
    
    for comp in components:
        group = list(comp)
        
        # 3. Build directed graph within each group (who is more upstream)
        # Use inclusion relationships or intersection size to determine direction
        DG = nx.DiGraph()
        DG.add_nodes_from(group)
        
        for i, m1 in enumerate(group):
            for m2 in group[i+1:]:
                v1 = I_selected[m1].fillna(0).astype(int)
                v2 = I_selected[m2].fillna(0).astype(int)
                
                # Inclusion: if v1 is a superset of v2, v1 should be above v2
                if (v1 >= v2).all() and (v1 > v2).any():
                    DG.add_edge(m1, m2)
                elif (v2 >= v1).all() and (v2 > v1).any():
                    DG.add_edge(m2, m1)
                else:
                    # Cross but not inclusion, use intersection size to determine direction
                    inter1 = (v1 & v2).sum()
                    inter2 = inter1  # symmetric
                    if inter1 > 0:
                        # The one with fewer total 1s goes above (more sparse mutation appears earlier)
                        if v1.sum() < v2.sum():
                            DG.add_edge(m1, m2)
                        else:
                            DG.add_edge(m2, m1)
        
        try:
            order = list(nx.topological_sort(DG))
            log.debug(f"Group of size {len(group)}: topological order determined")
        except nx.NetworkXUnfeasible:
            # If there are cyclic relationships (symmetric cross), sort by total 1 count ascending
            order = sorted(group, key=lambda m: I_selected[m].sum())
            log.debug(f"Group of size {len(group)}: cycle detected, sorted by mutation frequency")
        
        subtree_groups.append(order)
    
    log.info(f"Clustered {len(external_mutations)} mutations into {len(subtree_groups)} subtree groups")
    
    return subtree_groups




# -------------------------
# Main function for scaffold building
# -------------------------

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import copy
import logging
from tqdm import tqdm
from typing import Tuple, Dict, Any, Optional

def build_scaffold_tree(
    P_somatic: pd.DataFrame,
    V_somatic: pd.DataFrame,
    A_somatic: pd.DataFrame,
    C_somatic: pd.DataFrame,
    I_somatic: pd.DataFrame,
    df_reads_somatic: pd.DataFrame,
    df_features_new: pd.DataFrame,
    params: Dict[str, Any],
    is_filter_quality: str,
    outputpath_scaffold: str,
    sampleid: str,
    immune_mutations: list,
    df_celltype: Optional[str] = None,
    logger_obj=None
) -> Tuple[Any, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Main function for building scaffold tree from mutation data.
    
    Parameters
    ----------
    P_somatic : pd.DataFrame
        Probability matrix
    V_somatic : pd.DataFrame  
        Variant count matrix
    A_somatic : pd.DataFrame
        Alternative allele count matrix
    C_somatic : pd.DataFrame
        Coverage matrix
    I_somatic : pd.DataFrame
        Binary mutation matrix
    df_reads_somatic : pd.DataFrame
        Reads dataframe
    params : Dict[str, Any]
        Parameters dictionary
    is_filter_quality : str
        Whether to apply quality filtering ('yes' or 'no')
    outputpath_scaffold : str
        Output directory path
    sampleid : str
        Sample ID
    immune_mutations : list
        List of immune mutations
    df_celltype : pd.DataFrame, optional
        Celltype dataframe
    logger_obj : logging.Logger, optional
        Logger instance for logging messages. If None, uses global logger.
    
    Returns
    -------
    Tuple containing:
        - T_scaffold: Final scaffold tree
        - M_scaffold: Final imputed mutation matrix
        - df_flipping_spots: Dataframe of flipping spots
        - df_total_flipping_count: Dataframe of total flipping counts
        - final_cleaned_I_selected_withNA3: Final cleaned binary matrix with NA=3
        - final_cleaned_M_scaffold: Final cleaned mutation matrix
        - backbone_mutations: List of backbone mutations
        - checked_mutation_group: Validated mutation group dictionary
        - spots_to_split: List of spots that were split
        - group_mutations: List of mutations in groups
        - no_group_mutations: List of mutations not in any group
        - remained_mutations: List of remained mutations
        - conflict_mutations: List of conflict mutations
        - root_mutations: List of root mutations
    """
    
    # Set up logger
    log = logger_obj if logger_obj is not None else logging.getLogger(__name__)
    
    # ===== Helper function: Get ordered mutation list =====
    def get_ordered_mutations(df_values, df_features, min_cell_threshold=30, bin_size=5):
        """Get ordered mutation list from DataFrame (MCN -> MCF -> name)"""
        df_reordered, _ = reorder_columns_by_mutant_stats(
            df_values, df_features, 
            min_cell_threshold=min_cell_threshold, 
            bin_size=bin_size, 
            descending=True, 
            return_stats=True
        )
        return df_reordered.columns.tolist()
    
    # ===== Set random seed =====
    set_seed(params.get('random_seed', 42))
    
    # ------------------------------
    # Step 1: Check Celltype Data
    # ------------------------------
    
    if is_filter_quality == "yes":
        # ------------------------------
        # Step 2: Initial Filtration
        # ------------------------------
        log.info("Performing initial filtration...")
        kept_cells, kept_mutations, P_sub, V_sub, A_sub, C_sub, I_sub = initial_filter(
            P_somatic, V_somatic, A_somatic, C_somatic, I_somatic, params
        )
        
        # ------------------------------
        # Step 3: Coverage-Based Filtration
        # ------------------------------
        log.info("Performing coverage-based filtration...")
        scaffold_mutations, df_summary_filtration = coverage_filters(
            kept_mutations, df_reads_somatic, df_celltype, params, outputpath_scaffold, logger_obj=log
        )
        
        if not scaffold_mutations:
            raise ValueError("No mutations passed coverage filtration after removing germline variants")
        
        # coverage_filters returns scaffold_mutations already sorted (via kept_mutations order)
    else:
        log.info("Skipping quality filtering (is_filter_quality=no)")
        P_sub, V_sub, A_sub, C_sub, I_sub = P_somatic.copy(), V_somatic.copy(), A_somatic.copy(), C_somatic.copy(), I_somatic.copy()
        # ===== Sort scaffold_mutations =====
        scaffold_mutations = get_ordered_mutations(I_sub, df_features_new)
    
    # ------------------------------
    # Step 4: Create Final Scaffold Matrices
    # ------------------------------
    log.info("Creating final scaffold matrices...")
    I_scaffold = I_sub[scaffold_mutations].copy()
    P_scaffold = P_sub[scaffold_mutations].copy()
    V_scaffold = V_sub[scaffold_mutations].copy()
    A_scaffold = A_sub[scaffold_mutations].copy()
    C_scaffold = C_sub[scaffold_mutations].copy()
    df_reads_scaffold = df_reads_somatic.loc[I_scaffold.index.insert(0, 'bulk'), scaffold_mutations].copy()
    
    # ------------------------------
    # Step 5: Resolve Immune Mutations
    # ------------------------------
    log.info("Resolving immune mutations...")
    if len(immune_mutations)==0:
        I_resolved = I_scaffold.copy()
        P_resolved = P_scaffold.copy()
        V_resolved = V_scaffold.copy()
        A_resolved = A_scaffold.copy()
        C_resolved = C_scaffold.copy()
        df_reads_resolved = df_reads_scaffold.copy()
        spots_to_split = []
        I_somatic_resolved = I_somatic.copy()
        P_somatic_resolved = P_somatic.copy()
    else:
        I_resolved, P_resolved, V_resolved, A_resolved, C_resolved, df_reads_resolved, spots_to_split = resolved_spots_by_immune_mutations(
            I_scaffold, immune_mutations, P_scaffold, V_scaffold, A_scaffold, C_scaffold, df_reads_scaffold, p_threshold=0.5, logger_obj=log
        )
        I_somatic_resolved, P_somatic_resolved = split_spots_by_immune_mutations_scaffold(
            spots_to_split, [i for i in immune_mutations if i in I_somatic.columns], I_somatic, P_somatic
        )
    
    I_somatic_resolved_withNA3 = I_somatic_resolved.replace({np.nan: 3}).astype(int)
    I_somatic_resolved_withNA3.to_csv(os.path.join(outputpath_scaffold, "I_somatic_resolved_withNA3.txt"), sep="\t")
    final_cleaned_I_somatic_resolved = I_somatic_resolved.loc[:, (I_somatic_resolved == 1).any(axis=0)]
    final_cleaned_I_somatic_resolved = final_cleaned_I_somatic_resolved.loc[(final_cleaned_I_somatic_resolved == 1).any(axis=1)]
    final_cleaned_I_somatic_resolved_withNA3 = final_cleaned_I_somatic_resolved.replace({np.nan: 3}).astype(int)
    final_cleaned_I_somatic_resolved_withNA3.to_csv(os.path.join(outputpath_scaffold, "final_cleaned_I_somatic_resolved_withNA3.txt"), sep="\t")
    
    # ------------------------------
    # Step 6: Mutation Grouping with Leiden Algorithm
    # ------------------------------
    log.info("Performing mutation grouping with Leiden algorithm...")
    
    # Calculate adaptive resolution
    if 'resolution_of_graph' in params and params['resolution_of_graph'] is not None:
        resolution = params['resolution_of_graph']
        log.info(f"Using user-specified resolution: {resolution}")
    else:
        resolution, metrics = compute_adaptive_resolution(
            I_resolved,
            min_resolution=params.get('min_resolution', 0.5),
            max_resolution=params.get('max_resolution', 2.0),
            logger_obj=log
        )
        log.info(f"Adaptive resolution based on I_resolved: {resolution:.2f}")
        params['resolution_of_graph'] = resolution
    
    # Set cutoff_mcf_for_graph and cutoff_mcn_for_graph
    if is_filter_quality == "yes":
        cutoff_mcf_for_graph = 0.05
        cutoff_mcn_for_graph = 5
    else:
        cutoff_mcf_for_graph = 0
        cutoff_mcn_for_graph = 0
    
    # Calculate correlation for mutations
    clone_weights, pair_weights = get_correlation_graph_elements(
        I_resolved, 100, 42, cutoff_mcf_for_graph, cutoff_mcn_for_graph, logger_obj=log
    )
    
    # ===== Pass sorted scaffold_mutations =====
    mutation_group, partition, G_ig = leiden_mutation_groups(
        clone_weights, pair_weights, 
        outputpath_scaffold + "/" + sampleid + ".graph_for_mut_grouping.pdf", 
        params['resolution_of_graph'], 
        seed=get_seed(), 
        ordered_mutations=scaffold_mutations,
        logger_obj=log
    )
    
    # ===== Reorder mutation_group =====
    sorted_mutation_names = get_ordered_mutations(I_resolved, df_features_new)
    mutation_group = {mut: mutation_group[mut] for mut in sorted_mutation_names if mut in mutation_group}
    
    group_mutations = list(mutation_group.keys())
    no_group_mutations = [i for i in scaffold_mutations if i not in group_mutations]
    log.info(f"There are {len(no_group_mutations)} scaffold mutations lost after Leiden graph grouping...")
    
    I_selected_and_sorted, mut_df_sorted, group_to_muts, final_order = sort_I_hierarchical_freeze_ones_fixed(
        I_resolved, mutation_group, logger_obj=log
    )
    df_celltype_sub = df_celltype[df_celltype['barcode'].isin(I_selected_and_sorted.index)].copy()
    plot_heatmap_with_celltype_by_your_sorting(
        I_selected_and_sorted, df_celltype_sub, mutation_group, 
        list(mut_df_sorted['mutation']), 
        os.path.join(outputpath_scaffold, sampleid + ".heatmap_with_celltype_right_in_I_selected_and_sorted_after_graph_grouping.pdf"),
        logger_obj=log
    )
    
    # ------------------------------
    # Step 7: Select Backbone Mutations
    # ------------------------------
    log.info("Selecting backbone mutations...")
    backbones_of_group = select_founder_mutations(I_selected_and_sorted, mutation_group, logger_obj=log)
    
    # ===== Sort backbone_mutations =====
    backbone_mutations = sorted(list(backbones_of_group.values()))
    # Rebuild backbones_of_group maintaining order consistent with backbone_mutations
    backbone_to_group = {v: k for k, v in backbones_of_group.items()}
    backbones_of_group = {backbone_to_group[mut]: mut for mut in backbone_mutations}
    
    # Check if the group assigned for mutation is correct
    log.info("Validating mutation group assignments...")
    mutation_list_under_backbone_mutations = {
        backbone: [mutation for mutation, value in mutation_group.items() if value == mutation_group[backbone]] 
        for backbone in backbone_mutations
    }
    misgroup_mutations = []
    for check_mut in [i for i in list(mutation_group.keys()) if i not in backbone_mutations]:
        graph_key = next(key for key, mutations in mutation_list_under_backbone_mutations.items() if check_mut in mutations)
        best_backbone, intersection_counts = find_best_backbone_for_new_mutation_scaffold(
            mutation_list_under_backbone_mutations, I_somatic_resolved, I_somatic_resolved, check_mut, logger_obj=log
        )
        if graph_key != best_backbone:
            misgroup_mutations.append(check_mut)
    
    checked_mutation_group = {k: v for k, v in mutation_group.items() if k not in misgroup_mutations}
    
    # ===== Reorder checked_mutation_group =====
    checked_mutation_group = {mut: checked_mutation_group[mut] for mut in sorted_mutation_names if mut in checked_mutation_group}
    
    mutation_list_under_backbone_mutations = {
        backbone: [mutation for mutation, value in checked_mutation_group.items() if value == checked_mutation_group[backbone]] 
        for backbone in backbone_mutations
    }
    
    # ------------------------------
    # Step 8: Build Backbone Tree
    # ------------------------------
    log.info("Building backbone tree...")
    dir_backbone = os.path.join(outputpath_scaffold, "phylo_backbone_tree")
    if not os.path.exists(dir_backbone):
        os.makedirs(dir_backbone)
    
    T_B = build_backbone_tree(backbone_mutations)
    I_B = I_selected_and_sorted[backbone_mutations]
    I_B_withNA3 = I_B.replace({np.nan: 3}).astype(int)
    I_B_withNA3.to_csv(os.path.join(dir_backbone, "I_B_withNA3.txt"), sep="\t")
        
    M_B, cell_assignments = impute_backbone_clones(I_selected_and_sorted, backbone_mutations, checked_mutation_group, logger_obj=log)
    
    scp.ul.is_conflict_free_gusfield(M_B)
    WriteTfile(
        os.path.join(dir_backbone, "M_backbone_basedPivots.filtered_sites_inferred"), 
        M_B, M_B.index.tolist(), M_B.columns.tolist(), judge="yes", logger_obj=log
    )
    
    # Clean M_B: remove all zeros columns(muts) or rows(cells)
    final_cleaned_M_B = M_B.loc[:, (M_B != 0).any(axis=0)]
    final_cleaned_M_B = final_cleaned_M_B.loc[(final_cleaned_M_B != 0).any(axis=1)]
    WriteTfile(
        os.path.join(dir_backbone, "final_cleaned_M_B_basedPivots.filtered_sites_inferred"), 
        final_cleaned_M_B, final_cleaned_M_B.index.tolist(), final_cleaned_M_B.columns.tolist(), judge="yes", logger_obj=log
    )
    kept_rows = final_cleaned_M_B.index
    kept_cols = final_cleaned_M_B.columns
    final_cleaned_I_B_withNA3 = I_B_withNA3.loc[kept_rows, kept_cols]
    final_cleaned_I_B_withNA3.to_csv(os.path.join(dir_backbone, "final_cleaned_I_B_withNA3_for_circosPlot.txt"), sep="\t")
    
    # ------------------------------
    # Step 9: Sub-grouping within Backbone Groups
    # ------------------------------
    log.info("Performing sub-grouping within backbone groups...")
    
    # Prepare data for subsequent operations
    # ===== Use sorted checked_mutation_group keys =====
    sorted_group_mutations = [i for i in sorted_mutation_names if i in group_mutations]
    I_selected = I_selected_and_sorted[sorted_group_mutations]
    P_selected = P_resolved.loc[I_selected.index, I_selected.columns]
    V_selected = V_resolved.loc[I_selected.index, I_selected.columns]
    A_selected = A_resolved.loc[I_selected.index, I_selected.columns]
    C_selected = C_resolved.loc[I_selected.index, I_selected.columns]
    df_reads_selected = df_reads_resolved.loc[I_selected.index.insert(0, 'bulk'), I_selected.columns]
    
    T_current = copy.deepcopy(T_B)
    M_current = M_B.copy()
    M_current.insert(0, 'ROOT', 1)
    
    # NA weight settings for penalty calculation
    ω_NA = params.get('general_weight_NA', 0.001)
    fnfp_ratio = params.get('fnfp_ratio', 0.1)
    φ = params.get('phi', 1.0)
    
    # Refining positions for each mutation
    group_but_backbone_mutations = [i for i in list(checked_mutation_group.keys()) if i not in backbone_mutations]
    all_mutations = I_somatic.columns.tolist()
    
    # Identify hub clusters and use penalty function for mutation placement
    root_mutations = []
    sorted_I_selected, sorting_stats_of_I_selected = reorder_columns_by_mutant_stats(I_selected, df_features_new)
    sorted_P_selected, sorting_stats_of_P_selected = reorder_columns_by_mutant_stats(P_selected, df_features_new)
    
    results_of_subgrouping = perform_subgrouping_within_backbone_groups_and_build_initial_scaffold_tree(
        sorted_I_selected, 
        sorted_P_selected, 
        T_current, 
        M_current, 
        checked_mutation_group, 
        backbones_of_group, 
        df_reads_selected, 
        df_features_new, 
        outputpath_scaffold, 
        sampleid, 
        log, 
        params, 
        M_B, 
        root_mutations, 
        cutoff_mcf_for_graph, 
        cutoff_mcn_for_graph
    )
    complete_mutation_hierarchy, subgroup_backbone_mutations, subgroup_details, T_current, M_current, root_mutations, external_mutations, conflict_mutations = results_of_subgrouping
    log.info(f"There are {len(external_mutations)} scaffold mutations remained because of low-support...")
    
    # ------------------------------
    # Step 10: Calculate Penalty and Refine Placement
    # ------------------------------
    log.info("Refining placement with penalty calculation...")
    remained_mutations = []
    if len(external_mutations) > 0:
        sorted_external_mutations = get_ordered_mutations(I_selected, df_features_new)
        external_mutations_sorted = [mut for mut in sorted_external_mutations if mut in external_mutations]
        
        remained_mutations, conflict_mutations_adding, T_current, M_current = integrate_mutations_to_scaffold_within_group(
            sorted_attached_mutations=external_mutations_sorted,
            T_current=T_current,
            M_current=M_current,
            I_attached=I_selected,
            P_attached=P_selected,
            mutation_group=checked_mutation_group,
            ω_NA=ω_NA,
            fnfp_ratio=fnfp_ratio,
            φ=φ,
            logger_obj=log
        )
        conflict_mutations.extend(conflict_mutations_adding)
    else:
        log.info("No external scaffold mutations after initial building.")
    
    # ------------------------------
    # Step 11: Handle Misassigned Mutations
    # ------------------------------
    log.info("Processing misassigned mutations...")
    remained_mutations_by_misgroup = []
    if misgroup_mutations:
        subtree_groups_of_misgroups = cluster_external_mutations_by_intersection_scaffold(I_resolved, misgroup_mutations, logger_obj=log)
        
        remained_mutations_by_misgroup, conflict_mutations_by_misgroup, T_current, M_current, root_mutations_by_misgroup = process_misassigned_mutations_direct_to_root(
            subtree_groups=subtree_groups_of_misgroups,
            T_current=T_current,
            M_current=M_current,
            I_attached=I_selected,
            P_attached=P_selected,
            ω_NA=ω_NA,
            fnfp_ratio=fnfp_ratio,
            φ=φ,
            logger_obj=log
        )
        # Update external_mutations with those still in conflict
        remained_mutations.extend(remained_mutations_by_misgroup)
        conflict_mutations.extend(conflict_mutations_by_misgroup)
        root_mutations.extend(root_mutations_by_misgroup)
        
    else:
        log.info("No misassigned mutations to process")
    
    # ------------------------------
    # Step 12: Output Results
    # ------------------------------
    log.info("Outputting results ...")
    
    # Final scaffold tree and matrix    
    T_scaffold = copy.deepcopy(T_current)
    M_scaffold_initial = M_current.copy()
    M_scaffold_initial = M_scaffold_initial.drop(columns=['ROOT'], errors='ignore')
    mutations_on_T_scaffold = M_scaffold_initial.columns.to_series().apply(lambda x: x.split("|")).explode().unique().tolist()
    M_scaffold = split_merged_columns(M_scaffold_initial, mutations_on_T_scaffold)
    
    print_tree_logger(T_scaffold, logger_obj=log)
    log.info(f"Final scaffold matrix shape: {M_scaffold.shape}")
    
    # Create output directory
    dir_scaffold = os.path.join(outputpath_scaffold, "phylo_scaffold_tree")
    if not os.path.exists(dir_scaffold):
        os.makedirs(dir_scaffold)
    
    # Prepare binary matrix with NA=3
    I_selected_withNA3 = I_selected.replace({np.nan: 3}).astype(int)
    I_selected_withNA3.to_csv(os.path.join(dir_scaffold, "I_scaffold_withNA3.txt"), sep="\t")
    WriteTfile(
        os.path.join(dir_scaffold, "M_scaffold_basedPivots.filtered_sites_inferred"), 
        M_scaffold, M_scaffold.index.tolist(), M_scaffold.columns.tolist(), judge="yes", logger_obj=log
    )
    
    # Clean M_scaffold: remove all zeros columns(muts) or rows(cells)
    final_cleaned_M_scaffold = M_scaffold.loc[:, (M_scaffold != 0).any(axis=0)]
    final_cleaned_M_scaffold = final_cleaned_M_scaffold.loc[(final_cleaned_M_scaffold != 0).any(axis=1)]
    
    # Get kept rows and columns
    kept_rows = final_cleaned_M_scaffold.index
    kept_cols = final_cleaned_M_scaffold.columns
    
    # Extract from I_selected_withNA3
    final_cleaned_I_selected_withNA3 = I_selected_withNA3.loc[kept_rows, kept_cols]
    
    WriteTfile(
        os.path.join(dir_scaffold, "final_cleaned_M_scaffold_basedPivots.filtered_sites_inferred"), 
        final_cleaned_M_scaffold, final_cleaned_M_scaffold.index.tolist(), final_cleaned_M_scaffold.columns.tolist(), judge="yes", logger_obj=log
    )
    final_cleaned_I_selected_withNA3.to_csv(os.path.join(dir_scaffold, "final_cleaned_I_scaffold_withNA3_for_circosPlot.txt"), sep="\t")
    
    # ------------------------------
    # Step 13: Identify Flipping Spots
    # ------------------------------
    log.info("Identifying flipping spots ...")
    df_I_selected_withNA3_for_flipping = final_cleaned_I_selected_withNA3.copy()
    df_phylogeny = final_cleaned_M_scaffold.copy()
    
    # Get false_negative_flipping spots
    false_negative_flipping_spots = df_I_selected_withNA3_for_flipping.apply(
        lambda col: find_flipping_spots(col, df_phylogeny[col.name], condition_in_bin=0, condition_phylogeny=1, logger_obj=log)
    )
    
    # Get NAto1 spots
    NAto1_flipping_spots = df_I_selected_withNA3_for_flipping.apply(
        lambda col: find_flipping_spots(col, df_phylogeny[col.name], condition_in_bin=3, condition_phylogeny=1, logger_obj=log)
    )
    
    # Get false_positive_flipping spots
    false_positive_flipping_spots = df_I_selected_withNA3_for_flipping.apply(
        lambda col: find_flipping_spots(col, df_phylogeny[col.name], condition_in_bin=1, condition_phylogeny=0, logger_obj=log)
    )
    
    # Get NAto0 spots
    NAto0_flipping_spots = df_I_selected_withNA3_for_flipping.apply(
        lambda col: find_flipping_spots(col, df_phylogeny[col.name], condition_in_bin=3, condition_phylogeny=0, logger_obj=log)
    )
    
    # Process na list
    if false_negative_flipping_spots.empty:
        false_negative_flipping_spots = {col: [] for col in df_I_selected_withNA3_for_flipping.columns}
    
    if false_positive_flipping_spots.empty:
        false_positive_flipping_spots = {col: [] for col in df_I_selected_withNA3_for_flipping.columns}
    
    if NAto1_flipping_spots.empty:
        NAto1_flipping_spots = {col: [] for col in df_I_selected_withNA3_for_flipping.columns}
    
    if NAto0_flipping_spots.empty:
        NAto0_flipping_spots = {col: [] for col in df_I_selected_withNA3_for_flipping.columns}
    
    # Get flipping_spots dataframe
    df_flipping_spots = pd.DataFrame({
        'Mutation': df_I_selected_withNA3_for_flipping.columns,
        'false_negative_flipping_spots': [', '.join(false_negative_flipping_spots.get(col, [])) for col in df_I_selected_withNA3_for_flipping.columns],
        'false_positive_flipping_spots': [', '.join(false_positive_flipping_spots.get(col, [])) for col in df_I_selected_withNA3_for_flipping.columns],
        'NAto1_flipping_spots': [', '.join(NAto1_flipping_spots.get(col, [])) for col in df_I_selected_withNA3_for_flipping.columns],
        'NAto0_flipping_spots': [', '.join(NAto0_flipping_spots.get(col, [])) for col in df_I_selected_withNA3_for_flipping.columns]
    })
    df_flipping_spots.to_csv(os.path.join(dir_scaffold, "df_flipping_spots.txt"), sep="\t", index=False)
    
    # ------------------------------
    # Step 14: Calculate Total Flipping Counts
    # ------------------------------
    log.info("Calculating total flipping counts ...")
    total_FN_flipping = ((df_I_selected_withNA3_for_flipping == 0) & (df_phylogeny == 1)).sum().sum()
    total_FP_flipping = ((df_I_selected_withNA3_for_flipping == 1) & (df_phylogeny == 0)).sum().sum()
    total_NAto0 = ((df_I_selected_withNA3_for_flipping == 3) & (df_phylogeny == 0)).sum().sum()
    total_NAto1 = ((df_I_selected_withNA3_for_flipping == 3) & (df_phylogeny == 1)).sum().sum()
    
    log.info(f"Total False Negative flipping count: {total_FN_flipping}")
    log.info(f"Total False Positive flipping count: {total_FP_flipping}")
    log.info(f"Total NA to 0 flipping count: {total_NAto0}")
    log.info(f"Total NA to 1 flipping count: {total_NAto1}")
    
    df_total_flipping_count = pd.DataFrame({
        'total_flipping_False_Negative': [total_FN_flipping],
        'total_flipping_False_Positive': [total_FP_flipping],
        'total_flipping_NA_to_0': [total_NAto0],
        'total_flipping_NA_to_1': [total_NAto1]
    })
    df_total_flipping_count.to_csv(os.path.join(dir_scaffold, "df_total_flipping_count.txt"), sep="\t", index=False)
    
    log.info("Scaffold building completed successfully!")
    
    # Convert TreeNode to dictionary format and save as JSON file
    tree_dict = tree_to_dict(T_scaffold)
    with open(os.path.join(dir_scaffold, 'T_scaffold.json'), 'w') as f:
        json.dump(tree_dict, f, indent=4)
    
    # Save tree as text format
    T_scaffold.save_to_file(os.path.join(dir_scaffold, 'T_scaffold.txt'))
    
    return (
        T_scaffold, 
        M_scaffold, 
        df_flipping_spots, 
        df_total_flipping_count, 
        final_cleaned_I_selected_withNA3, 
        final_cleaned_M_scaffold, 
        backbone_mutations, 
        checked_mutation_group, 
        spots_to_split, 
        list(checked_mutation_group.keys()), 
        no_group_mutations,
        remained_mutations, 
        conflict_mutations, 
        root_mutations
    )



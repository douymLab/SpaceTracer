#!/usr/bin/env python3

"""
mutation_integrator.py

Integrate mutations into scaffold tree for PhyloSOLID (Methods Section 4).

Implements:
- 4.1 Dynamic programming-based phylogenetic likelihood evaluation for non-scaffold mutations
- 4.2 Bayesian posterior calculation for mutation classification (mosaic, heterozygous, reference homozygous)
- 4.3 Logistic regression-based classification using tree-based posterior probability and read-level features
- 4.4 Placement of high-confidence mosaic mutations onto the scaffold tree

Inputs:
- P: posterior probability matrix (cells x mutations)  
- G: genotype matrix (cells x mutations)  
- read_features: dataframe of read-level features (cells x mutations)  
- scaffold_tree: initial scaffold tree structure (TreeNode)
- params: dictionary of parameters with defaults

Outputs:
- full_tree: expanded scaffold tree with mosaic mutations integrated
- classified_mutations: list of classified mosaic mutations
- posterior_probabilities: posterior probabilities for each mutation
"""


import os
import numpy as np
import pandas as pd
import networkx as nx
import logging
import copy
import math
import random
from tqdm import tqdm
from copy import deepcopy
import itertools
from itertools import combinations
from collections import defaultdict, Counter
from typing import Set, List, Dict, Optional, Tuple, Any, Union
from scipy.special import logsumexp, comb, perm, beta
from anytree import Node, RenderTree
from src.germline_filter import pairwise_counts, jaccard_index, are_mutations_correlated, reorder_columns_by_mutant_stats
from src.scaffold_builder import TreeNode, tree_to_dict, print_tree_dict
from src.scaffold_builder import print_tree, add_new_mutation_to_tree_independent, split_merged_columns, WriteTfile, compute_bayesian_penalty_each_pos, compute_bayesian_penalty_each_chain_mut_by_pos, build_lineage_parent_dict_from_tree
from src.reproducibility import set_seed, deterministic_permutation, get_seed
import scphylo as scp
from scphylo.pl._helper import (
    _add_barplot,
    _add_chromplot,
    _clonal_cell_mutation_list,
    _get_tree,
    _newick_info2_mutation_list,
)

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

# Default parameters matching Methods Section 3
DEFAULT_PARAMS = {
    "posterior_threshold": 0.5,        # Phylogenetic mosaic threshold
    "maf_max_threshold": 0.3,           # Maximum mutation frequency threshold
    "maf_mean_threshold": 0.1,          # Mean mutation frequency threshold
    "coverage_threshold": 10,           # Read coverage threshold
    "logistic_regression_threshold": 0.5,  # Logistic regression decision threshold
    # Additional parameters can be added here as needed
}


# -------------------------
# Utils
# -------------------------

def count_list(mutation_series):
    # 替换 NaN 为 'NA'
    mutation_series = mutation_series.fillna('NA')
    # 计数
    return dict(Counter(mutation_series))

def log_sum_exp(log_probs):
    max_log_prob = np.max(log_probs)
    return max_log_prob + np.log(np.sum(np.exp(log_probs - max_log_prob)))

def logdiffexp(a, b):
       return logsumexp([a, b], b=[1, -1])


def divide_columns(row):
    if pd.notna(row['total']) and row['total'] != 0:
        return int(row['alt']) / int(row['total'])
    else:
        return np.nan

def reads2df(reads):
    df = pd.DataFrame(reads, columns=['data'])
    df['alt'] = df['data'].str.split('/').str[0]
    df['total'] = df['data'].str.split('/').str[1]
    df['af'] = df.apply(divide_columns, axis=1)
    df['is_mut'] = pd.to_numeric(df['alt']).apply(lambda x: 1 if x > 0 else 0)
    df['is_mut'] = df['is_mut'].fillna(0)
    return df    

def get_ternaryVec_NAto3(vector, df_reads_persite, cutoff):
    more_cutoff = np.where((vector > cutoff) & (vector <= 1), 1, np.where(vector == 3, 3, 0))
    mut_info = df_reads_persite['is_mut'].to_numpy()
    subclone = np.where((more_cutoff == 1) & (mut_info == 1), 1,
        np.where((more_cutoff == 0) & (mut_info == 1), 0,
            np.where((more_cutoff == 1) & (mut_info == 0), 0,
                np.where((more_cutoff == 3) | (mut_info == 3), 3, 0))))
    return subclone

def posterior2ter_NAto3_bothPosteriorMutallele(df_posterior, df_reads, cutoff):
    M_posterior = df_posterior.values
    M_reads = df_reads.iloc[1:,].values
    # Get the primarily binary martrix
    M_bin = M_posterior.copy()
    M_bin[np.isnan(M_bin)] = 3
    for i, v in enumerate(M_bin.T):
        M_bin[:,i] = get_ternaryVec_NAto3(v, reads2df(M_reads[:,i]), cutoff)
    df_bin = pd.DataFrame(M_bin).astype(int)
    df_bin.columns = df_posterior.columns
    df_bin.index = df_posterior.index
    df_bin.index.name = "cellIDxmutID"
    return df_bin

def count_conditions(I_attached_col, M_current_col):
    """
    Count the number of rows in I_attached and M_current meeting specific conditions:
    1. Number of rows where former is 1 and latter is 0
    2. Number of rows where former is 0 and latter is 1
    3. Number of rows where former is NaN and latter is 1
    4. Number of rows where former is NaN and latter is 0
    
    Parameters:
    I_attached_col: pandas Series, a column from I_attached
    M_current_col: pandas Series, a column from M_current
    
    Returns:
    dict: A dictionary containing four statistical values
    """
    # Convert to corresponding boolean values
    I_attached_col = I_attached_col.fillna(np.nan)
    M_current_col = M_current_col.fillna(np.nan)
    
    count_1_0 = ((I_attached_col == 1) & (M_current_col == 0)).sum()  # I_attached is 1 and M_current is 0
    count_0_1 = ((I_attached_col == 0) & (M_current_col == 1)).sum()  # I_attached is 0 and M_current is 1
    count_na_1 = ((I_attached_col.isna()) & (M_current_col == 1)).sum()  # I_attached is NaN and M_current is 1
    count_na_0 = ((I_attached_col.isna()) & (M_current_col == 0)).sum()  # I_attached is NaN and M_current is 0
    
    return {
        'count_fp_1_0': count_1_0,
        'count_fn_0_1': count_0_1,
        'count_na_1': count_na_1,
        'count_na_0': count_na_0
    }


def as_binary_mask(series, index):
    """
    Normalize a Series-like vector to an aligned boolean mask.
    """
    aligned = pd.Series(series, index=index) if not isinstance(series, pd.Series) else series.reindex(index)
    values = aligned.to_numpy(copy=False)
    
    if np.issubdtype(values.dtype, np.bool_):
        return pd.Series(values, index=index, copy=False)
    
    if np.issubdtype(values.dtype, np.integer):
        return pd.Series(values > 0, index=index)
    
    if np.issubdtype(values.dtype, np.floating):
        mask = np.isfinite(values) & (values > 0)
        return pd.Series(mask, index=index)
    
    numeric = pd.to_numeric(aligned, errors="coerce")
    numeric_values = numeric.to_numpy(dtype=np.float64, na_value=np.nan)
    mask = np.isfinite(numeric_values) & (numeric_values > 0)
    return pd.Series(mask, index=index)


def _coerce_bool_array(series, index):
    return series.reindex(index, fill_value=0).fillna(0).to_numpy(dtype=bool, copy=False)


def _coerce_float_array(series, index):
    return pd.to_numeric(series.reindex(index), errors="coerce").to_numpy(
        dtype=np.float64,
        na_value=np.nan,
    )


def _log_tree_snapshot(active_logger, label, tree):
    if not _VERBOSE_TREE_UPDATES:
        return
    active_logger.info(label)
    # print_tree(tree)


def _sanitize_matrix_for_conflict_check(matrix, active_logger, context_label, sample_limit=5):
    values = matrix.to_numpy(copy=False)
    try:
        finite_mask = np.isfinite(values)
        numeric_values = values.astype(np.float64, copy=False)
    except TypeError:
        numeric_values = matrix.apply(pd.to_numeric, errors="coerce").to_numpy(
            dtype=np.float64,
            copy=False,
            na_value=np.nan,
        )
        finite_mask = np.isfinite(numeric_values)
    
    if bool(finite_mask.all()):
        return matrix
    
    invalid_positions = np.argwhere(~finite_mask)
    samples = [
        f"{matrix.index[row]}::{matrix.columns[col]}={numeric_values[row, col]}"
        for row, col in invalid_positions[:sample_limit]
    ]
    active_logger.warning(
        "%s sanitizing matrix before conflict check: non_finite=%d sample=%s",
        context_label,
        int((~finite_mask).sum()),
        samples,
    )
    
    sanitized_values = np.nan_to_num(
        numeric_values,
        nan=0.0,
        posinf=1.0,
        neginf=0.0,
        copy=True,
    )
    binary_values = (sanitized_values > 0).astype(np.int8, copy=False)
    return pd.DataFrame(binary_values, index=matrix.index, columns=matrix.columns)


def _is_conflict_free_gusfield_safe(matrix, active_logger=None, context_label="ConflictCheck"):
    if active_logger is None:
        active_logger = logger
    try:
        return scp.ul.is_conflict_free_gusfield(matrix)
    except (pd.errors.IntCastingNaNError, ValueError, TypeError) as exc:
        active_logger.warning(
            "%s conflict check retry after non-binary/non-finite matrix error: %s",
            context_label,
            exc,
        )
        matrix_for_check = _sanitize_matrix_for_conflict_check(matrix, active_logger, context_label)
        return scp.ul.is_conflict_free_gusfield(matrix_for_check)


def _resolve_touched_columns_for_conflict_check(final_position, parent_dict, matrix_columns, new_mut):
    """Return the matrix columns whose updates could introduce new conflicts."""
    placement_type = final_position.get("placement_type")
    anchor = final_position.get("anchor")
    matrix_column_set = set(matrix_columns)
    touched_columns = []
    seen = set()
    
    if anchor is not None:
        for node in get_full_mutnode_chain_with_anchor(anchor, parent_dict):
            if node == "ROOT":
                continue
            resolved_node = node
            if placement_type == "on_node" and node == anchor and anchor != "ROOT":
                resolved_node = f"{anchor}|{new_mut}"
            if resolved_node in matrix_column_set and resolved_node not in seen:
                touched_columns.append(resolved_node)
                seen.add(resolved_node)
    
    if placement_type != "on_node" and new_mut in matrix_column_set and new_mut not in seen:
        touched_columns.append(new_mut)
    
    return touched_columns


def _is_conflict_free_local_update(matrix, touched_columns, active_logger=None, context_label="LocalConflictCheck", sample_limit=5):
    """
    Exact local conflict check using the four-gamete criterion.
    
    Only pairs involving touched columns can introduce a new conflict because
    untouched-vs-untouched pairs are unchanged by the latest placement step.
    """
    if active_logger is None:
        active_logger = logger
    
    if not touched_columns:
        return True
    
    missing = [col for col in touched_columns if col not in matrix.columns]
    if missing:
        active_logger.warning(
            "%s touched columns missing from matrix; falling back to full check. sample=%s",
            context_label,
            missing[:sample_limit],
        )
        return _is_conflict_free_gusfield_safe(matrix, active_logger, context_label)
    
    values = matrix.to_numpy(copy=False)
    if np.issubdtype(values.dtype, np.bool_):
        bool_values = values
    elif np.issubdtype(values.dtype, np.integer):
        bool_values = values > 0
    elif np.issubdtype(values.dtype, np.floating):
        bool_values = np.nan_to_num(values, nan=0.0, posinf=1.0, neginf=0.0, copy=True) > 0
    else:
        numeric_values = matrix.apply(pd.to_numeric, errors="coerce").to_numpy(
            dtype=np.float64,
            copy=False,
            na_value=np.nan,
        )
        bool_values = np.nan_to_num(numeric_values, nan=0.0, posinf=1.0, neginf=0.0, copy=True) > 0
    
    all_columns = matrix.columns
    touched_indices = all_columns.get_indexer(touched_columns)
    if np.any(touched_indices < 0):
        unresolved = [touched_columns[i] for i, idx in enumerate(touched_indices) if idx < 0]
        active_logger.warning(
            "%s unresolved touched indices; falling back to full check. sample=%s",
            context_label,
            unresolved[:sample_limit],
        )
        return _is_conflict_free_gusfield_safe(matrix, active_logger, context_label)
    
    for touched_name, touched_idx in zip(touched_columns, touched_indices):
        col = bool_values[:, touched_idx][:, None]
        has10 = np.any(col & ~bool_values, axis=0)
        has01 = np.any(~col & bool_values, axis=0)
        has11 = np.any(col & bool_values, axis=0)
        conflict_mask = has10 & has01 & has11
        conflict_mask[touched_idx] = False
        if bool(conflict_mask.any()):
            conflict_indices = np.flatnonzero(conflict_mask)[:sample_limit]
            conflict_columns = [all_columns[i] for i in conflict_indices]
            active_logger.warning(
                "%s detected local conflict: touched=%s conflicts=%s",
                context_label,
                touched_name,
                conflict_columns,
            )
            return False
    
    return True

# -------------------------
# Normlize input likelihood probability data
# -------------------------

def convert_to_log(df, epsilon=1e-10):
    """
    Convert the raw probability data to log format.
    Avoid the log-zero problem by taking the log(max(val, epsilon)).
    """
    # Use epsilon to prevent log(0) problems
    return np.log(np.maximum(df, epsilon))

def exp_normalize(log_prob):
    max_log_prob = log_prob.max()
    y = np.exp(log_prob - max_log_prob)
    return y / y.sum()

def normalize_columns(llmut_col, llunmut_col):
    # 转换为 log 格式
    log_llmut_col = llmut_col
    log_llunmut_col = llunmut_col
    paired_likelihoods_prod = np.vstack([log_llmut_col, log_llunmut_col]).T
    norm_llmut_prod, norm_llunmut_prod = np.log(np.array([exp_normalize(pair) for pair in paired_likelihoods_prod]).T)
    return pd.Series({"norm_llmut": norm_llmut_prod, "norm_llunmut": norm_llunmut_prod})

def apply_normalization(df, is_log=True):
    if is_log:
        process_df = df
    else:
        # Converts the original probability value to log format
        process_df = convert_to_log(df)
    
    norm_results = pd.DataFrame(index=df.index)
    num_columns = df.shape[1] // 2
    
    for i in range(num_columns):
        llmut_col = process_df.iloc[:, i]
        llunmut_col = process_df.iloc[:, num_columns + i]
        normalized = normalize_columns(llmut_col, llunmut_col)
        norm_results[f'norm_llmut_{i}'] = normalized['norm_llmut']
        norm_results[f'norm_llunmut_{i}'] = normalized['norm_llunmut']
    
    return norm_results


# -------------------------
# 4.1 Dynamic Programming-Based Phylogenetic Likelihood Evaluation
# &
# 4.2 Bayesian Posterior Calculation
# -------------------------

def intersect_is_self(v1, v0):
    # If True, v1 is a subset of v0
    v0_indices = [i for i, val in enumerate(v0) if val == 1]
    v1_indices = [i for i, val in enumerate(v1) if val == 1]
    v0_set = set(v0_indices)
    v1_set = set(v1_indices)
    if v1_set == v1_set.intersection(v0_set):
        return True
    else:
        return False

def get_allBranchSet_as_dict(df_phylogeny):
    # 强制每一列的元素转为字符串并拼接
    df = df_phylogeny.apply(lambda col: ''.join(map(str, col.astype(int))), axis=0)
    content_dict = defaultdict(list)
    for index, content in df.items():
        content_dict[content].append(index)
    final_dict = {content: '+'.join(str(indexes)) for content, indexes in content_dict.items()}
    return final_dict

def str2array(string_value):
    return np.array([int(char) for char in string_value])

def get_allBranchSet(M_tree):
    columns_to_delete = np.all(M_tree == 1, axis=0)
    if all(not item for item in columns_to_delete):
        M_tree_noRoot = M_tree[:, ~columns_to_delete]
        clusters_allBranchSet = [str2array(s) for s in get_allBranchSet_as_dict(pd.DataFrame(M_tree_noRoot)).keys()]
    else:
        clusters_allBranchSet = [str2array(s) for s in get_allBranchSet_as_dict(pd.DataFrame(M_tree)).keys()]
    clusters_allBranchSet = [cluster for cluster in clusters_allBranchSet if not np.all(cluster == 0)]
    return clusters_allBranchSet

def get_1stBranchSet(all_clusters):
    pivot_clusters = []
    for i, v1 in enumerate(all_clusters):
        is_pivot = True
        for j, v0 in enumerate(all_clusters):
            if i != j:
                if intersect_is_self(v1, v0):
                    if intersect_is_self(v0, v1):
                        continue
                    else:
                        is_pivot = False
                        break
        if is_pivot:
            pivot_clusters.append(v1)
    return pivot_clusters

def get_earlyBranchSet(no_pivot_clusters, pivot_clusters):
    early_clusters = []
    for pivot_cluster in pivot_clusters:
        early_clusters.append(pivot_cluster)
        subset_clusters = [cluster for cluster in no_pivot_clusters if intersect_is_self(cluster, pivot_cluster)]
        if len(subset_clusters)==0:
            continue
        else:
            # M_subset_clusters = np.vstack(subset_clusters).T
            # first_clusters = get_1stBranchSet(M_subset_clusters)
            first_subset_clusters = get_1stBranchSet(subset_clusters)
            if len(first_subset_clusters)==0:
                continue
            else:
                early_clusters = early_clusters + first_subset_clusters
    
    return early_clusters

def get_leafBranchSet(all_clusters):
    leaf_clusters = []
    for i, v1 in enumerate(all_clusters):
        is_leaf = True
        for j, v0 in enumerate(all_clusters):
            if i != j:
                if intersect_is_self(v0, v1):
                    is_leaf = False
        if is_leaf:
            leaf_clusters.append(v1)
    return leaf_clusters

def is_subset(cluster, node_cluster):
    # Gets the position of 1 in the cluster
    cluster_one_positions = np.where(cluster == 1)[0]
    # Determine if the node_cluster has at least one, but not all, 1s at these locations
    one_positions = np.sum(node_cluster[cluster_one_positions])
    # At least one 1 but not all 1
    return 0 < one_positions < len(cluster_one_positions)

def DP_calSomaticPosterior_withoutTree(site_llmut, site_llunmut, cell_num):
    log_sc_gt_likelihood_dict = {}
    log_sc_gt_likelihood_given_l_list = []
    log_sc_gt_posterior_given_l_list = []
    log_sc_gt_likelihood_dict[(0,0)]=np.log(1)
    for i in range(cell_num):
        log_sc_gt_likelihood_dict[(-1,i)] = -np.inf
    for i in range(1,int(cell_num)+1):
        log_sc_gt_likelihood_dict[(i,0)] = -np.inf
    for l in range(int(cell_num)+1):
        for n in range(1,int(cell_num)+1):
            # M_l_n = M_l_n-1 * P(Dn|Gn_is_not_mutate) + M_l-1_n-1 * P(Dn|Gn_is_mutated)
            log_sc_gt_likelihood_dict[(l,n)] = logsumexp([log_sc_gt_likelihood_dict[(l,n-1)] + site_llunmut[n-1], log_sc_gt_likelihood_dict[(l-1,n-1)] + site_llmut[n-1]])
        log_sc_gt_likelihood_given_l_list.append(log_sc_gt_likelihood_dict[(l,cell_num)])
    posterior_l_is_0_dp_path=exp_normalize(np.array(log_sc_gt_likelihood_given_l_list))[0]
    somatic_posterior_per_site = float(1 - posterior_l_is_0_dp_path)
    return somatic_posterior_per_site

def get_optimal_cell_combination(path_dict, cell_num, max_l):
    # Retrace the path to find mutant combinations
    l, n = max_l, cell_num
    mut_cell_indices = [0] * cell_num
    while (l, n) in path_dict and path_dict[(l, n)] is not None:
        prev_l, prev_n = path_dict[(l, n)]
        if prev_l == l - 1:  # mutation
            mut_cell_indices[n-1] = 1  # n-1 cell mutation
        l, n = prev_l, prev_n
    return mut_cell_indices

def DP_calSomaticPosterior_late(site_llmut, site_llunmut, leaf_cluster):
    ### Late: t least 1 cell
    non_existing = np.sum(site_llunmut[leaf_cluster==0])
    probs_llmut = site_llmut[leaf_cluster==1]
    probs_llunmut = site_llunmut[leaf_cluster==1]
    cell_num = len(probs_llmut)
    # The initializing state
    likelihood_dict = {}
    likelihood_list = []
    path_dict = {}  # To track the path
    # Initialize DP tables
    likelihood_dict[(0,0)] = np.log(1)
    path_dict[(0,0)] = None  # No previous state
    for i in range(cell_num):
        likelihood_dict[(-1,i)] = -np.inf
        path_dict[(-1,i)] = None
    for i in range(1, int(cell_num)+1):
        likelihood_dict[(i,0)] = -np.inf
        path_dict[(i,0)] = None
    # DP: n cells, l mutant cells
    for l in range(int(cell_num)+1):
        for n in range(1, int(cell_num)+1):
            # M_l_n = M_l_n-1 * P(Dn|Gn_is_not_mutate) + M_l-1_n-1 * P(Dn|Gn_is_mutated)
            if n == 1:  # Apply non_existing only at the first step
                option1 = likelihood_dict[(l, n-1)] + probs_llunmut[n-1] + non_existing
                option2 = likelihood_dict[(l-1, n-1)] + probs_llmut[n-1] + non_existing
            else:
                option1 = likelihood_dict[(l, n-1)] + probs_llunmut[n-1]
                option2 = likelihood_dict[(l-1, n-1)] + probs_llmut[n-1]
            likelihood_dict[(l, n)] = logsumexp([option1, option2])
            if option1 > option2:
                path_dict[(l,n)] = (l, n-1)
            else:
                path_dict[(l,n)] = (l-1, n-1)
        likelihood_list.append(likelihood_dict[(l,cell_num)])
    all_prob = logsumexp(likelihood_list)
    at_least_one_mutcell_prob = np.log(np.exp(all_prob) - np.exp(likelihood_list[0]))
    one_mutcell_prob = likelihood_list[1]
    # Find the combination of mutations with the greatest probability
    max_l = np.argmax([likelihood_dict[(l, cell_num)] for l in range(cell_num + 1)])
    mut_cells = get_optimal_cell_combination(path_dict, cell_num, max_l)
    return [at_least_one_mutcell_prob, likelihood_dict[(max_l, cell_num)], mut_cells, one_mutcell_prob]

def DP_calSomaticPosterior_early(site_llmut, site_llunmut, pivot_clusters):
    ### Early: at least 2 clusters but not all clusters
    cluster_llmut = np.array([np.sum(site_llmut[c == 1]) for c in pivot_clusters])
    cluster_llunmut = np.array([np.sum(site_llunmut[c == 1]) for c in pivot_clusters])
    cluster_num = len(cluster_llmut)
    # The initializing state
    likelihood_dict = {}
    likelihood_list = []
    path_dict = {}
    # Initialize DP tables
    likelihood_dict[(0,0)] = np.log(1)
    path_dict[(0,0)] = None  # No previous state
    for i in range(cluster_num):
        likelihood_dict[(-1, i)] = -np.inf
        path_dict[(-1,i)] = None
    for i in range(1, int(cluster_num)+1):
        likelihood_dict[(i, 0)] = -np.inf
        path_dict[(i,0)] = None
    # DP: n cells, l mutant cells
    for l in range(int(cluster_num)+1):
        for n in range(1, int(cluster_num)+1):
            # Calculate the two possible values for logsumexp
            option1 = likelihood_dict[(l,n-1)] + cluster_llunmut[n-1]
            option2 = likelihood_dict[(l-1,n-1)] + cluster_llmut[n-1]
            likelihood_dict[(l,n)] = logsumexp([option1, option2])
            # Choose the maximum of the two options
            if option1 > option2:
                path_dict[(l,n)] = (l, n-1)
            else:
                path_dict[(l,n)] = (l-1, n-1)
        likelihood_list.append(likelihood_dict[(l,cluster_num)])
    all_prob = logsumexp(likelihood_list)
    at_least_two_mutcluster_but_not_allcluster_prob = logdiffexp(all_prob, logsumexp([likelihood_list[0], likelihood_list[1], likelihood_list[-1]]))
    # Find the combination of mutations with the greatest probability
    max_l = range(2, cluster_num)[np.argmax([likelihood_dict[(l, cluster_num)] for l in range(2, cluster_num)])]  # at least 2 clusters but all clusters
    mut_clusters = get_optimal_cell_combination(path_dict, cluster_num, max_l)
    return [at_least_two_mutcluster_but_not_allcluster_prob, likelihood_dict[(max_l, cluster_num)], mut_clusters]

def DP_calSomaticPosterior_internal(site_llmut, site_llunmut, internal_subset):
    ### Internal: at least 2 clusters
    non_existing = np.sum(site_llunmut[sum(internal_subset)==0])
    cluster_llmut = np.array([np.sum(site_llmut[c==1]) for c in internal_subset])
    cluster_llunmut = np.array([np.sum(site_llunmut[c==1]) for c in internal_subset])
    cluster_num = len(cluster_llmut)
    # The initializing state
    likelihood_dict = {}
    path_dict = {}  # To track the path
    likelihood_list = [] 
    # Initial state
    likelihood_dict[(0,0)] = np.log(1)
    path_dict[(0,0)] = None
    for i in range(cluster_num):
        likelihood_dict[(-1,i)] = -np.inf
        path_dict[(-1,i)] = None
    for i in range(1, int(cluster_num)+1):
        likelihood_dict[(i,0)] = -np.inf
        path_dict[(i,0)] = None
    # DP: n cells, l mutant cells
    for l in range(int(cluster_num)+1):
        for n in range(1, int(cluster_num)+1):
            option1 = likelihood_dict[(l, n-1)] + cluster_llunmut[n-1] + non_existing
            option2 = likelihood_dict[(l-1, n-1)] + cluster_llmut[n-1] + non_existing
            likelihood_dict[(l,n)] = logsumexp([option1, option2])
            if option1 > option2:
                path_dict[(l, n)] = (l, n-1)
            else:
                path_dict[(l, n)] = (l-1, n-1)
        likelihood_list.append(likelihood_dict[(l,cluster_num)])
    all_prob = logsumexp(likelihood_list)
    at_least_two_mutcluster_prob = logdiffexp(all_prob, logsumexp([likelihood_list[0],likelihood_list[1]]))
    # Find the combination of mutations with the greatest probability
    max_l = range(2, cluster_num+1)[np.argmax([likelihood_dict[(l, cluster_num)] for l in range(2, cluster_num+1)])]  # at least 2 clusters but all clusters
    mut_clusters = get_optimal_cell_combination(path_dict, cluster_num, max_l)
    return [at_least_two_mutcluster_prob, likelihood_dict[(max_l, cluster_num)], mut_clusters]

def get_newSomaticPosterior(site_llmut, site_llunmut, node_clusters_dict):
    # node/cluster
    leaf_clusters = node_clusters_dict['leaf_clusters']
    internal_clusters = node_clusters_dict['internal_clusters']
    pivot_clusters = node_clusters_dict['pivot_clusters']
    node_clusters = leaf_clusters+internal_clusters+pivot_clusters
    node_num = np.sum(leaf_clusters) + len(internal_clusters)
    ### calculate posterior under conditions
    node_position_list = []
    branch_state_list = []
    cluster_list = []
    mut_posterior_list = []
    max_mut_posterior_list = []
    onecell_mut_posterior_list = []
    
    # Condition: late branches
    for idx,leaf_cluster in enumerate(leaf_clusters):
        node_position_list.append("leaf_clusters")
        branch_state_list.append("late")
        # dp
        late_results = DP_calSomaticPosterior_late(site_llmut, site_llunmut, leaf_cluster)
        # cluster
        late_cluster_predicted = np.zeros_like(leaf_cluster)
        indices_with_1 = np.where(leaf_cluster == 1)[0]
        late_cluster_predicted[indices_with_1] = late_results[2]
        cluster_list.append(late_cluster_predicted)
        # posterior
        mut_posterior_list.append(late_results[0])
        max_mut_posterior_list.append(late_results[1])
        onecell_mut_posterior_list.append(late_results[3])
    
    # Condition: internal branches
    for idx,internal_cluster in enumerate(internal_clusters):
        node_position_list.append("internal_clusters")
        branch_state_list.append("internal")
        # dp
        internal_subset = [node_cluster for node_cluster in node_clusters if is_subset(internal_cluster, node_cluster)]
        if len(internal_subset)>1:
            internal_results = DP_calSomaticPosterior_internal(site_llmut, site_llunmut, internal_subset)
            # cluster
            indices_with_1 = np.where(np.array(internal_results[2]) == 1)[0]
            internal_cluster_predicted = np.sum(np.array(internal_subset)[indices_with_1], axis=0)
            cluster_list.append(internal_cluster_predicted)
            # posterior
            mut_posterior_list.append(internal_results[0])
            max_mut_posterior_list.append(internal_results[1])
        elif len(internal_subset)==1:
            # cluster
            cluster_list.append(internal_subset[0])
            # posterior
            mut_posterior_list.append(np.sum(np.concatenate([site_llmut[internal_cluster==1], site_llunmut[internal_cluster==0]])))
            max_mut_posterior_list.append(np.sum(np.concatenate([site_llmut[internal_cluster==1], site_llunmut[internal_cluster==0]])))
    
    # Condition: early branches
    if len(pivot_clusters)>2:
        node_position_list.append("pivot_clusters")
        branch_state_list.append("early")
        # dp
        early_results = DP_calSomaticPosterior_early(site_llmut, site_llunmut, pivot_clusters)
        # cluster
        indices_with_1 = np.where(np.array(early_results[2]) == 1)[0]
        early_cluster_predicted = np.sum(np.array(pivot_clusters)[indices_with_1], axis=0)
        cluster_list.append(early_cluster_predicted)
        # posterior
        mut_posterior_list.append(early_results[0])
        max_mut_posterior_list.append(early_results[1])
    else:
        node_position_list.append("pivot_clusters")
        branch_state_list.append("none")
        cluster_list.append(np.zeros(len(site_llmut), dtype=int))
        mut_posterior_list.append(-np.inf)
        max_mut_posterior_list.append(-np.inf)
    
    # Collect final results
    df_allSP = pd.DataFrame({'node_position': node_position_list, 'cluster': cluster_list, 'branch_state': branch_state_list, 'somatic_posterior_conditional': max_mut_posterior_list})
    max_row = df_allSP.loc[df_allSP['somatic_posterior_conditional'].idxmax()]
    # Calculate posterior
    wrong_posterior = np.exp(np.sum(site_llmut))
    unmut_posterior = np.exp(np.sum(site_llunmut))
    mut_posterior = np.exp(logsumexp(mut_posterior_list))*(1/node_num)
    onecell_posterior = np.exp(logsumexp(onecell_mut_posterior_list))*(1/node_num)
    ### normalized SP
    total_prob = wrong_posterior + unmut_posterior + mut_posterior
    normalized_mut = mut_posterior / total_prob
    normalized_onecell = onecell_posterior / total_prob
    normalized_unmut = unmut_posterior / total_prob
    normalized_wrong = wrong_posterior / total_prob
    normalized_mut_max = np.exp(max_row['somatic_posterior_conditional'])*(1/node_num) / total_prob
    ### out frame
    df_outSP = max_row.copy()
    df_outSP = df_outSP.drop(['somatic_posterior_conditional'])
    df_outSP['somatic_posterior_per_site_max'] = normalized_mut_max
    df_outSP['somatic_posterior_per_site_onecell'] = normalized_onecell
    df_outSP['somatic_posterior_per_site'] = normalized_mut
    df_outSP['artifact_posterior_per_site'] = normalized_unmut
    df_outSP['germline_posterior_per_site'] = normalized_wrong
    
    return [df_allSP, df_outSP]

def all_newSomaticPosterior(df_llmut, df_llunmut, M_lowR):
    ### Step1: initiate 2 dataframe
    # Out somatic_posterior_per_site based on tree (normalization)
    df_newSP_out_init = pd.DataFrame(columns=['node_position', 'cluster', 'branch_state', 'somatic_posterior_per_site_max', 'somatic_posterior_per_site_onecell', 'somatic_posterior_per_site', 'artifact_posterior_per_site', 'germline_posterior_per_site'])
    df_newSP_out = df_newSP_out_init.dropna(axis=1, how='all')
    ### Step2: Calculation SP for each site
    M_llmut = df_llmut.values
    M_llunmut = df_llunmut.values
    withoutTree_posterior = []
    df_newSP_out = df_newSP_out.copy()
    all_nonnan_indices = []
    
    for i in tqdm(range(0, M_llunmut.shape[1])):
        new_mutid = df_llunmut.columns[i]
        site_llmut = M_llmut[:,i]
        site_llunmut = M_llunmut[:,i]
        ### Step3: Get all nodes/clusters we need
        # non-nan element index
        non_nan_indices_llmut = np.where(~np.isnan(site_llmut))[0]
        non_nan_indices_llunmut = np.where(~np.isnan(site_llunmut))[0]
        common_non_nan_indices = np.intersect1d(non_nan_indices_llmut, non_nan_indices_llunmut)
        all_nonnan_indices.append(common_non_nan_indices)
        M_nonnan = M_lowR[common_non_nan_indices,]
        site_llmut_nonnan = site_llmut[common_non_nan_indices]
        site_llunmut_nonnan = site_llunmut[common_non_nan_indices]
        ## Get low-resolution tree
        all_clusters = get_allBranchSet(M_nonnan)
        # full_leaf_clusters = get_leafBranchSet(all_clusters, pivot_clusters)
        pivot_clusters = get_1stBranchSet(all_clusters)
        no_pivot_clusters = [cluster for cluster in all_clusters if not any(np.array_equal(cluster, pivot_cluster) for pivot_cluster in pivot_clusters)]
        early_clusters = get_earlyBranchSet(no_pivot_clusters, pivot_clusters)
        if len(pivot_clusters)!=1:
            selected_clusters = pivot_clusters
        elif len(early_clusters)!=1:
            selected_clusters = early_clusters
        else:
            selected_clusters = all_clusters
        # # Fornmated low resolution tree 
        # M_pivot = np.array(selected_clusters).T
        ## Get all conditions clusters
        leaf_clusters = get_leafBranchSet(selected_clusters)
        internal_clusters = [cluster for cluster in selected_clusters if not any(np.array_equal(cluster, leaf_cluster) for leaf_cluster in leaf_clusters) and not any(np.array_equal(cluster, pivot_cluster) for pivot_cluster in pivot_clusters)]
        # Generate important clusters dictionary
        node_clusters_dict = {
            "leaf_clusters": leaf_clusters, 
            "internal_clusters": internal_clusters, 
            "pivot_clusters": pivot_clusters
        }
        ### Update somatic posterior per site
        ## prod results:
        each_newSP = get_newSomaticPosterior(site_llmut_nonnan, site_llunmut_nonnan, node_clusters_dict)
        # out posterior (mut, unmut, wrong)
        df_newSP_out = pd.concat([df_newSP_out, pd.DataFrame([each_newSP[1]])], ignore_index=True)
        withoutTree_posterior.append(DP_calSomaticPosterior_withoutTree(site_llmut_nonnan, site_llunmut_nonnan, len(site_llunmut_nonnan)))
    
    df_newSP_out.index = df_llunmut.columns
    df_newSP_out['nonnan_indices'] = all_nonnan_indices
    
    return [df_newSP_out, withoutTree_posterior]


# -------------------------
# Criteria of passing tree
# -------------------------

# Function to determine phylogeny_label based on the criteria
def determine_phylogeny_label_by_one_likelihood(row, pass_tree_cutoff, unpass_tree_cutoff):
    if row['mutant_cellnum'] == 1:
        return "cell_specific"
    elif row['mutant_cellnum'] == 0:
        return "absent"
    elif ((row['somatic_posterior_per_site'] >= pass_tree_cutoff and 
            row['somatic_posterior_per_site_onecell'] < unpass_tree_cutoff)):
        return "successful_pass"
    else:
        return "failed_pass"

def compare_elements_vectorized(val1, val2):
    """
    Vectorized comparison of two NumPy arrays and return the count for each flip type.
    """
    # Create a dictionary to store the flip count
    flip_counts = {'flipping_0_to_1': 0, 'flipping_1_to_0': 0,
                   'flipping_NA_to_0': 0, 'flipping_NA_to_1': 0}
    # Count the flip type
    flip_counts['flipping_NA_to_0'] = np.sum((val1 == 3) & (val2 == 0))
    flip_counts['flipping_NA_to_1'] = np.sum((val1 == 3) & (val2 == 1))
    flip_counts['flipping_1_to_0'] = np.sum((val1 == 1) & (val2 == 0))
    flip_counts['flipping_0_to_1'] = np.sum((val1 == 0) & (val2 == 1))
    return flip_counts




# -------------------------
# 4.3 Logistic Regression Classification
# -------------------------
# scdna_classifier.py
# scrna_classifier.py




# -------------------------
# 4.4 Placement of Potential Mosaic Mutations onto the Scaffold Tree
# -------------------------

def split_spots_by_immune_mutations(
    spots_to_split: list,  # List of spots to be split
    immune_mutations: list,  # List of immune mutations
    I_process: pd.DataFrame,  # I mutation matrix
    P_process: pd.DataFrame   # P matrix (posterior probability)
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Split the given spots_to_split into two versions based on immune mutations: 
    one with -immune suffix and the other with -non suffix.
    Also keep rows that do not need to be split.
    
    Parameters:
    ----------
    spots_to_split: list
        List of spots to be split (each spot is a row in the I_process matrix)
        
    immune_mutations: list
        List of immune mutations
    
    I_process: pd.DataFrame
        Mutation matrix (rows are spots, columns are mutations)
    
    P_process: pd.DataFrame
        Posterior probability matrix
    
    Returns:
    ----------
    I_process_resolved: pd.DataFrame
        Processed mutation matrix
        
    P_process_resolved: pd.DataFrame
        Processed posterior probability matrix
    """
    
    # Store split row data
    resolved_I_rows = []
    resolved_P_rows = []
    resolved_index = []
    
    # Iterate over each spot
    for spot_name, row in I_process.iterrows():
        # Check if this row needs to be split
        if spot_name in spots_to_split:
            # Get the original row for this spot
            I_row = row
            P_row = P_process.loc[spot_name]
            
            # Create -immune row, keep immune mutation columns, set others to NA
            immune_row = I_row.copy()
            immune_row[immune_mutations] = I_row[immune_mutations]  # Keep immune mutations
            immune_row[~I_row.index.isin(immune_mutations)] = np.nan  # Set non-immune mutations to NA
            
            P_immune_row = P_row.copy()
            P_immune_row[~P_row.index.isin(immune_mutations)] = np.nan  # Set non-immune mutations to NA
            
            # Create -non row, set immune mutation columns to NA, keep others as original
            non_immune_row = I_row.copy()
            non_immune_row[immune_mutations] = np.nan  # Set immune mutations to NA
            
            P_non_immune_row = P_row.copy()
            P_non_immune_row[immune_mutations] = np.nan  # Set immune mutations to NA
            
            # Add the split rows to the result lists
            resolved_I_rows.append(immune_row)
            resolved_P_rows.append(P_immune_row)
            resolved_index.append(f"{spot_name}-immune")
            
            resolved_I_rows.append(non_immune_row)
            resolved_P_rows.append(P_non_immune_row)
            resolved_index.append(f"{spot_name}-non")
        else:
            # Keep rows that don't need splitting as-is
            resolved_I_rows.append(row)
            resolved_P_rows.append(P_process.loc[spot_name])
            resolved_index.append(spot_name)
    
    # Combine the split rows into new DataFrames
    I_process_resolved = pd.DataFrame(resolved_I_rows, index=resolved_index, columns=I_process.columns)
    P_process_resolved = pd.DataFrame(resolved_P_rows, index=resolved_index, columns=P_process.columns)
    
    return I_process_resolved, P_process_resolved

def merge_mutations(M_current_each_mut, all_nodes_in_T_scaffold):
    # Function: Merge columns in M_current_each_mut according to all_nodes_in_T_scaffold
    # Create an empty list to store merged columns
    merged_columns = []
    
    # Iterate over all_nodes_in_T_scaffold
    for mutation in all_nodes_in_T_scaffold:
        # If the mutation consists of multiple mutations separated by '|'
        if '|' in mutation:
            # Split the mutation name to get multiple mutation columns
            mutations = mutation.split('|')
            # Get the data from these columns
            cols = [M_current_each_mut[col] for col in mutations if col in M_current_each_mut.columns]
            if cols:
                # Merge these columns using max (taking the maximum value, or modify as needed)
                merged_column = pd.concat(cols, axis=1).max(axis=1)
                merged_columns.append(merged_column)
            else:
                # If no matching columns are found (to avoid errors), fill with NaN
                merged_columns.append(pd.Series([pd.NA] * M_current_each_mut.shape[0], index=M_current_each_mut.index))
        else:
            # If it's a single mutation name, directly take the corresponding column
            if mutation in M_current_each_mut.columns:
                merged_columns.append(M_current_each_mut[mutation])
            else:
                # If the column is not found, add NaN
                merged_columns.append(pd.Series([pd.NA] * M_current_each_mut.shape[0], index=M_current_each_mut.index))
    
    # Combine the merged columns into a new DataFrame
    merged_df = pd.concat(merged_columns, axis=1)
    
    # Update column names to elements in all_nodes_in_T_scaffold
    merged_df.columns = all_nodes_in_T_scaffold
    
    return merged_df


def find_column_in_merged_columns(df, mutation_name):
    """
    Search for a mutation name in the DataFrame columns, supporting merged columns (separated by '|')
    
    Parameters:
    -----------
    df : pd.DataFrame
        The DataFrame to search
    mutation_name : str
        The mutation name to find
    
    Returns:
    --------
    str or None : The found column name, or None if not found
    """
    if mutation_name in df.columns:
        return mutation_name
    
    # Search in merged columns
    for col in df.columns:
        if '|' in col:
            muts_in_col = col.split('|')
            if mutation_name in muts_in_col:
                return col
    
    return None




# -------------------------
# find potiential positions that is intersected with new_mut based on whole tree for current tree
#   - intersection_nodes
#   - parent_dict
# -------------------------

import copy
from itertools import combinations
from typing import List, Dict, Any, Optional, Set, Tuple
import numpy as np


# ============================================================================
# Helper functions: Extract nodes and edges from TreeNode
# ============================================================================

def extract_nodes_edges_from_tree(tree):
    """
    Extract nodes and edges from a TreeNode object.
    
    Parameters
    ----------
    tree : TreeNode
        The tree to extract from
    
    Returns
    -------
    tuple : (nodes, edges)
        nodes: list of node dictionaries
        edges: list of (parent, child) tuples
    """
    nodes = []
    edges = []
    
    for node in tree.traverse():
        node_dict = {
            'name': node.name,
            'parent': node.parent.name if node.parent else None,
            'children': [child.name for child in node.children]
        }
        nodes.append(node_dict)
        if node.parent:
            edges.append((node.parent.name, node.name))
    
    return nodes, edges


def build_tree_parent_dict_from_nodes(nodes):
    """
    Build a parent-child dictionary from a nodes list.
    
    Parameters
    ----------
    nodes : list
        List of node dictionaries
    
    Returns
    -------
    dict
        Mapping of child name to parent name
    """
    parent_dict = {}
    for node in nodes:
        if node['parent'] is not None:
            parent_dict[node['name']] = node['parent']
    return parent_dict


def build_tree_parent_dict(tree):
    """
    Build a parent-child dictionary directly from the tree (legacy compatibility).
    
    Parameters
    ----------
    tree : TreeNode
        The tree to build from
    
    Returns
    -------
    dict
        Mapping of child name to parent name
    """
    parent_dict = {}
    for node in tree.traverse():
        for child in node.children:
            parent_dict[child.name] = node.name
    return parent_dict


def build_parent_dict_from_candidates(candidate_positions):
    """
    Build a parent-child dictionary from candidate positions.
    Used for subsequent tree construction and validation.
    
    Parameters
    ----------
    candidate_positions : list
        List of candidate position dictionaries
    
    Returns
    -------
    dict
        Mapping of child name to parent name
    """
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
        
        # on_node and new_leaf types do not need parent_dict updates
        # because their parent-child relationships are already in the tree structure
    
    return parent_dict


# ============================================================================
# Core candidate position generation functions (following original logic)
# ============================================================================

def _create_on_node_candidate_fast(base_tree_copy, node, new_mut):
    """
    Quickly create a candidate for placing on a node.
    Original logic: new mutation replaces the original node (becomes its parent).
    If the original node is ROOT, the new mutation becomes the new root.
    
    Parameters
    ----------
    base_tree_copy : TreeNode
        Deep copy of the base tree
    node : TreeNode
        The node to place on
    new_mut : str
        New mutation name
    
    Returns
    -------
    dict
        Candidate position dictionary
    """
    new_tree = copy.deepcopy(base_tree_copy)
    anchor_node = new_tree.find(node.name)
    
    # Replace node logic
    new_node = TreeNode(new_mut)
    for child in anchor_node.children:
        new_node.add_child(copy.deepcopy(child))
    
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
        "nodes": _extract_nodes_info(new_tree),
        "edges": _extract_edges_info(new_tree)
    }


def _create_new_leaf_candidate_fast(base_tree_copy, node, new_mut):
    """
    Quickly create a candidate for adding a new leaf.
    Original logic: new mutation becomes a child of the original node.
    
    Parameters
    ----------
    base_tree_copy : TreeNode
        Deep copy of the base tree
    node : TreeNode
        The parent node
    new_mut : str
        New mutation name
    
    Returns
    -------
    dict
        Candidate position dictionary
    """
    new_tree = copy.deepcopy(base_tree_copy)
    anchor_node = new_tree.find(node.name)
    new_leaf = TreeNode(new_mut)
    anchor_node.add_child(new_leaf)
    
    return {
        "placement_type": "new_leaf",
        "anchor": node.name,
        "meta": {},
        "nodes": _extract_nodes_info(new_tree),
        "edges": _extract_edges_info(new_tree)
    }


def _create_on_edge_candidate_fast(base_tree_copy, parent_node, child_node, new_mut):
    """
    Quickly create a candidate for placing on an edge.
    Original logic: insert a new node between parent and child.
    
    Parameters
    ----------
    base_tree_copy : TreeNode
        Deep copy of the base tree
    parent_node : TreeNode
        The parent node
    child_node : TreeNode
        The child node
    new_mut : str
        New mutation name
    
    Returns
    -------
    dict
        Candidate position dictionary
    """
    new_tree = copy.deepcopy(base_tree_copy)
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
        "nodes": _extract_nodes_info(new_tree),
        "edges": _extract_edges_info(new_tree)
    }


def _create_merge_candidate_fast(base_tree_copy, parent_node, children_combo, new_mut):
    """
    Quickly create a candidate for merging child nodes.
    Original logic: insert a new node between parent and multiple children.
    
    Parameters
    ----------
    base_tree_copy : TreeNode
        Deep copy of the base tree
    parent_node : TreeNode
        The parent node
    children_combo : list
        List of child nodes to merge
    new_mut : str
        New mutation name
    
    Returns
    -------
    dict
        Candidate position dictionary
    """
    new_tree = copy.deepcopy(base_tree_copy)
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
        "nodes": _extract_nodes_info(new_tree),
        "edges": _extract_edges_info(new_tree)
    }


def _extract_nodes_info(tree):
    """
    Extract node information from a tree.
    
    Parameters
    ----------
    tree : TreeNode
        The tree to extract from
    
    Returns
    -------
    list
        List of node dictionaries with name, parent, and children
    """
    return [{"name": n.name,
             "parent": n.parent.name if n.parent else None,
             "children": [c.name for c in n.children]} 
            for n in tree.traverse()]


def _extract_edges_info(tree):
    """
    Extract edge information from a tree.
    
    Parameters
    ----------
    tree : TreeNode
        The tree to extract from
    
    Returns
    -------
    list
        List of (parent, child) tuples
    """
    return [(n.parent.name, n.name) for n in tree.traverse() if n.parent]


# ============================================================================
# Find intersection nodes
# ============================================================================

def find_all_intersect_muts_from_tree_by_matrix(
    tree, matrix, target_mut, min_overlap=1,
    matrix_positive=None,
    matrix_positive_values=None,
    matrix_positive_col_to_idx=None,
    tree_nodes=None,
    node_to_mutations=None,
    logger_obj=None,
):
    """
    Return all tree nodes that intersect with target_mut in the matrix
    with at least min_overlap cells co-occurring.
    
    Parameters
    ----------
    tree : TreeNode
        The tree to search
    matrix : pd.DataFrame
        Mutation presence matrix
    target_mut : str
        Target mutation name
    min_overlap : int, default=1
        Minimum number of cells where both mutations are 1
    matrix_positive : pd.DataFrame, optional
        Pre-computed boolean matrix of presence
    matrix_positive_values : np.ndarray, optional
        Pre-computed numpy array of presence
    matrix_positive_col_to_idx : dict, optional
        Column name to index mapping
    tree_nodes : list, optional
        Pre-computed list of tree nodes
    node_to_mutations : dict, optional
        Pre-computed mapping of node names to mutation lists
    logger_obj : logging.Logger, optional
        Logger instance for logging messages. If None, uses global logger.
    
    Returns
    -------
    set
        Set of node names that intersect with target_mut
    """
    # Use provided logger or fall back to global logger
    log = logger_obj if logger_obj is not None else logging.getLogger(__name__)
    
    intersect_nodes = set()
    if target_mut not in matrix.columns:
        log.debug(f"Target mutation {target_mut} not found in matrix columns")
        return intersect_nodes
    
    if matrix_positive is None and matrix_positive_values is None:
        matrix_positive = matrix.eq(1).fillna(False)
    if matrix_positive_values is None:
        matrix_positive_values = matrix_positive.to_numpy(dtype=bool, copy=False)
    if matrix_positive_col_to_idx is None:
        source_columns = matrix_positive.columns if matrix_positive is not None else matrix.columns
        matrix_positive_col_to_idx = {col: idx for idx, col in enumerate(source_columns)}
    
    target_idx = matrix_positive_col_to_idx.get(target_mut)
    if target_idx is None:
        log.debug(f"Target mutation {target_mut} not found in matrix_positive columns")
        return intersect_nodes
    
    target_positive = matrix_positive_values[:, target_idx]
    if not bool(target_positive.any()):
        log.debug(f"Target mutation {target_mut} has no positive cells")
        return intersect_nodes
    
    overlap_counts = matrix_positive_values[target_positive].sum(axis=0, dtype=np.int64)
    nodes_to_scan = tree_nodes if tree_nodes is not None else tree.all_nodes()
    
    for node in nodes_to_scan:
        node_name = node.name if hasattr(node, "name") else node
        if node_name == "ROOT":
            continue
        
        mutations_on_node = node_to_mutations.get(node_name, ()) if node_to_mutations is not None else tuple(node_name.split("|"))
        for mut in mutations_on_node:
            if mut == target_mut:
                continue
            mut_idx = matrix_positive_col_to_idx.get(mut)
            if mut_idx is not None and int(overlap_counts[mut_idx]) >= min_overlap:
                intersect_nodes.add(node_name)
                break
    
    log.debug(f"Found {len(intersect_nodes)} intersection nodes for mutation {target_mut}")
    
    return intersect_nodes


def find_all_path_nodes(intersection_nodes, parent_dict, logger_obj=None):
    """
    Find all nodes that lie on paths connecting intersection nodes.
    
    Parameters
    ----------
    intersection_nodes : set or list
        Set/list of intersection node names
    parent_dict : dict
        Parent-child dictionary (tree_parent_dict)
    logger_obj : logging.Logger, optional
        Logger instance for logging messages. If None, uses global logger.
    
    Returns
    -------
    list
        Sorted list of all node names on paths between intersection nodes
    """
    # Use provided logger or fall back to global logger
    log = logger_obj if logger_obj is not None else logging.getLogger(__name__)
    
    all_path_nodes = set()
    all_path_nodes.add('ROOT')  # Always include ROOT
    
    # For each intersection node, find the path from ROOT to that node
    for node in intersection_nodes:
        path_to_root = get_path_to_root(node, parent_dict)
        all_path_nodes.update(path_to_root)
    
    # Find paths connecting different intersection nodes
    intersection_list = list(intersection_nodes)
    for i in range(len(intersection_list)):
        for j in range(i + 1, len(intersection_list)):
            path_between = get_path_between_nodes(intersection_list[i], intersection_list[j], parent_dict)
            all_path_nodes.update(path_between)
    
    # ============================================================
    # CRITICAL FIX: Return sorted list to ensure deterministic iteration order
    # ============================================================
    return sorted(all_path_nodes)


# ============================================================================
# Main function: Find candidate positions (following original logic)
# ============================================================================

def find_intersection_positions_within_tree_directly(
    T_current,
    new_mut: str,
    matrix,
    min_overlap: int = 1,
    intersection_nodes: Optional[Set] = None,
    tree_parent_dict: Optional[Dict] = None,
    node_lookup: Optional[Dict] = None,
    logger_obj=None,
):
    """
    Optimized version based on intersection analysis, directly finding relevant positions.
    Follows the original logic using TreeNode objects.
    """
    # Use provided logger or fall back to global logger
    log = logger_obj if logger_obj is not None else logging.getLogger(__name__)
    
    log.debug(f"Finding intersection positions for mutation: {new_mut}")
    
    # 1. Find all nodes that intersect with the target mutation
    if intersection_nodes is None:
        intersection_nodes = find_all_intersect_muts_from_tree_by_matrix(
            T_current, matrix, new_mut, min_overlap, logger_obj=log
        )
    
    if len(intersection_nodes) == 0:
        log.debug(f"No intersection nodes found for {new_mut}")
        return []
    
    # 2. Build tree parent-child dictionary
    if tree_parent_dict is None:
        tree_parent_dict = build_tree_parent_dict(T_current)
    
    # 3. Find all path nodes (returns sorted list for deterministic iteration)
    all_path_nodes = find_all_path_nodes(intersection_nodes, tree_parent_dict, logger_obj=log)
    
    # 4. Pre-create a deep copy of the base tree
    base_tree_copy = copy.deepcopy(T_current)
    
    # 5. Generate candidate positions only on relevant nodes
    candidate_positions = []
    
    # ============================================================
    # CRITICAL FIX: Iterate over sorted all_path_nodes to ensure deterministic order
    # ============================================================
    for node_name in all_path_nodes:
        if node_name == "ROOT":
            node = base_tree_copy.find(node_name)
            if node is None:
                continue
            candidate_positions.append(_create_on_node_candidate_fast(base_tree_copy, node, new_mut))
            continue
        
        node = base_tree_copy.find(node_name)
        if node is None:
            continue
        
        candidate_positions.append(_create_on_node_candidate_fast(base_tree_copy, node, new_mut))
        candidate_positions.append(_create_new_leaf_candidate_fast(base_tree_copy, node, new_mut))
        
        for child in node.children:
            if child.name in all_path_nodes:
                candidate_positions.append(_create_on_edge_candidate_fast(base_tree_copy, node, child, new_mut))
        
        if len(node.children) >= 2:
            path_children = [child for child in node.children if child.name in all_path_nodes]
            if len(path_children) >= 2:
                for r in range(2, min(4, len(path_children) + 1)):
                    for combo in combinations(path_children, r):
                        candidate_positions.append(_create_merge_candidate_fast(base_tree_copy, node, combo, new_mut))
    
    log.debug(f"Generated {len(candidate_positions)} candidate positions for {new_mut}")
    
    return candidate_positions


# ============================================================================
# Debug function: Print tree structure
# ============================================================================

def print_tree_from_nodes(nodes, edges, root_name="ROOT", logger_obj=None):
    """
    Print tree structure from nodes and edges using logger.
    
    Parameters
    ----------
    nodes : list
        List of node dictionaries
    edges : list
        List of (parent, child) tuples
    root_name : str, default="ROOT"
        Name of the root node
    logger_obj : logging.Logger, optional
        Logger instance for logging messages. If None, uses global logger.
    """
    # Use provided logger or fall back to global logger
    log = logger_obj if logger_obj is not None else logging.getLogger(__name__)
    
    def print_node(name, indent=0):
        node = next((n for n in nodes if n['name'] == name), None)
        if node is None:
            return
        log.info("  " * indent + f"└─ {name}")
        for child in node['children']:
            print_node(child, indent + 1)
    
    log.info(f"Tree structure (root: {root_name}):")
    print_node(root_name)


# # ============================================================================
# # Usage
# # ============================================================================
# # Call in the original way
# potential_positions = find_intersection_positions_within_tree_directly(
#     T_current=T_current,      # TreeNode object
#     new_mut=new_mut,          # New mutation name
#     matrix=I_attached,        # Mutation matrix
#     min_overlap=1,            # Minimum overlapping cells
#     intersection_nodes=intersection_nodes,  # Optional
#     tree_parent_dict=current_tree_state["tree_parent_dict"],  # Optional
#     node_lookup=current_tree_state["node_lookup"],  # Optional
#     logger_obj=cv_logger,     # CV-specific logger
# )

# # Print results (now using logger)
# for pos in potential_positions:
#     log.info(f"{pos['placement_type']} on {pos['anchor']}")


def get_path_to_root(node, tree_parent_dict):
    """Find the path from a node to ROOT"""
    path = []
    current = node
    while current in tree_parent_dict:
        path.append(current)
        current = tree_parent_dict[current]
        if current == 'ROOT':
            path.append('ROOT')
            break
    return path


def get_path_between_nodes(node1, node2, tree_parent_dict):
    """Find the path between two nodes"""
    # Find the path from node1 to ROOT
    path1 = get_path_to_root(node1, tree_parent_dict)
    # Find the path from node2 to ROOT  
    path2 = get_path_to_root(node2, tree_parent_dict)
    
    # Reverse paths so they start from ROOT
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
# Compute which clone new_mut most likely belongs to by repeating 100 times
# -------------------------

def compute_corr_cache_with_new_mut(I_attached, existing_muts, new_mut):
    """
    Compute correlation cache including the new mutation
    """
    # All mutations (existing + new)
    all_muts = existing_muts + [new_mut]
    
    corr_cache = {}
    
    # Compute correlations for all mutation pairs (including new with existing)
    for u, v in itertools.combinations(all_muts, 2):
        corr = are_mutations_correlated(I_attached, u, v)
        corr_cache[(u, v)] = corr
        corr_cache[(v, u)] = corr
    
    # Self-correlation is always True
    for m in all_muts:
        corr_cache[(m, m)] = True
    
    return corr_cache

def compute_new_mut_clone_affinity_correct(new_mut, mutation_clones_rescue, I_attached, n_shuffle=100, logger_obj=None):
    """
    Correct version: uses corr_cache that includes the new mutation
    """
    
    log = logger_obj if logger_obj is not None else logging.getLogger(__name__)
    
    # Get all existing mutations
    all_existing_muts = []
    for clone in mutation_clones_rescue.values():
        all_existing_muts.extend(clone)
    all_existing_muts = list(set(all_existing_muts))
    
    # Recompute corr_cache including the new mutation
    corr_cache = compute_corr_cache_with_new_mut(I_attached, all_existing_muts, new_mut)
    
    clone_affinity = {}
    detailed_scores = {}
    
    log.info(f"Calculating correlations with tree mutations for: {new_mut}")
    
    for clone_rep, clone_muts in mutation_clones_rescue.items():
        clone_key = tuple(sorted(clone_muts))
        
        # Check correlations between new mutation and clone members
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
        
        # Calculate relationship with clone representative mutation
        rep_correlation = 0
        if clone_rep in clone_muts:
            key1 = (new_mut, clone_rep)
            key2 = (clone_rep, new_mut)
            rep_correlation = 1 if (corr_cache.get(key1, False) or corr_cache.get(key2, False)) else 0
            counts = pairwise_counts(I_attached, new_mut, clone_rep)
        
        # Calculate weighted score (considering clone size)
        base_score = direct_ratio
        size_factor = min(1.0, len(clone_muts) / 10)  # Avoid overly large clones having too much influence
        final_score = base_score * (0.7 + 0.3 * size_factor)
        
        clone_affinity[clone_key] = final_score
        detailed_scores[clone_key] = {
            'direct_ratio': direct_ratio,
            'direct_correlations': direct_corr_count,
            'rep_correlation': rep_correlation,
            'clone_size': len(clone_muts)
        }
    
    return clone_affinity, detailed_scores

def select_max_affinity_clone(clone_affinity):
    # Filter to keep only entries with score > 0
    filtered_affinities = {k: v for k, v in clone_affinity.items() if v > 0}
    
    # If no entries with score > 0, return empty list
    if not filtered_affinities:
        return []
    
    # Find the maximum affinity score
    max_affinity = max(filtered_affinities.values())
    
    # Filter all clones with the maximum affinity score
    max_clones = [k for k, v in filtered_affinities.items() if v == max_affinity]
    
    return max_clones

def select_best_clone(detailed_scores):
    # Step 1: Find the maximum rep_correlation value
    max_rep_correlation = max([score['rep_correlation'] for score in detailed_scores.values()])
    
    # If max rep_correlation is 0, return empty list directly
    if max_rep_correlation == 0:
        return []
    
    # Step 2: Find all clones with the maximum rep_correlation
    max_rep_correlation_clones = [
        clone for clone, score in detailed_scores.items() if score['rep_correlation'] == max_rep_correlation
    ]
    
    # Step 3: If only one clone, return it directly
    if len(max_rep_correlation_clones) == 1:
        return max_rep_correlation_clones
    
    # Step 4: If multiple clones, compare their direct_ratio
    max_direct_ratio = max([detailed_scores[clone]['direct_ratio'] for clone in max_rep_correlation_clones])
    
    # If max direct_ratio is 0, return empty list
    if max_direct_ratio == 0:
        return []
    
    # Step 5: Return clones with the maximum direct_ratio
    best_clones = [
        clone for clone in max_rep_correlation_clones 
        if detailed_scores[clone]['direct_ratio'] == max_direct_ratio
    ]
    
    return best_clones




# -------------------------
# Calculate Bayesian penalty score
# -------------------------

def compute_bayesian_penalty_for_all_positions_consider_ROOT(
    new_mut, selected_positions, T_current, M_current, I_selected, P_selected, parent_dict, intersection_nodes, 
    ω_NA=0.001, fnfp_ratio=0.1, φ=1, logger_obj=None
):
    """
    Calculate penalties for all candidate positions (optimized version - consider_ROOT).
    
    Optimization strategies:
    1. Cache conflict masks (most time-consuming operation)
    2. Use NumPy vectorization instead of loops
    3. Cache M_current and I_selected column data
    4. Automatically decide whether to use parallel execution based on number of positions
    
    Parameters
    ----------
    new_mut : str
        New mutation to be processed
    selected_positions : list
        List of candidate positions, each containing placement_type, anchor, etc.
    T_current : TreeNode
        Current phylogenetic tree structure
    M_current : pd.DataFrame
        Current mutation matrix, index: cells, columns: mutations
    I_selected : pd.DataFrame
        Mutation presence matrix, index: cells, columns: mutations
    P_selected : pd.DataFrame
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
    df_penalty : pd.DataFrame
        DataFrame containing penalties for all candidate positions
    """
    # Use provided logger or fall back to global logger
    log = logger_obj if logger_obj is not None else logging.getLogger(__name__)
    
    import pandas as pd
    import numpy as np
    from joblib import Parallel, delayed
    
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
    
    # Cache M_current columns (convert to NumPy arrays)
    m_current_cache = {}
    for col in matrix_columns:
        if col != 'ROOT':
            m_current_cache[col] = M_current[col].reindex(matrix_index, fill_value=0).to_numpy(dtype=np.int8)
    
    # Cache I_selected columns (convert to NumPy arrays)
    i_selected_cache = {}
    for col in I_selected.columns:
        i_selected_cache[col] = I_selected[col].reindex(matrix_index, fill_value=np.nan).to_numpy(dtype=np.float64)
    
    # Cache conflict masks (most time-consuming operation)
    conflict_cache = {}
    
    def get_conflict_mask(anchor, sibling_nodes, exclude_nodes=None):
        """Get conflict mask with caching"""
        if exclude_nodes is None:
            exclude_nodes = []
        
        # Build cache key
        cache_key = (anchor, tuple(sorted(sibling_nodes)), tuple(sorted(exclude_nodes)))
        
        if cache_key in conflict_cache:
            return conflict_cache[cache_key]
        
        # Calculate conflict mask (using consider_ROOT version)
        lineage_parent = build_lineage_parent_dict_from_tree(T_current, anchor)
        lineage_conflict_nodes = get_all_conflict_nodes_outside_lineage(
            anchor, lineage_parent, matrix_columns, exclude_nodes=exclude_nodes
        )
        
        all_conflict_nodes = list(set(sibling_nodes + lineage_conflict_nodes))
        
        # Vectorized conflict mask construction
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
            chain_cache[anchor] = get_full_mutnode_chain_with_anchor(anchor, parent_dict)
        return chain_cache[anchor]
    
    # Cache node mutations
    node_mutations_cache = {}
    def get_node_mutations(node_name):
        if node_name not in node_mutations_cache:
            node_mutations_cache[node_name] = node_name.split("|")
        return node_mutations_cache[node_name]
    
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
        
        # ============================================================
        # Calculate penalty (using consider_ROOT version)
        # ============================================================
        
        full_mutnode_chain = get_full_chain(anchor)
        
        posterior_vec = P_selected[new_mut]
        input_binary_vec = I_selected[new_mut]
        
        new_mut_penalty, actual_na_flip_ratio, refined_ω_NA, φ_adjusted, weight_na_to_1, weight_na_to_0 = compute_dynamic_penalty(
            input_binary_vec, posterior_vec, imputed_vec, fnfp_ratio, ω_NA, φ,
            na_ratio, mut_ratio, placement_type, N_nodes
        )
        
        # Calculate chain penalty (vectorized optimization)
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
            
            # Get node vector
            if node in m_current_cache:
                node_array = m_current_cache[node]
            else:
                node_array = M_current[node].reindex(matrix_index, fill_value=0).to_numpy(dtype=np.float64)
            
            # Vectorized find cells to flip
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
                    weight_na_to_1,
                    weight_na_to_0,
                    fnfp_ratio
                )
                chain_penalty += mut_penalty
        
        total_chain_penalty = new_mut_penalty + chain_penalty
        
        log_N_nodes_penalty = np.log(N_nodes)
        
        # ============================================================
        # Key difference: consider_ROOT uses original φ (no dynamic adjustment)
        # ============================================================
        BIC_penalty = φ * np.log(N_nodes)  # Note: φ, not φ_adjusted
        
        root_penalty = 0
        if anchor == 'ROOT':
            root_penalty = np.log(N_nodes) * 0.5
        
        base_total_penalty = total_chain_penalty + log_N_nodes_penalty + BIC_penalty + merge_penalty + root_penalty
        
        # Use consider_ROOT version functions
        intersection_penalty = compute_intersection_based_penalty(
            new_mut, pos, intersection_nodes, M_current, I_selected, na_ratio, mut_ratio, actual_na_flip_ratio
        )
        
        hierarchy_penalty = compute_hierarchy_penalty(
            new_mut, pos, M_current, I_selected, parent_dict, na_ratio, mut_ratio, actual_na_flip_ratio
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
    
    if n_positions <= 20:
        # Serial execution (with cache optimization)
        results = [compute_single_position(idx, pos) for idx, pos in enumerate(selected_positions)]
    else:
        # Parallel execution (with cache optimization)
        n_jobs = -1
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


def apply_position_to_tree(
    new_mut, position, imputed_vec, T_current, M_current, I_selected, parent_dict
):
    """
    Apply the specific position to trees and matrices
    """
    M_updated = M_current.copy()
    T_updated = T_current.copy()
    
    anchor = position['anchor']
    placement_type = position['placement_type']
    
    # Get full mutation chain from anchor to ROOT
    full_mutnode_chain = get_full_mutnode_chain_with_anchor(anchor, parent_dict)
    
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
    
    return T_updated, M_updated




# Keep existing helper functions, as they are already defined in the first function
def get_all_conflict_nodes_outside_lineage(anchor, parent_dict, all_columns, exclude_nodes=None):
    """
    Get all potential conflict nodes outside the lineage
    
    Parameters:
    - anchor: Current anchor node
    - parent_dict: Parent dictionary
    - all_columns: All mutation nodes
    - exclude_nodes: List of nodes to exclude (e.g., merge_children)
    
    Returns:
    - conflict_nodes: List of all nodes outside the lineage
    """
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
    
    return conflict_nodes


def get_full_mutnode_chain_with_anchor(anchor, parent_dict):
    """
    Get the complete mutation chain from anchor to ROOT (including the anchor itself)
    """
    mutation_chain = [anchor]
    current = anchor
    
    while True:
        parent = parent_dict.get(current, None)
        if parent is None or parent == 'ROOT':
            break
        mutation_chain.append(parent)
        current = parent
    
    return mutation_chain


def compute_dynamic_penalty(input_binary_vec, posterior_vec, imputed_vec, fnfp_ratio, base_ω_NA, base_φ, na_ratio, mut_ratio, placement_type, N_nodes):
    """
    Dynamically adjust all penalty weights based on the actual NA→1 flip ratio
    """
    # Calculate the actual NA→1 flip ratio
    na_mask = input_binary_vec.isna()
    na_to_one_count = ((imputed_vec == 1) & na_mask).sum()
    total_na_count = na_mask.sum()
    
    if total_na_count > 0:
        actual_na_flip_ratio = na_to_one_count / total_na_count
    else:
        actual_na_flip_ratio = 0
    
    # Calculate w_fn value
    w_fn = fnfp_ratio * 1  # default w_fn = 0.1
    
    # More refined dynamic adjustment: make w_na approach w_fn when flip ratio is high
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
    
    # Ensure it does not exceed w_fn (since NA flips are inherently less reliable than FN)
    refined_ω_NA = min(refined_ω_NA, w_fn * 0.95)
    
    # Dynamically adjust BIC penalty based on NA characteristics and actual flip ratio
    φ_adjusted = compute_dynamic_bic_penalty(
        base_φ, na_ratio, mut_ratio, actual_na_flip_ratio, placement_type, N_nodes
    )
    
    # Compute penalty using adjusted weights
    penalty, weight_na_to_1, weight_na_to_0 = compute_bayesian_penalty_each_pos(
        input_binary_vec, posterior_vec, imputed_vec, fnfp_ratio, refined_ω_NA
    )
    
    return penalty, actual_na_flip_ratio, refined_ω_NA, φ_adjusted, weight_na_to_1, weight_na_to_0


def compute_dynamic_bic_penalty(base_φ, na_ratio, mut_ratio, actual_na_flip_ratio, placement_type, N_nodes):
    """
    Dynamically adjust BIC penalty based on NA characteristics and actual flip ratio
    """
    # Base adjustment: high NA + low mutation rate -> reduce BIC penalty
    if na_ratio > 0.7 and mut_ratio < 0.1:
        base_adjustment = 0.3
    elif na_ratio > 0.5 and mut_ratio < 0.2:
        base_adjustment = 0.6
    else:
        base_adjustment = 1.0
    
    # Further adjustment based on actual flip ratio
    if actual_na_flip_ratio > 0.4:
        # If actual flip ratio is very high, this placement may be unreliable
        # For new node creation, be more lenient (reduce penalty)
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
    
    # Node count consideration: more nodes -> more lenient penalty for new node creation
    node_adjustment = 1.0
    if N_nodes > 20 and placement_type in ['new_leaf', 'on_edge']:
        node_adjustment = 0.9
    elif N_nodes > 50 and placement_type in ['new_leaf', 'on_edge']:
        node_adjustment = 0.8
    
    φ_adjusted = base_φ * base_adjustment * flip_adjustment * node_adjustment
    
    return max(φ_adjusted, 0.1)  # Ensure it doesn't go too low


def compute_intersection_based_penalty(new_mut, position, intersection_nodes, M_current, I_selected, na_ratio, mut_ratio, actual_na_flip_ratio):
    """
    Adjust penalty based on intersection patterns and actual NA flips
    """
    extra_penalty = 0
    placement_type = position['placement_type']
    anchor = position['anchor']
    
    # Get mutant cells for new_mut
    mut_cells = set(I_selected.index[I_selected[new_mut].fillna(0) == 1])
    
    if placement_type == 'on_node':
        anchor_cells = set(M_current.index[M_current[anchor] == 1])
        
        # Case 1: If new_mut strongly intersects with only one node, but is placed on a different early node
        if len(intersection_nodes) == 1:
            sole_intersection = list(intersection_nodes)[0]
            if anchor != sole_intersection and anchor in [n for n in M_current.columns if n != 'ROOT']:
                # Adjust penalty based on actual flip ratio
                flip_based_multiplier = 1.0 + actual_na_flip_ratio  # More flips -> heavier penalty
                extra_penalty += np.log(len(M_current.columns)) * 0.8 * flip_based_multiplier
        
        # Case 2: If anchor covers far more cells than new_mut's mutant cells
        if len(anchor_cells) > len(mut_cells) * 3:
            # Adjust FN penalty based on actual flip ratio
            fn_penalty = np.log(len(anchor_cells) - len(mut_cells)) * 0.3
            extra_penalty += fn_penalty * (1.0 + actual_na_flip_ratio)
            
        # Case 3: High NA mutation placed on a non-intersecting node with high actual flip ratio
        if na_ratio > 0.7 and anchor not in intersection_nodes and actual_na_flip_ratio > 0.3:
            extra_penalty += np.log(len(M_current.columns)) * 0.5
    
    elif placement_type in ['new_leaf', 'on_edge']:
        # For high NA + low mutation rate with low actual flip ratio, reduce penalty for creating new nodes
        if na_ratio > 0.7 and mut_ratio < 0.1 and actual_na_flip_ratio < 0.2:
            bonus = -np.log(len(M_current.columns)) * 0.4 * (1.0 - actual_na_flip_ratio)
            extra_penalty += bonus
            
        # For mutations intersecting with multiple nodes, reduce penalty for creating new nodes
        if len(intersection_nodes) >= 2:
            bonus = -np.log(len(M_current.columns)) * 0.3
            extra_penalty += bonus
    
    return extra_penalty


def compute_hierarchy_penalty(new_mut, position, M_current, I_selected, parent_dict, 
                            na_ratio, mut_ratio, actual_na_flip_ratio):
    """
    Compute hierarchy合理性 penalty, considering actual NA flips
    """
    hierarchy_penalty = 0
    placement_type = position['placement_type']
    anchor = position['anchor']
    
    mut_cells = set(I_selected.index[I_selected[new_mut].fillna(0) == 1])
    
    if placement_type == 'on_node':
        anchor_cells = set(M_current.index[M_current[anchor] == 1])
        
        # Check if this anchor has children
        anchor_children = find_children_of_node(anchor, M_current.columns, parent_dict)
        
        if anchor_children:
            # If anchor has children, and new_mut has a specific pattern with children
            children_intersection = False
            for child in anchor_children:
                child_cells = set(M_current.index[M_current[child] == 1])
                if mut_cells.issubset(child_cells) and len(mut_cells) < len(child_cells) * 0.8:
                    children_intersection = True
                    break
            
            if children_intersection:
                # Adjust hierarchy不合理 penalty based on actual flip ratio
                flip_multiplier = 1.0 + actual_na_flip_ratio
                hierarchy_penalty += np.log(len(anchor_children) + 1) * 0.4 * flip_multiplier
    
    elif placement_type == 'new_leaf':
        parent = anchor
        parent_cells = set(M_current.index[M_current[parent] == 1])
        
        # If new_mut's cells are a proper subset of parent cells, this is a reasonable subclone
        if mut_cells.issubset(parent_cells) and len(mut_cells) < len(parent_cells) * 0.8:
            # Give a reward, with magnitude depending on actual flip ratio (fewer flips -> more reward)
            bonus_multiplier = 1.0 - actual_na_flip_ratio  # Fewer flips -> more reward
            hierarchy_penalty -= np.log(len(parent_cells) - len(mut_cells)) * 0.2 * bonus_multiplier
            
        # For high NA mutations with low actual flip ratio, further reward reasonable subclone placement
        if na_ratio > 0.7 and mut_cells.issubset(parent_cells) and actual_na_flip_ratio < 0.2:
            bonus_multiplier = 1.0 - actual_na_flip_ratio
            hierarchy_penalty -= np.log(len(parent_cells)) * 0.3 * bonus_multiplier
    
    return hierarchy_penalty


def find_children_of_node(node, all_columns, parent_dict):
    """Find direct children of a node"""
    children = []
    for col in all_columns:
        if col == node or col == 'ROOT':
            continue
        if parent_dict.get(col) == node:
            children.append(col)
    return children



# -------------------------
# Calculate fp_ratio and fn_ratio
# -------------------------

##### Calculate across the entire tree
def calculate_fp_fn_ratios_across_tree(M_checkpoint, I_attached):
    """
    Calculate fp_ratio and fn_ratio for each mutation
    
    Parameters:
    M_checkpoint: Imputed genotype matrix (binary, 0/1)
    I_attached: Original input genotype matrix (0, 1, NaN)
    
    Returns:
    DataFrame: Contains fp_ratio and fn_ratio for each mutation
    dict: Contains list of other mutations with fp>0 for each mutation
    """
    mutations = M_checkpoint.columns
    results = []
    fp_mutations_dict = {}  # Store other mutations with fp>0 for each mutation
    
    for mut in mutations:
        # 1. Calculate fp_ratio
        mutant_cells = M_checkpoint[M_checkpoint[mut] == 1].index
        fp_ratios = []
        fp_positive_muts = []  # Store other mutations with fp>0 for current mutation
        
        for other_mut in mutations:
            if other_mut == mut:
                continue
                
            # Calculate FP for other mutations within the current mutation's clone
            other_original = I_attached.loc[mutant_cells, other_mut]
            other_imputed = M_checkpoint.loc[mutant_cells, other_mut]
            # FP: original is 1 but imputed is 0
            fp_count = ((other_original == 1) & (other_imputed == 0)).sum()
            
            # Calculate FP for other mutations across all cells
            total_other_original = I_attached.loc[:, other_mut]
            total_other_imputed = M_checkpoint.loc[:, other_mut]
            total_fp_count = ((total_other_original == 1) & (total_other_imputed == 0)).sum()
            
            if total_fp_count > 0:
                fp_ratio = fp_count / total_fp_count
                fp_ratios.append(fp_ratio)
                
                # If fp_ratio > 0, record this other_mut
                if fp_ratio > 0:
                    fp_positive_muts.append(other_mut)
            else:
                fp_ratios.append(0)
        
        # Store other mutations with fp>0 for the current mutation
        fp_mutations_dict[mut] = fp_positive_muts
        
        avg_fp_ratio = np.mean(fp_ratios) if fp_ratios else 0
        
        # 2. Calculate fn_ratio (unchanged)
        mut_original = I_attached.loc[mutant_cells, mut]
        mut_imputed = M_checkpoint.loc[mutant_cells, mut]
        
        # Cells where this mutation is 0 in the original data
        zero_cells = mut_original[mut_original == 0].index
        
        # Calculate cells with intersections with other mutations
        I_attached_sub_mut = I_attached.loc[mutant_cells, :]
        row_counts = I_attached_sub_mut.eq(1).sum(axis=1)
        intersect_cells = mut_original[(mut_original == 1) & (row_counts > 1)]
        
        total_zeros = len(zero_cells)
        total_intersect = len(intersect_cells)
        
        if (total_zeros + total_intersect) > 0:
            fn_ratio = total_zeros / (total_zeros + total_intersect)
        else:
            fn_ratio = 0
        
        results.append({
            'identifier': mut,
            'fp_ratio': avg_fp_ratio,
            'fn_ratio': fn_ratio
        })
    
    return pd.DataFrame(results), fp_mutations_dict

# # Usage
# ratios_df, fp_mutations_dict = calculate_fp_fn_ratios_across_tree(M_checkpoint, I_attached)
# logger.info(ratios_df)
# logger.info(fp_mutations_dict)


##### Calculate within each subclone in the tree
def calculate_fp_ratios_within_subclone(M_checkpoint, I_attached, mutation_clones_for_subclone):
    """
    Calculate fp_ratio for each mutation within each subclone, and record other mutations with FP>0
    
    Parameters:
    M_checkpoint: Imputed genotype matrix (binary, 0/1)
    I_attached: Original input genotype matrix (0, 1, NaN)
    mutation_clones_for_subclone: Dictionary where keys are subclone representative mutations, 
                                  values are lists of all mutations in that subclone
    
    Returns:
    DataFrame: Contains fp_ratio for each mutation within each subclone
    dict: List of other mutations with FP>0 for each mutation
    """
    results_for_out_subclone_muts = []
    results_for_in_subclone_muts = []
    # New: Store other mutations with FP>0 for each mutation
    fp_positive_mutations_dict_for_out_subclone_muts = {}
    fp_positive_mutations_dict_for_in_subclone_muts = {}
    
    # Iterate through each subclone
    for subclone_rep, subclone_mutations in mutation_clones_for_subclone.items():
        # Get all cells for the current subclone (cells with at least one mutation from this subclone)
        subclone_cells = set()
        for mut in subclone_mutations:
            mutant_cells = M_checkpoint[M_checkpoint[mut] == 1].index
            subclone_cells.update(mutant_cells)
        
        subclone_cells = list(subclone_cells)
        
        # First calculate fp count and fp ratio for all mutations
        for mut in subclone_mutations:
            # Get mutant cells for the current mutation
            mutant_cells = M_checkpoint[M_checkpoint[mut] == 1].index
            
            # Get other mutations in the subclone (excluding itself)
            other_mutations_out_subclone = [m for m in M_checkpoint.columns if m != mut and m not in subclone_mutations]
            
            fp_ratios = []
            # New: Store other mutations with FP>0 for the current mutation
            fp_positive_for_current_mut = []
            
            # Only calculate FP for other mutations within the same subclone
            for other_mut in other_mutations_out_subclone:
                # Calculate FP for other mutations within the current mutation's clone
                other_original = I_attached.loc[mutant_cells, other_mut]
                other_imputed = M_checkpoint.loc[mutant_cells, other_mut]
                
                # FP: original is 1 but imputed is 0
                fp_count = ((other_original == 1) & (other_imputed == 0)).sum()
                
                # Calculate FP for other mutations across all cells in the current subclone
                subclone_other_original = I_attached.loc[subclone_cells, other_mut]
                subclone_other_imputed = M_checkpoint.loc[subclone_cells, other_mut]
                subclone_fp_count = ((subclone_other_original == 1) & (subclone_other_imputed == 0)).sum()
                
                if subclone_fp_count > 0:
                    fp_ratio = fp_count / subclone_fp_count
                    fp_ratios.append(fp_ratio)
                    # New: If FP ratio > 0, record this other_mut
                    if fp_ratio > 0:
                        fp_positive_for_current_mut.append(other_mut)
                else:
                    fp_ratios.append(0)
            
            # Calculate average fp_ratio
            avg_fp_ratio = np.mean(fp_ratios) if fp_ratios else 0
            
            results_for_out_subclone_muts.append({
                'identifier': mut,
                'subclone_representative': subclone_rep,
                'subclone_size': len(subclone_mutations),
                'fp_ratio_within_subclone_for_out_subclone_muts': avg_fp_ratio,
                'count_fp_mutations_for_out_subclone_muts': len(fp_positive_for_current_mut),
                'ratio_fp_mutations_for_out_subclone_muts': len(fp_positive_for_current_mut)/len(M_checkpoint.columns)
            })
            
            # New: Store other mutations with FP>0 for the current mutation
            fp_positive_mutations_dict_for_out_subclone_muts[mut] = fp_positive_for_current_mut
        
        # If subclone has only one mutation, cannot calculate fp_ratio, set to 0
        if len(subclone_mutations) <= 1:
            for mut in subclone_mutations:
                results_for_in_subclone_muts.append({
                    'identifier': mut,
                    'fp_ratio_within_subclone_for_in_subclone_muts': 0,
                    'count_fp_mutations_for_in_subclone_muts': 0,
                    'ratio_fp_mutations_for_in_subclone_muts': 0
                })
                # New: For single-mutation subclones, the list of other mutations with FP>0 is empty
                fp_positive_mutations_dict_for_in_subclone_muts[mut] = []        
            
        else:
            
            # Iterate through each mutation in the subclone
            for mut in subclone_mutations:
                # Get mutant cells for the current mutation
                mutant_cells = M_checkpoint[M_checkpoint[mut] == 1].index
                
                # Get other mutations in the subclone (excluding itself)
                other_mutations_in_subclone = [m for m in subclone_mutations if m != mut]
                
                fp_ratios = []
                # New: Store other mutations with FP>0 for the current mutation
                fp_positive_for_current_mut = []
                
                # Only calculate FP for other mutations within the same subclone
                for other_mut in other_mutations_in_subclone:
                    # Calculate FP for other mutations within the current mutation's clone
                    other_original = I_attached.loc[mutant_cells, other_mut]
                    other_imputed = M_checkpoint.loc[mutant_cells, other_mut]
                    
                    # FP: original is 1 but imputed is 0
                    fp_count = ((other_original == 1) & (other_imputed == 0)).sum()
                    
                    # Calculate FP for other mutations across all cells in the current subclone
                    subclone_other_original = I_attached.loc[subclone_cells, other_mut]
                    subclone_other_imputed = M_checkpoint.loc[subclone_cells, other_mut]
                    subclone_fp_count = ((subclone_other_original == 1) & (subclone_other_imputed == 0)).sum()
                    
                    if subclone_fp_count > 0:
                        fp_ratio = fp_count / subclone_fp_count
                        fp_ratios.append(fp_ratio)
                        # New: If FP ratio > 0, record this other_mut
                        if fp_ratio > 0:
                            fp_positive_for_current_mut.append(other_mut)
                    else:
                        fp_ratios.append(0)
                
                # Calculate average fp_ratio
                avg_fp_ratio = np.mean(fp_ratios) if fp_ratios else 0
                
                results_for_in_subclone_muts.append({
                    'identifier': mut,
                    'fp_ratio_within_subclone_for_in_subclone_muts': avg_fp_ratio,
                    'count_fp_mutations_for_in_subclone_muts': len(fp_positive_for_current_mut),
                    'ratio_fp_mutations_for_in_subclone_muts': len(fp_positive_for_current_mut)/len(M_checkpoint.columns)
                })
                
                # New: Store other mutations with FP>0 for the current mutation
                fp_positive_mutations_dict_for_in_subclone_muts[mut] = fp_positive_for_current_mut
    
    df_results = pd.merge(
        pd.DataFrame(results_for_out_subclone_muts),
        pd.DataFrame(results_for_in_subclone_muts),
        on="identifier",
        how='outer')
    
    return df_results, fp_positive_mutations_dict_for_out_subclone_muts, fp_positive_mutations_dict_for_in_subclone_muts

# Usage
# fp_ratios_df, fp_positive_mutations_dict_for_out_subclone_muts, fp_positive_mutations_dict_for_in_subclone_muts = calculate_fp_ratios_within_subclone(M_checkpoint, I_attached, mutation_clones_for_subclone)
# print(ratios_df)


def calculate_fp_ratios_persite_within_subclone(M_checkpoint, I_attached, mutation_clones_for_subclone):
    """
    Calculate fp_ratio_persite for each mutation within each subclone
    (Number of FP divided by total number of 1s in the original data for that mutation within the subclone)
    
    Parameters:
    M_checkpoint: Imputed genotype matrix (binary, 0/1)
    I_attached: Original input genotype matrix (0, 1, NaN)
    mutation_clones_for_subclone: Dictionary where keys are subclone representative mutations,
                                  values are lists of all mutations in that subclone
    
    Returns:
    DataFrame: Contains fp_count_persite and fp_ratio_persite for each mutation within each subclone
    """
    results = []
    
    # Iterate through each subclone
    for subclone_rep, subclone_mutations in mutation_clones_for_subclone.items():
        # Get all cells for the current subclone (cells with at least one mutation from this subclone)
        subclone_cells = set()
        for mut in subclone_mutations:
            mutant_cells = M_checkpoint[M_checkpoint[mut] == 1].index
            subclone_cells.update(mutant_cells)
        subclone_cells = list(subclone_cells)
        
        # Iterate through each mutation in the subclone
        for mut in subclone_mutations:
            # Calculate the number of 1s in the original data for the current mutation within the subclone
            original_data_in_subclone = I_attached.loc[subclone_cells, mut]
            total_ones_in_subclone = (original_data_in_subclone == 1).sum()
            
            # If the mutation has no 1s in the original data within the subclone, set fp_ratio to 0
            if total_ones_in_subclone == 0:
                results.append({
                    'identifier': mut,
                    'subclone_representative': subclone_rep,
                    'fp_count_persite': 0,
                    'fp_ratio_persite': 0,
                    'total_ones_in_subclone': 0,
                    'subclone_size': len(subclone_mutations)
                })
                continue
            
            # Calculate the FP count for the current mutation itself
            # Across all subclone cells: original is 1 but imputed is 0
            original_mut_in_subclone = I_attached.loc[subclone_cells, mut]
            imputed_mut_in_subclone = M_checkpoint.loc[subclone_cells, mut]
            
            fp_count = ((original_mut_in_subclone == 1) & (imputed_mut_in_subclone == 0)).sum()
            
            # Calculate fp_ratio
            fp_ratio = fp_count / total_ones_in_subclone
            
            results.append({
                'identifier': mut,
                'subclone_representative': subclone_rep,
                'fp_count_persite': fp_count,
                'fp_ratio_persite': fp_ratio,
                'total_ones_in_subclone': total_ones_in_subclone,
                'subclone_size': len(subclone_mutations)
            })
    
    return pd.DataFrame(results)


def calculate_fp_ratio_per_mutation_with_fp_mutations_dict(M_checkpoint, I_attached):
    """
    Calculate fp_ratio_per_mutation for each mutation across all cells, and output the list of other mutations with fp>0 for each mutation
    
    Parameters:
    M_checkpoint: Imputed genotype matrix (binary, 0/1)
    I_attached: Original input genotype matrix (0, 1, NaN)
    
    Returns:
    tuple: (DataFrame, dict)
        DataFrame: Contains fp_count and fp_ratio_per_mutation for each mutation
        dict: Contains list of other mutations with fp>0 for each mutation
    """
    results = []
    fp_mutations_dict = {}  # Store list of other mutations with fp>0 for each mutation
    
    # Iterate through each mutation
    for mut in I_attached.columns:
        # Calculate the number of 1s in the original data for this mutation across all cells
        original_data = I_attached[mut]
        total_ones = (original_data == 1).sum()
        total_zeros = (original_data == 0).sum()
        
        # If the mutation has no 1s in the original data, set fp_ratio to 0
        if total_ones == 0:
            results.append({
                'identifier': mut,
                'fp_cells_count': 0,
                'fp_cells_ratio_per_mutation': 0,
                'total_ones_cells': 0,
                'total_zeros_cells': 0,
                'total_coverage_cells': 0
            })
            fp_mutations_dict[mut] = []  # Empty list
            continue
        
        # Calculate the FP count for this mutation: original is 1 but imputed is 0
        fp_count = ((original_data == 1) & (M_checkpoint[mut] == 0)).sum()
        
        # Calculate fp_ratio_per_mutation
        fp_ratio = fp_count / total_ones
        
        # Find other mutations with FP>0 for this mutation
        fp_positive_muts = []
        for other_mut in I_attached.columns:
            if other_mut == mut:
                continue
                
            # Find cells where this mutation is originally 1 but imputed as 0
            fp_cells_mask = (original_data == 1) & (M_checkpoint[mut] == 0)
            if fp_cells_mask.any():
                # In these FP cells, check the status of other_mut
                fp_cells = original_data[fp_cells_mask].index
                other_in_fp_cells = I_attached.loc[fp_cells, other_mut]
                
                # If other_mut has at least one 1 in these FP cells, record it
                if (other_in_fp_cells == 1).any():
                    fp_positive_muts.append(other_mut)
        
        results.append({
            'identifier': mut,
            'fp_cells_count': fp_count,
            'fp_cells_ratio_per_mutation': fp_ratio,
            'total_ones_cells': total_ones,
            'total_zeros_cells': total_zeros,
            'total_coverage_cells': total_ones + total_zeros
        })
        
        fp_mutations_dict[mut] = fp_positive_muts
    
    return pd.DataFrame(results), fp_mutations_dict

# # Usage
# df_fp_ratio, fp_mutations_dict = calculate_fp_ratio_per_mutation_with_fp_mutations_dict(M_checkpoint, I_attached)

def calculate_fp_ratio_per_cell(M_checkpoint, I_attached):
    """
    Calculate fp_ratio_per_cell for each cell across all mutations
    (Number of FP divided by total number of 1s in the original data for that cell across all mutations)
    
    Parameters:
    M_checkpoint: Imputed genotype matrix (binary, 0/1)
    I_attached: Original input genotype matrix (0, 1, NaN)
    
    Returns:
    DataFrame: Contains fp_count and fp_ratio_per_cell for each cell
    """
    results = []
    
    # Iterate through each cell
    for cell in I_attached.index:
        # Calculate the number of 1s in the original data for this cell across all mutations
        original_data = I_attached.loc[cell]
        total_ones = (original_data == 1).sum()
        total_zeros = (original_data == 0).sum()
        
        # If the cell has no 1s in the original data, set fp_ratio to 0
        if total_ones == 0:
            results.append({
                'cell_id': cell,
                'fp_muts_count': 0,
                'fp_muts_ratio_per_cell': 0,
                'total_ones_muts': 0,
                'total_zeros_muts': 0,
                'total_coverage_muts': 0
            })
            continue
        
        # Calculate the FP count for this cell: original is 1 but imputed is 0
        fp_count = ((original_data == 1) & (M_checkpoint.loc[cell] == 0)).sum()
        
        # Calculate fp_ratio_per_cell
        fp_ratio = fp_count / total_ones
        
        results.append({
            'cell_id': cell,
            'fp_muts_count': fp_count,
            'fp_muts_ratio_per_cell': fp_ratio,
            'total_ones_muts': total_ones,
            'total_zeros_muts': total_zeros,
            'total_coverage_muts': total_ones+total_zeros
        })
    
    return pd.DataFrame(results)

def calculate_comprehensive_fp_metrics(M_checkpoint, I_attached):
    """
    Calculate all comprehensive FP metrics
    
    Parameters:
    M_checkpoint: Imputed genotype matrix (binary, 0/1)
    I_attached: Original input genotype matrix (0, 1, NaN)
    
    Returns:
    tuple: (per_mutation_df, per_cell_df, overall_metrics)
    """
    # Calculate per-mutation FP ratio
    per_mutation_df, fp_mutations_dict = calculate_fp_ratio_per_mutation_with_fp_mutations_dict(M_checkpoint, I_attached)
    
    # Calculate per-cell FP ratio
    per_cell_df = calculate_fp_ratio_per_cell(M_checkpoint, I_attached)
    
    # Calculate overall metrics
    total_fp = per_mutation_df['fp_cells_count'].sum()
    total_original_ones = per_mutation_df['total_ones_cells'].sum()
    overall_fp_ratio = total_fp / total_original_ones if total_original_ones > 0 else 0
    
    overall_metrics = {
        'total_fp_count': total_fp,
        'total_original_ones': total_original_ones,
        'overall_fp_ratio': overall_fp_ratio,
        'mean_fp_ratio_per_mutation': per_mutation_df['fp_cells_ratio_per_mutation'].mean(),
        'median_fp_ratio_per_mutation': per_mutation_df['fp_cells_ratio_per_mutation'].median(),
        'mean_fp_ratio_per_cell': per_cell_df['fp_muts_ratio_per_cell'].mean(),
        'median_fp_ratio_per_cell': per_cell_df['fp_muts_ratio_per_cell'].median()
    }
    
    return per_mutation_df, per_cell_df, overall_metrics, fp_mutations_dict

# # Calculate all metrics
# per_mutation_df, per_cell_df, overall_metrics = calculate_comprehensive_fp_metrics(M_full, I_attached)

# Individual functions can be used separately
# per_mutation_df = calculate_fp_ratio_per_mutation(M_checkpoint, I_attached)
# per_cell_df = calculate_fp_ratio_per_cell(M_checkpoint, I_attached)




def get_all_daughter_mutations(node):
    """Best implementation: Get all descendant mutation names of a node"""
    if not node.children:
        return []
    
    daughter_mutations = []
    for child in node.children:
        # Use existing traverse method, concise and efficient
        for descendant in child.traverse():
            daughter_mutations.append(descendant.name)
    
    return daughter_mutations

def find_ordered_branch_groups_for_rehanged_mutations_with_keys_as_earlist(tree_root, mutation_list):
    """
    Optimized version: Use DFS to find all related mutation groups, return dictionary with earliest mutation as key
    """
    target_mutations = set(mutation_list)
    all_groups = []
    
    def dfs(node, current_branch):
        current_branch.append(node.name)
        
        # Check how many target mutations are on the current branch
        branch_mutations = [m for m in current_branch if m in target_mutations]
        
        # If the current branch has at least 2 target mutations
        if len(branch_mutations) >= 2:
            # Check if this group is already contained in other groups
            is_new_group = True
            for existing_group in all_groups:
                if set(branch_mutations).issubset(set(existing_group)):
                    is_new_group = False
                    break
            
            if is_new_group:
                all_groups.append(branch_mutations.copy())
        
        # Recursively search child nodes
        if hasattr(node, 'children') and node.children:
            for child in node.children:
                dfs(child, current_branch)
        
        # Backtrack
        current_branch.pop()
    
    # Start DFS from root node
    dfs(tree_root, [])
    
    # Filter: keep only maximal continuous groups
    final_groups = []
    for group in all_groups:
        is_maximal = True
        for other_group in all_groups:
            if group != other_group and set(group).issubset(set(other_group)):
                is_maximal = False
                break
        if is_maximal:
            final_groups.append(group)
    
    # Build dictionary: key is the first mutation of each group, value is the entire group
    result_dict = {}
    for group in final_groups:
        earliest_mutation = group[0]  # The first element is the earliest mutation
        result_dict[earliest_mutation] = group
    
    return result_dict

# # Usage with optimized version
# related_mutations_dict = find_ordered_branch_groups_for_rehanged_mutations_with_keys_as_earlist(
#     T_current, 
#     rehanged_mutations_by_fpratio_within_subclone_but_backbone
# )

# logger.info("Dictionary of related mutation groups with earliest mutation as key:")
# for earliest_mutation, group in related_mutations_dict.items():
#     logger.info(f"Key: {earliest_mutation}")
#     logger.info(f"Value: {group}")
#     logger.info("")


def find_ordered_branch_groups_for_rehanged_mutations_with_keys_as_earlist_relaxed(tree_root, mutation_list, logger_obj=None):
    """
    Relaxed version: Allow a single mutation to form a group
    """
    
    log = logger_obj if logger_obj is not None else logging.getLogger(__name__)
    
    log.info("=== Relaxed version: Single mutation allowed ===")
    log.info(f"Input mutation_list: {mutation_list}")
    
    target_mutations = set(mutation_list)
    all_groups = []
    
    def dfs(node, current_branch):
        current_branch.append(node.name)
        
        # Check how many target mutations are on the current branch
        branch_target_mutations = [m for m in current_branch if m in target_mutations]
        
        # Relaxed condition: record if there is at least 1 target mutation
        if len(branch_target_mutations) >= 1:
            log.info(f"Current node: {node.name}, branch target mutations: {branch_target_mutations}")
            
            # Check if this group is already contained in other groups
            is_new_group = True
            for existing_group in all_groups:
                if set(branch_target_mutations).issubset(set(existing_group)):
                    is_new_group = False
                    break
            
            if is_new_group:
                all_groups.append(branch_target_mutations.copy())
                log.info(f"  ➕ Added new group: {branch_target_mutations}")
        
        # Recursively search child nodes
        if hasattr(node, 'children') and node.children:
            for child in node.children:
                dfs(child, current_branch)
        
        # Backtrack
        current_branch.pop()
    
    # Start DFS from root node
    log.info("Starting DFS traversal...")
    dfs(tree_root, [])
    
    log.info(f"All groups found: {all_groups}")
    
    # Build dictionary: key is the first mutation of each group, value is the entire group
    result_dict = {}
    for group in all_groups:
        earliest_mutation = group[0]  # The first element is the earliest mutation
        result_dict[earliest_mutation] = group
    
    log.info(f"Final result dictionary: {result_dict}")
    return result_dict




##### Reasonably remove some mutations from TreeNode and Matrix format tree for re-attachment
def remove_mutations_from_tree_and_matrix(root: TreeNode, M_current: pd.DataFrame, rehanged_mutations: list):
    import copy
    T_removed = root.copy()
    M_removed = M_current.copy()
    
    # Iterate through each node, process rehanged_mutations
    for node in T_removed.all_nodes():
        if node.name == "ROOT":
            continue
        muts = node.name.split("|")
        remaining_muts = [m for m in muts if m not in rehanged_mutations]
        
        if len(remaining_muts) == 0:
            # All mutations in the node are removed -> delete the node
            parent = node.parent
            if parent is None:
                raise ValueError("Can not remove ROOT")
            for c in list(node.children):
                parent.add_child(c)
            parent.remove_child(node)
            # Delete corresponding matrix column
            if node.name in M_removed.columns:
                M_removed = M_removed.drop(columns=[node.name])
        else:
            # Node still has remaining mutations -> update node name
            new_node_name = "|".join(remaining_muts)
            if node.name != new_node_name:
                # Update tree node name
                old_name = node.name
                node.name = new_node_name
                # Update matrix column name
                if old_name in M_removed.columns:
                    M_removed = M_removed.rename(columns={old_name: new_node_name})
    
    # Ensure ROOT column exists
    if 'ROOT' not in M_removed.columns:
        M_removed.insert(0, 'ROOT', 1)
    
    # Ensure column order matches tree node order
    final_columns = ['ROOT'] + T_removed.all_names_no_root()
    # Sort matrix columns, fill missing columns with all zeros
    for col in final_columns:
        if col not in M_removed.columns and col != 'ROOT':
            M_removed[col] = 0
    M_removed = M_removed[final_columns]
        
    return T_removed, M_removed


def calculate_flip_counts_per_site(df1, df2):
    # Make sure both dataframes have the same shape
    if df1.shape != df2.shape:
        raise ValueError("The shapes of the two dataframes do not match.")
    # Get column names
    column_names = df1.columns
    # Define a function to calculate the flip count of a column
    def count_flips(column):
        val1 = df1[column].values
        val2 = df2[column].values
        # Calculate the flip count using vectorized operations
        return compare_elements_vectorized(val1, val2)
    # Apply the count flips function to each column
    flip_counts_results = df1.columns.to_series().apply(count_flips)
    # Convert the result to a dataframe
    flip_counts_df = pd.DataFrame(list(flip_counts_results))
    flip_counts_df.columns = ['flipping_False_Negative', 'flipping_False_Positive', 'flipping_NA_to_0', 'flipping_NA_to_1']
    flip_counts_df.index = df1.columns
    return flip_counts_df


# Assuming T_current is the current tree structure, new_mut is the current mutation
def generate_new_leaf_on_root(T_current: TreeNode, new_mut: str):
    # Deep copy the tree to avoid modifying the original
    new_tree = deepcopy(T_current)
    
    # Find the root node (assuming the root node name is "ROOT")
    root_node = new_tree.find("ROOT")
    
    # Create a new leaf node
    new_leaf = TreeNode(new_mut)
    
    # Add the new leaf node to the root node's children list
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
    
    return new_leaf_position


import networkx as nx
def cluster_external_mutations_by_intersection(I_selected, external_mutations, min_shared=1):
    """
    Cluster external_mutations into subtree groups based on intersection, and provide a reasonable addition order within each group.
    
    Parameters:
    - I_selected: DataFrame, mutation × cell matrix
    - external_mutations: list[str], mutations to process
    - min_shared: int, two mutations are considered to have intersection if shared cells >= min_shared
    
    Returns:
    - list[list[str]]: Each subtree (mutation group), order within group can be used directly with add_new_mutation_to_tree_independent()
    """
    # 1. Build intersection graph
    G = nx.Graph()
    G.add_nodes_from(external_mutations)
    
    for i, m1 in enumerate(external_mutations):
        for m2 in external_mutations[i+1:]:
            if m1 not in I_selected.columns or m2 not in I_selected.columns:
                continue
            v1 = I_selected[m1].fillna(0).astype(int)
            v2 = I_selected[m2].fillna(0).astype(int)
            inter = (v1 & v2).sum()
            if inter >= min_shared:
                G.add_edge(m1, m2)
    
    # 2. Find connected components (each component is a subtree)
    subtree_groups = []
    for comp in nx.connected_components(G):
        group = list(comp)
        
        # 3. Build a directed graph within each group (which mutation is higher up)
        # Use inclusion relationship or intersection size to determine direction
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
                    # Intersecting but not inclusive, determine direction using intersection size
                    inter1 = (v1 & v2).sum()
                    inter2 = inter1  # symmetric
                    if inter1 > 0:
                        # Put the one with fewer total 1s above (more sparse mutation appears earlier)
                        if v1.sum() < v2.sum():
                            DG.add_edge(m1, m2)
                        else:
                            DG.add_edge(m2, m1)
        
        try:
            order = list(nx.topological_sort(DG))
        except nx.NetworkXUnfeasible:
            # If there are cyclic relationships (symmetric intersection), sort by total 1 count ascending
            order = sorted(group, key=lambda m: I_selected[m].sum())
        
        subtree_groups.append(order)
    
    return subtree_groups




##### Partition barcode clones based on mutation clone from TreeNode format tree
def get_mutation_clone_and_backbone_mut_as_keys_by_first_level_with_frequency(root: TreeNode, df: pd.DataFrame, logger_obj=None):
    """
    Partition clones by the first level of ROOT, fully splitting compound mutations into individual mutations.
    For compound mutation roots, select the mutation with the highest frequency in the dataframe as the key.
    
    Parameters:
        root: Root node of the phylogenetic tree
        df: Mutation dataframe, rows are cells, columns are mutations, values are 0/1/NaN
    
    Returns:
        Dict[str, List[str]]: 
            key: Root mutation with the highest frequency in the dataframe
            value: All individual mutations contained within the clone
    """
    
    log = logger_obj if logger_obj is not None else logging.getLogger(__name__)
    
    def split_compound_mutations(mutation_name):
        """Split compound mutation name into a list of individual mutations"""
        if '|' in mutation_name:
            return mutation_name.split('|')
        else:
            return [mutation_name]
    
    def get_mutation_frequency(mutation, dataframe):
        """Calculate the frequency of a mutation in the dataframe (proportion of cells with value 1)"""
        if mutation not in dataframe.columns:
            return 0
        col_data = dataframe[mutation]             # Keep NaN values
        if len(col_data) == 0:
            return 0
        return (col_data == 1).sum() / len(col_data)
    
    clone_dict = {}
    for child in root.children:
        # Get all node names in the subtree and split compound mutations
        all_mutations = []
        for node in child.traverse():
            all_mutations.extend(split_compound_mutations(node.name))
        
        # Process compound mutations in the root node
        root_mutations = split_compound_mutations(child.name)
        
        # If the root has only one mutation, use it directly
        if len(root_mutations) == 1:
            clone_key = root_mutations[0]
        else:
            # For compound mutations, select the one with the highest frequency in the dataframe
            mutation_frequencies = []
            for mutation in root_mutations:
                freq = get_mutation_frequency(mutation, df)
                mutation_frequencies.append((mutation, freq))
            
            # Sort by frequency and select the highest one
            mutation_frequencies.sort(key=lambda x: x[1], reverse=True)
            clone_key = mutation_frequencies[0][0]
            
            log.info(f"Compound mutation {child.name} split into: {root_mutations}")
            log.info(f"Frequency statistics: {mutation_frequencies}")
            log.info(f"Selected as key: {clone_key}\n")
        
        clone_dict[clone_key] = all_mutations
    
    return clone_dict


def find_best_backbone_node(mutation_list_under_backbone_nodes, M_current, intersection_counts_under_backbone_node):
    """
    Simplified version: only return the best backbone node
    """
    if not intersection_counts_under_backbone_node:
        return None
    
    max_count = max(intersection_counts_under_backbone_node.values())
    max_nodes = [node for node, count in intersection_counts_under_backbone_node.items() 
                if count == max_count]
    
    if len(max_nodes) == 1:
        return max_nodes[0]
    
    # Tie handling: find the actual existing column
    best_node = None
    best_normalized = -1
    
    for node in max_nodes:
        # ===== Key fix: find the actual column name in M_current =====
        actual_col = find_column_in_merged_columns(M_current, node)
        if actual_col is None:
            actual_col = node
        
        if actual_col not in M_current.columns:
            continue
        
        backbone_cells = M_current[M_current[actual_col] == 1].index
        backbone_cell_count = len(backbone_cells)
        
        if backbone_cell_count > 0:
            normalized = intersection_counts_under_backbone_node[node] / backbone_cell_count
            if normalized > best_normalized:
                best_normalized = normalized
                best_node = node
        else:
            if best_normalized < 0:
                best_node = node
                best_normalized = 0
    
    return best_node

def calculate_intersection_counts_under_backbone_nodes(
    mutation_list_under_backbone_nodes,
    M_current,
    I_attached,
    new_mut,
    i_attached_positive=None,
    i_attached_positive_values=None,
    i_attached_positive_col_to_idx=None,
    backbone_mutation_indices=None,
    backbone_positive_masks=None,
):
    """
    Calculate co-occurrence counts between new_mut and mutation lists under each backbone node
    
    Parameters:
    mutation_list_under_backbone_nodes: dict, mapping from backbone node to mutation list
    M_current: DataFrame, cell states on backbone nodes (0/1)
    I_attached: DataFrame, cell states on other mutations (0/1/NaN)
    new_mut: str, the new mutation to check
    
    Returns:
    dict: co-occurrence count for each backbone node
    """
    
    # Check if new_mut exists in I_attached
    if new_mut not in I_attached.columns:
        raise ValueError(f"Mutation {new_mut} not found in I_attached dataframe")
    
    intersection_counts = {}
    if i_attached_positive is None:
        i_attached_positive = I_attached.reindex(index=M_current.index).eq(1).fillna(False)
    if i_attached_positive_values is None:
        i_attached_positive_values = i_attached_positive.to_numpy(dtype=bool, copy=False)
    if i_attached_positive_col_to_idx is None:
        i_attached_positive_col_to_idx = {
            col: idx for idx, col in enumerate(i_attached_positive.columns)
        }
    
    new_mut_idx = i_attached_positive_col_to_idx.get(new_mut)
    if new_mut_idx is None:
        raise ValueError(f"Mutation {new_mut} not found in I_attached dataframe")
    
    new_mut_positive_mask = i_attached_positive_values[:, new_mut_idx]
    
    for backbone_node, mutation_list in mutation_list_under_backbone_nodes.items():
        if backbone_node not in M_current.columns:
            intersection_counts[backbone_node] = 0
            continue
        
        if backbone_positive_masks is not None and backbone_node in backbone_positive_masks:
            backbone_positive_mask = backbone_positive_masks[backbone_node]
        else:
            backbone_positive_mask = _coerce_bool_array(M_current[backbone_node], M_current.index)
        
        selected_rows = backbone_positive_mask & new_mut_positive_mask
        if not bool(selected_rows.any()):
            intersection_counts[backbone_node] = 0
            continue
        
        if backbone_mutation_indices is not None and backbone_node in backbone_mutation_indices:
            valid_indices = backbone_mutation_indices[backbone_node]
        else:
            valid_indices = [
                i_attached_positive_col_to_idx[mutation]
                for mutation in mutation_list
                if mutation in i_attached_positive_col_to_idx
            ]
        
        if not valid_indices:
            intersection_counts[backbone_node] = 0
            continue
        
        intersection_counts[backbone_node] = int(i_attached_positive_values[selected_rows][:, valid_indices].sum())
    
    return intersection_counts

# # Usage
# new_mut = 'chr8_91087030_T_G'
# result = calculate_intersection_counts_under_backbone_nodes(mutation_list_under_backbone_nodes, M_current, I_attached, new_mut)
# logger.info(result)

def find_best_backbone_for_new_mutation(
    mutation_list_under_backbone_nodes,
    M_current,
    I_attached,
    new_mut,
    i_attached_positive=None,
    i_attached_positive_values=None,
    i_attached_positive_col_to_idx=None,
    backbone_mutation_indices=None,
    backbone_positive_masks=None,
):
    """
    Full function: calculate co-occurrence counts and find the best backbone node
    """
    # Step 1: Calculate co-occurrence counts
    intersection_counts = calculate_intersection_counts_under_backbone_nodes(
        mutation_list_under_backbone_nodes,
        M_current,
        I_attached,
        new_mut,
        i_attached_positive=i_attached_positive,
        i_attached_positive_values=i_attached_positive_values,
        i_attached_positive_col_to_idx=i_attached_positive_col_to_idx,
        backbone_mutation_indices=backbone_mutation_indices,
        backbone_positive_masks=backbone_positive_masks,
    )
    
    # Step 2: Find the best backbone node
    best_backbone = find_best_backbone_node(
        mutation_list_under_backbone_nodes, M_current, intersection_counts
    )
    
    return best_backbone, intersection_counts

# # Usage with full function
# best_backbone, intersection_counts = find_best_backbone_for_new_mutation(
#     mutation_list_under_backbone_nodes, M_current, I_attached, new_mut
# )
# logger.info(f"Best backbone node: {best_backbone}")
# logger.info(f"Co-occurrence counts for all backbones: {intersection_counts}")


def assign_clone_labels(M_full: pd.DataFrame, mutation_clones: dict) -> pd.DataFrame:
    """
    Assign clone labels to each barcode based on the mutations it contains.
    Returns a new DataFrame containing barcodes and clone labels.
    """
    # Initialize an empty dataframe for storing labels and colors
    result_df = pd.DataFrame(columns=["label", "color", "backbone_mutation"])
    
    # Iterate through each mutation group (i.e., each clone)
    for clone_idx, (mutation_group, mutations) in enumerate(mutation_clones.items(), start=1):
        # Assign labels and colors to barcodes containing this mutation group
        for barcode in M_full.index:
            # If this barcode contains any mutation from the current group
            if M_full.loc[barcode, mutations].sum() > 0:
                # Use pd.concat instead of append
                new_row = pd.DataFrame({"label": [barcode], "color": [f'C{clone_idx}'], "backbone_mutation": mutation_group})
                result_df = pd.concat([result_df, new_row], ignore_index=True)
    
    return result_df


def get_first_level_backbone_nodes(root: TreeNode) -> List[str]:
    """
    Get all first-level backbone nodes
    """
    first_level_nodes = []
    for child in root.children:
        first_level_nodes.append(child.name)
    return first_level_nodes


def _build_attach_tree_state(root: TreeNode) -> Dict[str, Any]:
    """
    Build a one-time tree state cache for the Step6 attachment phase to reduce repeated traversal within a single mutation.
    """
    traversed_nodes = list(root.traverse())
    node_lookup = {}
    tree_parent_dict = {}
    node_to_mutations = {}
    
    for node in traversed_nodes:
        node_lookup[node.name] = node
        node_to_mutations[node.name] = tuple(node.name.split("|")) if node.name != "ROOT" else ("ROOT",)
        for child in node.children:
            tree_parent_dict[child.name] = node.name
    
    node_list_under_backbone_nodes = {}
    mutation_list_under_backbone_nodes = {}
    for child in root.children:
        node_names = []
        mutation_names = []
        for node in child.traverse():
            node_names.append(node.name)
            mutation_names.extend(node_to_mutations[node.name])
        node_list_under_backbone_nodes[child.name] = node_names
        mutation_list_under_backbone_nodes[child.name] = mutation_names
    
    return {
        "tree_nodes": traversed_nodes,
        "node_lookup": node_lookup,
        "tree_parent_dict": tree_parent_dict,
        "node_to_mutations": node_to_mutations,
        "node_list_under_backbone_nodes": node_list_under_backbone_nodes,
        "mutation_list_under_backbone_nodes": mutation_list_under_backbone_nodes,
    }




# -------------------------
# Calculate intersection/FN_flip per mutation for parent muts, then identify mutations to rehang
# -------------------------

import pandas as pd
import numpy as np

def find_mutation_column(mutation, columns):
    """Find the column in M_current that contains the specified mutation"""
    for col in columns:
        if '|' in col:
            if mutation in col.split('|'):
                return col
        else:
            if mutation == col:
                return col
    return None

def get_all_parents_for_mutation(T_current, mutation):
    """Get all parent mutations for a given mutation (all ancestors from direct parent up to ROOT)"""
    parent_dict = build_lineage_parent_dict_from_tree(T_current, mutation)
    
    all_parents = []
    current = mutation
    
    # Traverse up the parent chain until reaching ROOT
    while current in parent_dict:
        parent = parent_dict[current]
        all_parents.append(parent)
        current = parent
    
    return all_parents

def calculate_intersection_and_inter_vs_fn_flipping_ratio_per_mutation(T_current, M_current, I_attached, logger_obj=None):
    """
    Correct version: each mutation returns only one row of statistics.
    
    Parameters
    ----------
    T_current : TreeNode
        Current phylogenetic tree
    M_current : pd.DataFrame
        Current mutation matrix
    I_attached : pd.DataFrame
        Attached mutation information
    logger_obj : logging.Logger, optional
        Logger instance for logging messages. If None, uses global logger.
    
    Returns
    -------
    pd.DataFrame
        Statistics for each mutation
    """
    # Use provided logger or fall back to global logger
    log = logger_obj if logger_obj is not None else logging.getLogger(__name__)
    
    results = []
    mutations = I_attached.columns.tolist()
    
    # Align the indices of the two DataFrames
    common_index = I_attached.index.intersection(M_current.index)
    I_aligned = I_attached.loc[common_index]
    M_aligned = M_current.loc[common_index]
    
    log.debug(f"Calculating intersection statistics for {len(mutations)} mutations")
    
    for mutation in mutations:
        # Find the corresponding column of the current mutation in M_current
        m_col = find_mutation_column(mutation, M_aligned.columns)
        
        # Get all parent mutations
        all_parents = get_all_parents_for_mutation(T_current, m_col)
        
        # 2. If the node is found, get all its direct children
        all_children = []
        if T_current.find(m_col):
            children = T_current.find(m_col).children
            for child in children:
                all_children.append(child.name)
        else:
            log.debug(f"Node not found: {m_col}")
        
        if m_col is None:
            log.debug(f"Warning: mutation {mutation} not found in M_current")
            # If not found, fill with NaN
            result = {
                'mutation': mutation,
                'total_retained_cells': np.nan,
                'total_parent_intersection': np.nan,
                'total_parent_FNflipping': np.nan,
                'parent_retention_ratio': np.nan,
                'parent_count': np.nan,
                'mutation_clone_size': np.nan,
                'intersection_cell_count_on_mutation': np.nan,
                'intersection_cell_ratio_on_mutation': np.nan,
                'intersection_cells_on_mutation_parents': None,
                'intersection_cells_on_mutation_children': None
            }
            results.append(result)
            continue
        
        # Find cells that remain as 1
        retained_mask = (I_aligned[mutation] == 1) & (M_aligned[m_col] == 1)
        retained_cells = I_aligned.index[retained_mask]
        total_retained = len(retained_cells)
        
        if len(retained_cells) == 0 or not all_parents:
            # No retained cells or no parent mutations
            result = {
                'mutation': mutation,
                'total_retained_cells': len(retained_cells),
                'total_parent_intersection': 0,
                'total_parent_FNflipping': 0,
                'parent_retention_ratio': 1,
                'parent_count': len(all_parents),
                'mutation_clone_size': total_retained,
                'intersection_cell_count_on_mutation': total_retained,
                'intersection_cell_ratio_on_mutation': 1,
                'intersection_cells_on_mutation_parents': None,
                'intersection_cells_on_mutation_children': None
            }
            results.append(result)
            continue
        
        # Calculate overall statistics for all parent mutations
        total_parent_intersection = 0
        total_parent_FNflipping = 0
        intersection_cells_on_mutation_parents = []
        intersection_cells_on_mutation_children = []
        
        for parent in all_parents:
            parent_m_col = find_mutation_column(parent, M_aligned.columns)
            if parent_m_col is None:
                continue
            
            for cell in retained_cells:
                if M_aligned.loc[cell, parent_m_col] == 1:
                    if parent in I_aligned.columns and pd.notna(I_aligned.loc[cell, parent]):
                        if I_aligned.loc[cell, parent] == 1:
                            total_parent_intersection += 1
                            intersection_cells_on_mutation_parents.append(cell)
                        else:
                            total_parent_FNflipping += 1
        
        for child in all_children:
            child_m_col = find_mutation_column(child, M_aligned.columns)
            if child_m_col is None:
                continue
            
            for cell in retained_cells:
                if M_aligned.loc[cell, child_m_col] == 1:
                    if child in I_aligned.columns and pd.notna(I_aligned.loc[cell, child]):
                        if I_aligned.loc[cell, child] == 1:
                            intersection_cells_on_mutation_children.append(cell)
        
        unique_intersection_cells_on_mutation_parents = list(set(intersection_cells_on_mutation_parents))
        unique_intersection_cells_on_mutation_children = list(set(intersection_cells_on_mutation_children))
        intersection_cell_ratio_on_mutation = len(unique_intersection_cells_on_mutation_parents) / total_retained if total_retained > 0 else 0
        
        # Calculate overall ratio: intersection / (intersection + FN)
        total_parent_events = total_parent_intersection + total_parent_FNflipping
        parent_ratio = total_parent_intersection / total_parent_events if total_parent_events > 0 else 0.0
        
        result = {
            'mutation': mutation,
            'total_retained_cells': len(retained_cells),
            'total_parent_intersection': total_parent_intersection,
            'total_parent_FNflipping': total_parent_FNflipping,
            'parent_retention_ratio': parent_ratio,
            'parent_count': len(all_parents),
            'mutation_clone_size': total_retained,
            'intersection_cell_count_on_mutation': len(unique_intersection_cells_on_mutation_parents),
            'intersection_cell_ratio_on_mutation': intersection_cell_ratio_on_mutation,
            'intersection_cells_on_mutation_parents': intersection_cells_on_mutation_parents,
            'intersection_cells_on_mutation_children': intersection_cells_on_mutation_children
        }
        results.append(result)
    
    if results:
        log.debug(f"Calculated statistics for {len(results)} mutations")
        return pd.DataFrame(results)
    else:
        log.warning("No results calculated")
        return pd.DataFrame()


def process_matrices_by_removed_some_mutations_from_tree(M_current, I_attached):
    """
    Process two matrices:
    1. Extract the corresponding columns from I_attached that are in M_current (excluding ROOT)
    2. In the extracted I_attached data, set cells to 0 in M_current for rows that are all 0 or NaN
    
    Parameters:
    M_current: DataFrame, current matrix
    I_attached: DataFrame, attached matrix
    
    Returns:
    tuple: (I_attached_removed_outgroup, M_current_modified)
    """
    # Step 1: Extract columns from I_attached that correspond to M_current (excluding ROOT)
    m_columns = [col for col in M_current.columns if col != "ROOT"]
    split_columns = [col.split("|")[0] if "|" in col else col for col in m_columns]
    
    # Extract corresponding columns from I_attached
    I_attached_removed_outgroup = I_attached[split_columns]
    
    # Step 2: Find rows in I_attached_removed_outgroup that are all 0 or NaN
    rows_to_zero = I_attached_removed_outgroup[
        (I_attached_removed_outgroup == 0) | (I_attached_removed_outgroup.isna())
    ].all(axis=1)
    
    # Get indices of rows to modify
    zero_rows_indices = rows_to_zero[rows_to_zero].index
    
    # Set corresponding rows in M_current (except ROOT column) to 0
    M_current_modified = M_current.copy()
    M_current_modified.loc[zero_rows_indices, M_current_modified.columns != "ROOT"] = 0
    
    return I_attached_removed_outgroup, M_current_modified


def process_conflicting_cells_stay_outgroup(M_current_split_and_noROOT, I_attached_only_outgroup, logger_obj=None):
    """
    Process conflicting cells in two dataframes, return two processed dataframes
    
    Parameters:
    M_current_split_and_noROOT: First dataframe
    I_attached_only_outgroup: Second dataframe
    
    Returns:
    M_processed: Processed first dataframe
    I_processed: Processed second dataframe
    """
    
    log = logger_obj if logger_obj is not None else logging.getLogger(__name__)
    
    # Copy original dataframes to avoid modifying the originals
    M_processed = M_current_split_and_noROOT.copy()
    I_processed = I_attached_only_outgroup.copy()
    
    # Ensure both dataframes have the same row indices
    common_cells = M_processed.index.intersection(I_processed.index)
    M_processed = M_processed.loc[common_cells]
    I_processed = I_processed.loc[common_cells]
    
    # Find cells that have 1s in both dataframes
    conflicting_cells = []
    
    for cell in common_cells:
        # Get the row for this cell in both dataframes
        M_row = M_processed.loc[cell]
        I_row = I_processed.loc[cell]
        
        # Check if both dataframes have at least one 1
        # Note: I_attached_only_outgroup may have NaN values, handle accordingly
        M_has_ones = (M_row == 1).any()
        I_has_ones = (I_row.fillna(0) == 1).any()  # Treat NaN as 0
        
        if M_has_ones and I_has_ones:
            conflicting_cells.append(cell)
    
    log.info(f"Found {len(conflicting_cells)} cells with 1s in both dataframes")
    
    # Process each conflicting cell
    for cell in conflicting_cells:
        M_row = M_processed.loc[cell]
        I_row = I_processed.loc[cell]
        
        # Count number of 1s (for I dataframe, treat NaN as 0)
        M_ones_count = (M_row == 1).sum()
        I_ones_count = (I_row.fillna(0) == 1).sum()
        
        # Process according to rules
        if M_ones_count > I_ones_count:
            # Case 1: M has more 1s - set the entire row in I to 0
            I_processed.loc[cell] = 0
            log.debug(f"  -> M has more 1s, keeping M unchanged, setting I row to 0")
            
        elif I_ones_count > M_ones_count:
            # Case 2: I has more 1s - set the entire row in M to 0
            M_processed.loc[cell] = 0
            log.debug(f"  -> I has more 1s, setting M row to 0, keeping I unchanged")
            
        else:
            # Case 3: Equal number of 1s - set the entire row in M to 0, keep I unchanged
            M_processed.loc[cell] = 0
            log.debug(f"  -> Equal number of 1s, setting M row to 0, keeping I unchanged")
    
    return M_processed, I_processed

# # Usage
# M_current_split_and_noROOT_processed, I_attached_only_outgroup_processed = process_conflicting_cells(
#     M_current_split_and_noROOT, 
#     I_attached_only_outgroup
# )

# logger.info("Processed M dataframe:")
# logger.info(M_current_split_and_noROOT_processed)
# logger.info("\nProcessed I dataframe:")
# logger.info(I_attached_only_outgroup_processed)


def process_conflicting_cells_stay_maintree(M_current_split_and_noROOT, I_attached_only_outgroup, logger_obj=None):
    """
    Process conflicting cells in two dataframes
    
    Parameters:
    M_current_split_and_noROOT: First dataframe
    I_attached_only_outgroup: Second dataframe
    
    Returns:
    M_current_split_and_noROOT_processed_unconflict: Processed first dataframe
    I_attached_only_outgroup_processed: Processed second dataframe
    """
    
    log = logger_obj if logger_obj is not None else logging.getLogger(__name__)
    
    # Copy original dataframes to avoid modifying the originals
    M_processed = M_current_split_and_noROOT.copy()
    I_processed = I_attached_only_outgroup.copy()
    
    # 1. Find rows that have at least one 1 in both dataframes
    # For M dataframe, find rows with at least one 1
    M_has_ones = (M_processed == 1).any(axis=1)
    # For I dataframe, find rows with at least one 1 (handle NaN values)
    I_has_ones = (I_processed == 1).fillna(False)
    I_has_ones = I_has_ones.any(axis=1)
    
    # Rows with 1s in both dataframes
    conflicting_cells = M_has_ones & I_has_ones
    conflicting_indices = conflicting_cells[conflicting_cells].index
    
    log.info(f"Found {len(conflicting_indices)} conflicting cells")
    
    # 2. Process each conflicting cell
    for cell in conflicting_indices:
        # Count number of 1s in each dataframe
        M_ones_count = (M_processed.loc[cell] == 1).sum()
        I_ones_count = (I_processed.loc[cell] == 1).sum()
        
        log.debug(f"Cell {cell}: M has {M_ones_count} 1s, I has {I_ones_count} 1s")
        
        # Process according to rules
        if M_ones_count > I_ones_count:
            # M has more 1s: set all 1s in I row to 0
            I_row = I_processed.loc[cell]
            I_row[I_row == 1] = 0
            I_processed.loc[cell] = I_row
            log.debug(f"  -> M has more 1s: keeping M unchanged, setting I row 1s to 0")
            
        elif I_ones_count > M_ones_count:
            # I has more 1s: set all 1s in M row to 0
            M_row = M_processed.loc[cell]
            M_row[M_row == 1] = 0
            M_processed.loc[cell] = M_row
            log.debug(f"  -> I has more 1s: setting M row 1s to 0, keeping I unchanged")
            
        else:
            # Equal number of 1s: keep M unchanged, set all 1s in I row to 0
            I_row = I_processed.loc[cell]
            I_row[I_row == 1] = 0
            I_processed.loc[cell] = I_row
            log.debug(f"  -> Equal number of 1s: keeping M unchanged, setting I row 1s to 0")
    
    return M_processed, I_processed

# # Usage
# M_current_split_and_noROOT_processed, I_attached_only_outgroup_processed = process_conflicting_cells(
#     M_current_split_and_noROOT, I_attached_only_outgroup
# )


def process_conflicting_cells_allto_outgroup(M_current_split_and_noROOT, I_attached_only_outgroup, logger_obj=None):
    """
    Efficient version: set all conflicting cells in M to 0.
    
    Parameters
    ----------
    M_current_split_and_noROOT : pd.DataFrame
        First dataframe (M matrix)
    I_attached_only_outgroup : pd.DataFrame
        Second dataframe (I matrix)
    logger_obj : logging.Logger, optional
        Logger instance for logging messages. If None, uses global logger.
    
    Returns
    -------
    tuple : (M_processed, I_processed)
        Processed dataframes
    """
    # Use provided logger or fall back to global logger
    log = logger_obj if logger_obj is not None else logging.getLogger(__name__)
    
    M_processed = M_current_split_and_noROOT.copy()
    I_processed = I_attached_only_outgroup.copy()
    
    # Find cells with 1s in both dataframes
    M_has_ones = (M_processed == 1).any(axis=1)
    I_has_ones = (I_processed == 1).fillna(False).any(axis=1)
    conflicting_indices = (M_has_ones & I_has_ones)
    
    log.info(f"Found {conflicting_indices.sum()} cells with 1s in both dataframes")
    
    # Set conflicting rows in M to 0
    M_processed.loc[conflicting_indices] = 0
    
    return M_processed, I_processed

# # Usage
# M_current_split_and_noROOT_processed, I_attached_only_outgroup_processed = process_conflicting_cells_allto_outgroup(
#     M_current_split_and_noROOT, I_attached_only_outgroup
# )


# View detailed information of conflicting cells
def show_conflict_details(I_attached_removed_outgroup, I_attached_only_outgroup, n_examples=3, logger_obj=None):
    """
    Show detailed information about conflicting cells.
    
    Parameters
    ----------
    I_attached_removed_outgroup : pd.DataFrame
        First dataframe
    I_attached_only_outgroup : pd.DataFrame
        Second dataframe
    n_examples : int, default=3
        Number of examples to show
    logger_obj : logging.Logger, optional
        Logger instance for logging messages. If None, uses global logger.
    """
    # Use provided logger or fall back to global logger
    log = logger_obj if logger_obj is not None else logging.getLogger(__name__)
    
    common_cells = I_attached_removed_outgroup.index.intersection(I_attached_only_outgroup.index)
    
    # Convert data to integer type
    removed_int = I_attached_removed_outgroup.loc[common_cells].fillna(0).astype(int)
    only_int = I_attached_only_outgroup.loc[common_cells].fillna(0).astype(int)
    
    # Find cells with 1s in both dataframes
    has_1_in_removed = (removed_int.sum(axis=1) > 0)
    has_1_in_only = (only_int.sum(axis=1) > 0)
    conflict_cells = has_1_in_removed & has_1_in_only
    
    conflict_cells_index = conflict_cells[conflict_cells].index
    
    log.info(f"\n=== Conflicting Cells Details ===")
    log.info(f"Total conflicting cells found: {len(conflict_cells_index)}")
    
    if len(conflict_cells_index) > 0:
        # Count number of 1s
        count_removed = removed_int.loc[conflict_cells_index].sum(axis=1)
        count_only = only_int.loc[conflict_cells_index].sum(axis=1)
        
        # Show first few examples
        for i, cell in enumerate(conflict_cells_index[:n_examples]):
            log.info(f"\nConflict cell {i+1}: {cell}")
            log.info(f"  Number of 1s in first dataframe: {count_removed[cell]}")
            log.info(f"  Number of 1s in second dataframe: {count_only[cell]}")
            log.info(f"  1 positions in first dataframe: {removed_int.loc[cell][removed_int.loc[cell] == 1].index.tolist()}")
            log.info(f"  1 positions in second dataframe: {only_int.loc[cell][only_int.loc[cell] == 1].index.tolist()}")

# Usage
# show_conflict_details(M_current_split_and_noROOT, I_attached_only_outgroup, n_examples=54)


def create_column_mapping(source_columns, target_columns):
    """
    Create a mapping from source column names to target column names
    """
    column_mapping = {}
    
    # Create normalized versions of target column names (sorted by mutation)
    target_normalized = {}
    for col in target_columns:
        if '|' in col:
            mutations = sorted(col.split('|'))
            normalized = '|'.join(mutations)
            target_normalized[normalized] = col
        else:
            target_normalized[col] = col
    
    # Map source column names to target column names
    for col in source_columns:
        if '|' in col:
            # For merged columns, normalize and look up corresponding target column name
            mutations = sorted(col.split('|'))
            normalized = '|'.join(mutations)
            if normalized in target_normalized:
                column_mapping[col] = target_normalized[normalized]
            else:
                # If no match found, keep the original column name
                column_mapping[col] = col
        else:
            # For non-merged columns, look up directly
            if col in target_normalized:
                column_mapping[col] = target_normalized[col]
            else:
                column_mapping[col] = col
    
    return column_mapping

# # Usage
# source_cols = merge_duplicate_columns(M_current_split_and_noROOT_processed).columns
# target_cols = M_current_noROOT.columns

# column_mapping = create_column_mapping(source_cols, target_cols)

# # Rename columns
# M_current_split_renamed = merge_duplicate_columns(M_current_split_and_noROOT_processed).rename(columns=column_mapping)

# # Now the column names of both dataframes should match
# logger.info("M_current_noROOT columns:", M_current_noROOT.columns.tolist()[:5])
# logger.info("Renamed columns:", M_current_split_renamed.columns.tolist()[:5])



# -------------------------
# Calculate intersection/FN_flip per mutation for parent muts, then identify cells that can be deleted when no suitable mutation is found
# -------------------------

def calculate_intersection_and_flipping_to_1_count_per_cell(M_for_fp_ratio_and_fn_ratio_duplicatecells, I_attached):
    """
    Calculate two values for each cell:
    1. Number of mutations that are 1 in both dataframes (intersection)
    2. Number of mutations that flip from 0 or NA to 1 from I_attached to M_current (flipping to 1)
    
    Parameters:
    M_for_fp_ratio_and_fn_ratio_duplicatecells: Current subclone mutation matrix (1/0)
    I_attached: Reference mutation matrix (may contain 0, 1, NaN)
    
    Returns:
    DataFrame: Contains both count values for each cell
    """
    
    # Ensure the indices and columns of both DataFrames are aligned
    common_cells = M_for_fp_ratio_and_fn_ratio_duplicatecells.index.intersection(I_attached.index)
    common_mutations = M_for_fp_ratio_and_fn_ratio_duplicatecells.columns.intersection(I_attached.columns)
    
    # Align data
    M_aligned = M_for_fp_ratio_and_fn_ratio_duplicatecells.loc[common_cells, common_mutations]
    I_aligned = I_attached.loc[common_cells, common_mutations]
    
    results = []
    
    for cell in common_cells:
        # Get the row data for the current cell in both dataframes
        M_row = M_aligned.loc[cell]
        I_row = I_aligned.loc[cell]
        
        # 1. Calculate intersection: number of mutations that are 1 in both dataframes
        intersection_count = ((M_row == 1) & (I_row == 1)).sum()
        
        # 2. Calculate flipping to 1: number of mutations that flip from I_attached (0 or NA) to M_current (1)
        # Condition: I_attached is 0 or NaN, and M_current is 1
        flipping_count_fn = ((I_row == 0) & (M_row == 1)).sum()
        flipping_count_NAto1 = ((I_row.isna()) & (M_row == 1)).sum()
        flipping_count = flipping_count_fn + flipping_count_NAto1
        
        results.append({
            'cell': cell,
            'intersection_count': intersection_count,
            'flipping_count_fn': flipping_count_fn,
            'flipping_count_NAto1': flipping_count_NAto1,
            'flipping_to_1_count': flipping_count
        })
    
    # Create result DataFrame
    result_df = pd.DataFrame(results).set_index('cell')
    
    return result_df

# # Usage
# result = calculate_intersection_and_flipping_to_1_count_per_cell(M_for_fp_ratio_and_fn_ratio_duplicatecells, I_attached)
# logger.info(result)




# -------------------------
# 4.5 Main Integration Function
# -------------------------

##### Dynamic progremming
def run_dp_pass_tree(
    data,
    df_features_new,
    M_scaffold,
    outputpath_full,
    scaffold_mutations,
    p_thresh=0.5,
    pass_tree_cutoff=0.9,
    unpass_tree_cutoff=0.1,
    is_log_value_for_likelihoods=True,
    logger_obj=None
):
    """
    Run DP-based posterior computation and mutation filtering through scaffold tree.
    
    Parameters
    ----------
    data : dict
        Dictionary containing P, ll_mut, ll_unmut, df_reads, features, etc.
    M_scaffold : pd.DataFrame
        Scaffold matrix (binary 0/1).
    outputpath_full : str
        Output path.
    scaffold_mutations : list
        List of initial scaffold mutations.
    p_thresh : float
        Threshold for binarizing posterior.
    pass_tree_cutoff : float
        Posterior threshold for "pass tree" mutations.
    unpass_tree_cutoff : float
        Posterior threshold for "fail tree" mutations.
    is_log_value_for_likelihoods : bool
        Whether to treat likelihoods as log values for normalization.
    logger_obj : logging.Logger, optional
        Logger instance for logging messages. If None, uses global logger.
    
    Returns
    -------
    df_combined : pd.DataFrame
        Integrated feature table with phylogeny_label and flipping_counts.
    attached_mutations : list
        Mutations that can be attached to the tree (successful_pass + cell_specific).
    """
    # Use provided logger or fall back to global logger
    log = logger_obj if logger_obj is not None else logging.getLogger(__name__)
    
    # Step 0: Align data indices
    orderd_index = M_scaffold.index.tolist()
    in_posterior = data['P'].reindex(index=orderd_index)
    in_llmut = data['ll_mut'].reindex(index=orderd_index)
    in_llunmut = data['ll_unmut'].reindex(index=orderd_index)
    in_reads = data['df_reads'].reindex(index=[['bulk'] + orderd_index])
    
    # Step 1: Normalize likelihoods
    log.info("Normalizing likelihoods...")
    combined_df = pd.concat([in_llmut, in_llunmut], axis=1)
    normalized_result = apply_normalization(combined_df, is_log=is_log_value_for_likelihoods)
    normalized_llmut = normalized_result.filter(like='norm_llmut')
    normalized_llmut.columns = in_llmut.columns
    normalized_llunmut = normalized_result.filter(like='norm_llunmut')
    normalized_llunmut.columns = in_llunmut.columns
    
    # Step 2: Calculate posterior
    log.info("Calculating somatic posterior probabilities...")
    df_newSP_out, withoutTree_posterior = all_newSomaticPosterior(
        normalized_llmut, normalized_llunmut, M_scaffold.values.astype(int)
    )
    
    # Step 3: Integrate features
    log.info("Integrating features...")
    out_features = df_features_new.T.drop(['somatic_posterior_persite'], axis=1)
    out_features['somatic_posterior_per_site_old'] = df_features_new.T['somatic_posterior_persite']
    out_features['withoutTree_posterior'] = withoutTree_posterior
    out_features = pd.concat([out_features, df_newSP_out], axis=1)
    
    # Step 4: Add phylogeny label
    for col in ['somatic_posterior_per_site', 'somatic_posterior_per_site_onecell']:
        out_features[col] = out_features[col].astype(float)
    out_features['mutant_cellnum'] = out_features['mutant_cellnum'].astype(int)
    out_features['phylogeny_label'] = out_features.apply(
        determine_phylogeny_label_by_one_likelihood,
        axis=1,
        args=(pass_tree_cutoff, unpass_tree_cutoff,),
    )
    
    # Step 5: Calculate flipping counts
    log.info("Calculating flipping counts...")
    df_binary_withNA3 = posterior2ter_NAto3_bothPosteriorMutallele(in_posterior, in_reads, p_thresh)
    
    flipping_counts_list = []
    for i in range(len(out_features['cluster'])):
        init_bin = df_binary_withNA3.iloc[:, i].tolist()
        nonnan_idx = np.array(out_features['nonnan_indices'].iloc[i].tolist())
        init_bin_nonnan = np.array([init_bin[j] for j in nonnan_idx])
        max_cluster = out_features['cluster'].iloc[i]
        flipping_counts_per_site = compare_elements_vectorized(init_bin_nonnan, max_cluster)
        flipping_counts_list.append(flipping_counts_per_site)
    
    df_flip_counts = pd.DataFrame(flipping_counts_list)
    df_flip_counts.columns = [f'tree_{flip_type}' for flip_type in df_flip_counts.columns.tolist()]
    df_flip_counts.index = out_features.index
    assert df_flip_counts.index.equals(out_features.index), "Indexes of df_flip_counts and out_features do not match."
    
    # Step 6: Merge
    log.info("Merging results...")
    df_combined = pd.concat([out_features, df_flip_counts], axis=1)
    df_combined['cluster'] = df_combined['cluster'].apply(lambda x: ','.join(map(str, x)))
    df_combined['nonnan_indices'] = df_combined['nonnan_indices'].apply(lambda x: ','.join(map(str, x)))
    
    output_file = os.path.join(outputpath_full, "out_features.somatic_posterior_basedTree.txt")
    df_combined.to_csv(output_file, sep="\t")
    log.info(f"Results saved to: {output_file}")
    
    # Step 7: Mutation classification
    passtree_mutations = list(df_combined[df_combined['phylogeny_label'] == 'successful_pass'].index)
    onecell_mutations = list(df_combined[df_combined['phylogeny_label'] == 'cell_specific'].index)
    
    log.info(f"Pass tree mutations: {len(passtree_mutations)}, Cell-specific mutations: {len(onecell_mutations)}")
    
    return df_combined, passtree_mutations, onecell_mutations




##### integrate_mutations_to_scaffold
def attach_mutations_to_current_tree(
    sorted_attached_mutations, T_current, M_current, I_attached, P_attached, 
    ω_NA, fnfp_ratio, φ, logger_obj=None, root_mutations=None, max_retries=None
):
    """
    Process external mutations and integrate them into the phylogenetic tree (with rollback and retry support).
    
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
        Maximum number of candidate positions to try. None means try all.
    
    Returns
    -------
    tuple : (external_mutations, conflict_mutations, T_current, M_current, root_mutations)
    """
    # Use provided logger or fall back to global logger
    log = logger_obj if logger_obj is not None else logging.getLogger(__name__)
    
    if root_mutations is None:
        root_mutations = []
    
    external_mutations = []
    conflict_mutations = []
    
    matrix_index = M_current.index
    i_attached_positive = I_attached.reindex(index=matrix_index).eq(1).fillna(False)
    i_attached_positive_values = i_attached_positive.to_numpy(dtype=bool, copy=False)
    i_attached_positive_col_to_idx = {
        col: idx for idx, col in enumerate(i_attached_positive.columns)
    }
    tree_state = None
    tree_state_dirty = True
    
    def ensure_tree_state():
        nonlocal tree_state, tree_state_dirty
        if tree_state is None or tree_state_dirty:
            tree_state = _build_attach_tree_state(T_current)
            backbone_positive_masks = {}
            backbone_mutation_indices = {}
            for backbone_node, mutation_list in tree_state["mutation_list_under_backbone_nodes"].items():
                if backbone_node in M_current.columns:
                    backbone_positive_masks[backbone_node] = _coerce_bool_array(
                        M_current[backbone_node],
                        matrix_index,
                    )
                backbone_mutation_indices[backbone_node] = [
                    i_attached_positive_col_to_idx[mutation]
                    for mutation in mutation_list
                    if mutation in i_attached_positive_col_to_idx
                ]
            tree_state["backbone_positive_masks"] = backbone_positive_masks
            tree_state["backbone_mutation_indices"] = backbone_mutation_indices
            tree_state_dirty = False
        return tree_state
    
    for new_mut in tqdm(sorted_attached_mutations, desc="Processing mutations", unit="mutation"):
        log.info(f"Processing mutation: {new_mut}")
        current_tree_state = ensure_tree_state()
        
        # Determine which backbone clone new_mut belongs to
        best_backbone, intersection_counts = find_best_backbone_for_new_mutation(
            current_tree_state["mutation_list_under_backbone_nodes"],
            M_current,
            I_attached,
            new_mut,
            i_attached_positive=i_attached_positive,
            i_attached_positive_values=i_attached_positive_values,
            i_attached_positive_col_to_idx=i_attached_positive_col_to_idx,
            backbone_mutation_indices=current_tree_state["backbone_mutation_indices"],
            backbone_positive_masks=current_tree_state["backbone_positive_masks"],
        )
        assigned_nodes = current_tree_state["node_list_under_backbone_nodes"][best_backbone]
        
        # Find intersection nodes
        intersection_nodes = find_all_intersect_muts_from_tree_by_matrix(
            T_current,
            I_attached,
            new_mut,
            matrix_positive=i_attached_positive,
            matrix_positive_values=i_attached_positive_values,
            matrix_positive_col_to_idx=i_attached_positive_col_to_idx,
            tree_nodes=current_tree_state["tree_nodes"],
            node_to_mutations=current_tree_state["node_to_mutations"],
        )
        if len(intersection_nodes) == 0:
            external_mutations.append(new_mut)
            log.info(f"Mutation {new_mut} added to external_mutations (no intersection found)")
            continue
        
        # Get candidate positions using optimized method
        potential_positions = find_intersection_positions_within_tree_directly(
            T_current,
            new_mut,
            I_attached,
            min_overlap=1,
            intersection_nodes=intersection_nodes,
            tree_parent_dict=current_tree_state["tree_parent_dict"],
            node_lookup=current_tree_state["node_lookup"],
        )
        parent_dict = build_parent_dict_from_candidates(potential_positions)
        selected_positions = [p for p in potential_positions if p['anchor'] in assigned_nodes]
        
        # Check if any candidate positions were found
        if len(selected_positions) == 0:
            selected_positions = [p for p in potential_positions if p['anchor'] != 'ROOT']
        
        # ---- Backup current state ----
        M_backup = M_current.copy()
        T_backup = T_current.copy()
        
        # ---- Calculate penalties for all candidate positions ----
        df_penalty = compute_bayesian_penalty_for_all_positions_consider_ROOT(
            new_mut, selected_positions, T_current, M_current, I_attached, P_attached, 
            parent_dict, intersection_nodes, ω_NA=ω_NA, fnfp_ratio=fnfp_ratio, φ=φ,
            logger_obj=log
        )
        
        if df_penalty.empty:
            external_mutations.append(new_mut)
            log.warning(f"Mutation {new_mut} added to external_mutations (no valid penalty scores)")
            continue
        
        # Filter out positions where imputed_vec is all zeros
        df_valid = df_penalty[df_penalty['imputed_vec'].apply(lambda x: x.sum() > 0)]
        if df_valid.empty:
            external_mutations.append(new_mut)
            log.warning(f"Mutation {new_mut} added to external_mutations (all imputed_vec are zero)")
            continue
        
        df_sorted = df_valid.sort_values('total_penalty')
        
        # Determine number of candidates to try
        if max_retries is None:
            candidates_to_try = df_sorted
        else:
            candidates_to_try = df_sorted.head(max_retries)
        
        # ---- Iterate through candidate positions ----
        placed = False
        for attempt, (idx, row) in enumerate(candidates_to_try.iterrows()):
            
            # Restore backup
            M_current = M_backup.copy()
            T_current = T_backup.copy()
            
            try:
                # Apply this position to tree and matrix
                T_current, M_current = apply_position_to_tree(
                    new_mut, row['position'], row['imputed_vec'], 
                    T_current, M_current, I_attached, parent_dict
                )
                
                # Check for conflicts
                if scp.ul.is_conflict_free_gusfield(M_current):
                    placed = True
                    tree_state_dirty = True
                    break
                else:
                    log.warning(f"✗ Position {row['position_index']} caused conflict, trying next candidate")
                    
            except Exception as e:
                log.error(f"Error placing mutation at position {row['position_index']}: {e}")
                continue
        
        # ---- Handle result ----
        if not placed:
            M_current = M_backup.copy()
            T_current = T_backup.copy()
            conflict_mutations.append(new_mut)
            log.warning(f"Mutation {new_mut} added to conflict_mutations (all {len(candidates_to_try)} candidates failed)")
    
    return external_mutations, conflict_mutations, T_current, M_current, root_mutations


def process_rescue_mutations(
    sorted_rescue_mutations, T_current, M_current, I_attached, P_attached, 
    mutation_clones_rescue, ω_NA, fnfp_ratio, φ, logger_obj=None, root_mutations=None, max_retries=None
):
    """
    Process rescue mutations and integrate them into the phylogenetic tree with clone affinity analysis.
    
    Parameters
    ----------
    sorted_rescue_mutations : list
        Sorted list of rescue mutations to process
    T_current : TreeNode
        Current phylogenetic tree structure
    M_current : pd.DataFrame
        Current mutation matrix
    I_attached : pd.DataFrame
        Attached mutation information
    P_attached : pd.DataFrame
        Attached probability information
    mutation_clones_rescue : dict
        Mutation clones for rescue
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
        Maximum number of candidate positions to try. None means try all.
    
    Returns
    -------
    tuple : (external_mutations, conflict_mutations, T_current, M_current, root_mutations)
    """
    # Use provided logger or fall back to global logger
    log = logger_obj if logger_obj is not None else logging.getLogger(__name__)
    
    if root_mutations is None:
        root_mutations = []
    
    external_mutations = []
    conflict_mutations = []
    
    for new_mut in tqdm(sorted_rescue_mutations, desc="Processing rescue mutations", unit="mutation"):
        log.info(f"Processing rescue mutation: {new_mut}")
        
        # Find intersection nodes
        intersection_nodes = find_all_intersect_muts_from_tree_by_matrix(T_current, I_attached, new_mut)
        if len(intersection_nodes) == 0:
            external_mutations.append(new_mut)
            log.info(f"Mutation {new_mut} added to external_mutations (no intersection found)")
            continue
        
        # Get candidate positions using optimized method
        potential_positions = find_intersection_positions_within_tree_directly(
            T_current,
            new_mut,
            I_attached,
            min_overlap=1,
            intersection_nodes=intersection_nodes,
        )
        parent_dict = build_parent_dict_from_candidates(potential_positions)
        
        # Check if any candidate positions were found
        if len(potential_positions) == 0:
            external_mutations.append(new_mut)
            log.info(f"Mutation {new_mut} added to external_mutations (no candidate positions found despite having intersection nodes)")
            continue
        
        # Select which clone to place under based on intersection
        clone_affinity, detailed_scores = compute_new_mut_clone_affinity_correct(
            new_mut, mutation_clones_rescue, I_attached, n_shuffle=100
        )
        assigned_clone = select_best_clone(detailed_scores)
        
        if len(assigned_clone) == 0:
            external_mutations.append(new_mut)
            continue
        
        assigned_clone_muts = []
        for clone in assigned_clone:
            assigned_clone_muts.extend(list(clone))
        
        # Filter candidate positions by clone affinity
        selected_positions = [
            position for position in potential_positions 
            if position['anchor'] in assigned_clone_muts
        ]
        
        # ---- Backup current state ----
        M_backup = M_current.copy()
        T_backup = T_current.copy()
        
        # ---- Calculate penalties for all candidate positions ----
        df_penalty = compute_bayesian_penalty_for_all_positions_consider_ROOT(
            new_mut, selected_positions, T_current, M_current, I_attached, P_attached, 
            parent_dict, intersection_nodes, ω_NA=ω_NA, fnfp_ratio=fnfp_ratio, φ=φ,
            logger_obj=log
        )
        
        if df_penalty.empty:
            external_mutations.append(new_mut)
            log.warning(f"Mutation {new_mut} added to external_mutations (no valid penalty scores)")
            continue
        
        df_valid = df_penalty[df_penalty['imputed_vec'].apply(lambda x: x.sum() > 0)]
        if df_valid.empty:
            external_mutations.append(new_mut)
            log.warning(f"Mutation {new_mut} added to external_mutations (all imputed_vec are zero)")
            continue
        
        df_sorted = df_valid.sort_values('total_penalty')
        
        if max_retries is None:
            candidates_to_try = df_sorted
        else:
            candidates_to_try = df_sorted.head(max_retries)
        
        # ---- Iterate through candidate positions ----
        placed = False
        for attempt, (idx, row) in enumerate(candidates_to_try.iterrows()):
            
            M_current = M_backup.copy()
            T_current = T_backup.copy()
            
            try:
                T_current, M_current = apply_position_to_tree(
                    new_mut, row['position'], row['imputed_vec'], 
                    T_current, M_current, I_attached, parent_dict
                )
                
                if scp.ul.is_conflict_free_gusfield(M_current):
                    placed = True
                    break
                else:
                    log.warning(f"✗ Position {row['position_index']} caused conflict, trying next candidate")
                    
            except Exception as e:
                log.error(f"Error placing mutation at position {row['position_index']}: {e}")
                continue
        
        if placed:
            log.info(f"Updated tree! Current M_current is conflict-free, shape: {M_current.shape}")
        else:
            M_current = M_backup.copy()
            T_current = T_backup.copy()
            conflict_mutations.append(new_mut)
            log.warning(f"Mutation {new_mut} added to conflict_mutations (all candidates failed)")
    
    log.info(f"Processing complete. External: {len(external_mutations)}, Conflict: {len(conflict_mutations)}")
    
    return external_mutations, conflict_mutations, T_current, M_current, root_mutations




###### All the remaining unprocessed items will be added to the new node of the ROOT directory.
def process_external_mutations_by_subtree_groups(
    subtree_groups, T_current, M_current, I_attached, P_attached, 
    ω_NA, fnfp_ratio, φ, logger_obj=None, root_mutations=None, max_retries=None
):
    """
    Process external mutations through subtree groups, supporting both multi-mutation groups
    and single-mutation groups separately (with rollback and retry support).
    
    Parameters
    ----------
    subtree_groups : list of lists
        List of subtree groups, each containing related mutations
    T_current : TreeNode
        Current phylogenetic tree structure
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
        Maximum number of candidate positions to try. None means try all.
    
    Returns
    -------
    tuple : (remained_mutations, conflict_mutations, T_current, M_current, root_mutations)
    """
    # Use provided logger or fall back to global logger
    log = logger_obj if logger_obj is not None else logging.getLogger(__name__)
    
    if root_mutations is None:
        root_mutations = []
    
    remained_mutations = []
    conflict_mutations = []
    
    matrix_index = M_current.index
    i_attached_positive = I_attached.reindex(index=matrix_index).eq(1).fillna(False)
    i_attached_positive_values = i_attached_positive.to_numpy(dtype=bool, copy=False)
    i_attached_positive_col_to_idx = {
        col: idx for idx, col in enumerate(i_attached_positive.columns)
    }
    tree_state = None
    tree_state_dirty = True
    
    def ensure_tree_state():
        nonlocal tree_state, tree_state_dirty
        if tree_state is None or tree_state_dirty:
            tree_state = _build_attach_tree_state(T_current)
            backbone_positive_masks = {}
            backbone_mutation_indices = {}
            for backbone_node, mutation_list in tree_state["mutation_list_under_backbone_nodes"].items():
                if backbone_node in M_current.columns:
                    backbone_positive_masks[backbone_node] = _coerce_bool_array(
                        M_current[backbone_node],
                        matrix_index,
                    )
                backbone_mutation_indices[backbone_node] = [
                    i_attached_positive_col_to_idx[mutation]
                    for mutation in mutation_list
                    if mutation in i_attached_positive_col_to_idx
                ]
            tree_state["backbone_positive_masks"] = backbone_positive_masks
            tree_state["backbone_mutation_indices"] = backbone_mutation_indices
            tree_state_dirty = False
        return tree_state
    
    # Separate subtree groups into multi-mutation and singleton groups
    multi_mut_subtree_groups = [g for g in subtree_groups if len(g) > 1]
    singleton_subtree_groups = [g for g in subtree_groups if len(g) == 1]
    log.info(f"Found {len(multi_mut_subtree_groups)} multi-mutation groups and {len(singleton_subtree_groups)} singleton groups")
    
    ##### Process multi-mutation groups (length > 1)
    for group_idx, group in enumerate(tqdm(multi_mut_subtree_groups, desc="Processing multiple subtrees")):
        log.info(f"Building subtree {group_idx+1}/{len(multi_mut_subtree_groups)} for mutations: {group}")
        
        # Sort by number of 1s in I_attached (descending)
        sorted_group = sorted(group, key=lambda subtree_mut: I_attached[subtree_mut].sum(), reverse=True)
        
        # Add mutations one by one, first attached to ROOT
        reattached_mutations = []
        
        for idx, subtree_mut in enumerate(tqdm(sorted_group, desc="Processing mutations in group")):
            log.info(f"Processing mutation {idx+1}/{len(sorted_group)}: {subtree_mut}")
            
            if idx == 0:
                # First mutation attached directly to ROOT
                T_rollback = copy.deepcopy(T_current)
                M_rollback = M_current.copy()
                final_position = generate_new_leaf_on_root(T_current, subtree_mut)
                T_current = add_new_mutation_to_tree_independent(subtree_mut, T_current, final_position)
                M_current[subtree_mut] = I_attached[subtree_mut].fillna(0).astype(int)
                tree_state_dirty = True
                
                # Use incremental conflict check
                touched_columns = [subtree_mut] if subtree_mut in M_current.columns else []
                if not _is_conflict_free_local_update(
                    M_current,
                    touched_columns,
                    log,
                    f"Steps_SubtreeGroup_First[{subtree_mut}]"
                ):
                    log.warning(f"First mutation {subtree_mut} caused conflict, rolling back")
                    T_current = copy.deepcopy(T_rollback)
                    M_current = M_rollback.copy()
                    reattached_mutations.append(subtree_mut)
                    continue
            
            else:
                current_tree_state = ensure_tree_state()
                
                # Determine which backbone clone subtree_mut belongs to
                best_backbone, intersection_counts = find_best_backbone_for_new_mutation(
                    current_tree_state["mutation_list_under_backbone_nodes"],
                    M_current,
                    I_attached,
                    subtree_mut,
                    i_attached_positive=i_attached_positive,
                    i_attached_positive_values=i_attached_positive_values,
                    i_attached_positive_col_to_idx=i_attached_positive_col_to_idx,
                    backbone_mutation_indices=current_tree_state["backbone_mutation_indices"],
                    backbone_positive_masks=current_tree_state["backbone_positive_masks"],
                )
                assigned_nodes = current_tree_state["node_list_under_backbone_nodes"][best_backbone]
                
                # Find intersection nodes
                intersection_nodes = find_all_intersect_muts_from_tree_by_matrix(
                    T_current,
                    I_attached,
                    subtree_mut,
                    matrix_positive=i_attached_positive,
                    matrix_positive_values=i_attached_positive_values,
                    matrix_positive_col_to_idx=i_attached_positive_col_to_idx,
                    tree_nodes=current_tree_state["tree_nodes"],
                    node_to_mutations=current_tree_state["node_to_mutations"],
                )
                if len(intersection_nodes) == 0:
                    reattached_mutations.append(subtree_mut)
                    log.info(f"Mutation {subtree_mut} added to reattached_mutations (no intersection found)")
                    continue
                
                # Get candidate positions using optimized method
                potential_positions = find_intersection_positions_within_tree_directly(
                    T_current,
                    subtree_mut,
                    I_attached,
                    min_overlap=1,
                    intersection_nodes=intersection_nodes,
                    tree_parent_dict=current_tree_state["tree_parent_dict"],
                    node_lookup=current_tree_state["node_lookup"],
                )
                parent_dict = build_parent_dict_from_candidates(potential_positions)
                selected_positions = [p for p in potential_positions if p['anchor'] in assigned_nodes]
                
                # Check if any candidate positions were found
                if len(selected_positions) == 0:
                    reattached_mutations.append(subtree_mut)
                    log.info(f"Mutation {subtree_mut} added to reattached_mutations (no candidate positions found despite having intersection nodes)")
                    continue
                
                # ---- Backup current state ----
                M_backup = M_current.copy()
                T_backup = T_current.copy()
                
                # ---- Calculate penalties for all candidate positions ----
                df_penalty = compute_bayesian_penalty_for_all_positions_consider_ROOT(
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
                
                if max_retries is None:
                    candidates_to_try = df_sorted
                else:
                    candidates_to_try = df_sorted.head(max_retries)
                
                # ---- Iterate through candidate positions ----
                placed = False
                for attempt, (idx_row, row) in enumerate(candidates_to_try.iterrows()):
                    
                    M_current = M_backup.copy()
                    T_current = T_backup.copy()
                    
                    try:
                        T_current, M_current = apply_position_to_tree(
                            subtree_mut, row['position'], row['imputed_vec'], 
                            T_current, M_current, I_attached, parent_dict
                        )
                        
                        # Use incremental conflict check
                        touched_columns = _resolve_touched_columns_for_conflict_check(row['position'], parent_dict, M_current.columns, subtree_mut)
                        if _is_conflict_free_local_update(
                            M_current,
                            touched_columns,
                            log,
                            f"Steps_SubtreeGroup[{subtree_mut}]_attempt{attempt}"
                        ):
                            placed = True
                            tree_state_dirty = True
                            log.info(f"✓ Mutation {subtree_mut} placed at position {row['position_index']} (attempt {attempt+1})")
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
        
        # Process reattached mutations (those not placed in the first round)
        second_reattached_mutations = []
        if len(reattached_mutations) > 0:            
            log.info(f"Processing {len(reattached_mutations)} reattached mutations for group {group_idx+1}")
            
            sorted_reattached_mutations = [i for i in I_attached.columns if i in reattached_mutations]
            for subtree_mut in tqdm(sorted_reattached_mutations, desc="Processing re-attached mutations"):
                log.info(f"Processing re-attached mutation: {subtree_mut}")                
                current_tree_state = ensure_tree_state()

                # Determine which backbone clone subtree_mut belongs to
                best_backbone, intersection_counts = find_best_backbone_for_new_mutation(
                    current_tree_state["mutation_list_under_backbone_nodes"],
                    M_current,
                    I_attached,
                    subtree_mut,
                    i_attached_positive=i_attached_positive,
                    i_attached_positive_values=i_attached_positive_values,
                    i_attached_positive_col_to_idx=i_attached_positive_col_to_idx,
                    backbone_mutation_indices=current_tree_state["backbone_mutation_indices"],
                    backbone_positive_masks=current_tree_state["backbone_positive_masks"],
                )
                assigned_nodes = current_tree_state["node_list_under_backbone_nodes"][best_backbone]
                
                # Find intersection nodes
                intersection_nodes = find_all_intersect_muts_from_tree_by_matrix(
                    T_current,
                    I_attached,
                    subtree_mut,
                    matrix_positive=i_attached_positive,
                    matrix_positive_values=i_attached_positive_values,
                    matrix_positive_col_to_idx=i_attached_positive_col_to_idx,
                    tree_nodes=current_tree_state["tree_nodes"],
                    node_to_mutations=current_tree_state["node_to_mutations"],
                )
                if len(intersection_nodes) == 0:
                    second_reattached_mutations.append(subtree_mut)
                    log.info(f"Mutation {subtree_mut} added to second_reattached_mutations (no intersection found)")
                    continue
                
                # Get candidate positions using optimized method
                potential_positions = find_intersection_positions_within_tree_directly(
                    T_current,
                    subtree_mut,
                    I_attached,
                    min_overlap=1,
                    intersection_nodes=intersection_nodes,
                    tree_parent_dict=current_tree_state["tree_parent_dict"],
                    node_lookup=current_tree_state["node_lookup"],
                )
                parent_dict = build_parent_dict_from_candidates(potential_positions)
                selected_positions = [p for p in potential_positions if p['anchor'] in assigned_nodes]
                
                if len(selected_positions) == 0:
                    second_reattached_mutations.append(subtree_mut)
                    log.info(f"Mutation {subtree_mut} added to second_reattached_mutations (no candidate positions found despite having intersection nodes)")
                    continue
                
                # ---- Backup current state ----
                M_backup = M_current.copy()
                T_backup = T_current.copy()
                
                # ---- Calculate penalties for all candidate positions ----
                df_penalty = compute_bayesian_penalty_for_all_positions_consider_ROOT(
                    subtree_mut, selected_positions, T_current, M_current, I_attached, P_attached, 
                    parent_dict, intersection_nodes, ω_NA=ω_NA, fnfp_ratio=fnfp_ratio, φ=φ,
                    logger_obj=log
                )
                
                if df_penalty.empty:
                    second_reattached_mutations.append(subtree_mut)
                    log.warning(f"Mutation {subtree_mut} added to second_reattached_mutations (no valid penalty scores)")
                    continue
                
                df_valid = df_penalty[df_penalty['imputed_vec'].apply(lambda x: x.sum() > 0)]
                if df_valid.empty:
                    second_reattached_mutations.append(subtree_mut)
                    log.warning(f"Mutation {subtree_mut} added to second_reattached_mutations (all imputed_vec are zero)")
                    continue
                
                df_sorted = df_valid.sort_values('total_penalty')
                
                if max_retries is None:
                    candidates_to_try = df_sorted
                else:
                    candidates_to_try = df_sorted.head(max_retries)
                
                # ---- Iterate through candidate positions ----
                placed = False
                for attempt, (idx_row, row) in enumerate(candidates_to_try.iterrows()):
                    
                    M_current = M_backup.copy()
                    T_current = T_backup.copy()
                    
                    try:
                        T_current, M_current = apply_position_to_tree(
                            subtree_mut, row['position'], row['imputed_vec'], 
                            T_current, M_current, I_attached, parent_dict
                        )
                        
                        # Use incremental conflict check
                        touched_columns = _resolve_touched_columns_for_conflict_check(row['position'], parent_dict, M_current.columns, subtree_mut)
                        if _is_conflict_free_local_update(
                            M_current,
                            touched_columns,
                            log,
                            f"Steps_SubtreeReattach[{subtree_mut}]_attempt{attempt}"
                        ):
                            log.info(f"✓ Mutation {subtree_mut} successfully placed (score={row['total_penalty']:.4f})")
                            placed = True
                            tree_state_dirty = True
                            break
                        else:
                            log.warning(f"✗ Position {row['position_index']} caused conflict, trying next candidate")
                            
                    except Exception as e:
                        log.error(f"Error placing mutation at position {row['position_index']}: {e}")
                        continue
                
                if not placed:
                    M_current = M_backup.copy()
                    T_current = T_backup.copy()
                    second_reattached_mutations.append(subtree_mut)
                    log.warning(f"Mutation {subtree_mut} added to second_reattached_mutations (all candidates failed)")
        
        # Record mutations still not processed
        if second_reattached_mutations:
            log.warning(f"Group {group_idx+1} has {len(second_reattached_mutations)} mutations still remaining: {second_reattached_mutations}")
            remained_mutations.extend(second_reattached_mutations)
    
    ##### Process singleton mutation groups (length = 1)
    for group in tqdm(singleton_subtree_groups, desc="Processing singleton subtrees"):
        subtree_mut = group[0]
        log.info(f"Attaching singleton mutation directly to ROOT: {subtree_mut}")
        
        T_rollback = copy.deepcopy(T_current)
        M_rollback = M_current.copy()            
        
        final_position = generate_new_leaf_on_root(T_current, subtree_mut)
        T_current = add_new_mutation_to_tree_independent(subtree_mut, T_current, final_position)
        M_current[subtree_mut] = I_attached[subtree_mut].fillna(0).astype(int)
        tree_state_dirty = True
                
        # Use incremental conflict check
        touched_columns = [subtree_mut] if subtree_mut in M_current.columns else []
        if not _is_conflict_free_local_update(
            M_current,
            touched_columns,
            log,
            f"Steps_SubtreeSingleton[{subtree_mut}]"
        ):
            log.warning(f"Conflict detected after adding {subtree_mut}, rolling back")
            T_current = copy.deepcopy(T_rollback)
            M_current = M_rollback.copy()
            tree_state_dirty = True
            remained_mutations.append(subtree_mut)
            log.info(f"Mutation {subtree_mut} added to remained_mutations due to conflict")
            continue
    
    log.info("All external mutations have been processed successfully.")
    
    if not _is_conflict_free_gusfield_safe(M_current, log, "Steps_SubtreeGroups[final]"):
        log.error("⚠️ Final matrix has conflicts! But continuing as per configuration.")
    
    return remained_mutations, conflict_mutations, T_current, M_current, root_mutations





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

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
_VERBOSE_TREE_UPDATES = os.environ.get("PHYLOSOLID_VERBOSE_TREES", "").lower() in {"1", "true", "yes", "on"}

# Default parameters matching Methods Section 3
DEFAULT_PARAMS = {
    "posterior_threshold": 0.5,        # Phylogenetic mosaic threshold
    "maf_max_threshold": 0.3,           # 最大突变频率阈值
    "maf_mean_threshold": 0.1,          # 平均突变频率阈值
    "coverage_threshold": 10,           # 读取的覆盖度阈值
    "logistic_regression_threshold": 0.5,  # Logistic回归判定阈值
    # 其他必要的参数可以继续添加
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
    计算 I_attached 和 M_current 中对应行的条件个数：
    1. 前者为1且后者为0的个数
    2. 前者为0且后者为1的个数
    3. 前者为NaN且后者为1的个数
    4. 前者为NaN且后者为0的个数
    
    Parameters:
    I_attached_col: pandas Series, 来自 I_attached 的某一列
    M_current_col: pandas Series, 来自 M_current 的某一列
    
    Returns:
    dict: 包含四个统计值的字典
    """
    # 转换为对应的布尔值
    I_attached_col = I_attached_col.fillna(np.nan)
    M_current_col = M_current_col.fillna(np.nan)
    
    count_1_0 = ((I_attached_col == 1) & (M_current_col == 0)).sum()  # I_attached 为 1 且 M_current 为 0
    count_0_1 = ((I_attached_col == 0) & (M_current_col == 1)).sum()  # I_attached 为 0 且 M_current 为 1
    count_na_1 = ((I_attached_col.isna()) & (M_current_col == 1)).sum()  # I_attached 为 NaN 且 M_current 为 1
    count_na_0 = ((I_attached_col.isna()) & (M_current_col == 0)).sum()  # I_attached 为 NaN 且 M_current 为 0
    
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
    spots_to_split: list,  # 需要拆分的spot列表
    immune_mutations: list,  # 免疫突变的列表
    I_process: pd.DataFrame,  # I 突变矩阵
    P_process: pd.DataFrame   # P 矩阵 (posterior probability)
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    拆分给定的 spots_to_split，基于免疫突变拆分成两个版本：一个是 -immune 后缀，另一个是 -non 后缀。
    同时，保留不需要拆分的行。
    
    参数:
    ----------
    spots_to_split: list
        需要拆分的spot列表（每个spot是矩阵I_process中的行）
        
    immune_mutations: list
        免疫突变的列表
    
    I_process: pd.DataFrame
        突变矩阵 (行是spot，列是突变)
    
    P_process: pd.DataFrame
        Posterior probability 矩阵
    
    返回值:
    ----------
    I_process_resolved: pd.DataFrame
        处理后的突变矩阵
        
    P_process_resolved: pd.DataFrame
        处理后的 Posterior probability 矩阵
    """
    
    # 存储拆分后的行数据
    resolved_I_rows = []
    resolved_P_rows = []
    resolved_index = []
    
    # 遍历每个spot
    for spot_name, row in I_process.iterrows():
        # 如果该行需要拆分
        if spot_name in spots_to_split:
            # 获取该spot的原始行
            I_row = row
            P_row = P_process.loc[spot_name]
            
            # 创建 -immune 行，保留免疫突变的列，其他列为NA
            immune_row = I_row.copy()
            immune_row[immune_mutations] = I_row[immune_mutations]  # 免疫突变保留
            immune_row[~I_row.index.isin(immune_mutations)] = np.nan  # 非免疫突变置为NA
            
            P_immune_row = P_row.copy()
            P_immune_row[~P_row.index.isin(immune_mutations)] = np.nan  # 非免疫突变置为NA
            
            # 创建 -non 行，免疫突变列置为NA，其他列为原始值
            non_immune_row = I_row.copy()
            non_immune_row[immune_mutations] = np.nan  # 免疫突变置为NA
            
            P_non_immune_row = P_row.copy()
            P_non_immune_row[immune_mutations] = np.nan  # 免疫突变置为NA
            
            # 将拆分后的行分别加入到结果列表中
            resolved_I_rows.append(immune_row)
            resolved_P_rows.append(P_immune_row)
            resolved_index.append(f"{spot_name}-immune")
            
            resolved_I_rows.append(non_immune_row)
            resolved_P_rows.append(P_non_immune_row)
            resolved_index.append(f"{spot_name}-non")
        else:
            # 不需要拆分的行直接保留
            resolved_I_rows.append(row)
            resolved_P_rows.append(P_process.loc[spot_name])
            resolved_index.append(spot_name)
    
    # 将拆分后的行合并成新的DataFrame
    I_process_resolved = pd.DataFrame(resolved_I_rows, index=resolved_index, columns=I_process.columns)
    P_process_resolved = pd.DataFrame(resolved_P_rows, index=resolved_index, columns=P_process.columns)
    
    return I_process_resolved, P_process_resolved

def merge_mutations(M_current_each_mut, all_nodes_in_T_scaffold):
    # 函数：按照 all_nodes_in_T_scaffold 合并 M_current_each_mut 中的列
    # 创建一个空的列表用于存储合并后的列
    merged_columns = []
    
    # 遍历 all_nodes_in_T_scaffold
    for mutation in all_nodes_in_T_scaffold:
        # 如果突变是由 | 分隔的多个突变
        if '|' in mutation:
            # 分割突变名称，获取多个突变列
            mutations = mutation.split('|')
            # 获取这些列对应的数据
            cols = [M_current_each_mut[col] for col in mutations if col in M_current_each_mut.columns]
            if cols:
                # 合并这些列，使用 max 来合并（取最大的值，或者根据需求修改）
                merged_column = pd.concat(cols, axis=1).max(axis=1)
                merged_columns.append(merged_column)
            else:
                # 如果没有匹配的列（避免出错），用 NaN 填充
                merged_columns.append(pd.Series([pd.NA] * M_current_each_mut.shape[0], index=M_current_each_mut.index))
        else:
            # 如果是单一突变名称，直接取对应列
            if mutation in M_current_each_mut.columns:
                merged_columns.append(M_current_each_mut[mutation])
            else:
                # 如果没有找到该列，添加 NaN
                merged_columns.append(pd.Series([pd.NA] * M_current_each_mut.shape[0], index=M_current_each_mut.index))
    
    # 将合并后的列组合成新的 DataFrame
    merged_df = pd.concat(merged_columns, axis=1)
    
    # 更新列名为 all_nodes_in_T_scaffold 中的元素
    merged_df.columns = all_nodes_in_T_scaffold
    
    return merged_df


def find_column_in_merged_columns(df, mutation_name):
    """
    在DataFrame的列中查找突变名，支持合并列（用'|'分隔）
    
    Parameters:
    -----------
    df : pd.DataFrame
        要搜索的DataFrame
    mutation_name : str
        要查找的突变名
    
    Returns:
    --------
    str or None : 找到的列名，如果没找到返回None
    """
    if mutation_name in df.columns:
        return mutation_name
    
    # 在合并列中查找
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

# 配置日志
logger = logging.getLogger(__name__)

def find_intersection_positions_within_tree_directly(T_current: TreeNode, new_mut: str, matrix, min_overlap=1):
    """
    基于交集分析的优化版本，直接找到相关位置
    """
    # 1. 找到所有与目标突变有交集的节点
    intersection_nodes = find_all_intersect_muts_from_tree_by_matrix(
        T_current, matrix, new_mut, min_overlap
    )
    
    
    if len(intersection_nodes) == 0:
        logger.debug(f"No intersection nodes found for {new_mut}")
        return []  # 没有交集，返回空列表
    
    # 2. 构建树的父子关系字典
    tree_parent_dict = build_tree_parent_dict(T_current)
    
    # 3. 找到所有相关路径上的节点
    all_path_nodes = find_all_path_nodes(intersection_nodes, tree_parent_dict)
    
    
    # 4. 预先创建基础树的深拷贝
    base_tree_copy = deepcopy(T_current)
    
    # 5. 只在这些相关节点上生成候选位置
    candidate_positions = []
    
    for node_name in all_path_nodes:
        if node_name == "ROOT":
            # 保留 ROOT 上的 on_node 类型的 position
            node = base_tree_copy.find(node_name)
            if node is None:
                logger.warning(f"Node {node_name} not found in tree")
                continue
            
            # 生成 ROOT 上的 on_node 类型的候选位置
            candidate_positions.append(_create_on_node_candidate_fast(base_tree_copy, node, new_mut))
            
            continue  # ROOT 的其他类型位置仍然跳过
        
        node = base_tree_copy.find(node_name)
        if node is None:
            logger.warning(f"Node {node_name} not found in tree")
            continue
        
        # --- 1) 放在 node 上 ---
        candidate_positions.append(_create_on_node_candidate_fast(base_tree_copy, node, new_mut))
        
        # --- 2) 新 leaf ---
        candidate_positions.append(_create_new_leaf_candidate_fast(base_tree_copy, node, new_mut))
        
        # --- 3) 放在每条 edge 上 ---
        for child in node.children:
            if child.name in all_path_nodes:  # 只考虑路径上的子节点
                candidate_positions.append(_create_on_edge_candidate_fast(base_tree_copy, node, child, new_mut))
        
        # --- 4) 新 parent merge 多个子节点 ---
        if len(node.children) >= 2:
            # 只考虑路径上的子节点组合
            path_children = [child for child in node.children if child.name in all_path_nodes]
            if len(path_children) >= 2:
                # 限制组合数量避免爆炸
                for r in range(2, min(4, len(path_children) + 1)):
                    for combo in combinations(path_children, r):
                        candidate_positions.append(_create_merge_candidate_fast(base_tree_copy, node, combo, new_mut))
    
    return candidate_positions


def build_tree_parent_dict(tree):
    """直接从树构建父子关系字典"""
    parent_dict = {}
    for node in tree.traverse():
        for child in node.children:
            parent_dict[child.name] = node.name
    return parent_dict


def _create_on_node_candidate_fast(base_tree_copy, node, new_mut):
    """快速创建放在节点上的候选"""
    new_tree = deepcopy(base_tree_copy)
    anchor_node = new_tree.find(node.name)
    
    # 替换节点逻辑
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
        "nodes": _extract_nodes_info(new_tree),
        "edges": _extract_edges_info(new_tree)
    }


def _create_new_leaf_candidate_fast(base_tree_copy, node, new_mut):
    """快速创建新叶子的候选"""
    new_tree = deepcopy(base_tree_copy)
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
    """快速创建放在边上的候选"""
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
        "nodes": _extract_nodes_info(new_tree),
        "edges": _extract_edges_info(new_tree)
    }


def _create_merge_candidate_fast(base_tree_copy, parent_node, children_combo, new_mut):
    """快速创建合并子节点的候选"""
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
        "nodes": _extract_nodes_info(new_tree),
        "edges": _extract_edges_info(new_tree)
    }


def _extract_nodes_info(tree):
    """提取节点信息"""
    return [{"name": n.name,
             "parent": n.parent.name if n.parent else None,
             "children": [c.name for c in n.children]} 
            for n in tree.traverse()]


def _extract_edges_info(tree):
    """提取边信息"""
    return [(n.parent.name, n.name) for n in tree.traverse() if n.parent]


def build_parent_dict_from_candidates(candidate_positions):
    """从 candidate_positions 中构建父子关系字典"""
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


# def find_all_intersect_muts_from_tree_by_matrix(tree, matrix, target_mut, min_overlap=1):
#     """
#     返回所有与 target_mut 在 matrix 中有至少 min_overlap 个细胞共同出现的突变集合。
#     """
#     intersect_muts = set()
    
#     # 遍历所有节点（除 ROOT）
#     for node in tree.all_nodes():
#         if node.name == "ROOT":
#             continue
        
#         for mut in node.name.split("|"):
#             if mut == target_mut:
#                 continue  # 跳过自身
            
#             if mut not in matrix.columns or target_mut not in matrix.columns:
#                 continue  # 确保突变存在于矩阵中
            
#             # 计算两突变的共现数量
#             mut_vec = matrix[mut].fillna(0)
#             target_vec = matrix[target_mut].fillna(0)
#             N11 = ((mut_vec == 1) & (target_vec == 1)).sum()
            
#             # 共现数量足够才认为有交集
#             if N11 >= min_overlap:
#                 intersect_muts.add(mut)
    
#     return intersect_muts

def find_all_intersect_muts_from_tree_by_matrix(tree, matrix, target_mut, min_overlap=1):
    """
    返回所有与 target_mut 在 matrix 中有至少 min_overlap 个细胞共同出现的树节点集合。
    """
    intersect_nodes = set()
    
    # 遍历所有节点（除 ROOT）
    for node in tree.all_nodes():
        if node.name == "ROOT":
            continue
        
        has_intersection = False
        
        # 检查节点中的每个突变是否与目标突变有交集
        for mut in node.name.split("|"):
            if mut == target_mut:
                continue  # 跳过自身
            
            if mut not in matrix.columns or target_mut not in matrix.columns:
                continue  # 确保突变存在于矩阵中
            
            # 计算两突变的共现数量
            mut_vec = matrix[mut].fillna(0)
            target_vec = matrix[target_mut].fillna(0)
            N11 = ((mut_vec == 1) & (target_vec == 1)).sum()
            
            # 如果任何一个突变有足够的共现，就标记这个节点
            if N11 >= min_overlap:
                has_intersection = True
                break  # 只要节点中有一个突变有交集就够了
        
        # 如果节点中有任何突变与目标突变有交集，就添加整个节点
        if has_intersection:
            intersect_nodes.add(node.name)
    
    return intersect_nodes

def find_all_path_nodes(intersection_nodes, tree_parent_dict):
    """
    找到所有在连接交集节点的路径上的节点
    """
    all_path_nodes = set()
    all_path_nodes.add('ROOT')  # 总是包含 ROOT
    
    # 对于每个交集节点，找到从 ROOT 到该节点的路径
    for node in intersection_nodes:
        path_to_root = get_path_to_root(node, tree_parent_dict)
        all_path_nodes.update(path_to_root)
    
    # 找到连接不同交集节点的路径
    intersection_list = list(intersection_nodes)
    for i in range(len(intersection_list)):
        for j in range(i + 1, len(intersection_list)):
            path_between = get_path_between_nodes(intersection_list[i], intersection_list[j], tree_parent_dict)
            all_path_nodes.update(path_between)
    
    return all_path_nodes


def get_path_to_root(node, tree_parent_dict):
    """找到从节点到 ROOT 的路径"""
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
    """找到两个节点之间的路径"""
    # 找到从 node1 到 ROOT 的路径
    path1 = get_path_to_root(node1, tree_parent_dict)
    # 找到从 node2 到 ROOT 的路径  
    path2 = get_path_to_root(node2, tree_parent_dict)
    
    # 反转路径，使其从 ROOT 开始
    path1_from_root = list(reversed(path1))
    path2_from_root = list(reversed(path2))
    
    # 找到最近公共祖先 (LCA)
    lca = None
    for i in range(min(len(path1_from_root), len(path2_from_root))):
        if path1_from_root[i] == path2_from_root[i]:
            lca = path1_from_root[i]
        else:
            break
    
    if lca is None:
        return []
    
    # 构建完整路径: node1 -> LCA -> node2
    lca_index1 = path1_from_root.index(lca)
    lca_index2 = path2_from_root.index(lca)
    
    path_node1_to_lca = path1_from_root[lca_index1:]
    path_lca_to_node2 = path2_from_root[lca_index2:]
    
    # 合并路径，去掉重复的 LCA
    full_path = path_node1_to_lca + path_lca_to_node2[1:]
    return full_path




# -------------------------
# 用重复 100 遍的方法计算 new_mut 最应该属于哪一个 clone
# -------------------------

def compute_corr_cache_with_new_mut(I_attached, existing_muts, new_mut):
    """
    计算包含新突变的相关性缓存
    """
    # 所有突变（现有突变 + 新突变）
    all_muts = existing_muts + [new_mut]
    
    corr_cache = {}
    
    # 计算所有突变对的相关性（包括新突变与现有突变）
    for u, v in itertools.combinations(all_muts, 2):
        corr = are_mutations_correlated(I_attached, u, v)
        corr_cache[(u, v)] = corr
        corr_cache[(v, u)] = corr
    
    # 自相关设为True
    for m in all_muts:
        corr_cache[(m, m)] = True
    
    return corr_cache

def compute_new_mut_clone_affinity_correct(new_mut, mutation_clones_rescue, I_attached, n_shuffle=100):
    """
    正确版本：使用包含新突变的corr_cache
    """
    # 获取所有现有突变
    all_existing_muts = []
    for clone in mutation_clones_rescue.values():
        all_existing_muts.extend(clone)
    all_existing_muts = list(set(all_existing_muts))
    
    # 重新计算包含新突变的corr_cache
    # print("重新计算包含新突变的相关性缓存...")
    corr_cache = compute_corr_cache_with_new_mut(I_attached, all_existing_muts, new_mut)
    
    clone_affinity = {}
    detailed_scores = {}
    
    print(f"calculate the correlation with tree mutations for : {new_mut}")
    
    for clone_rep, clone_muts in mutation_clones_rescue.items():
        clone_key = tuple(sorted(clone_muts))
        
        # 检查新突变与克隆成员的相关性
        direct_corr_count = 0
        correlation_details = []
        
        for existing_mut in clone_muts:
            key1 = (new_mut, existing_mut)
            key2 = (existing_mut, new_mut)
            
            is_corr = corr_cache.get(key1, False) or corr_cache.get(key2, False)
            
            if is_corr:
                direct_corr_count += 1
                # 获取详细的count数据用于调试
                counts = pairwise_counts(I_attached, new_mut, existing_mut)
                correlation_details.append(f"{existing_mut}: N11={counts['N11']}, J={jaccard_index(I_attached, new_mut, existing_mut):.3f}")
            else:
                counts = pairwise_counts(I_attached, new_mut, existing_mut)
                correlation_details.append(f"{existing_mut}: N11={counts['N11']} (not corr)")
        
        direct_ratio = direct_corr_count / len(clone_muts) if clone_muts else 0
        
        # print(f"\n克隆 {clone_rep}: {len(clone_muts)}个突变")
        # print(f"  直接关联: {direct_corr_count}/{len(clone_muts)} = {direct_ratio:.3f}")
        
        # # 显示前几个相关性详情
        # for detail in correlation_details[:3]:
        #     print(f"    {detail}")
        # if len(correlation_details) > 3:
        #     print(f"    ... 还有 {len(correlation_details) - 3} 个突变")
        
        # 计算与克隆代表突变的关系
        rep_correlation = 0
        if clone_rep in clone_muts:
            key1 = (new_mut, clone_rep)
            key2 = (clone_rep, new_mut)
            rep_correlation = 1 if (corr_cache.get(key1, False) or corr_cache.get(key2, False)) else 0
            counts = pairwise_counts(I_attached, new_mut, clone_rep)
            # print(f"  与代表突变关联: {rep_correlation} (N11={counts['N11']}, J={jaccard_index(I_attached, new_mut, clone_rep):.3f})")
        
        # 计算加权分数（考虑克隆大小）
        base_score = direct_ratio
        size_factor = min(1.0, len(clone_muts) / 10)  # 避免过大克隆过度影响
        final_score = base_score * (0.7 + 0.3 * size_factor)
        
        clone_affinity[clone_key] = final_score
        detailed_scores[clone_key] = {
            'direct_ratio': direct_ratio,
            'direct_correlations': direct_corr_count,
            'rep_correlation': rep_correlation,
            'clone_size': len(clone_muts)
        }
        
        # print(f"  最终分数: {final_score:.3f}")
    
    return clone_affinity, detailed_scores

# # 使用示例
# new_mutation = "chr16_67944430_G_C"

# clone_affinity, detailed_scores = compute_new_mut_clone_affinity_correct(
#     new_mutation, 
#     mutation_clones_rescue, 
#     I_attached,
#     n_shuffle=100
# )

def select_max_affinity_clone(clone_affinity):
    # 筛选出大于0的键值对
    filtered_affinities = {k: v for k, v in clone_affinity.items() if v > 0}
    
    # 如果没有大于0的项，返回空列表
    if not filtered_affinities:
        return []
    
    # 找到最大的 affinity score
    max_affinity = max(filtered_affinities.values())
    
    # 筛选出所有具有最大 affinity score 的 clone
    max_clones = [k for k, v in filtered_affinities.items() if v == max_affinity]
    
    return max_clones

def select_best_clone(detailed_scores):
    # 第一步: 找到 rep_correlation 最大的值
    max_rep_correlation = max([score['rep_correlation'] for score in detailed_scores.values()])
    
    # 如果 rep_correlation 最大值是 0, 直接返回空列表
    if max_rep_correlation == 0:
        return []
    
    # 第二步: 找到所有 rep_correlation 为最大值的 clone
    max_rep_correlation_clones = [
        clone for clone, score in detailed_scores.items() if score['rep_correlation'] == max_rep_correlation
    ]
    
    # 第三步: 如果只有一个 clone, 直接返回它
    if len(max_rep_correlation_clones) == 1:
        return max_rep_correlation_clones
    
    # 第四步: 如果有多个 clone, 比较它们的 direct_ratio
    max_direct_ratio = max([detailed_scores[clone]['direct_ratio'] for clone in max_rep_correlation_clones])
    
    # 如果 direct_ratio 最大值是 0, 返回空列表
    if max_direct_ratio == 0:
        return []
    
    # 第五步: 返回 direct_ratio 最大值的 clone
    best_clones = [
        clone for clone in max_rep_correlation_clones 
        if detailed_scores[clone]['direct_ratio'] == max_direct_ratio
    ]
    
    return best_clones

# best_clones = select_best_clone(detailed_scores)




# -------------------------
# 计算 bayesian 罚分
# -------------------------

def compute_bayesian_penalty_for_all_positions_consider_ROOT(
    new_mut, selected_positions, T_current, M_current, I_selected, P_selected, parent_dict, intersection_nodes, 
    ω_NA=0.001, fnfp_ratio=0.1, φ=1
):
    """
    计算所有候选位置的罚分，不修改 M_current
    
    Returns:
    --------
    df_penalty : pandas.DataFrame
        包含所有候选位置罚分的DataFrame
    """
    import pandas as pd
    import numpy as np
    
    results = []
    
    if len(selected_positions) == 0:
        logger.warning(f"No selected positions for mutation {new_mut}")
        return pd.DataFrame()
    
    new_mut_bin_vector = I_selected[new_mut].replace({pd.NA: np.nan}).fillna(0).astype(int)
        
    input_binary_vec_full = I_selected[new_mut].replace({pd.NA: np.nan})
    na_ratio = input_binary_vec_full.isna().mean()
    mut_ratio = input_binary_vec_full.fillna(0).mean()
    N_nodes_beforeT = len(T_current.all_nodes())
    
    for idx, pos in enumerate(selected_positions):
        placement_type = pos['placement_type']
        anchor = pos['anchor']
        
        imputed_vec = pd.Series(0, index=M_current.index)
        merge_penalty = 0
        
        # ============================================================
        # 根据放置类型计算 imputed vector
        # ============================================================
        
        if placement_type == 'on_node':
            parent = anchor
            vec_parent = M_current[parent] if parent != 'ROOT' else pd.Series(1, index=M_current.index)
            
            sibling_nodes = [n['name'] for n in pos['nodes'] if n['parent'] == parent and n['name'] != new_mut]
            lineage_conflict_nodes = get_all_conflict_nodes_outside_lineage(
                parent, build_lineage_parent_dict_from_tree(T_current, anchor), M_current.columns
            )
            all_conflict_nodes = list(set(sibling_nodes + lineage_conflict_nodes))
            
            vec_conflicts = pd.Series(0, index=M_current.index)
            for conflict in all_conflict_nodes:
                conflict_series = M_current[conflict].reindex(M_current.index, fill_value=0)
                vec_conflicts = vec_conflicts | conflict_series
            
            new_mut_cleaned = new_mut_bin_vector & ~vec_conflicts
            anchor_series = M_current[anchor].reindex(new_mut_cleaned.index, fill_value=0)
            imputed_vec = (anchor_series | new_mut_cleaned).astype(int)
            N_nodes = N_nodes_beforeT + 1
            
        elif placement_type == 'new_leaf':
            parent = anchor
            vec_parent = M_current[parent] if parent != 'ROOT' else pd.Series(1, index=M_current.index)
            
            sibling_nodes = [n['name'] for n in pos['nodes'] if n['parent'] == parent and n['name'] != new_mut]
            lineage_conflict_nodes = get_all_conflict_nodes_outside_lineage(
                parent, build_lineage_parent_dict_from_tree(T_current, anchor), M_current.columns
            )
            all_conflict_nodes = list(set(sibling_nodes + lineage_conflict_nodes))
            
            vec_conflicts = pd.Series(0, index=M_current.index)
            for conflict in all_conflict_nodes:
                conflict_series = M_current[conflict].reindex(M_current.index, fill_value=0)
                vec_conflicts = vec_conflicts | conflict_series
            
            new_mut_cleaned = new_mut_bin_vector & ~vec_conflicts
            imputed_vec = new_mut_cleaned.astype(int)
            N_nodes = N_nodes_beforeT + 2
            
        elif placement_type == 'on_edge':
            parent = anchor
            child = pos['meta']['child']
            vec_parent = M_current[parent] if parent != 'ROOT' else pd.Series(1, index=M_current.index)
            vec_child = M_current[child]
            
            sibling_nodes = [n['name'] for n in pos['nodes'] if n['parent']==parent and n['name'] not in [child,new_mut]]
            lineage_conflict_nodes = get_all_conflict_nodes_outside_lineage(
                parent, build_lineage_parent_dict_from_tree(T_current, anchor), M_current.columns
            )
            all_conflict_nodes = list(set(sibling_nodes + lineage_conflict_nodes))
            
            vec_conflicts = pd.Series(0, index=M_current.index)
            for conflict in all_conflict_nodes:
                conflict_series = M_current[conflict].reindex(M_current.index, fill_value=0)
                vec_conflicts = vec_conflicts | conflict_series
            
            new_mut_cleaned = new_mut_bin_vector & ~vec_conflicts
            vec_child_aligned = vec_child.reindex(new_mut_cleaned.index, fill_value=0)
            imputed_vec = (vec_child_aligned | new_mut_cleaned).astype(int)
            N_nodes = N_nodes_beforeT + 2
            
        elif placement_type == 'new_parent_merge':
            parent = anchor
            merge_children = pos['meta']['merge_children']
            vec_parent = M_current[parent] if parent != 'ROOT' else pd.Series(1, index=M_current.index)
            
            vec_children = pd.Series(0, index=M_current.index)
            for c in merge_children:
                child_series = M_current[c].reindex(M_current.index, fill_value=0)
                vec_children = vec_children | child_series
            
            sibling_nodes = [n['name'] for n in pos['nodes'] if n['parent']==parent and n['name'] not in merge_children+[new_mut]]
            lineage_conflict_nodes = get_all_conflict_nodes_outside_lineage(
                parent, build_lineage_parent_dict_from_tree(T_current, anchor), M_current.columns, exclude_nodes=merge_children
            )
            all_conflict_nodes = list(set(sibling_nodes + lineage_conflict_nodes))
            
            vec_conflicts = pd.Series(0, index=M_current.index)
            for conflict in all_conflict_nodes:
                conflict_series = M_current[conflict].reindex(M_current.index, fill_value=0)
                vec_conflicts = vec_conflicts | conflict_series
            
            new_mut_cleaned = new_mut_bin_vector & ~vec_conflicts
            vec_children_aligned = vec_children.reindex(new_mut_cleaned.index, fill_value=0)
            imputed_vec = (vec_children_aligned | new_mut_cleaned).astype(int)
            merge_penalty = np.log(len(merge_children)) * 0.5
            N_nodes = N_nodes_beforeT + 2
            
        else:
            raise ValueError(f"Unknown placement_type: {placement_type}")
        
        # ============================================================
        # 计算罚分
        # ============================================================
        
        full_mutnode_chain = get_full_mutnode_chain_with_anchor(anchor, parent_dict)
        
        posterior_vec = P_selected[new_mut]
        input_binary_vec = I_selected[new_mut]
        
        new_mut_penalty, actual_na_flip_ratio, refined_ω_NA, φ_adjusted, weight_na_to_1, weight_na_to_0 = compute_dynamic_penalty(
            input_binary_vec, posterior_vec, imputed_vec, fnfp_ratio, ω_NA, φ,
            na_ratio, mut_ratio, placement_type, N_nodes
        )
        
        chain_penalty = 0
        chain_mutations_count = 0
        
        for node in full_mutnode_chain:
            if node == 'ROOT':
                continue
            
            mutations_on_node = node.split("|")
            
            for mutation in mutations_on_node:
                if mutation == new_mut:
                    continue
                    
                mut_input_binary_vec = I_selected[mutation].replace({pd.NA: np.nan})
                mut_posterior_vec = P_selected[mutation]
                
                mut_new_vec = M_current[node].copy()
                cells_should_be_1 = imputed_vec[imputed_vec == 1].index
                cells_to_flip = []
                for cell in cells_should_be_1:
                    if mut_new_vec[cell] == 0 or pd.isna(mut_new_vec[cell]):
                        cells_to_flip.append(cell)
                        mut_new_vec[cell] = 1
                
                mut_penalty = compute_bayesian_penalty_each_chain_mut_by_pos(
                    mut_input_binary_vec[cells_to_flip], mut_posterior_vec[cells_to_flip], 
                    mut_new_vec[cells_to_flip], weight_na_to_1, weight_na_to_0, fnfp_ratio
                )
                
                chain_penalty += mut_penalty
                chain_mutations_count += 1
        
        total_chain_penalty = new_mut_penalty + chain_penalty
        
        log_N_nodes_penalty = np.log(N_nodes)
        BIC_penalty = φ * np.log(N_nodes)
        root_penalty = 0
        if anchor == 'ROOT':
            root_penalty = np.log(N_nodes) * 0.5
        
        base_total_penalty = total_chain_penalty + log_N_nodes_penalty + BIC_penalty + merge_penalty + root_penalty
        
        intersection_penalty = compute_intersection_based_penalty(
            new_mut, pos, intersection_nodes, M_current, I_selected, na_ratio, mut_ratio, actual_na_flip_ratio
        )
        
        hierarchy_penalty = compute_hierarchy_penalty(
            new_mut, pos, M_current, I_selected, parent_dict, na_ratio, mut_ratio, actual_na_flip_ratio
        )
        
        total_penalty = base_total_penalty + intersection_penalty + hierarchy_penalty
        
        results.append({
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
        })
    
    df_penalty = pd.DataFrame(results)
    return df_penalty


def apply_position_to_tree(
    new_mut, position, imputed_vec, T_current, M_current, I_selected, parent_dict
):
    """
    将特定位置应用到树和矩阵上
    """
    M_updated = M_current.copy()
    T_updated = T_current.copy()
    
    anchor = position['anchor']
    
    full_mutnode_chain = get_full_mutnode_chain_with_anchor(anchor, parent_dict)
    
    if not imputed_vec.index.equals(M_updated.index):
        imputed_vec = imputed_vec.reindex(M_updated.index, fill_value=0)
    
    cells_with_final_one = imputed_vec[imputed_vec == 1].index.tolist()
    if len(cells_with_final_one) > 0:
        for cell in cells_with_final_one:
            for mutation in full_mutnode_chain:
                if M_updated.loc[cell, mutation] == 0:
                    M_updated.loc[cell, mutation] = 1
    
    M_updated[new_mut] = imputed_vec
    T_updated = add_new_mutation_to_tree_independent(new_mut, T_updated, position)
    
    return T_updated, M_updated


def compute_bayesian_penalty_for_positions_consider_ROOT(
    new_mut, selected_positions, T_current, M_current, I_selected, P_selected, parent_dict, intersection_nodes, 
    ω_NA=0.001, fnfp_ratio=0.1, φ=1
):
    import copy
    import pandas as pd
    import numpy as np
    results = []
    
    # 如果没有候选位置，直接返回 None
    if len(selected_positions) == 0:
        logger.warning(f"No selected positions for mutation {new_mut}")
        return None, None, pd.DataFrame(), M_current
    
    new_mut_bin_vector = I_selected[new_mut].replace({pd.NA: np.nan}).fillna(0).astype(int)
        
    # 计算突变的特征用于动态调整惩罚
    input_binary_vec_full = I_selected[new_mut].replace({pd.NA: np.nan})
    na_ratio = input_binary_vec_full.isna().mean()
    mut_ratio = input_binary_vec_full.fillna(0).mean()
    N_nodes_beforeT = len(T_current.all_nodes())
    
    for idx, pos in enumerate(selected_positions):
        placement_type = pos['placement_type']
        anchor = pos['anchor']
        # N_nodes_beforeT = len([node['name'] for node in pos['nodes']])
        
        # 默认 imputed vector
        imputed_vec = pd.Series(0, index=M_current.index)
        merge_penalty = 0
        
        # -------------------------
        # 修正：同时考虑直系sibling和lineage之外的冲突
        # -------------------------
        if placement_type == 'on_node':
            parent = anchor
            vec_parent = M_current[parent] if parent != 'ROOT' else pd.Series(1, index=M_current.index)
            
            # 获取直系sibling冲突节点
            sibling_nodes = [n['name'] for n in pos['nodes'] if n['parent'] == parent and n['name'] != new_mut]
            
            # 获取lineage之外的所有冲突节点
            lineage_conflict_nodes = get_all_conflict_nodes_outside_lineage(parent, build_lineage_parent_dict_from_tree(T_current, anchor), M_current.columns)
            
            # 合并所有冲突节点（去重）
            all_conflict_nodes = list(set(sibling_nodes + lineage_conflict_nodes))
            
            # 构建冲突向量
            vec_conflicts = pd.Series(0, index=M_current.index)
            for conflict in all_conflict_nodes:
                vec_conflicts |= M_current[conflict]
            
            # 正确的逻辑：现有节点的向量与清理后的new_mut合并
            new_mut_cleaned = new_mut_bin_vector & ~vec_conflicts
            imputed_vec = (M_current[anchor] | new_mut_cleaned).astype(int)
            N_nodes = N_nodes_beforeT + 1
            
        elif placement_type == 'new_leaf':
            parent = anchor
            vec_parent = M_current[parent] if parent != 'ROOT' else pd.Series(1, index=M_current.index)
            
            # 获取直系sibling冲突节点
            sibling_nodes = [n['name'] for n in pos['nodes'] if n['parent'] == parent and n['name'] != new_mut]
            
            # 获取lineage之外的所有冲突节点
            lineage_conflict_nodes = get_all_conflict_nodes_outside_lineage(parent, build_lineage_parent_dict_from_tree(T_current, anchor), M_current.columns)
            
            # 合并所有冲突节点（去重）
            all_conflict_nodes = list(set(sibling_nodes + lineage_conflict_nodes))
            
            # 构建冲突向量
            vec_conflicts = pd.Series(0, index=M_current.index)
            for conflict in all_conflict_nodes:
                vec_conflicts |= M_current[conflict]
            
            # 正确的逻辑：先排除所有冲突，再与parent交集
            new_mut_cleaned = new_mut_bin_vector & ~vec_conflicts
            imputed_vec = new_mut_cleaned.astype(int)
            N_nodes = N_nodes_beforeT + 2
            
        elif placement_type == 'on_edge':
            parent = anchor
            child = pos['meta']['child']
            vec_parent = M_current[parent] if parent != 'ROOT' else pd.Series(1, index=M_current.index)
            vec_child = M_current[child]
            
            # 获取直系sibling冲突节点
            sibling_nodes = [n['name'] for n in pos['nodes'] if n['parent']==parent and n['name'] not in [child,new_mut]]
            
            # 获取lineage之外的所有冲突节点
            lineage_conflict_nodes = get_all_conflict_nodes_outside_lineage(parent, build_lineage_parent_dict_from_tree(T_current, anchor), M_current.columns)
            
            # 合并所有冲突节点（去重）
            all_conflict_nodes = list(set(sibling_nodes + lineage_conflict_nodes))
            
            # 构建冲突向量
            vec_conflicts = pd.Series(0, index=M_current.index)
            for conflict in all_conflict_nodes:
                vec_conflicts |= M_current[conflict]
            
            # 正确的逻辑：child ∪ (清理后的new_mut ∩ parent)
            new_mut_cleaned = new_mut_bin_vector & ~vec_conflicts
            imputed_vec = (vec_child | new_mut_cleaned).astype(int)
            N_nodes = N_nodes_beforeT + 2
            
        elif placement_type == 'new_parent_merge':
            parent = anchor
            merge_children = pos['meta']['merge_children']
            vec_parent = M_current[parent] if parent != 'ROOT' else pd.Series(1, index=M_current.index)
            
            # 构建children的联合向量
            vec_children = pd.Series(0, index=M_current.index)
            for c in merge_children:
                vec_children |= M_current[c]
            
            # 获取直系sibling冲突节点
            sibling_nodes = [n['name'] for n in pos['nodes'] if n['parent']==parent and n['name'] not in merge_children+[new_mut]]
            
            # 获取lineage之外的所有冲突节点（排除merge_children）
            lineage_conflict_nodes = get_all_conflict_nodes_outside_lineage(parent, build_lineage_parent_dict_from_tree(T_current, anchor), M_current.columns, exclude_nodes=merge_children)
            
            # 合并所有冲突节点（去重）
            all_conflict_nodes = list(set(sibling_nodes + lineage_conflict_nodes))
            
            # 构建冲突向量
            vec_conflicts = pd.Series(0, index=M_current.index)
            for conflict in all_conflict_nodes:
                vec_conflicts |= M_current[conflict]
            
            # 正确的逻辑：children ∪ (清理后的new_mut ∩ parent)
            new_mut_cleaned = new_mut_bin_vector & ~vec_conflicts
            imputed_vec = (vec_children | new_mut_cleaned).astype(int)
            merge_penalty = np.log(len(merge_children)) * 0.5
            N_nodes = N_nodes_beforeT + 2
            
        else:
            raise ValueError(f"Unknown placement_type: {placement_type}")
        
        # -------------------------
        # 计算完整突变链的总罚分
        # -------------------------
        # 获取从锚点到ROOT的完整突变链
        full_mutnode_chain = get_full_mutnode_chain_with_anchor(anchor, parent_dict)
        
        # 计算new_mut本身的罚分
        posterior_vec = P_selected[new_mut]
        input_binary_vec = I_selected[new_mut]
        
        new_mut_penalty, actual_na_flip_ratio, refined_ω_NA, φ_adjusted, weight_na_to_1, weight_na_to_0 = compute_dynamic_penalty(
            input_binary_vec, posterior_vec, imputed_vec, fnfp_ratio, ω_NA, φ,
            na_ratio, mut_ratio, placement_type, N_nodes
        )
        
        # 计算完整突变链中其他突变的罚分
        chain_penalty = 0
        chain_na_flip_ratio = 0
        chain_mutations_count = 0
        
        for node in full_mutnode_chain:
            if node == 'ROOT':
                continue
            
            mutations_on_node = node.split("|")
            
            for mutation in mutations_on_node:
                if mutation == new_mut:  # 跳过new_mut本身，因为已经计算过了
                    continue
                    
                # 获取该突变的原始数据
                mut_input_binary_vec = I_selected[mutation].replace({pd.NA: np.nan})
                mut_posterior_vec = P_selected[mutation]
                
                # 计算该突变在new_mut放置后的新向量
                # 对于split_full_mutmutation_chain_mutations中的突变，如果它们在final_imputed_vec为1的细胞中当前为0或NA，需要翻转为1
                mut_new_vec = M_current[node].copy()
                cells_should_be_1 = imputed_vec[imputed_vec == 1].index
                # 找出要计算的额外发生翻转的 index
                cells_to_flip = []
                for cell in cells_should_be_1:
                    if mut_new_vec[cell] == 0 or pd.isna(mut_new_vec[cell]):
                        cells_to_flip.append(cell)
                        mut_new_vec[cell] = 1
                
                # 计算该突变的罚分
                mut_penalty = compute_bayesian_penalty_each_chain_mut_by_pos(
                    mut_input_binary_vec[cells_to_flip], mut_posterior_vec[cells_to_flip], mut_new_vec[cells_to_flip], weight_na_to_1, weight_na_to_0, fnfp_ratio
                )
                
                chain_penalty += mut_penalty
                chain_mutations_count += 1
        
        # 总罚分 = new_mut罚分 + 完整突变链罚分
        total_chain_penalty = new_mut_penalty + chain_penalty
        
        log_N_nodes_penalty = np.log(N_nodes)
        # BIC_penalty = φ_adjusted * np.log(N_nodes)  # 使用动态调整后的φ
        BIC_penalty = φ * np.log(N_nodes)  # 使用动态调整后的φ
        root_penalty = 0
        if anchor == 'ROOT':
            root_penalty = np.log(N_nodes) * 0.5
        
        # 基础总罚分
        base_total_penalty = total_chain_penalty + log_N_nodes_penalty + BIC_penalty + merge_penalty + root_penalty
        
        # -------------------------
        # 添加基于intersection模式和层级的精细调整
        # -------------------------
        intersection_penalty = compute_intersection_based_penalty(
            new_mut, pos, intersection_nodes, M_current, I_selected, na_ratio, mut_ratio, actual_na_flip_ratio
        )
        
        hierarchy_penalty = compute_hierarchy_penalty(
            new_mut, pos, M_current, I_selected, parent_dict, na_ratio, mut_ratio, actual_na_flip_ratio
        )
        
        total_penalty = base_total_penalty + intersection_penalty + hierarchy_penalty
        
        results.append({
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
        })
    
    # 初始化df_penalty避免作用域问题
    df_penalty = pd.DataFrame()
    
    try:
        # -------------------------
        # 安全地选出最小 penalty 的位置
        # -------------------------
        if len(results) == 0:
            logger.warning(f"No valid results for mutation {new_mut}")
            return None, None, pd.DataFrame(), M_current
        
        df_penalty = pd.DataFrame(results)
        
        # 检查 DataFrame 是否为空
        if df_penalty.empty:
            logger.warning(f"Empty penalty DataFrame for mutation {new_mut}")
            return None, None, df_penalty, M_current
        
        # 检查 total_penalty 列是否存在
        if 'total_penalty' not in df_penalty.columns:
            logger.warning(f"total_penalty column missing for mutation {new_mut}. Columns: {df_penalty.columns.tolist()}")
            return None, None, df_penalty, M_current
        
        # 检查 total_penalty 列是否有有效值
        if df_penalty['total_penalty'].isna().all():
            logger.warning(f"All total_penalty values are NaN for mutation {new_mut}")
            return None, None, df_penalty, M_current
        
        # 查找最小 penalty 并检查 imputed_vec 是否全为 0
        valid_results = []
        for idx, row in df_penalty.iterrows():
            if row['imputed_vec'].sum() == 0:
                continue
            valid_results.append(row)
        
        if not valid_results:
            logger.warning(f"No valid results (with non-zero imputed_vec) for mutation {new_mut}")
            return None, None, df_penalty, M_current
        
        # 选择最小 penalty 的有效位置
        df_valid_penalty = pd.DataFrame(valid_results)
        min_idx = df_valid_penalty['base_total_penalty'].idxmin()
        min_row = df_valid_penalty.loc[min_idx]
        
        final_position = min_row['position']
        final_imputed_vec = min_row['imputed_vec'].copy()
        
        # -------------------------
        # 修正：简化的M_current更新逻辑
        # -------------------------
        anchor = final_position['anchor']
        
        # 获取从锚点到ROOT的完整突变链
        full_mutnode_chain = get_full_mutnode_chain_with_anchor(anchor, parent_dict)
        
        # 找出所有在final_imputed_vec中为1的细胞
        cells_with_final_one = final_imputed_vec[final_imputed_vec == 1].index.tolist()
        
        if len(cells_with_final_one) > 0:
            # 对于每个在final_imputed_vec中为1的细胞，确保其父节点链也为1
            for cell in cells_with_final_one:
                for mutation in full_mutnode_chain:
                    # 只有当该突变当前为0时才填充为1
                    if M_current.loc[cell, mutation] == 0:
                        M_current.loc[cell, mutation] = 1
        
        return final_position, final_imputed_vec.astype(int), df_penalty, M_current
    
    except Exception as e:
        logger.warning(f"Error selecting minimum penalty for mutation {new_mut}: {e}")
        return None, None, df_penalty, M_current


# 保留原有的辅助函数，因为它们已经在第一个函数中定义过了
def get_all_conflict_nodes_outside_lineage(anchor, parent_dict, all_columns, exclude_nodes=None):
    """
    获取lineage之外的所有可能冲突节点
    
    Parameters:
    - anchor: 当前锚点节点
    - parent_dict: 父节点字典
    - all_columns: 所有突变节点
    - exclude_nodes: 要排除的节点列表（如merge_children）
    
    Returns:
    - conflict_nodes: lineage之外的所有节点列表
    """
    if exclude_nodes is None:
        exclude_nodes = []
    
    # 获取当前lineage的所有节点（从anchor到ROOT的路径）
    lineage_nodes = set()
    current = anchor
    while current is not None:
        lineage_nodes.add(current)
        parent = parent_dict.get(current, None)
        if parent is None or parent == 'ROOT':
            break
        current = parent
    lineage_nodes.add('ROOT')
    
    # 获取所有不在当前lineage中的节点
    conflict_nodes = []
    for node in all_columns:
        if node == 'ROOT':
            continue
        if node not in lineage_nodes and node not in exclude_nodes:
            conflict_nodes.append(node)
    
    return conflict_nodes


def get_full_mutnode_chain_with_anchor(anchor, parent_dict):
    """
    获取从锚点到ROOT的完整突变链（包括锚点本身）
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
    基于实际NA→1翻转比例动态调整所有惩罚权重
    """
    # 计算实际的NA→1翻转比例
    na_mask = input_binary_vec.isna()
    na_to_one_count = ((imputed_vec == 1) & na_mask).sum()
    total_na_count = na_mask.sum()
    
    if total_na_count > 0:
        actual_na_flip_ratio = na_to_one_count / total_na_count
    else:
        actual_na_flip_ratio = 0
    
    # 计算w_fn的值
    w_fn = fnfp_ratio * 1  # 默认w_fn = 0.1
    
    # 更精细的动态调整：让w_na在翻转比例高时接近w_fn
    if actual_na_flip_ratio > 0.8:  # 超过80%的NA被翻转为1
        refined_ω_NA = w_fn * 0.9  # 几乎等于w_fn
    elif actual_na_flip_ratio > 0.6:  # 超过60%的NA被翻转为1
        refined_ω_NA = w_fn * 0.7  # 达到w_fn的70%
    elif actual_na_flip_ratio > 0.4:  # 超过40%的NA被翻转为1
        refined_ω_NA = w_fn * 0.5  # 达到w_fn的50%
    elif actual_na_flip_ratio > 0.2:  # 超过20%的NA被翻转为1
        refined_ω_NA = base_ω_NA * 4.0  # 大幅增加惩罚
    elif actual_na_flip_ratio > 0.1:  # 超过10%的NA被翻转为1
        refined_ω_NA = base_ω_NA * 2.0
    elif actual_na_flip_ratio > 0.05:  # 超过5%的NA被翻转为1
        refined_ω_NA = base_ω_NA * 1.5
    else:
        refined_ω_NA = base_ω_NA  # 翻转很少，保持原惩罚
    
    # 确保不会超过w_fn（因为NA翻转本质上比FN更不可靠）
    refined_ω_NA = min(refined_ω_NA, w_fn * 0.95)
    
    # 根据NA特征和实际翻转比例动态调整BIC惩罚
    φ_adjusted = compute_dynamic_bic_penalty(
        base_φ, na_ratio, mut_ratio, actual_na_flip_ratio, placement_type, N_nodes
    )
    
    # 使用调整后的权重计算惩罚
    penalty, weight_na_to_1, weight_na_to_0 = compute_bayesian_penalty_each_pos(
        input_binary_vec, posterior_vec, imputed_vec, fnfp_ratio, refined_ω_NA
    )
    
    return penalty, actual_na_flip_ratio, refined_ω_NA, φ_adjusted, weight_na_to_1, weight_na_to_0


def compute_dynamic_bic_penalty(base_φ, na_ratio, mut_ratio, actual_na_flip_ratio, placement_type, N_nodes):
    """
    基于NA特征和实际翻转比例动态调整BIC惩罚
    """
    # 基础调整：高NA小突变降低BIC惩罚
    if na_ratio > 0.7 and mut_ratio < 0.1:
        base_adjustment = 0.3
    elif na_ratio > 0.5 and mut_ratio < 0.2:
        base_adjustment = 0.6
    else:
        base_adjustment = 1.0
    
    # 基于实际翻转的进一步调整
    if actual_na_flip_ratio > 0.4:
        # 如果实际翻转比例很高，说明这个放置可能不可靠
        # 对于新节点创建，给予更多宽容（降低惩罚）
        # 对于现有节点放置，增加惩罚（因为可能造成大量不可靠翻转）
        if placement_type in ['new_leaf', 'on_edge']:
            flip_adjustment = 0.8  # 新节点创建，进一步降低惩罚
        else:
            flip_adjustment = 1.2  # 现有节点放置，增加惩罚
    elif actual_na_flip_ratio > 0.2:
        if placement_type in ['new_leaf', 'on_edge']:
            flip_adjustment = 0.9
        else:
            flip_adjustment = 1.1
    else:
        flip_adjustment = 1.0
    
    # 节点数量考虑：节点越多，对新节点创建的惩罚应该相对更宽容
    node_adjustment = 1.0
    if N_nodes > 20 and placement_type in ['new_leaf', 'on_edge']:
        node_adjustment = 0.9
    elif N_nodes > 50 and placement_type in ['new_leaf', 'on_edge']:
        node_adjustment = 0.8
    
    φ_adjusted = base_φ * base_adjustment * flip_adjustment * node_adjustment
    
    return max(φ_adjusted, 0.1)  # 确保不会降得太低


def compute_intersection_based_penalty(new_mut, position, intersection_nodes, M_current, I_selected, na_ratio, mut_ratio, actual_na_flip_ratio):
    """
    根据intersection模式和实际NA翻转调整罚分
    """
    extra_penalty = 0
    placement_type = position['placement_type']
    anchor = position['anchor']
    
    # 获取new_mut的mutant cells
    mut_cells = set(I_selected.index[I_selected[new_mut].fillna(0) == 1])
    
    if placement_type == 'on_node':
        anchor_cells = set(M_current.index[M_current[anchor] == 1])
        
        # 情况1：如果new_mut只与一个节点强相交，但被放在不同的early node
        if len(intersection_nodes) == 1:
            sole_intersection = list(intersection_nodes)[0]
            if anchor != sole_intersection and anchor in [n for n in M_current.columns if n != 'ROOT']:
                # 结合实际翻转比例调整惩罚
                flip_based_multiplier = 1.0 + actual_na_flip_ratio  # 翻转越多，惩罚越重
                extra_penalty += np.log(len(M_current.columns)) * 0.8 * flip_based_multiplier
        
        # 情况2：如果anchor涵盖的细胞远多于new_mut的mutant cells
        if len(anchor_cells) > len(mut_cells) * 3:
            # 结合实际翻转比例调整FN惩罚
            fn_penalty = np.log(len(anchor_cells) - len(mut_cells)) * 0.3
            extra_penalty += fn_penalty * (1.0 + actual_na_flip_ratio)
            
        # 情况3：高NA突变被放在非相交节点上且实际翻转比例高
        if na_ratio > 0.7 and anchor not in intersection_nodes and actual_na_flip_ratio > 0.3:
            extra_penalty += np.log(len(M_current.columns)) * 0.5
    
    elif placement_type in ['new_leaf', 'on_edge']:
        # 对于高NA小突变且实际翻转比例低的，减少开新节点的惩罚
        if na_ratio > 0.7 and mut_ratio < 0.1 and actual_na_flip_ratio < 0.2:
            bonus = -np.log(len(M_current.columns)) * 0.4 * (1.0 - actual_na_flip_ratio)
            extra_penalty += bonus
            
        # 对于与多个节点相交的突变，减少开新节点的惩罚
        if len(intersection_nodes) >= 2:
            bonus = -np.log(len(M_current.columns)) * 0.3
            extra_penalty += bonus
    
    return extra_penalty


def compute_hierarchy_penalty(new_mut, position, M_current, I_selected, parent_dict, 
                            na_ratio, mut_ratio, actual_na_flip_ratio):
    """
    计算层级合理性惩罚，考虑实际NA翻转
    """
    hierarchy_penalty = 0
    placement_type = position['placement_type']
    anchor = position['anchor']
    
    mut_cells = set(I_selected.index[I_selected[new_mut].fillna(0) == 1])
    
    if placement_type == 'on_node':
        anchor_cells = set(M_current.index[M_current[anchor] == 1])
        
        # 检查这个anchor是否有children
        anchor_children = find_children_of_node(anchor, M_current.columns, parent_dict)
        
        if anchor_children:
            # 如果anchor有children，且new_mut与children有特定模式
            children_intersection = False
            for child in anchor_children:
                child_cells = set(M_current.index[M_current[child] == 1])
                if mut_cells.issubset(child_cells) and len(mut_cells) < len(child_cells) * 0.8:
                    children_intersection = True
                    break
            
            if children_intersection:
                # 结合实际翻转比例调整层级不合理惩罚
                flip_multiplier = 1.0 + actual_na_flip_ratio
                hierarchy_penalty += np.log(len(anchor_children) + 1) * 0.4 * flip_multiplier
    
    elif placement_type == 'new_leaf':
        parent = anchor
        parent_cells = set(M_current.index[M_current[parent] == 1])
        
        # 如果new_mut的细胞是parent细胞的真子集，这是一个合理的子克隆
        if mut_cells.issubset(parent_cells) and len(mut_cells) < len(parent_cells) * 0.8:
            # 给予奖励，奖励幅度取决于实际翻转比例（翻转越少奖励越多）
            bonus_multiplier = 1.0 - actual_na_flip_ratio  # 翻转越少，奖励越多
            hierarchy_penalty -= np.log(len(parent_cells) - len(mut_cells)) * 0.2 * bonus_multiplier
            
        # 对于高NA突变且实际翻转比例低的，进一步奖励合理的子克隆放置
        if na_ratio > 0.7 and mut_cells.issubset(parent_cells) and actual_na_flip_ratio < 0.2:
            bonus_multiplier = 1.0 - actual_na_flip_ratio
            hierarchy_penalty -= np.log(len(parent_cells)) * 0.3 * bonus_multiplier
    
    return hierarchy_penalty


def find_children_of_node(node, all_columns, parent_dict):
    """找到节点的直接children"""
    children = []
    for col in all_columns:
        if col == node or col == 'ROOT':
            continue
        if parent_dict.get(col) == node:
            children.append(col)
    return children



# -------------------------
# 计算 fp_ratio 和 fn_ratio
# -------------------------

##### 在整个树上计算
def calculate_fp_fn_ratios_across_tree(M_checkpoint, I_attached):
    """
    计算每个突变的fp_ratio和fn_ratio
    
    Parameters:
    M_checkpoint: imputed之后的基因型矩阵 (binary, 0/1)
    I_attached: 原始输入的基因型矩阵 (0, 1, NaN)
    
    Returns:
    DataFrame: 包含每个突变的fp_ratio和fn_ratio
    dict: 包含每个突变对应的fp>0的其他突变列表
    """
    mutations = M_checkpoint.columns
    results = []
    fp_mutations_dict = {}  # 存储每个突变对应的fp>0的其他突变
    
    for mut in mutations:
        # 1. 计算fp_ratio
        mutant_cells = M_checkpoint[M_checkpoint[mut] == 1].index
        fp_ratios = []
        fp_positive_muts = []  # 存储当前突变fp>0的其他突变
        
        for other_mut in mutations:
            if other_mut == mut:
                continue
                
            # 在当前突变的clone中计算其他突变的FP
            other_original = I_attached.loc[mutant_cells, other_mut]
            other_imputed = M_checkpoint.loc[mutant_cells, other_mut]
            # FP: 原始为1但imputed为0
            fp_count = ((other_original == 1) & (other_imputed == 0)).sum()
            
            # 在所有 cells 中计算其他突变的FP
            total_other_original = I_attached.loc[:, other_mut]
            total_other_imputed = M_checkpoint.loc[:, other_mut]
            total_fp_count = ((total_other_original == 1) & (total_other_imputed == 0)).sum()
            
            if total_fp_count > 0:
                fp_ratio = fp_count / total_fp_count
                fp_ratios.append(fp_ratio)
                
                # 如果fp_ratio > 0，记录这个other_mut
                if fp_ratio > 0:
                    fp_positive_muts.append(other_mut)
            else:
                fp_ratios.append(0)
        
        # 将当前突变的fp>0的其他突变存入字典
        fp_mutations_dict[mut] = fp_positive_muts
        
        avg_fp_ratio = np.mean(fp_ratios) if fp_ratios else 0
        
        # 2. 计算fn_ratio (保持不变)
        mut_original = I_attached.loc[mutant_cells, mut]
        mut_imputed = M_checkpoint.loc[mutant_cells, mut]
        
        # 该突变在原始数据中为0的细胞
        zero_cells = mut_original[mut_original == 0].index
        
        # 计算与其他突变有intersection的细胞数
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

# # 使用
# ratios_df, fp_mutations_dict = calculate_fp_fn_ratios_across_tree(M_checkpoint, I_attached)
# print(ratios_df)
# print(fp_mutations_dict)


##### 在树中的每一个 subclone 内部计算
def calculate_fp_ratios_within_subclone(M_checkpoint, I_attached, mutation_clones_for_subclone):
    """
    在每个subclone内部计算每个突变的fp_ratio，并记录FP>0的其他突变
    
    Parameters:
    M_checkpoint: imputed之后的基因型矩阵 (binary, 0/1)
    I_attached: 原始输入的基因型矩阵 (0, 1, NaN)
    mutation_clones_for_subclone: 字典，键为subclone代表突变，值为该subclone包含的所有突变列表
    
    Returns:
    DataFrame: 包含每个突变在每个subclone中的fp_ratio
    dict: 每个突变对应的FP>0的其他突变列表
    """
    results_for_out_subclone_muts = []
    results_for_in_subclone_muts = []
    # 新增：存储每个突变对应的FP>0的其他突变
    fp_positive_mutations_dict_for_out_subclone_muts = {}
    fp_positive_mutations_dict_for_in_subclone_muts = {}
    
    # 遍历每个subclone
    for subclone_rep, subclone_mutations in mutation_clones_for_subclone.items():
        # 获取当前 subclone 的所有 cells（这些 cells 至少有一个该 subclone 的突变）
        subclone_cells = set()
        for mut in subclone_mutations:
            mutant_cells = M_checkpoint[M_checkpoint[mut] == 1].index
            subclone_cells.update(mutant_cells)
        
        subclone_cells = list(subclone_cells)
        
        # 先算所有突变的 fp count 和 fp ratio
        for mut in subclone_mutations:
            # 获取当前突变的mutant cells
            mutant_cells = M_checkpoint[M_checkpoint[mut] == 1].index
            
            # 获取当前突变在subclone中的其他突变（排除自身）
            other_mutations_out_subclone = [m for m in M_checkpoint.columns if m != mut and m not in subclone_mutations]
            
            fp_ratios = []
            # 新增：存储当前突变对应的FP>0的其他突变
            fp_positive_for_current_mut = []
            
            # 只计算与当前突变在同一subclone中的其他突变的FP
            for other_mut in other_mutations_out_subclone:
                # 在当前突变的clone中计算其他突变的FP
                other_original = I_attached.loc[mutant_cells, other_mut]
                other_imputed = M_checkpoint.loc[mutant_cells, other_mut]
                
                # FP: 原始为1但imputed为0
                fp_count = ((other_original == 1) & (other_imputed == 0)).sum()
                
                # 在当前subclone的所有cells中计算其他突变的FP
                subclone_other_original = I_attached.loc[subclone_cells, other_mut]
                subclone_other_imputed = M_checkpoint.loc[subclone_cells, other_mut]
                subclone_fp_count = ((subclone_other_original == 1) & (subclone_other_imputed == 0)).sum()
                
                if subclone_fp_count > 0:
                    fp_ratio = fp_count / subclone_fp_count
                    fp_ratios.append(fp_ratio)
                    # 新增：如果FP ratio > 0，记录这个other_mut
                    if fp_ratio > 0:
                        fp_positive_for_current_mut.append(other_mut)
                else:
                    fp_ratios.append(0)
            
            # 计算平均fp_ratio
            avg_fp_ratio = np.mean(fp_ratios) if fp_ratios else 0
            
            results_for_out_subclone_muts.append({
                'identifier': mut,
                'subclone_representative': subclone_rep,
                'subclone_size': len(subclone_mutations),
                'fp_ratio_within_subclone_for_out_subclone_muts': avg_fp_ratio,
                'count_fp_mutations_for_out_subclone_muts': len(fp_positive_for_current_mut),
                'ratio_fp_mutations_for_out_subclone_muts': len(fp_positive_for_current_mut)/len(M_checkpoint.columns)
            })
            
            # 新增：将当前突变的FP>0的其他突变存入字典
            fp_positive_mutations_dict_for_out_subclone_muts[mut] = fp_positive_for_current_mut
        
        # 如果subclone只有一个突变，无法计算fp_ratio，设为0
        if len(subclone_mutations) <= 1:
            for mut in subclone_mutations:
                results_for_in_subclone_muts.append({
                    'identifier': mut,
                    'fp_ratio_within_subclone_for_in_subclone_muts': 0,
                    'count_fp_mutations_for_in_subclone_muts': 0,
                    'ratio_fp_mutations_for_in_subclone_muts': 0
                })
                # 新增：对于单个突变的subclone，FP>0的其他突变列表为空
                fp_positive_mutations_dict_for_in_subclone_muts[mut] = []        
            
        else:
            
            # 遍历subclone中的每个突变
            for mut in subclone_mutations:
                # 获取当前突变的mutant cells
                mutant_cells = M_checkpoint[M_checkpoint[mut] == 1].index
                
                # 获取当前突变在subclone中的其他突变（排除自身）
                other_mutations_in_subclone = [m for m in subclone_mutations if m != mut]
                
                fp_ratios = []
                # 新增：存储当前突变对应的FP>0的其他突变
                fp_positive_for_current_mut = []
                
                # 只计算与当前突变在同一subclone中的其他突变的FP
                for other_mut in other_mutations_in_subclone:
                    # 在当前突变的clone中计算其他突变的FP
                    other_original = I_attached.loc[mutant_cells, other_mut]
                    other_imputed = M_checkpoint.loc[mutant_cells, other_mut]
                    
                    # FP: 原始为1但imputed为0
                    fp_count = ((other_original == 1) & (other_imputed == 0)).sum()
                    
                    # 在当前subclone的所有cells中计算其他突变的FP
                    subclone_other_original = I_attached.loc[subclone_cells, other_mut]
                    subclone_other_imputed = M_checkpoint.loc[subclone_cells, other_mut]
                    subclone_fp_count = ((subclone_other_original == 1) & (subclone_other_imputed == 0)).sum()
                    
                    if subclone_fp_count > 0:
                        fp_ratio = fp_count / subclone_fp_count
                        fp_ratios.append(fp_ratio)
                        # 新增：如果FP ratio > 0，记录这个other_mut
                        if fp_ratio > 0:
                            fp_positive_for_current_mut.append(other_mut)
                    else:
                        fp_ratios.append(0)
                
                # 计算平均fp_ratio
                avg_fp_ratio = np.mean(fp_ratios) if fp_ratios else 0
                
                results_for_in_subclone_muts.append({
                    'identifier': mut,
                    'fp_ratio_within_subclone_for_in_subclone_muts': avg_fp_ratio,
                    'count_fp_mutations_for_in_subclone_muts': len(fp_positive_for_current_mut),
                    'ratio_fp_mutations_for_in_subclone_muts': len(fp_positive_for_current_mut)/len(M_checkpoint.columns)
                })
                
                # 新增：将当前突变的FP>0的其他突变存入字典
                fp_positive_mutations_dict_for_in_subclone_muts[mut] = fp_positive_for_current_mut
    
    df_results = pd.merge(
        pd.DataFrame(results_for_out_subclone_muts),
        pd.DataFrame(results_for_in_subclone_muts),
        on="identifier",
        how='outer')
    
    return df_results, fp_positive_mutations_dict_for_out_subclone_muts, fp_positive_mutations_dict_for_in_subclone_muts

# 使用示例
# fp_ratios_df, fp_positive_mutations_dict_for_out_subclone_muts, fp_positive_mutations_dict_for_in_subclone_muts = calculate_fp_ratios_within_subclone(M_checkpoint, I_attached, mutation_clones_for_subclone)
# print(ratios_df)


def calculate_fp_ratios_persite_within_subclone(M_checkpoint, I_attached, mutation_clones_for_subclone):
    """
    在每个subclone内部计算每个突变自身的fp_ratio_persite
    (FP个数除以该突变在subclone中的原始数据所有1的个数)
    
    Parameters:
    M_checkpoint: imputed之后的基因型矩阵 (binary, 0/1)
    I_attached: 原始输入的基因型矩阵 (0, 1, NaN)
    mutation_clones_for_subclone: 字典，键为subclone代表突变，值为该subclone包含的所有突变列表
    
    Returns:
    DataFrame: 包含每个突变在每个subclone中的fp_count_persite和fp_ratio_persite
    """
    results = []
    
    # 遍历每个subclone
    for subclone_rep, subclone_mutations in mutation_clones_for_subclone.items():
        # 获取当前subclone的所有cells（这些cells至少有一个该subclone的突变）
        subclone_cells = set()
        for mut in subclone_mutations:
            mutant_cells = M_checkpoint[M_checkpoint[mut] == 1].index
            subclone_cells.update(mutant_cells)
        subclone_cells = list(subclone_cells)
        
        # 遍历subclone中的每个突变
        for mut in subclone_mutations:
            # 计算当前突变在subclone中的原始数据所有1的个数
            original_data_in_subclone = I_attached.loc[subclone_cells, mut]
            total_ones_in_subclone = (original_data_in_subclone == 1).sum()
            
            # 如果该突变在subclone的原始数据中没有1，则fp_ratio设为0
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
            
            # 计算当前突变自身的FP数量
            # 在所有subclone细胞中：原始为1但imputed为0
            original_mut_in_subclone = I_attached.loc[subclone_cells, mut]
            imputed_mut_in_subclone = M_checkpoint.loc[subclone_cells, mut]
            
            fp_count = ((original_mut_in_subclone == 1) & (imputed_mut_in_subclone == 0)).sum()
            
            # 计算fp_ratio
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
    计算每个突变在所有细胞中的fp_ratio_per_mutation，并输出每个突变对应的fp>0的其他突变列表
    
    Parameters:
    M_checkpoint: imputed之后的基因型矩阵 (binary, 0/1)
    I_attached: 原始输入的基因型矩阵 (0, 1, NaN)
    
    Returns:
    tuple: (DataFrame, dict)
        DataFrame: 包含每个突变的fp_count和fp_ratio_per_mutation
        dict: 包含每个突变对应的fp>0的其他突变列表
    """
    results = []
    fp_mutations_dict = {}  # 存储每个突变对应的fp>0的其他突变列表
    
    # 遍历每个突变
    for mut in I_attached.columns:
        # 计算该突变在所有细胞中的原始数据所有1的个数
        original_data = I_attached[mut]
        total_ones = (original_data == 1).sum()
        total_zeros = (original_data == 0).sum()
        
        # 如果该突变在原始数据中没有1，则fp_ratio设为0
        if total_ones == 0:
            results.append({
                'identifier': mut,
                'fp_cells_count': 0,
                'fp_cells_ratio_per_mutation': 0,
                'total_ones_cells': 0,
                'total_zeros_cells': 0,
                'total_coverage_cells': 0
            })
            fp_mutations_dict[mut] = []  # 空列表
            continue
        
        # 计算该突变的FP数量：原始为1但imputed为0
        fp_count = ((original_data == 1) & (M_checkpoint[mut] == 0)).sum()
        
        # 计算fp_ratio_per_mutation
        fp_ratio = fp_count / total_ones
        
        # 找出该突变FP>0对应的其他突变
        fp_positive_muts = []
        for other_mut in I_attached.columns:
            if other_mut == mut:
                continue
                
            # 找出该突变原始为1但imputed为0的细胞
            fp_cells_mask = (original_data == 1) & (M_checkpoint[mut] == 0)
            if fp_cells_mask.any():
                # 在这些FP细胞中，检查other_mut的情况
                fp_cells = original_data[fp_cells_mask].index
                other_in_fp_cells = I_attached.loc[fp_cells, other_mut]
                
                # 如果other_mut在这些FP细胞中有至少一个1，则记录
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

# # 使用示例
# df_fp_ratio, fp_mutations_dict = calculate_fp_ratio_per_mutation_with_fp_mutations_dict(M_checkpoint, I_attached)

def calculate_fp_ratio_per_cell(M_checkpoint, I_attached):
    """
    计算每个细胞在所有突变中的fp_ratio_per_cell
    (FP个数除以该细胞在所有突变中的原始数据所有1的个数)
    
    Parameters:
    M_checkpoint: imputed之后的基因型矩阵 (binary, 0/1)
    I_attached: 原始输入的基因型矩阵 (0, 1, NaN)
    
    Returns:
    DataFrame: 包含每个细胞的fp_count和fp_ratio_per_cell
    """
    results = []
    
    # 遍历每个细胞
    for cell in I_attached.index:
        # 计算该细胞在所有突变中的原始数据所有1的个数
        original_data = I_attached.loc[cell]
        total_ones = (original_data == 1).sum()
        total_zeros = (original_data == 0).sum()
        
        # 如果该细胞在原始数据中没有1，则fp_ratio设为0
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
        
        # 计算该细胞的FP数量：原始为1但imputed为0
        fp_count = ((original_data == 1) & (M_checkpoint.loc[cell] == 0)).sum()
        
        # 计算fp_ratio_per_cell
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
    综合计算所有FP指标
    
    Parameters:
    M_checkpoint: imputed之后的基因型矩阵 (binary, 0/1)
    I_attached: 原始输入的基因型矩阵 (0, 1, NaN)
    
    Returns:
    tuple: (per_mutation_df, per_cell_df, overall_metrics)
    """
    # 计算每突变的FP比率
    per_mutation_df, fp_mutations_dict = calculate_fp_ratio_per_mutation_with_fp_mutations_dict(M_checkpoint, I_attached)
    
    # 计算每细胞的FP比率
    per_cell_df = calculate_fp_ratio_per_cell(M_checkpoint, I_attached)
    
    # 计算总体指标
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

# # 计算所有指标
# per_mutation_df, per_cell_df, overall_metrics = calculate_comprehensive_fp_metrics(M_full, I_attached)

# 可以单独使用某个函数
# per_mutation_df = calculate_fp_ratio_per_mutation(M_checkpoint, I_attached)
# per_cell_df = calculate_fp_ratio_per_cell(M_checkpoint, I_attached)




def get_all_daughter_mutations(node):
    """最佳实现：获取节点的所有后代突变名称"""
    if not node.children:
        return []
    
    daughter_mutations = []
    for child in node.children:
        # 使用现有的traverse方法，简洁高效
        for descendant in child.traverse():
            daughter_mutations.append(descendant.name)
    
    return daughter_mutations

def find_ordered_branch_groups_for_rehanged_mutations_with_keys_as_earlist(tree_root, mutation_list):
    """
    优化版本：使用DFS找出所有相关突变组，并以最早突变为key返回字典
    """
    target_mutations = set(mutation_list)
    all_groups = []
    
    def dfs(node, current_branch):
        current_branch.append(node.name)
        
        # 检查当前分支上有多少目标突变
        branch_mutations = [m for m in current_branch if m in target_mutations]
        
        # 如果当前分支有至少2个目标突变
        if len(branch_mutations) >= 2:
            # 检查这个组是否已经被包含在其他组中
            is_new_group = True
            for existing_group in all_groups:
                if set(branch_mutations).issubset(set(existing_group)):
                    is_new_group = False
                    break
            
            if is_new_group:
                all_groups.append(branch_mutations.copy())
        
        # 递归搜索子节点
        if hasattr(node, 'children') and node.children:
            for child in node.children:
                dfs(child, current_branch)
        
        # 回溯
        current_branch.pop()
    
    # 从根节点开始DFS
    dfs(tree_root, [])
    
    # 过滤：只保留最大的连续组
    final_groups = []
    for group in all_groups:
        is_maximal = True
        for other_group in all_groups:
            if group != other_group and set(group).issubset(set(other_group)):
                is_maximal = False
                break
        if is_maximal:
            final_groups.append(group)
    
    # 构建字典：key是每个组的第一个突变，value是整个组
    result_dict = {}
    for group in final_groups:
        earliest_mutation = group[0]  # 第一个元素就是最早的突变
        result_dict[earliest_mutation] = group
    
    return result_dict

# # 使用优化版本
# related_mutations_dict = find_ordered_branch_groups_for_rehanged_mutations_with_keys_as_earlist(
#     T_current, 
#     rehanged_mutations_by_fpratio_within_subclone_but_backbone
# )

# print("以最早突变为key的关联突变组字典:")
# for earliest_mutation, group in related_mutations_dict.items():
#     print(f"Key: {earliest_mutation}")
#     print(f"Value: {group}")
#     print()


def find_ordered_branch_groups_for_rehanged_mutations_with_keys_as_earlist_relaxed(tree_root, mutation_list):
    """
    放宽版本：允许单个突变也形成组
    """
    print("=== 放宽版本：允许单个突变 ===")
    print(f"输入的 mutation_list: {mutation_list}")
    
    target_mutations = set(mutation_list)
    all_groups = []
    
    def dfs(node, current_branch):
        current_branch.append(node.name)
        
        # 检查当前分支上有多少目标突变
        branch_target_mutations = [m for m in current_branch if m in target_mutations]
        
        # 放宽条件：只要有至少1个目标突变就记录
        if len(branch_target_mutations) >= 1:
            print(f"当前节点: {node.name}, 分支目标突变: {branch_target_mutations}")
            
            # 检查这个组是否已经被包含在其他组中
            is_new_group = True
            for existing_group in all_groups:
                if set(branch_target_mutations).issubset(set(existing_group)):
                    is_new_group = False
                    break
            
            if is_new_group:
                all_groups.append(branch_target_mutations.copy())
                print(f"  ➕ 添加新组: {branch_target_mutations}")
        
        # 递归搜索子节点
        if hasattr(node, 'children') and node.children:
            for child in node.children:
                dfs(child, current_branch)
        
        # 回溯
        current_branch.pop()
    
    # 从根节点开始DFS
    print("开始DFS遍历...")
    dfs(tree_root, [])
    
    print(f"所有找到的组: {all_groups}")
    
    # 构建字典：key是每个组的第一个突变，value是整个组
    result_dict = {}
    for group in all_groups:
        earliest_mutation = group[0]  # 第一个元素就是最早的突变
        result_dict[earliest_mutation] = group
    
    print(f"最终结果字典: {result_dict}")
    return result_dict




##### 将 TreeNode 和 Matrix 格式的树合理删除一些突变以备重新挂树
def remove_mutations_from_tree_and_matrix(root: TreeNode, M_current: pd.DataFrame, rehanged_mutations: list):
    import copy
    T_removed = root.copy()
    M_removed = M_current.copy()
    
    # 遍历每个节点，处理 rehanged_mutations
    for node in T_removed.all_nodes():
        if node.name == "ROOT":
            continue
        muts = node.name.split("|")
        remaining_muts = [m for m in muts if m not in rehanged_mutations]
        
        if len(remaining_muts) == 0:
            # 节点所有 mutation 都被移除 → 删除节点
            parent = node.parent
            if parent is None:
                raise ValueError("Can not remove ROOT")
            for c in list(node.children):
                parent.add_child(c)
            parent.remove_child(node)
            # 对应矩阵列删掉
            if node.name in M_removed.columns:
                M_removed = M_removed.drop(columns=[node.name])
        else:
            # 节点还有剩余 mutation → 更新 node 名
            new_node_name = "|".join(remaining_muts)
            if node.name != new_node_name:
                # 更新树节点名
                old_name = node.name
                node.name = new_node_name
                # 更新矩阵列名
                if old_name in M_removed.columns:
                    M_removed = M_removed.rename(columns={old_name: new_node_name})
    
    # 最终保证 ROOT 列存在
    if 'ROOT' not in M_removed.columns:
        M_removed.insert(0, 'ROOT', 1)
    
    # 保证列顺序和树节点顺序一致
    final_columns = ['ROOT'] + T_removed.all_names_no_root()
    # 对矩阵列进行排序，如果缺失列则补全全0列
    for col in final_columns:
        if col not in M_removed.columns and col != 'ROOT':
            M_removed[col] = 0
    M_removed = M_removed[final_columns]
        
    return T_removed, M_removed


def calculate_flip_counts_per_site(df1, df2):
    # Make sure both data boxes have the same shape
    if df1.shape != df2.shape:
        raise ValueError("The shapes of the two dataframes do not match.")
    # Get column name
    column_names = df1.columns
    # Define a function to calculate the flip count of a column
    def count_flips(column):
        val1 = df1[column].values
        val2 = df2[column].values
        # Calculate the flip count using vectorization operations
        return compare_elements_vectorized(val1, val2)
    # Apply the count flips function to each column
    flip_counts_results = df1.columns.to_series().apply(count_flips)
    # Convert the result to a data frame
    flip_counts_df = pd.DataFrame(list(flip_counts_results))
    flip_counts_df.columns = ['flipping_False_Negative', 'flipping_False_Positive', 'flipping_NA_to_0', 'flipping_NA_to_1']
    flip_counts_df.index = df1.columns
    return flip_counts_df


# 假设 T_current 是当前树结构，new_mut 是当前突变
def generate_new_leaf_on_root(T_current: TreeNode, new_mut: str):
    # 深拷贝树，避免修改原树
    new_tree = deepcopy(T_current)
    
    # 找到根节点（假设根节点名字是 "ROOT"）
    root_node = new_tree.find("ROOT")
    
    # 创建新的叶子节点
    new_leaf = TreeNode(new_mut)
    
    # 将新的叶子节点添加到根节点的子节点列表中
    root_node.add_child(new_leaf)
    
    # 生成相应的候选位置
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
    基于 intersection 将 external_mutations 聚成若干个子树组，并给出组内合理添加顺序。
    
    参数:
    - I_selected: DataFrame，突变×细胞的矩阵
    - external_mutations: list[str]，待处理突变
    - min_shared: int，两突变之间共同为1的细胞数量 >= min_shared 才认为有 intersection
    
    返回:
    - list[list[str]]: 每个子树（mutation group），组内顺序可直接用于 add_new_mutation_to_tree_independent()
    """
    # 1. 构建 intersection graph
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
    
    # 2. 找出连通子图（每个子图就是一个 subtree）
    subtree_groups = []
    for comp in nx.connected_components(G):
        group = list(comp)
        
        # 3. 在每个组内部建立 directed graph (谁更靠上)
        # 用包含关系或 intersection 大小决定方向
        DG = nx.DiGraph()
        DG.add_nodes_from(group)
        for i, m1 in enumerate(group):
            for m2 in group[i+1:]:
                v1 = I_selected[m1].fillna(0).astype(int)
                v2 = I_selected[m2].fillna(0).astype(int)
                
                # inclusion: 若 v1 是 v2 的超集，则 v1 应在 v2 之上
                if (v1 >= v2).all() and (v1 > v2).any():
                    DG.add_edge(m1, m2)
                elif (v2 >= v1).all() and (v2 > v1).any():
                    DG.add_edge(m2, m1)
                else:
                    # 交叉但非包含，用 intersection size 判断方向
                    inter1 = (v1 & v2).sum()
                    inter2 = inter1  # symmetric
                    if inter1 > 0:
                        # 用总 1 数少的放上面（更稀疏的 mutation 更早出现）
                        if v1.sum() < v2.sum():
                            DG.add_edge(m1, m2)
                        else:
                            DG.add_edge(m2, m1)
        
        try:
            order = list(nx.topological_sort(DG))
        except nx.NetworkXUnfeasible:
            # 如果有循环关系（对称交叉），就按 1 数量从小到大排序
            order = sorted(group, key=lambda m: I_selected[m].sum())
        
        subtree_groups.append(order)
    
    return subtree_groups




##### 根据树 TreeNode 格式的 mutation clone 划分 barcode clone

def get_mutation_clone_and_backbone_mut_as_keys_by_first_level_with_frequency(root: TreeNode, df: pd.DataFrame) -> Dict[str, List[str]]:
    """
    按 ROOT 的第一层级划分 clone，将复合突变完全拆分为独立突变。
    对于复合突变的根，选择在数据框中出现频率最高的突变作为键。
    
    参数:
        root: 系统发育树的根节点
        df: 突变数据框，行是细胞，列是突变，值为0/1/NaN
    
    返回:
        Dict[str, List[str]]: 
            key: 在数据框中出现频率最高的根突变
            value: clone 内包含的所有独立 mutation
    """
    def split_compound_mutations(mutation_name):
        """将复合突变名称拆分为独立突变列表"""
        if '|' in mutation_name:
            return mutation_name.split('|')
        else:
            return [mutation_name]
    
    def get_mutation_frequency(mutation, dataframe):
        """计算突变在数据框中的出现频率（值为1的比例）"""
        if mutation not in dataframe.columns:
            return 0
        # col_data = dataframe[mutation].dropna()  # 去掉NaN值
        col_data = dataframe[mutation]             # 保留NaN值
        if len(col_data) == 0:
            return 0
        return (col_data == 1).sum() / len(col_data)
    
    clone_dict = {}
    for child in root.children:
        # 获取子树中所有节点的名称并拆分复合突变
        all_mutations = []
        for node in child.traverse():
            all_mutations.extend(split_compound_mutations(node.name))
        
        # 处理根节点的复合突变
        root_mutations = split_compound_mutations(child.name)
        
        # 如果根节点只有一个突变，直接使用
        if len(root_mutations) == 1:
            clone_key = root_mutations[0]
        else:
            # 对于复合突变，选择在数据框中出现频率最高的那个
            mutation_frequencies = []
            for mutation in root_mutations:
                freq = get_mutation_frequency(mutation, df)
                mutation_frequencies.append((mutation, freq))
            
            # 按频率排序，选择频率最高的
            mutation_frequencies.sort(key=lambda x: x[1], reverse=True)
            clone_key = mutation_frequencies[0][0]
            
            print(f"复合突变 {child.name} 拆分为: {root_mutations}")
            print(f"频率统计: {mutation_frequencies}")
            print(f"选择作为键: {clone_key}\n")
        
        clone_dict[clone_key] = all_mutations
    
    return clone_dict

def get_node_clone_and_backbone_node_as_keys_by_first_level(root: TreeNode) -> Dict[str, List[str]]:
    """
    按 ROOT 的第一层级划分 clone，直接使用节点名称作为键，不拆分复合突变。
    
    参数:
        root: 系统发育树的根节点
    
    返回:
        Dict[str, List[str]]: 
            key: 第一层级子节点的名称
            value: clone 内包含的所有节点名称
    """
    clone_dict = {}
    for child in root.children:
        # 获取子树中所有节点的名称
        all_nodes = []
        for node in child.traverse():
            all_nodes.append(node.name)
        
        # 直接使用子节点名称作为键
        clone_key = child.name
        clone_dict[clone_key] = all_nodes
    
    return clone_dict

def get_mutation_clone_and_backbone_node_as_keys_by_first_level(root: TreeNode) -> Dict[str, List[str]]:
    """
    按 ROOT 的第一层级划分 clone，使用节点名称作为键，但值中的突变拆分为独立突变。
    """
    clone_dict = {}
    for child in root.children:
        # 获取子树中所有节点的名称并拆分复合突变
        all_mutations = []
        for node in child.traverse():
            if '|' in node.name:
                all_mutations.extend(node.name.split('|'))
            else:
                all_mutations.append(node.name)
        
        # 使用原始节点名称作为键
        clone_dict[child.name] = all_mutations
    
    return clone_dict


def calculate_intersection_counts_under_backbone_nodes(
    mutation_list_under_backbone_nodes, M_current, I_attached, new_mut
):
    """
    计算new_mut与每个backbone node下mutation list的共现数量
    """
    if new_mut not in I_attached.columns:
        raise ValueError(f"突变 {new_mut} 不在 I_attached 数据框中")
    
    intersection_counts = {}
    
    for backbone_node, mutation_list in mutation_list_under_backbone_nodes.items():
        # ===== 关键修复：找到backbone_node在M_current中的实际列名 =====
        actual_col = find_column_in_merged_columns(M_current, backbone_node)
        
        if actual_col is None:
            # 如果找不到，尝试用backbone_node本身
            actual_col = backbone_node
        
        # 检查列是否存在于M_current中
        if actual_col not in M_current.columns:
            print(f"Warning: column '{actual_col}' not found in M_current, skipping...")
            intersection_counts[backbone_node] = 0
            continue
        
        # 1. 找到当前backbone node为1的细胞
        backbone_cells = M_current[M_current[actual_col] == 1].index
        
        # 2. 过滤出同时在I_attached中的细胞
        common_cells = backbone_cells.intersection(I_attached.index)
        
        if len(common_cells) == 0:
            intersection_counts[backbone_node] = 0
            continue
        
        # 3. 获取这些细胞在new_mut上的状态
        new_mut_status = I_attached.loc[common_cells, new_mut]
        
        # 4. 过滤出new_mut为1的细胞
        new_mut_positive_cells = new_mut_status[new_mut_status == 1].index
        
        if len(new_mut_positive_cells) == 0:
            intersection_counts[backbone_node] = 0
            continue
        
        # 5. 对于每个mutation in mutation_list，计算与new_mut同时为1的数量
        total_intersection = 0
        
        for mutation in mutation_list:
            if mutation in I_attached.columns:
                mut_status = I_attached.loc[new_mut_positive_cells, mutation]
                intersection_count = (mut_status == 1).sum()
                total_intersection += intersection_count
        
        intersection_counts[backbone_node] = total_intersection
    
    return intersection_counts


# # 使用示例
# new_mut = 'chr8_91087030_T_G'
# result = calculate_intersection_counts_under_backbone_nodes(mutation_list_under_backbone_nodes, M_current, I_attached, new_mut)
# print(result)

def find_best_backbone_for_new_mutation(
    mutation_list_under_backbone_nodes, M_current, I_attached, new_mut
):
    """
    完整函数：计算共现数量并找到最佳backbone node
    """
    # 第一步：计算共现数量
    intersection_counts = calculate_intersection_counts_under_backbone_nodes(
        mutation_list_under_backbone_nodes, M_current, I_attached, new_mut
    )
    
    # 第二步：找到最佳backbone node
    best_backbone = find_best_backbone_node(
        mutation_list_under_backbone_nodes, M_current, intersection_counts
    )
    
    return best_backbone, intersection_counts


def find_best_backbone_node(mutation_list_under_backbone_nodes, M_current, intersection_counts_under_backbone_node):
    """
    简化版本：只返回最佳backbone node
    """
    if not intersection_counts_under_backbone_node:
        return None
    
    max_count = max(intersection_counts_under_backbone_node.values())
    max_nodes = [node for node, count in intersection_counts_under_backbone_node.items() 
                if count == max_count]
    
    if len(max_nodes) == 1:
        return max_nodes[0]
    
    # 平局处理：找实际存在的列
    best_node = None
    best_normalized = -1
    
    for node in max_nodes:
        # ===== 关键修复：找到节点在M_current中的实际列名 =====
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

# # 使用完整函数
# best_backbone, intersection_counts = find_best_backbone_for_new_mutation(
#     mutation_list_under_backbone_nodes, M_current, I_attached, new_mut
# )
# print(f"最佳backbone node: {best_backbone}")
# print(f"所有backbone的共现数量: {intersection_counts}")


def assign_clone_labels(M_full: pd.DataFrame, mutation_clones: dict) -> pd.DataFrame:
    """
    为每个 barcode 根据它包含的 mutation 分配 clone label。
    返回一个包含 barcode 和 clone label 的新的 DataFrame。
    """
    # 初始化一个空的数据框，用于存储 label 和 color
    result_df = pd.DataFrame(columns=["label", "color", "backbone_mutation"])
    
    # 遍历每个 mutation group（即每个 clone）
    for clone_idx, (mutation_group, mutations) in enumerate(mutation_clones.items(), start=1):
        # 为包含这个 mutation group 的 barcodes 分配 label 和 color
        for barcode in M_full.index:
            # 如果该 barcode 包含当前 group 中的任何一个 mutation
            if M_full.loc[barcode, mutations].sum() > 0:
                # 使用 pd.concat 替代 append
                new_row = pd.DataFrame({"label": [barcode], "color": [f'C{clone_idx}'], "backbone_mutation": mutation_group})
                result_df = pd.concat([result_df, new_row], ignore_index=True)
    
    return result_df


def get_first_level_backbone_nodes(root: TreeNode) -> List[str]:
    """
    获取第一层级的所有backbone nodes
    """
    first_level_nodes = []
    for child in root.children:
        first_level_nodes.append(child.name)
    return first_level_nodes




# -------------------------
# 计算 parent muts 中的 intersection/FN_flip per mutations，最后找到一些要重挂的突变
# -------------------------

import pandas as pd
import numpy as np

def find_mutation_column(mutation, columns):
    """在M_current的列中查找包含指定突变的列"""
    for col in columns:
        if '|' in col:
            if mutation in col.split('|'):
                return col
        else:
            if mutation == col:
                return col
    return None

def get_all_parents_for_mutation(T_current, mutation):
    """获取某个突变的所有父突变（从直接父突变到ROOT路径上的所有祖先）"""
    parent_dict = build_lineage_parent_dict_from_tree(T_current, mutation)
    
    all_parents = []
    current = mutation
    
    # 沿着父链向上追溯，直到ROOT
    while current in parent_dict:
        parent = parent_dict[current]
        all_parents.append(parent)
        current = parent
    
    return all_parents

def calculate_intersection_and_inter_vs_fn_flipping_ratio_per_mutation(T_current, M_current, I_attached):
    """正确版本：每个突变只返回一行统计结果"""
    
    results = []
    mutations = I_attached.columns.tolist()
    
    # 对齐两个DataFrame的索引
    common_index = I_attached.index.intersection(M_current.index)
    I_aligned = I_attached.loc[common_index]
    M_aligned = M_current.loc[common_index]
    
    for mutation in mutations:
        # 找到当前突变在 M_current 中的对应列
        m_col = find_mutation_column(mutation, M_aligned.columns)
        
        # 获取所有父突变
        all_parents = get_all_parents_for_mutation(T_current, m_col)
        
        # 2. 如果找到了节点，获取它的所有直接子节点
        all_children = []
        if T_current.find(m_col):
            children = T_current.find(m_col).children
            for child in children:
                all_children.append(child.name)
        else:
            print(f"未找到节点: {m_col}")
        
        if m_col is None:
            print(f"警告: 在 M_current 中找不到突变 {mutation}")
            # 如果找不到，填充NaN
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
        
        # 找出保持为1的cells
        retained_mask = (I_aligned[mutation] == 1) & (M_aligned[m_col] == 1)
        retained_cells = I_aligned.index[retained_mask]
        total_retained = len(retained_cells)
        
        if len(retained_cells) == 0 or not all_parents:
            # 没有保持的cells或没有父突变
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
        
        # 计算所有父突变的总体统计
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
        intersection_cell_ratio_on_mutation = len(unique_intersection_cells_on_mutation_parents)/total_retained
        # 计算总体比例：intersection / (intersection + FN)
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
        return pd.DataFrame(results)
    else:
        return pd.DataFrame()


def process_matrices_by_removed_some_mutations_from_tree(M_current, I_attached):
    """
    处理两个矩阵：
    1. 从 I_attached 中提取 M_current 中对应的列（ROOT 除外）
    2. 在 I_attached 提取的数据中，将全为0或NaN的行在 M_current 中对应的cells设为0
    
    参数:
    M_current: DataFrame, 当前矩阵
    I_attached: DataFrame, 附加矩阵
    
    返回:
    tuple: (I_attached_removed_outgroup, M_current_modified)
    """
    # 第一步：从 I_attached 中提取 M_current 中对应的列（ROOT 除外）
    m_columns = [col for col in M_current.columns if col != "ROOT"]
    split_columns = [col.split("|")[0] if "|" in col else col for col in m_columns]
    
    # 从 I_attached 中提取对应的列
    I_attached_removed_outgroup = I_attached[split_columns]
    
    # 第二步：找出在 I_attached_removed_outgroup 中全为 0 或 NaN 的行
    rows_to_zero = I_attached_removed_outgroup[
        (I_attached_removed_outgroup == 0) | (I_attached_removed_outgroup.isna())
    ].all(axis=1)
    
    # 获取需要修改的行的索引
    zero_rows_indices = rows_to_zero[rows_to_zero].index
    
    # 在 M_current 中将对应的行（除了 ROOT 列）设为 0
    M_current_modified = M_current.copy()
    M_current_modified.loc[zero_rows_indices, M_current_modified.columns != "ROOT"] = 0
    
    return I_attached_removed_outgroup, M_current_modified


def process_conflicting_cells_stay_outgroup(M_current_split_and_noROOT, I_attached_only_outgroup):
    """
    处理两个数据框中存在冲突的cells，返回两个处理后的数据框
    
    参数:
    M_current_split_and_noROOT: 第一个数据框
    I_attached_only_outgroup: 第二个数据框
    
    返回:
    M_processed: 处理后的第一个数据框
    I_processed: 处理后的第二个数据框
    """
    
    # 复制原始数据框，避免修改原数据
    M_processed = M_current_split_and_noROOT.copy()
    I_processed = I_attached_only_outgroup.copy()
    
    # 确保两个数据框的行索引相同
    common_cells = M_processed.index.intersection(I_processed.index)
    M_processed = M_processed.loc[common_cells]
    I_processed = I_processed.loc[common_cells]
    
    # 找到在两个数据框中都有1的cells
    conflicting_cells = []
    
    for cell in common_cells:
        # 获取该cell在两个数据框中的行
        M_row = M_processed.loc[cell]
        I_row = I_processed.loc[cell]
        
        # 检查在两个数据框中是否都有至少一个1
        # 注意：I_attached_only_outgroup中可能有NaN值，需要处理
        M_has_ones = (M_row == 1).any()
        I_has_ones = (I_row.fillna(0) == 1).any()  # 将NaN视为0
        
        if M_has_ones and I_has_ones:
            conflicting_cells.append(cell)
    
    print(f"找到 {len(conflicting_cells)} 个在两个数据框中都有1的cells")
    
    # 处理每个存在冲突的cell
    for cell in conflicting_cells:
        M_row = M_processed.loc[cell]
        I_row = I_processed.loc[cell]
        
        # 计算1的数量（对于I数据框，将NaN视为0）
        M_ones_count = (M_row == 1).sum()
        I_ones_count = (I_row.fillna(0) == 1).sum()
        
        # print(f"Cell {cell}: M中有{M_ones_count}个1, I中有{I_ones_count}个1")
        
        # 根据规则处理
        if M_ones_count > I_ones_count:
            # 情况1: M中1更多 - 在I中将该cell整行变为0
            I_processed.loc[cell] = 0
            # print(f"  -> M中1更多，保持M不变，将I中该行变为0")
            
        elif I_ones_count > M_ones_count:
            # 情况2: I中1更多 - 在M中将该cell整行变为0
            M_processed.loc[cell] = 0
            # print(f"  -> I中1更多，将M中该行变为0，保持I不变")
            
        else:
            # 情况3: 1的数量相等 - 在M中将该cell整行变为0，保持I不变
            M_processed.loc[cell] = 0
            # print(f"  -> 1的数量相等，将M中该行变为0，保持I不变")
    
    return M_processed, I_processed

# # 使用函数
# M_current_split_and_noROOT_processed, I_attached_only_outgroup_processed = process_conflicting_cells(
#     M_current_split_and_noROOT, 
#     I_attached_only_outgroup
# )

# print("处理后的M数据框:")
# print(M_current_split_and_noROOT_processed)
# print("\n处理后的I数据框:")
# print(I_attached_only_outgroup_processed)


def process_conflicting_cells_stay_maintree(M_current_split_and_noROOT, I_attached_only_outgroup):
    """
    处理两个数据框中存在冲突的cells
    
    参数:
    M_current_split_and_noROOT: 第一个数据框
    I_attached_only_outgroup: 第二个数据框
    
    返回:
    M_current_split_and_noROOT_processed_unconflict: 处理后的第一个数据框
    I_attached_only_outgroup_processed: 处理后的第二个数据框
    """
    
    # 复制原始数据框，避免修改原数据
    M_processed = M_current_split_and_noROOT.copy()
    I_processed = I_attached_only_outgroup.copy()
    
    # 1. 找出在两个数据框中都有至少一个1的行
    # 对于M数据框，找出有至少一个1的行
    M_has_ones = (M_processed == 1).any(axis=1)
    # 对于I数据框，找出有至少一个1的行（注意处理NaN值）
    I_has_ones = (I_processed == 1).fillna(False)
    I_has_ones = I_has_ones.any(axis=1)
    
    # 两个数据框都有1的行
    conflicting_cells = M_has_ones & I_has_ones
    conflicting_indices = conflicting_cells[conflicting_cells].index
    
    print(f"找到 {len(conflicting_indices)} 个存在冲突的cells")
    
    # 2. 处理每个冲突的cell
    for cell in conflicting_indices:
        # 计算每个数据框中1的数量
        M_ones_count = (M_processed.loc[cell] == 1).sum()
        I_ones_count = (I_processed.loc[cell] == 1).sum()
        
        print(f"Cell {cell}: M中有{M_ones_count}个1, I中有{I_ones_count}个1")
        
        # 根据规则处理
        if M_ones_count > I_ones_count:
            # M中1多：I中该行所有1变成0
            I_row = I_processed.loc[cell]
            I_row[I_row == 1] = 0
            I_processed.loc[cell] = I_row
            print(f"  -> M中1多：保持M不变，I中该行1变为0")
            
        elif I_ones_count > M_ones_count:
            # I中1多：M中该行所有1变成0
            M_row = M_processed.loc[cell]
            M_row[M_row == 1] = 0
            M_processed.loc[cell] = M_row
            print(f"  -> I中1多：M中该行1变为0，保持I不变")
            
        else:
            # 1的数量相等：M保持不变，I中该行所有1变成0
            I_row = I_processed.loc[cell]
            I_row[I_row == 1] = 0
            I_processed.loc[cell] = I_row
            print(f"  -> 1的数量相等：保持M不变，I中该行1变为0")
    
    return M_processed, I_processed

# # 使用示例
# M_current_split_and_noROOT_processed, I_attached_only_outgroup_processed = process_conflicting_cells(
#     M_current_split_and_noROOT, I_attached_only_outgroup
# )


def process_conflicting_cells_allto_outgroup(M_current_split_and_noROOT, I_attached_only_outgroup):
    """
    高效版本的处理函数
    """
    
    M_processed = M_current_split_and_noROOT.copy()
    I_processed = I_attached_only_outgroup.copy()
    
    # 找出在两个数据框中都有至少一个1的行
    M_has_ones = (M_processed == 1).any(axis=1)
    I_has_ones = (I_processed == 1).fillna(False).any(axis=1)
    conflicting_indices = (M_has_ones & I_has_ones)
    
    print(f"找到 {conflicting_indices.sum()} 个在两个数据框中都存在1的cells")
    
    # 在M数据框中将这些冲突的行全部设为0
    M_processed.loc[conflicting_indices] = 0
    
    return M_processed, I_processed

# # 使用示例
# M_current_split_and_noROOT_processed, I_attached_only_outgroup_processed = process_conflicting_cells_allto_outgroup(
#     M_current_split_and_noROOT, I_attached_only_outgroup
# )


# 查看冲突cells的详细信息
def show_conflict_details(I_attached_removed_outgroup, I_attached_only_outgroup, n_examples=3):
    """显示冲突cells的详细信息"""
    
    common_cells = I_attached_removed_outgroup.index.intersection(I_attached_only_outgroup.index)
    
    # 将数据转换为整数类型
    removed_int = I_attached_removed_outgroup.loc[common_cells].fillna(0).astype(int)
    only_int = I_attached_only_outgroup.loc[common_cells].fillna(0).astype(int)
    
    # 找出在两个数据框中都至少有一个1的cells
    has_1_in_removed = (removed_int.sum(axis=1) > 0)
    has_1_in_only = (only_int.sum(axis=1) > 0)
    conflict_cells = has_1_in_removed & has_1_in_only
    
    conflict_cells_index = conflict_cells[conflict_cells].index
    
    print(f"\n=== 冲突cells详细信息 ===")
    print(f"总共发现 {len(conflict_cells_index)} 个冲突cells")
    
    if len(conflict_cells_index) > 0:
        # 计算1的数量
        count_removed = removed_int.loc[conflict_cells_index].sum(axis=1)
        count_only = only_int.loc[conflict_cells_index].sum(axis=1)
        
        # 显示前几个示例
        for i, cell in enumerate(conflict_cells_index[:n_examples]):
            print(f"\n冲突cell {i+1}: {cell}")
            print(f"  在第一个数据框中的1数量: {count_removed[cell]}")
            print(f"  在第二个数据框中的1数量: {count_only[cell]}")
            print(f"  第一个数据框中的1位置: {removed_int.loc[cell][removed_int.loc[cell] == 1].index.tolist()}")
            print(f"  第二个数据框中的1位置: {only_int.loc[cell][only_int.loc[cell] == 1].index.tolist()}")

# # 显示冲突cells的详细信息
# show_conflict_details(M_current_split_and_noROOT, I_attached_only_outgroup, n_examples=54)


def create_column_mapping(source_columns, target_columns):
    """
    创建从源列名到目标列名的映射
    """
    column_mapping = {}
    
    # 创建目标列名的规范化版本（按突变排序）
    target_normalized = {}
    for col in target_columns:
        if '|' in col:
            mutations = sorted(col.split('|'))
            normalized = '|'.join(mutations)
            target_normalized[normalized] = col
        else:
            target_normalized[col] = col
    
    # 映射源列名到目标列名
    for col in source_columns:
        if '|' in col:
            # 对于合并列，规范化后查找对应的目标列名
            mutations = sorted(col.split('|'))
            normalized = '|'.join(mutations)
            if normalized in target_normalized:
                column_mapping[col] = target_normalized[normalized]
            else:
                # 如果找不到对应，保持原列名
                column_mapping[col] = col
        else:
            # 对于非合并列，直接查找
            if col in target_normalized:
                column_mapping[col] = target_normalized[col]
            else:
                column_mapping[col] = col
    
    return column_mapping

# # 使用示例
# source_cols = merge_duplicate_columns(M_current_split_and_noROOT_processed).columns
# target_cols = M_current_noROOT.columns

# column_mapping = create_column_mapping(source_cols, target_cols)

# # 重命名列
# M_current_split_renamed = merge_duplicate_columns(M_current_split_and_noROOT_processed).rename(columns=column_mapping)

# # 现在两个数据框的列名应该一致了
# print("M_current_noROOT columns:", M_current_noROOT.columns.tolist()[:5])
# print("Renamed columns:", M_current_split_renamed.columns.tolist()[:5])



# -------------------------
# 计算 parent muts 中的 intersection/FN_flip per mutations，最后找到一些找不到合适的突变的可以删除的 cells
# -------------------------

def calculate_intersection_and_flipping_to_1_count_per_cell(M_for_fp_ratio_and_fn_ratio_duplicatecells, I_attached):
    """
    计算每个cell的两个值：
    1. 在两个数据框中同时为1的突变数量（intersection）
    2. 从I_attached到M_current从0或NA变成1的突变数量（flipping to 1）
    
    参数:
    M_for_fp_ratio_and_fn_ratio_duplicatecells: 当前subclone的突变矩阵 (1/0)
    I_attached: 参考突变矩阵 (可能包含0, 1, NaN)
    
    返回:
    DataFrame: 包含每个cell的两个计数值
    """
    
    # 确保两个DataFrame的索引和列对齐
    common_cells = M_for_fp_ratio_and_fn_ratio_duplicatecells.index.intersection(I_attached.index)
    common_mutations = M_for_fp_ratio_and_fn_ratio_duplicatecells.columns.intersection(I_attached.columns)
    
    # 对齐数据
    M_aligned = M_for_fp_ratio_and_fn_ratio_duplicatecells.loc[common_cells, common_mutations]
    I_aligned = I_attached.loc[common_cells, common_mutations]
    
    results = []
    
    for cell in common_cells:
        # 获取当前cell在两个数据框中的行数据
        M_row = M_aligned.loc[cell]
        I_row = I_aligned.loc[cell]
        
        # 1. 计算intersection: 两个数据框中同时为1的突变数量
        intersection_count = ((M_row == 1) & (I_row == 1)).sum()
        
        # 2. 计算flipping to 1: 从I_attached(0或NA)到M_current(1)的突变数量
        # 条件：I_attached中是0或NaN，且M_current中是1
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
    
    # 创建结果DataFrame
    result_df = pd.DataFrame(results).set_index('cell')
    
    return result_df

# # 使用示例
# result = calculate_intersection_and_flipping_to_1_count_per_cell(M_for_fp_ratio_and_fn_ratio_duplicatecells, I_attached)
# print(result)











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
    is_log_value_for_likelihoods=True
):
    """
    Run DP-based posterior computation and mutation filtering through scaffold tree.
    
    Parameters
    ----------
    data : dict
        包含 P, ll_mut, ll_unmut, df_reads, features 等关键输入。
    M_scaffold : pd.DataFrame
        树的 scaffold matrix (binary 0/1)。
    outputpath_full : str
        输出路径。
    scaffold_mutations : list
        初始 scaffold 的突变列表。
    p_thresh : float
        二值化 posterior 的阈值。
    pass_tree_cutoff : float
        “通过树” 突变的后验阈值。
    unpass_tree_cutoff : float
        “未通过树” 突变的后验阈值。
    is_log_value_for_likelihoods : bool
        是否将似然视为 log 值进行归一化。
    
    Returns
    -------
    df_combined : pd.DataFrame
        整合后的特征表，包含 phylogeny_label 和 flipping_counts。
    attached_mutations : list
        可挂载到树上的突变列表（successful_pass + cell_specific）。
    """
    
    # Step 0: 对齐数据索引
    orderd_index = M_scaffold.index.tolist()
    in_posterior = data['P'].reindex(index=orderd_index)
    in_llmut = data['ll_mut'].reindex(index=orderd_index)
    in_llunmut = data['ll_unmut'].reindex(index=orderd_index)
    in_reads = data['df_reads'].reindex(index=[['bulk'] + orderd_index])
    
    # Step 1: normalize likelihoods
    combined_df = pd.concat([in_llmut, in_llunmut], axis=1)
    normalized_result = apply_normalization(combined_df, is_log=is_log_value_for_likelihoods)
    normalized_llmut = normalized_result.filter(like='norm_llmut')
    normalized_llmut.columns = in_llmut.columns
    normalized_llunmut = normalized_result.filter(like='norm_llunmut')
    normalized_llunmut.columns = in_llunmut.columns
    
    # Step 2: calculate posterior
    df_newSP_out, withoutTree_posterior = all_newSomaticPosterior(
        normalized_llmut, normalized_llunmut, M_scaffold.values.astype(int)
    )
    
    # Step 3: features整合
    out_features = df_features_new.T.drop(['somatic_posterior_persite'], axis=1)
    out_features['somatic_posterior_per_site_old'] = df_features_new.T['somatic_posterior_persite']
    out_features['withoutTree_posterior'] = withoutTree_posterior
    out_features = pd.concat([out_features, df_newSP_out], axis=1)
    
    # Step 4: 添加 phylogeny label
    for col in ['somatic_posterior_per_site', 'somatic_posterior_per_site_onecell']:
        out_features[col] = out_features[col].astype(float)
    out_features['mutant_cellnum'] = out_features['mutant_cellnum'].astype(int)
    out_features['phylogeny_label'] = out_features.apply(
        determine_phylogeny_label_by_one_likelihood,
        axis=1,
        args=(pass_tree_cutoff, unpass_tree_cutoff,),
    )
    
    # Step 5: flipping counts
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
    
    # Step 6: 合并
    df_combined = pd.concat([out_features, df_flip_counts], axis=1)
    df_combined['cluster'] = df_combined['cluster'].apply(lambda x: ','.join(map(str, x)))
    df_combined['nonnan_indices'] = df_combined['nonnan_indices'].apply(lambda x: ','.join(map(str, x)))
    
    output_file = os.path.join(outputpath_full, "out_features.somatic_posterior_basedTree.txt")
    df_combined.to_csv(output_file, sep="\t")
    
    # Step 7: mutation分类
    passtree_mutations = list(df_combined[df_combined['phylogeny_label'] == 'successful_pass'].index)
    onecell_mutations = list(df_combined[df_combined['phylogeny_label'] == 'cell_specific'].index)
    
    return df_combined, passtree_mutations, onecell_mutations



##### integrate_mutations_to_scaffold
def attach_mutations_to_current_tree(
    sorted_attached_mutations, T_current, M_current, I_attached, P_attached, 
    ω_NA, fnfp_ratio, φ, logger, root_mutations=None, max_retries=None
):
    """
    处理外部突变并将其整合到进化树中（支持回滚和重试）
    
    Parameters:
    -----------
    sorted_attached_mutations : list
        排序的待处理突变列表
    T_current : dict
        当前进化树结构
    M_current : pandas.DataFrame
        当前突变矩阵
    I_attached : 
        附加的突变信息
    P_attached :
        附加的概率信息
    ω_NA : float
        NA值的权重参数
    fnfp_ratio : float
        假阴性假阳性比率
    φ : float
        贝叶斯罚分参数
    logger : logging.Logger
        日志记录器
    root_mutations : list, optional
        根突变列表
    max_retries : int, optional
        最大尝试候选位置数，None表示尝试全部
    
    Returns:
    --------
    tuple : (external_mutations, conflict_mutations, T_current, M_current, root_mutations)
    """
    
    if root_mutations is None:
        root_mutations = []
    
    external_mutations = []
    conflict_mutations = []
    
    for new_mut in tqdm(sorted_attached_mutations, desc="Processing mutations", unit="mutation"):
        
        # 确定 new_mut 应该属于哪一个 backbone clone
        mutation_list_under_backbone_nodes = get_mutation_clone_and_backbone_node_as_keys_by_first_level(T_current)
        node_list_under_backbone_nodes = get_node_clone_and_backbone_node_as_keys_by_first_level(T_current)
        best_backbone, intersection_counts = find_best_backbone_for_new_mutation(
            mutation_list_under_backbone_nodes, M_current, I_attached, new_mut
        )
        assigned_nodes = node_list_under_backbone_nodes[best_backbone]
        
        # 找到交集节点
        intersection_nodes = find_all_intersect_muts_from_tree_by_matrix(T_current, I_attached, new_mut)
        if len(intersection_nodes) == 0:
            external_mutations.append(new_mut)
            continue
        
        # 使用优化方法获取候选位置
        potential_positions = find_intersection_positions_within_tree_directly(
            T_current, new_mut, I_attached, min_overlap=1
        )
        parent_dict = build_parent_dict_from_candidates(potential_positions)
        selected_positions = [p for p in potential_positions if p['anchor'] in assigned_nodes]
        
        # 检查是否找到候选位置
        if len(selected_positions) == 0:
            selected_positions = [p for p in potential_positions if p['anchor'] != 'ROOT']
        
        # ---- 备份当前状态 ----
        M_backup = M_current.copy()
        T_backup = T_current.copy()
        
        # ---- 计算所有候选位置的罚分 ----
        df_penalty = compute_bayesian_penalty_for_all_positions_consider_ROOT(
            new_mut, selected_positions, T_current, M_current, I_attached, P_attached, 
            parent_dict, intersection_nodes, ω_NA=ω_NA, fnfp_ratio=fnfp_ratio, φ=φ
        )
        
        if df_penalty.empty:
            external_mutations.append(new_mut)
            logger.warning(f"Mutation {new_mut} added to external_mutations (no valid penalty scores)")
            continue
        
        # 过滤掉 imputed_vec 全为 0 的位置
        df_valid = df_penalty[df_penalty['imputed_vec'].apply(lambda x: x.sum() > 0)]
        if df_valid.empty:
            external_mutations.append(new_mut)
            logger.warning(f"Mutation {new_mut} added to external_mutations (all imputed_vec are zero)")
            continue
        
        df_sorted = df_valid.sort_values('total_penalty')
        
        # 决定尝试的候选数
        if max_retries is None:
            candidates_to_try = df_sorted
        else:
            candidates_to_try = df_sorted.head(max_retries)
        
        # ---- 循环尝试候选位置 ----
        placed = False
        for attempt, (idx, row) in enumerate(candidates_to_try.iterrows()):
            
            # 恢复备份
            M_current = M_backup.copy()
            T_current = T_backup.copy()
            
            try:
                # 应用该位置到树和矩阵
                T_current, M_current = apply_position_to_tree(
                    new_mut, row['position'], row['imputed_vec'], 
                    T_current, M_current, I_attached, parent_dict
                )
                
                # 检查冲突
                if scp.ul.is_conflict_free_gusfield(M_current):
                    placed = True
                    break
                else:
                    logger.warning(f"✗ Position {row['position_index']} caused conflict, trying next candidate")
                    
            except Exception as e:
                logger.error(f"Error placing mutation at position {row['position_index']}: {e}")
                continue
        
        # ---- 处理结果 ----
        if not placed:
            M_current = M_backup.copy()
            T_current = T_backup.copy()
            conflict_mutations.append(new_mut)
            logger.warning(f"Mutation {new_mut} added to conflict_mutations (all {len(candidates_to_try)} candidates failed)")
    
    
    return external_mutations, conflict_mutations, T_current, M_current, root_mutations


def process_rescue_mutations(
    sorted_rescue_mutations, T_current, M_current, I_attached, P_attached, 
    mutation_clones_rescue, ω_NA, fnfp_ratio, φ, logger, root_mutations=None, max_retries=None
):
    """
    处理救援突变并将其整合到进化树中，包含克隆亲和性分析（支持回滚和重试）
    
    Returns:
    --------
    tuple : (external_mutations, conflict_mutations, T_current, M_current, root_mutations)
    """
    
    if root_mutations is None:
        root_mutations = []
    
    external_mutations = []
    conflict_mutations = []
    
    for new_mut in tqdm(sorted_rescue_mutations, desc="Processing rescue mutations", unit="mutation"):
        
        # 找到交集节点
        intersection_nodes = find_all_intersect_muts_from_tree_by_matrix(T_current, I_attached, new_mut)
        if len(intersection_nodes) == 0:
            external_mutations.append(new_mut)
            continue
        
        # 使用优化方法获取候选位置
        potential_positions = find_intersection_positions_within_tree_directly(
            T_current, new_mut, I_attached, min_overlap=1
        )
        parent_dict = build_parent_dict_from_candidates(potential_positions)
        
        if len(potential_positions) == 0:
            external_mutations.append(new_mut)
            continue
        
        # 基于 intersection 的情况先选出应该放在哪一个 clone 下
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
        
        # 根据克隆亲和性筛选候选位置
        selected_positions = [
            position for position in potential_positions 
            if position['anchor'] in assigned_clone_muts
        ]
        
        # ---- 备份当前状态 ----
        M_backup = M_current.copy()
        T_backup = T_current.copy()
        
        # ---- 计算所有候选位置的罚分 ----
        df_penalty = compute_bayesian_penalty_for_all_positions_consider_ROOT(
            new_mut, selected_positions, T_current, M_current, I_attached, P_attached, 
            parent_dict, intersection_nodes, ω_NA=ω_NA, fnfp_ratio=fnfp_ratio, φ=φ
        )
        
        if df_penalty.empty:
            external_mutations.append(new_mut)
            logger.warning(f"Mutation {new_mut} added to external_mutations (no valid penalty scores)")
            continue
        
        df_valid = df_penalty[df_penalty['imputed_vec'].apply(lambda x: x.sum() > 0)]
        if df_valid.empty:
            external_mutations.append(new_mut)
            logger.warning(f"Mutation {new_mut} added to external_mutations (all imputed_vec are zero)")
            continue
        
        df_sorted = df_valid.sort_values('total_penalty')
        
        if max_retries is None:
            candidates_to_try = df_sorted
        else:
            candidates_to_try = df_sorted.head(max_retries)
        
        # ---- 循环尝试候选位置 ----
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
                    logger.warning(f"✗ Position {row['position_index']} caused conflict, trying next candidate")
                    
            except Exception as e:
                logger.error(f"Error placing mutation at position {row['position_index']}: {e}")
                continue
        
        if placed:
            # print_tree(T_current)
            logger.info(f"Updated tree! Current M_current is conflict-free, shape: {M_current.shape}")
        else:
            M_current = M_backup.copy()
            T_current = T_backup.copy()
            conflict_mutations.append(new_mut)
            logger.warning(f"Mutation {new_mut} added to conflict_mutations (all candidates failed)")
    
    logger.info(f"Processing complete. External: {len(external_mutations)}, Conflict: {len(conflict_mutations)}")
    
    return external_mutations, conflict_mutations, T_current, M_current, root_mutations




###### 全部走完依旧未处理的再加到 ROOT 的新节点中
def process_external_mutations_by_subtree_groups(
    subtree_groups, T_current, M_current, I_attached, P_attached, 
    ω_NA, fnfp_ratio, φ, logger, root_mutations=None, max_retries=None
):
    """
    通过子树组处理外部突变，支持多突变组和单突变组的分别处理（支持回滚和重试）
    
    Parameters:
    -----------
    subtree_groups : list of lists
        子树组列表，每个组包含相关的突变
    T_current : dict
        当前进化树结构
    M_current : pandas.DataFrame
        当前突变矩阵
    I_attached : 
        附加的突变信息矩阵
    P_attached :
        附加的概率信息矩阵
    ω_NA : float
        NA值的权重参数
    fnfp_ratio : float
        假阴性假阳性比率
    φ : float
        贝叶斯罚分参数
    logger : logging.Logger
        日志记录器
    root_mutations : list, optional
        根突变列表，如果为None则自动创建
    max_retries : int, optional
        最大尝试候选位置数，None表示尝试全部
    
    Returns:
    --------
    tuple : (remained_mutations, conflict_mutations, T_current, M_current, root_mutations)
        剩余未处理的突变列表、冲突突变列表、更新后的树、更新后的矩阵、根突变列表
    """
    
    if root_mutations is None:
        root_mutations = []
    
    remained_mutations = []
    conflict_mutations = []
    
    # 分离子树组和单元素组
    multi_mut_subtree_groups = [g for g in subtree_groups if len(g) > 1]
    singleton_subtree_groups = [g for g in subtree_groups if len(g) == 1]
    
    
    ##### 处理长度 >1 的子树组
    for group_idx, group in enumerate(tqdm(multi_mut_subtree_groups, desc="Processing multiple subtrees")):
        
        # 根据 I_attached 中每个 mutation 的 1 的个数排序（降序）
        sorted_group = sorted(group, key=lambda subtree_mut: I_attached[subtree_mut].sum(), reverse=True)
        
        # 按顺序一个一个加到树上，先挂到 ROOT 下
        reattached_mutations = []
        
        for idx, subtree_mut in enumerate(tqdm(sorted_group, desc="Processing mutations in group")):
            T_rollback = copy.deepcopy(T_current)
            M_rollback = M_current.copy()
            
            
            if idx == 0:
                # 第一个 mutation 直接挂到 ROOT
                final_position = generate_new_leaf_on_root(T_current, subtree_mut)
                T_current = add_new_mutation_to_tree_independent(subtree_mut, T_current, final_position)
                M_current[subtree_mut] = I_attached[subtree_mut].fillna(0).astype(int)
                
                # 检查冲突
                if not scp.ul.is_conflict_free_gusfield(M_current):
                    logger.warning(f"First mutation {subtree_mut} caused conflict, rolling back")
                    T_current = copy.deepcopy(T_rollback)
                    M_current = M_rollback.copy()
                    reattached_mutations.append(subtree_mut)
                    continue
            
            else:
                # 确定 subtree_mut 应该属于哪一个 backbone clone
                mutation_list_under_backbone_nodes = get_mutation_clone_and_backbone_node_as_keys_by_first_level(T_current)
                node_list_under_backbone_nodes = get_node_clone_and_backbone_node_as_keys_by_first_level(T_current)
                best_backbone, intersection_counts = find_best_backbone_for_new_mutation(
                    mutation_list_under_backbone_nodes, M_current, I_attached, subtree_mut
                )
                assigned_nodes = node_list_under_backbone_nodes[best_backbone]
                
                # 找到交集节点
                intersection_nodes = find_all_intersect_muts_from_tree_by_matrix(T_current, I_attached, subtree_mut)
                if len(intersection_nodes) == 0:
                    reattached_mutations.append(subtree_mut)
                    continue
                
                # 使用优化方法获取候选位置
                potential_positions = find_intersection_positions_within_tree_directly(
                    T_current, subtree_mut, I_attached, min_overlap=1
                )
                parent_dict = build_parent_dict_from_candidates(potential_positions)
                selected_positions = [p for p in potential_positions if p['anchor'] in assigned_nodes]
                
                # 检查是否找到候选位置
                if len(selected_positions) == 0:
                    reattached_mutations.append(subtree_mut)
                    continue
                
                # ---- 备份当前状态 ----
                M_backup = M_current.copy()
                T_backup = T_current.copy()
                
                # ---- 计算所有候选位置的罚分 ----
                df_penalty = compute_bayesian_penalty_for_all_positions_consider_ROOT(
                    subtree_mut, selected_positions, T_current, M_current, I_attached, P_attached, 
                    parent_dict, intersection_nodes, ω_NA=ω_NA, fnfp_ratio=fnfp_ratio, φ=φ
                )
                
                if df_penalty.empty:
                    reattached_mutations.append(subtree_mut)
                    logger.warning(f"Mutation {subtree_mut} added to reattached_mutations (no valid penalty scores)")
                    continue
                
                df_valid = df_penalty[df_penalty['imputed_vec'].apply(lambda x: x.sum() > 0)]
                if df_valid.empty:
                    reattached_mutations.append(subtree_mut)
                    logger.warning(f"Mutation {subtree_mut} added to reattached_mutations (all imputed_vec are zero)")
                    continue
                
                df_sorted = df_valid.sort_values('total_penalty')
                
                if max_retries is None:
                    candidates_to_try = df_sorted
                else:
                    candidates_to_try = df_sorted.head(max_retries)
                
                # ---- 循环尝试候选位置 ----
                placed = False
                for attempt, (idx_row, row) in enumerate(candidates_to_try.iterrows()):
                    
                    M_current = M_backup.copy()
                    T_current = T_backup.copy()
                    
                    try:
                        T_current, M_current = apply_position_to_tree(
                            subtree_mut, row['position'], row['imputed_vec'], 
                            T_current, M_current, I_attached, parent_dict
                        )
                        
                        if scp.ul.is_conflict_free_gusfield(M_current):
                            placed = True
                            break
                        else:
                            logger.warning(f"✗ Position {row['position_index']} caused conflict, trying next candidate")
                            
                    except Exception as e:
                        logger.error(f"Error placing mutation at position {row['position_index']}: {e}")
                        continue
                
                if not placed:
                    M_current = M_backup.copy()
                    T_current = T_backup.copy()
                    reattached_mutations.append(subtree_mut)
                    logger.warning(f"Mutation {subtree_mut} added to reattached_mutations (all candidates failed)")
        
        # 处理重新挂载的突变（第一轮未处理的）
        second_reattached_mutations = []
        if len(reattached_mutations) > 0:            
            
            sorted_reattached_mutations = [i for i in I_attached.columns if i in reattached_mutations]
            for subtree_mut in tqdm(sorted_reattached_mutations, desc="Processing re-attached mutations"):
                T_rollback = copy.deepcopy(T_current)
                M_rollback = M_current.copy()
                
                
                # 确定 subtree_mut 应该属于哪一个 backbone clone
                mutation_list_under_backbone_nodes = get_mutation_clone_and_backbone_node_as_keys_by_first_level(T_current)
                node_list_under_backbone_nodes = get_node_clone_and_backbone_node_as_keys_by_first_level(T_current)
                best_backbone, intersection_counts = find_best_backbone_for_new_mutation(
                    mutation_list_under_backbone_nodes, M_current, I_attached, subtree_mut
                )
                assigned_nodes = node_list_under_backbone_nodes[best_backbone]
                
                # 找到交集节点
                intersection_nodes = find_all_intersect_muts_from_tree_by_matrix(T_current, I_attached, subtree_mut)
                if len(intersection_nodes) == 0:
                    second_reattached_mutations.append(subtree_mut)
                    continue
                
                # 使用优化方法获取候选位置
                potential_positions = find_intersection_positions_within_tree_directly(
                    T_current, subtree_mut, I_attached, min_overlap=1
                )
                parent_dict = build_parent_dict_from_candidates(potential_positions)
                selected_positions = [p for p in potential_positions if p['anchor'] in assigned_nodes]
                
                if len(selected_positions) == 0:
                    second_reattached_mutations.append(subtree_mut)
                    continue
                
                # ---- 备份当前状态 ----
                M_backup = M_current.copy()
                T_backup = T_current.copy()
                
                # ---- 计算所有候选位置的罚分 ----
                df_penalty = compute_bayesian_penalty_for_all_positions_consider_ROOT(
                    subtree_mut, selected_positions, T_current, M_current, I_attached, P_attached, 
                    parent_dict, intersection_nodes, ω_NA=ω_NA, fnfp_ratio=fnfp_ratio, φ=φ
                )
                
                if df_penalty.empty:
                    second_reattached_mutations.append(subtree_mut)
                    logger.warning(f"Mutation {subtree_mut} added to second_reattached_mutations (no valid penalty scores)")
                    continue
                
                df_valid = df_penalty[df_penalty['imputed_vec'].apply(lambda x: x.sum() > 0)]
                if df_valid.empty:
                    second_reattached_mutations.append(subtree_mut)
                    logger.warning(f"Mutation {subtree_mut} added to second_reattached_mutations (all imputed_vec are zero)")
                    continue
                
                df_sorted = df_valid.sort_values('total_penalty')
                
                if max_retries is None:
                    candidates_to_try = df_sorted
                else:
                    candidates_to_try = df_sorted.head(max_retries)
                
                # ---- 循环尝试候选位置 ----
                placed = False
                for attempt, (idx_row, row) in enumerate(candidates_to_try.iterrows()):
                    
                    M_current = M_backup.copy()
                    T_current = T_backup.copy()
                    
                    try:
                        T_current, M_current = apply_position_to_tree(
                            subtree_mut, row['position'], row['imputed_vec'], 
                            T_current, M_current, I_attached, parent_dict
                        )
                        
                        if scp.ul.is_conflict_free_gusfield(M_current):
                            logger.info(f"✓ Mutation {subtree_mut} successfully placed (score={row['total_penalty']:.4f})")
                            placed = True
                            break
                        else:
                            logger.warning(f"✗ Position {row['position_index']} caused conflict, trying next candidate")
                            
                    except Exception as e:
                        logger.error(f"Error placing mutation at position {row['position_index']}: {e}")
                        continue
                
                if not placed:
                    M_current = M_backup.copy()
                    T_current = T_backup.copy()
                    second_reattached_mutations.append(subtree_mut)
                    logger.warning(f"Mutation {subtree_mut} added to second_reattached_mutations (all candidates failed)")
        
        # 记录仍然未处理的突变
        if second_reattached_mutations:
            logger.warning(f"Group {group_idx+1} has {len(second_reattached_mutations)} mutations still remaining: {second_reattached_mutations}")
            remained_mutations.extend(second_reattached_mutations)
    
    ##### 处理长度 =1 的单元素组
    for group in tqdm(singleton_subtree_groups, desc="Processing singleton subtrees"):
        T_rollback = copy.deepcopy(T_current)
        M_rollback = M_current.copy()            
        
        subtree_mut = group[0]
        
        final_position = generate_new_leaf_on_root(T_current, subtree_mut)
        T_current = add_new_mutation_to_tree_independent(subtree_mut, T_current, final_position)
        M_current[subtree_mut] = I_attached[subtree_mut].fillna(0).astype(int)
        
        # print_tree(T_current)
        
        if not scp.ul.is_conflict_free_gusfield(M_current):
            logger.warning(f"Conflict detected after adding {subtree_mut}, rolling back")
            
            # 回滚操作：从矩阵中移除这个突变
            T_current = copy.deepcopy(T_rollback)
            M_current = M_rollback.copy()
            
            # 把这个突变放到remained_mutations
            remained_mutations.append(subtree_mut)
            continue
    
    
    return remained_mutations, conflict_mutations, T_current, M_current, root_mutations





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




def get_random_chromosome_position(snp_name):
    """
    根据snp_name解析或生成随机的染色体和位置信息
    """
    parts = snp_name.split('_')
    
    # 如果格式正确，包含下划线分割的部分
    if len(parts) >= 2:
        chromosome = parts[0]
        position = parts[1]
        return chromosome, position
    else:
        # 生成随机的染色体和位置
        chromosome = str(random.randint(1, 22))  # 1-22号染色体
        position = str(random.randint(100000, 999999))  # 6位数字位置
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
    
    logger.info("Applying initial filtration (Section 3.1)")
    
    # Cell filtering: keep cells with at least one supported mutation (posterior > threshold)
    cells_with_V = (V > 0).any(axis=1)
    cells_with_P = (P > p_thresh).any(axis=1)
    kept_cells_mask = cells_with_V & cells_with_P
    kept_cells = P.index[kept_cells_mask].tolist()
    
    if len(kept_cells) == 0:
        raise ValueError("No cells passed initial filtration")
        
    logger.info(f"Kept {len(kept_cells)} cells after initial filtration")
    
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
        
    logger.info(f"Kept {len(kept_mutations)} mutations after initial filtration")
    
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

def filter_scaffold_muts_by_na_proportion(filtered_sites, df_reads, df_celltype, na_prop_thresh=0.9):
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
    Returns
    -------
    scaffold_mutations : list
        High-confidence shared mutations for building scaffold phylogeny
    NA_prop : pd.DataFrame
        NA proportion (mutation × cell type)
    """
    # --- 1. 去掉 bulk 行 ---
    reads = df_reads.drop(index='bulk', errors='ignore')
    # --- 2. 获取所有 cell type ---
    cell_types = df_celltype['cell_type'].unique()
    # --- 3. 构建 coverage matrix（有覆盖记 1，无覆盖或 NA 记 0） ---
    def has_coverage(val):
        if pd.isna(val):
            return 0
        try:
            _, total = val.split('/')
            return 1 if int(total) > 0 else 0
        except:
            return 0
    coverage_matrix = reads.applymap(has_coverage)
    # --- 4. 计算 NA proportion ---
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
    # --- 5. 计算 cutoff ---
    if na_prop_thresh is not None:
        theta = pd.Series(na_prop_thresh, index=cell_types)
    else:
        theta = df_NA_prop.quantile(0.75, axis=0)
    # --- 6. 判断 informative ---
    informative = df_NA_prop.lt(theta, axis=1)
    # --- 7. 筛选 high-confidence scaffold mutations ---
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
    print("Step1 (coverage-based) mutations:", len(scaffold_mutations))
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


def coverage_filters(kept_mutations, df_reads, df_celltype, params, outputpath):
    """
    高置信 scaffold mutation 过滤（基于 coverage）
    Step1: cross-cell-type coverage (NA proportion)
    Step2: CV filter (with median safety)
    两步并行作用于原始输入，最终结果取并集
    
    Parameters
    ----------
    kept_mutations : list
        候选突变
    df_reads : pd.DataFrame
        Rows = cells (第一行为 'bulk'), columns = mutations, values = 'mut/total' 或 NaN
    df_celltype : pd.DataFrame
        DataFrame 包含 'barcode' 和 'cell_type' 列
    params : dict
        包含 'na_prop_thresh_global' 和 'cv_thresh'
    outputpath : str, optional
        保存 summary csv 的路径
    
    Returns
    -------
    final_scaffold_mutations : list
        高置信 shared scaffold mutations
    summary_df : pd.DataFrame
        每个 mutation 的 median / CV / mean / std / pass_filter
    df_NA_prop : pd.DataFrame
        Step1 的 NA proportion
    """
    import numpy as np
    import pandas as pd
    
    if params is None:
        params = DEFAULT_PARAMS
    
    na_prop_thresh = params["na_prop_thresh_global"]
    cv_thresh = params["cv_thresh"]
    
    logger.info("Applying coverage-based filtration (Section 3.2)")
    
    # --- Step1: coverage-based filter ---
    step1_mutations, df_NA_prop = filter_scaffold_muts_by_na_proportion(
        kept_mutations, df_reads, df_celltype, na_prop_thresh
    )
    logger.info("Section 3.2.1) Selection of ubiquitously expressed regions across cell (types)")
    print("=====> Step1 (coverage-based) mutations:", len(step1_mutations))
    
    # --- Step2: CV filter ---
    # Convert the reads to total read counts, treating NA as 0
    reads_matrix_withoutNAcells = df_reads.drop(index='bulk', errors='ignore').applymap(get_total_reads_withoutNAcells)
    reads_matrix_withNAcells = df_reads.drop(index='bulk', errors='ignore').applymap(get_total_reads_withNAcells)
    reads_matrix_withNAcells = reads_matrix_withNAcells.applymap(lambda v: 1 if (not pd.isna(v) and v >= 1) else (0 if not pd.isna(v) else np.nan))
    step2_mutations = []
    median_dict = {}
    cv_dict = {}
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
        if cv <= cv_thresh:
            step2_mutations.append(mut)
    
    logger.info("Section 3.2.2) Selection of regions with relatively uniform read coverage")
    print("=====> Step2 (CV filter) mutations:", len(step2_mutations))
    
    # --- Union ---
    final_scaffold_mutations = list(set(step1_mutations) | set(step2_mutations))
    print("=====> Final scaffold mutations (union):", len(final_scaffold_mutations))
    
    # --- Summary ---
    df_cv_stats = pd.DataFrame({
        "cov_median": pd.Series(median_dict),
        "cov_CV": pd.Series(cv_dict),
        "cov_mean": pd.Series(mean_dict),
        "cov_std": pd.Series(std_dict)
    })
    df_cv_stats["pass_CV"] = df_cv_stats.index.isin(step2_mutations)
    df_cv_stats["pass_NA"] = df_cv_stats.index.isin(step1_mutations)
    df_cv_stats["pass_cov"] = df_cv_stats.index.isin(final_scaffold_mutations)
    
    df_summary = pd.concat([df_cv_stats, df_NA_prop], axis=1)
    
    # --- Save ---
    if outputpath is not None:
        os.makedirs(outputpath, exist_ok=True)
        df_summary.to_csv(os.path.join(outputpath, "Summary_df_in_scaffold_filtration.csv"))
    return final_scaffold_mutations, df_summary


def calculate_cv_for_subgrouping(df_reads_resolved, mutations, cv_threshold=6, logger=None):
    """
    Calculate CV values for mutations to filter low-CV mutations for subgrouping.
    
    Parameters:
    -----------
    df_reads : pd.DataFrame
        Reads dataframe containing read count information
    mutations : list
        List of mutations to calculate CV for
    cv_thresh : float
        CV threshold for filtering (default: 0.5)
    logger : logging.Logger
        Logger instance for logging messages
        
    Returns:
    --------
    Tuple containing:
        - low_cv_mutations: List of mutations with CV <= threshold
        - cv_stats: DataFrame with CV statistics for all mutations
    """
    
    # Convert the reads to total read counts, treating NA as 0
    reads_matrix_withoutNAcells = df_reads_resolved.drop(index='bulk', errors='ignore').applymap(get_total_reads_withoutNAcells)
    reads_matrix_withNAcells = df_reads_resolved.drop(index='bulk', errors='ignore').applymap(get_total_reads_withNAcells)
    reads_matrix_withNAcells = reads_matrix_withNAcells.applymap(lambda v: 1 if (not pd.isna(v) and v >= 1) else (0 if not pd.isna(v) else np.nan))
    low_cv_mutations = []
    median_dict = {}
    cv_dict = {}
    mean_dict = {}
    std_dict = {}
    for mut in mutations:
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
        if cv <= cv_threshold:
            low_cv_mutations.append(mut)
    
    print("=====> In this subgroup, low_cv_mutations:", len(low_cv_mutations))
    
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
            shuffled = list(np.random.permutation(other_muts))
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

def compute_clone_and_pair_weights(muts, corr_cache, n_shuffle=100):
    """
    Parameters
    ----------
    muts : list of str
        所有 mutation ID
    corr_cache : dict
        {(mut1, mut2): True/False}  两两 mutation 是否 correlated
    n_shuffle : int
        每个 mutation 的 shuffle 次数
    
    Returns
    -------
    clone_weights : dict
        {tuple(mut_ids): weight}  每个 clone 的全局权重，包括单节点 clone
    pair_weights : dict
        {tuple(m1,m2): weight}  每个 mutation pair 的权重（只考虑长度≥2的 clone）
    """
    clone_weights = defaultdict(float)  # 全局 clone 权重累加
    
    for ref in muts:
        other_muts = [m for m in muts if m != ref]
        ref_clone_counter = defaultdict(int)  # 记录当前 reference 下每个 clone 的计数
        
        for _ in range(n_shuffle):
            shuffled = list(np.random.permutation(other_muts))
            remaining = [ref] + shuffled.copy()
            
            while remaining:
                curr_ref = remaining[0]
                current_clone = [curr_ref]
                next_remaining = []
                
                for m in remaining[1:]:
                    key1 = (curr_ref, m)
                    key2 = (m, curr_ref)
                    # 检查 corr_cache，避免 KeyError
                    is_corr = corr_cache.get(key1, corr_cache.get(key2, False))
                    if is_corr:
                        current_clone.append(m)
                    else:
                        next_remaining.append(m)
                
                # 当前 shuffle 的 clone 计数
                ref_clone_counter[tuple(sorted(current_clone))] += 1
                remaining = next_remaining
        
        # Step: normalize by n_shuffle → clone proportion for this reference
        for clone, count in ref_clone_counter.items():
            clone_weights[clone] += count / n_shuffle  # 累加到全局
    
    # Step: 计算 pair 权重，只考虑长度≥2的 clone
    pair_weights = defaultdict(float)
    for clone, weight in clone_weights.items():
        if len(clone) > 1:
            for m1, m2 in itertools.combinations(sorted(clone), 2):
                pair_weights[(m1, m2)] += weight
    
    return clone_weights, pair_weights

from typing import List, Tuple, Dict
from collections import defaultdict
import itertools
import numpy as np
import pandas as pd

def get_correlation_graph_elements(I_S: pd.DataFrame, n_shuffle: int = 100, seed: int = 42, cutoff_mcf_for_graph: float = 0.05, cutoff_mcn_for_graph: int = 5) -> Tuple[Dict[Tuple[str], float], Dict[Tuple[str,str], float]]:
    """
    Compute clone_weights and pair_weights, filtering out low-support singleton mutations.
    
    Args:
        I_S: binary matrix (cells x mutations), values 0/1/NA
        n_shuffle: number of random permutations per reference mutation
        seed: random seed
        min_frac: minimum mutant cell fraction for singleton mutations
        min_cells: minimum mutant cell number for singleton mutations
    
    Returns:
        clone_weights: dict mapping clone tuples to weight
        pair_weights: dict mapping mutation pairs to weight
    """
    np.random.seed(seed)
    muts = list(I_S.columns)
    n_mut = len(muts)
    
    # Step 1: precompute pairwise correlation cache
    corr_cache = {}
    for u, v in itertools.combinations(muts, 2):
        corr = are_mutations_correlated(I_S, u, v)
        corr_cache[(u, v)] = corr
        corr_cache[(v, u)] = corr
    
    for m in muts:
        corr_cache[(m, m)] = True
    
    # Step 2: compute clone weights and pair weights
    clone_weights, pair_weights = compute_clone_and_pair_weights(muts, corr_cache, n_shuffle=n_shuffle)
    
    # Step 3: 计算每个 mutation 的突变 fraction 和突变细胞数
    mutant_cell_fraction = {mut: I_S[mut].mean(skipna=True) for mut in muts}
    mutant_cell_number = {mut: I_S[mut].sum(skipna=True) for mut in muts}
    
    # Step 4: 删除低支持 singleton mutations 对应的 clone
    count = 0
    for mut in muts:
        if clone_weights.get((mut,), 0) == n_mut:  # singleton clone
            frac = mutant_cell_fraction.get(mut, 0)
            num = mutant_cell_number.get(mut, 0)
            if frac <= cutoff_mcf_for_graph or num <= cutoff_mcn_for_graph:  # scDNA 要再重新定义, 这里可能就不卡了
                count +=1
                print(f"Filter out singleton low-support mutation: {mut}, frac={frac:.3f}, num={num}")
                clone_weights.pop((mut,), None)  # 删除这个 clone
    
    print(f"The number of filtered singleton, low-support mutations is: {count}")
    return clone_weights, pair_weights

import matplotlib.pyplot as plt
import networkx as nx
def plot_mutation_graph(G_ig, mutation_group, pdf_file, figsize=(8,8), edge_scale=0.2, seed=42):
    """
    可视化 mutation graph，节点颜色表示群组，边宽表示权重。
    
    Parameters
    ----------
    G_ig : igraph Graph
        已构建的 igraph 图
    mutation_group : dict
        {mutation_id: group_id} 每个 mutation 对应的群组
    figsize : tuple
        图像大小
    edge_scale : float
        边权重放大系数
    seed : int
        布局随机种子
    """
    # 1. 转换为 NetworkX 图
    G_nx = nx.Graph()
    for v in G_ig.vs:
        G_nx.add_node(v['name'])
    for e in G_ig.es:
        m1 = G_ig.vs[e.source]['name']
        m2 = G_ig.vs[e.target]['name']
        G_nx.add_edge(m1, m2, weight=float(e['weight']))
    
    # 2. 节点颜色
    groups = [mutation_group[n] for n in G_nx.nodes()]
    unique_groups = list(set(groups))
    color_map = plt.cm.get_cmap('tab20', len(unique_groups))
    node_colors = [color_map(g) for g in groups]
    
    # 3. 边宽
    edges = G_nx.edges()
    edge_weights = [G_nx[u][v]['weight'] for u,v in edges]
    edge_widths = [w*edge_scale for w in edge_weights]
    
    # 4. 布局
    pos = nx.spring_layout(G_nx, seed=seed, k=0.5)  # k 控制节点间距
    
    # 5. 绘制
    plt.figure(figsize=figsize)
    nx.draw_networkx_nodes(G_nx, pos, node_color=node_colors, node_size=200)  # 节点小一点
    nx.draw_networkx_edges(G_nx, pos, width=edge_widths, alpha=0.7)
    nx.draw_networkx_labels(G_nx, pos, font_size=10, font_color='black')
    plt.title("Mutation Graph with Leiden Groups")
    plt.axis('off')
    plt.margins(x=0.2, y=0.2)        # 给四周加 margin
    plt.tight_layout(pad=2.0)        # 额外空白
    plt.savefig(pdf_file, dpi=300)
    plt.close()


import igraph as ig
import leidenalg
def leiden_mutation_groups(clone_weights, pair_weights, pdf_file, resolution=1, seed=42):
    """
    根据 clone_weights 和 pair_weights 构建加权共现图，并使用 Leiden 算法划分 mutation group。
    
    Parameters
    ----------
    clone_weights : dict
        {tuple(mutations): weight}  每个 clone 的全局权重，包括单节点 clone
    pair_weights : dict
        {tuple(m1,m2): weight}  每个 mutation pair 的权重（只考虑长度>=2的 clone）
    resolution : float
        Leiden 算法分辨率参数
    seed : int
        随机种子
    
    Returns
    -------
    mutation_group : dict
        {mutation_id: group_id} 每个 mutation 对应的群组
    partition : leidenalg VertexPartition
        Leiden 算法返回的 partition 对象（可用于可视化等）
    G_ig : igraph Graph
        构建的 igraph 图
    """
    # 1. 收集所有 mutation（包括孤立节点）
    all_mutations = set()
    for clone in clone_weights.keys():
        all_mutations.update(clone)
    
    # 2. 构建 igraph 图
    G_ig = ig.Graph()
    G_ig.add_vertices(list(all_mutations))  # 所有 mutation 作为节点
    
    # 添加边（只考虑长度>=2的 clone）
    for (m1, m2), w in pair_weights.items():
        G_ig.add_edge(m1, m2, weight=float(w))
    
    # 3. 运行 Leiden 算法
    partition = leidenalg.find_partition(
        G_ig,
        leidenalg.RBConfigurationVertexPartition,
        weights='weight',
        resolution_parameter=resolution,
        seed=seed
    )
    
    # 4. 输出 mutation -> group 字典
    mutation_group = {}
    for idx, community in enumerate(partition):
        for v in community:
            mutation_group[G_ig.vs[v]['name']] = idx
    
    # 5. 绘图
    plot_mutation_graph(G_ig, mutation_group, pdf_file)
    
    return mutation_group, partition, G_ig


##### 在每一个 graph 划分的 groups 中找到 hub group
def detect_hub_clusters(G_ig, mutation_group):
    """
    基于加权度中心性检测hub群体
    """
    # 1. 构建群体级别的图（就是你例子中的方法）
    cluster_graph = build_cluster_graph(G_ig, mutation_group)
    
    # 2. 计算每个群体的加权度
    cluster_degrees = {}
    for cluster_id in set(mutation_group.values()):
        weighted_degree = 0
        for edge in cluster_graph.es:
            source_cluster = cluster_graph.vs[edge.source]['name']
            target_cluster = cluster_graph.vs[edge.target]['name']
            
            if source_cluster == cluster_id or target_cluster == cluster_id:
                weighted_degree += edge['weight']
        
        cluster_degrees[cluster_id] = weighted_degree
    
    # 3. 识别hub群体（阈值可调整）
    hub_threshold = np.percentile(list(cluster_degrees.values()), 75)
    hub_clusters = [cluster_id for cluster_id, degree in cluster_degrees.items() 
                   if degree > hub_threshold]
    
    return hub_clusters, cluster_degrees

def build_cluster_graph(G_ig, mutation_group):
    """
    构建群体级别的加权图（就是你图示的方法）
    """
    clusters = set(mutation_group.values())
    cluster_graph = ig.Graph()
    cluster_graph.add_vertices(list(clusters))
    
    # 计算群体间的连接权重
    inter_cluster_weights = {}
    for edge in G_ig.es:
        source_mut = G_ig.vs[edge.source]['name']
        target_mut = G_ig.vs[edge.target]['name']
        
        source_cluster = mutation_group[source_mut]
        target_cluster = mutation_group[target_mut]
        
        if source_cluster != target_cluster:
            pair = tuple(sorted([source_cluster, target_cluster]))
            inter_cluster_weights[pair] = inter_cluster_weights.get(pair, 0) + edge['weight']
    
    # 添加边
    for (cluster1, cluster2), weight in inter_cluster_weights.items():
        cluster_graph.add_edge(cluster1, cluster2, weight=weight)
    
    return cluster_graph




##### 根据 immune mutations 拆分 spots
def resolved_spots_by_immune_mutations(I_scaffold, immune_mutations, P_scaffold, V_scaffold, A_scaffold, C_scaffold, df_reads_scaffold, p_threshold=0.5):
    """
    根据immune_mutations拆分I_scaffold矩阵中的spot，并同时修改P_scaffold、C_scaffold、V_scaffold和A_scaffold中的相关突变。 (pandas DataFrame, 行是spot，列是mutation； 需要根据I_scaffold的拆分同时处理其他几个矩阵)
    
    参数:
    I_scaffold: 已经根据immune_mutations拆分 (Binary indicator matrix)
    immune_mutations: list, immune mutation的ID列表
    P_scaffold: Posterior probability matrix
    V_scaffold: Variant/Mutant allele frequency matrix
    A_scaffold: Alternative allele count matrix
    C_scaffold: Coverage matrix
    df_reads_scaffold: Reads info (alt/total) matrix
    p_threshold: float, 判断是否为突变的阈值，默认为 0.5

    返回:
    I_resolved: 处理后的I_scaffold矩阵 (Binary indicator matrix)
    P_resolved: 处理后的P_scaffold矩阵 (Posterior probability matrix)
    V_resolved: 处理后的V_scaffold矩阵 (Variant/Mutant allele frequency matrix)
    A_resolved: 处理后的A_scaffold矩阵 (Alternative allele count matrix)
    C_resolved: 处理后的C_scaffold矩阵 (Coverage matrix)
    df_reads_resolved: df_reads_scaffold (Reads info (alt/total) matrix)
    """
    
    # 找出在I_scaffold列中实际存在的immune mutations
    actual_immune_mutations = [mut for mut in immune_mutations if mut in I_scaffold.columns]
    print(f"找到 {len(actual_immune_mutations)} 个immune mutations在数据中")
    
    # 找出non-immune mutations (在I_scaffold中但不在immune_mutations中的列)
    non_immune_mutations = [col for col in I_scaffold.columns if col not in immune_mutations]
    print(f"找到 {len(non_immune_mutations)} 个non-immune mutations")
    
    # 存储结果的列表
    resolved_rows = []
    resolved_index = []
    
    P_resolved_rows = []
    C_resolved_rows = []
    V_resolved_rows = []
    A_resolved_rows = []
    df_reads_resolved_rows = []
    
    count = 0
    spots_to_split = []  # 用于记录需要拆分的spot的ID
    
    # 遍历每一行(每个spot)
    for spot_name, row in I_scaffold.iterrows():
        # 将<NA>转换为NaN以便处理
        I_row_clean = row.replace('<NA>', np.nan)
        
        # 对每个突变位点应用 p_threshold 判断
        I_row_clean = I_row_clean.apply(lambda x: 1 if pd.notna(x) and x > p_threshold else np.nan)
        
        # 检查是否有immune mutations为1
        immune_has_1 = False
        if actual_immune_mutations:
            immune_values = I_row_clean[actual_immune_mutations]
            immune_has_1 = any(immune_values == 1)
        
        # 检查是否有non-immune mutations为1
        non_immune_has_1 = False
        if non_immune_mutations:
            non_immune_values = I_row_clean[non_immune_mutations]
            non_immune_has_1 = any(non_immune_values == 1)
        
        # 如果同时有immune和non-immune mutations为1，需要拆分
        if immune_has_1 and non_immune_has_1:
            count += 1
            spots_to_split.append(spot_name)  # 记录需要拆分的spot
            
            # 创建immune版本的行
            immune_row = row.copy()
            # 将non-immune mutations设为NaN
            immune_row[non_immune_mutations] = np.nan
            resolved_rows.append(immune_row)
            resolved_index.append(f"{spot_name}-immune")
            
            # 创建non-immune版本的行
            non_immune_row = row.copy()
            # 将immune mutations设为NaN
            non_immune_row[actual_immune_mutations] = np.nan
            resolved_rows.append(non_immune_row)
            resolved_index.append(f"{spot_name}-non")
            
            # 同时处理 P_scaffold
            immune_P_row = P_scaffold.loc[spot_name].copy()
            immune_P_row[non_immune_mutations] = np.nan  # 将non-immune列设为NaN
            P_resolved_rows.append(immune_P_row)
            non_immune_P_row = P_scaffold.loc[spot_name].copy()
            non_immune_P_row[actual_immune_mutations] = np.nan  # 将immune列设为NaN
            P_resolved_rows.append(non_immune_P_row)
            
            # 更新V_scaffold
            immune_M_row = V_scaffold.loc[spot_name].copy()
            immune_M_row[non_immune_mutations] = np.nan
            V_resolved_rows.append(immune_M_row)
            non_immune_M_row = V_scaffold.loc[spot_name].copy()
            non_immune_M_row[actual_immune_mutations] = np.nan
            V_resolved_rows.append(non_immune_M_row)
            
            # 更新A_scaffold
            immune_A_row = A_scaffold.loc[spot_name].copy()
            immune_A_row[non_immune_mutations] = np.nan
            A_resolved_rows.append(immune_A_row)
            non_immune_A_row = A_scaffold.loc[spot_name].copy()
            non_immune_A_row[actual_immune_mutations] = np.nan
            A_resolved_rows.append(non_immune_A_row)
            
            # 更新C_scaffold
            immune_C_row = C_scaffold.loc[spot_name].copy()
            immune_C_row[non_immune_mutations] = np.nan
            C_resolved_rows.append(immune_C_row)
            non_immune_C_row = C_scaffold.loc[spot_name].copy()
            non_immune_C_row[actual_immune_mutations] = np.nan
            C_resolved_rows.append(non_immune_C_row)
            
            # 更新df_reads_scaffold
            immune_df_reads_row = df_reads_scaffold.loc[spot_name].copy()
            immune_df_reads_row[non_immune_mutations] = np.nan
            df_reads_resolved_rows.append(immune_df_reads_row)
            non_immune_df_reads_row = df_reads_scaffold.loc[spot_name].copy()
            non_immune_df_reads_row[actual_immune_mutations] = np.nan
            df_reads_resolved_rows.append(non_immune_df_reads_row)
                
        else:
            # 不需要拆分，直接保留原行
            resolved_rows.append(row)
            resolved_index.append(spot_name)
            
            # 同时处理 P_scaffold
            P_resolved_rows.append(P_scaffold.loc[spot_name])
            V_resolved_rows.append(V_scaffold.loc[spot_name])
            A_resolved_rows.append(A_scaffold.loc[spot_name])
            C_resolved_rows.append(C_scaffold.loc[spot_name])
            df_reads_resolved_rows.append(df_reads_scaffold.loc[spot_name])
    
    # 创建新的DataFrame
    I_resolved = pd.DataFrame(resolved_rows, index=resolved_index, columns=I_scaffold.columns)
    P_resolved = pd.DataFrame(P_resolved_rows, index=resolved_index, columns=P_scaffold.columns)
    V_resolved = pd.DataFrame(V_resolved_rows, index=resolved_index, columns=V_scaffold.columns)
    A_resolved = pd.DataFrame(A_resolved_rows, index=resolved_index, columns=A_scaffold.columns)
    C_resolved = pd.DataFrame(C_resolved_rows, index=resolved_index, columns=C_scaffold.columns)
    # df_reads 要把 bulk 行加上
    df_reads_resolved = pd.DataFrame(df_reads_resolved_rows, index=resolved_index, columns=df_reads_scaffold.columns)
    bulk_row = df_reads_scaffold.loc['bulk']
    df_reads_resolved_with_bulk = pd.concat([pd.DataFrame(bulk_row).T, df_reads_resolved])
    df_reads_resolved_with_bulk.index = ['bulk'] + df_reads_resolved.index.tolist()
    
    print(f"处理完成: 原始spot数 {len(I_scaffold)}, 处理后spot数 {len(I_resolved)}")
    print(f"一共有 {count} 个需要被拆分的 spots")
    # print(f"需要拆分的spots的ID: {spots_to_split}")  # 输出拆分的spot ID
    
    return I_resolved, P_resolved, V_resolved, A_resolved, C_resolved, df_reads_resolved_with_bulk, spots_to_split


# I_resolved, P_resolved, V_resolved, A_resolved, C_resolved = resolved_spots_by_immune_mutations(I_scaffold, immune_mutations, P_scaffold, V_scaffold, A_scaffold, C_scaffold, p_threshold=0.5)


def split_spots_by_immune_mutations_scaffold(
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


def sort_I_hierarchical_freeze_ones_fixed(I_S, mutation_group):
    """
    对突变矩阵 I_S 进行排序：
    - 列顺序：先按 group 升序，组内按 MCF (Mutation Cell Fraction) 降序
    - 行顺序：逐列“冻结”1，把包含 1 的 cell 提前，剩下的继续往后排序
    
    Parameters
    ----------
    I_S : pd.DataFrame
        突变矩阵，行是 cell，列是 mutation，值为 {0, 1, NA}
    mutation_group : dict
        {mutation_id: group_id}，每个突变所属的分组信息
    
    Returns
    -------
    I_sorted : pd.DataFrame
        已排序后的矩阵
    mut_df_sorted : pd.DataFrame
        突变信息表，包含 mutation, group, mcf，并已排序
    group_to_muts : dict
        每个 group 下对应的 mutation 列表
    final_order : list
        行（cell）的最终顺序
    """
    
    # Step 1: 保留 I_S 中存在且在 mutation_group 里的突变
    selected_muts = [m for m in I_S.columns if m in mutation_group]
    if len(selected_muts) == 0:
        raise ValueError("在 I_S.columns 中未找到任何 mutation_group 的 keys。")
    
    # Step 2: 计算每个突变的 group 和 MCF (n1 / (n1+n0)，NA 不计入)
    #   - n1: 出现 1 的 cell 数
    #   - n0: 出现 0 的 cell 数
    #   - mcf: 突变在 cell 中的“支持度”
    mut_info = []
    for mut in selected_muts:
        grp = mutation_group[mut]
        col = I_S[mut]
        n1 = (col == 1).fillna(False).sum()
        n0 = (col == 0).fillna(False).sum()
        mcf = n1 / (n1 + n0) if (n1 + n0) > 0 else -1.0
        mut_info.append({"mutation": mut, "group": grp, "mcf": mcf})
    mut_df = pd.DataFrame(mut_info)
    
    # Step 3: 按 group 升序、组内 mcf 降序对突变排序
    mut_df_sorted = mut_df.sort_values(
        by=["group", "mcf"], ascending=[True, False]
    ).reset_index(drop=True)
    sorted_mutations = mut_df_sorted["mutation"].tolist()
    
    # Step 4: 行排序（逐列冻结 1）
    #   - 遍历突变列
    #   - 每一列把值为 1 的 cell 提前（冻结），其余继续下一列判断
    remaining = list(I_S.index)   # 初始剩余 cell（保持原始顺序）
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
        # 冻结 1（直接排到前面）
        final_order.extend(ones)
        # 剩余 cell 进入下一轮排序
        remaining = zeros + nas
    
    # Step 5: 把剩下的 cell 接到后面
    final_order.extend(remaining)
    
    # Step 6: 生成排序后的 DataFrame
    I_sorted = I_S.loc[final_order, sorted_mutations].copy()
    
    # Step 7: 验证排序是否正确（包括 NA 匹配）
    orig_sub = I_S.loc[final_order, sorted_mutations]
    mask_eq = orig_sub.eq(I_sorted)
    mask_both_na = orig_sub.isna() & I_sorted.isna()
    mask_ok = mask_eq.fillna(False) | mask_both_na
    if not mask_ok.all().all():
        bad = (~mask_ok).sum().sum()
        raise RuntimeError(f"排序后有 {bad} 个 cell-mut 值不匹配（不应该发生）。")
    
    # Step 8: 输出每个 group 对应的突变列表
    group_to_muts = mut_df_sorted.groupby("group")["mutation"].apply(list).to_dict()
    
    return I_sorted, mut_df_sorted, group_to_muts, final_order

# I_selected_and_sorted, mut_df_sorted, group_to_muts, final_order = sort_I_hierarchical_freeze_ones_fixed(I_S, mutation_group)

def map_group_to_backbone_mutations(mutation_group, group_to_muts, backbone_mutations):
    """
    根据 mutation_group 中的组号，将 group_to_muts 中的组号替换为对应的 backbone_mutations 中的突变值。
    
    参数:
    mutation_group (dict): 每个突变到其对应组号的映射，格式为 {mutation: group}
    group_to_muts (dict): 每个组号到突变列表的映射，格式为 {group: [mutations]}
    backbone_mutations (list): backbone_mutations 列表，包含突变标识符
    
    返回:
    dict: 以 backbone_mutations 中的突变为键，每个组对应的突变列表为值，格式为 {backbone_mutation: [mutations]}
    """
    # Step 1: 创建一个映射，从 mutation 到其对应的组号
    mutation_to_group = {mutation: group for mutation, group in mutation_group.items()}
    
    # Step 2: 使用这个映射，将 group_to_muts 中的 keys 替换为对应的 backbone_mutations
    group_to_muts_with_backbone = {}
    
    # 遍历每个组号，重新命名为 backbone_mutations 中的突变值
    for group, mutations in group_to_muts.items():
        # 找到与组号对应的 backbone mutation
        new_key = [mutation for mutation in backbone_mutations if mutation_to_group.get(mutation) == group]
        if new_key:  # 确保找到了对应的 backbone mutation
            group_to_muts_with_backbone[new_key[0]] = mutations
    
    return group_to_muts_with_backbone

# # 使用示例
# group_to_muts_with_backbone = map_group_to_backbone_mutations(mutation_group, group_to_muts, backbone_mutations)




def compute_corr_cache_with_new_mut_scaffold(I_attached, existing_muts, new_mut):
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

def compute_new_mut_clone_affinity_correct_scaffold(new_mut, mutation_clones_rescue, I_attached, n_shuffle=100):
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
    corr_cache = compute_corr_cache_with_new_mut_scaffold(I_attached, all_existing_muts, new_mut)
    
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

def select_best_clone_scaffold(detailed_scores):
    """
    选择最佳的clone scaffold
    
    Args:
        detailed_scores: clone的详细评分字典，key是包含突变的元组
    """
    # 第一步: 找到 rep_correlation 最大的值
    max_rep_correlation = max([score['rep_correlation'] for score in detailed_scores.values()])
    
    if max_rep_correlation == 0:
        return []
    
    # 第二步: 找到所有 rep_correlation 为最大值的 clone
    max_rep_correlation_clones = [
        clone for clone, score in detailed_scores.items() if score['rep_correlation'] == max_rep_correlation
    ]
    
    if len(max_rep_correlation_clones) == 1:
        return max_rep_correlation_clones
    
    # 第四步: 比较 direct_ratio
    max_direct_ratio = max([detailed_scores[clone]['direct_ratio'] for clone in max_rep_correlation_clones])
    
    if max_direct_ratio == 0:
        return []
    
    # 第五步: 找到 direct_ratio 最大值的 clone
    best_clones = [
        clone for clone in max_rep_correlation_clones 
        if detailed_scores[clone]['direct_ratio'] == max_direct_ratio
    ]
    
    if len(best_clones) == 1:
        return best_clones
    
    # 第六步: 如果还有多个，选择突变数量最多的（直接用len(clone)）
    max_mutation_count = max([len(clone) for clone in best_clones])
    best_clones_with_most_mutations = [
        clone for clone in best_clones 
        if len(clone) == max_mutation_count
    ]
    
    if len(best_clones_with_most_mutations) == 1:
        return best_clones_with_most_mutations
    
    # 第七步: 如果突变数量也相同，按第一个突变名称排序取第一个
    return [sorted(best_clones_with_most_mutations, key=lambda x: x[0])[0]]

def add_new_mutation_to_clone(mutation_group_with_non_group_mutations, assigned_clone, new_mut):
    """
    根据 assigned_clone 中的突变，将 new_mut 添加到相应的组中。
    
    参数:
    mutation_group_with_non_group_mutations (dict): 包含突变及其组号的字典。
    assigned_clone (list of tuple): 每个元素是一个包含突变的元组，表示属于同一组的突变。
    new_mut (str): 新的突变，应该被添加到 assigned_clone 中的某个组。
    
    返回:
    dict: 更新后的 mutation_group_with_non_group_mutations 字典。
    """
    # Step 1: 找到 new_mut 所在的组号
    group_of_new_mut = None
    for mutation in assigned_clone[0]:  # 假设每个 tuple 中只有一组突变
        if mutation in mutation_group_with_non_group_mutations:
            group_of_new_mut = mutation_group_with_non_group_mutations[mutation]
            break  # 找到组号，退出循环
    
    if group_of_new_mut is not None:
        # Step 2: 将 new_mut 添加到该组
        mutation_group_with_non_group_mutations[new_mut] = group_of_new_mut
    else:
        print(f"Error: Could not find group for mutation {new_mut}")
    
    return mutation_group_with_non_group_mutations


# def process_new_mut(new_mut, group_to_muts_with_backbone, I_somatic_resolved, mutation_group_with_non_group_mutations):
#     """
#     处理一个新的突变，计算它的克隆亲和力，选择最佳克隆，并更新 mutation_group。
#     """
#     clone_affinity, detailed_scores = compute_new_mut_clone_affinity_correct_scaffold(
#         new_mut, 
#         group_to_muts_with_backbone, 
#         I_somatic_resolved,
#         n_shuffle=100
#     )
#     assigned_clone = select_best_clone_scaffold(detailed_scores)
    
#     if len(assigned_clone) > 0:
#         # 更新 mutation_group_with_non_group_mutations
#         updated_group = add_new_mutation_to_clone(mutation_group_with_non_group_mutations, assigned_clone, new_mut)
#         return updated_group
#     return mutation_group_with_non_group_mutations

# def parallel_process_mutations(I_somatic, mutation_clones_rescue, I_somatic_resolved, mutation_group_with_non_group_mutations, group_mutations):
#     """
#     并行处理每个新的突变。
#     """
#     # Step 1: 确定需要处理的突变
#     mutations_to_process = [i for i in I_somatic.columns if i not in group_mutations]
    
#     # Step 2: 使用线程池并行处理突变
#     with ThreadPoolExecutor() as executor:
#         futures = []
        
#         # 提交任务到线程池
#         for new_mut in mutations_to_process:
#             futures.append(executor.submit(process_new_mut, new_mut, mutation_clones_rescue, I_somatic_resolved, mutation_group_with_non_group_mutations))
        
#         # 等待任务完成，并收集结果
#         for future in as_completed(futures):
#             mutation_group_with_non_group_mutations = future.result()
    
#     return mutation_group_with_non_group_mutations

# # # 使用示例
# # mutation_group_with_non_group_mutations = parallel_process_mutations(I_somatic, mutation_clones_rescue, I_somatic_resolved, mutation_group_with_non_group_mutations, group_mutations)

# # # 打印最终结果
# # print(mutation_group_with_non_group_mutations)




def plot_heatmap_with_celltype_by_your_sorting(I_raw, df_celltype, mutation_group, your_sorting_muts, pdf_file):
    """
    绘制突变矩阵的热图（带 cell type 注释条 + 行列统计条形图 + legend）。
    
    Parameters
    ----------
    I_raw : pd.DataFrame
        原始突变矩阵 (cell x mutation)，元素为 {0,1,NA}
    df_celltype : pd.DataFrame
        每个 cell 的 cell_type 信息（必须包含 "cell_type" 列）
    mutation_group : dict
        {mutation_id: group_id}，指定每个 mutation 的分组
    your_sorting_muts : list
        预先定义好的 mutation 顺序（列顺序）
    pdf_file : str
        保存图像的 PDF 文件路径
    """
    
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap, BoundaryNorm
    from matplotlib.patches import Patch
    
    # -------------------
    # Step 1: 确定 cell 顺序（保持原始顺序，并筛选需要的 cell）
    # -------------------
    ordered_leaf_nodes = list(I_raw.index)
    cells_to_use = [c for c in ordered_leaf_nodes if c in I_raw.index]
    I_sorted_cells = I_raw.loc[cells_to_use]
    
    # -------------------
    # Step 2: 确定 mutation 顺序（使用外部提供的排序）
    # -------------------
    mutation_order = your_sorting_muts
    I_sorted = I_sorted_cells[mutation_order]
    
    # -------------------
    # Step 3: 转换矩阵数值（NA → 2，用于单独上色）
    # -------------------
    I_numeric = I_sorted.fillna(np.nan).apply(pd.to_numeric, errors="coerce")
    I_plot = I_numeric.copy()
    I_plot = I_plot.where(~I_plot.isna(), 2)  # 把 NA 填为 2
    
    # -------------------
    # Step 4: 计算行/列突变数量统计（用于条形图）
    # -------------------
    row_sums = I_numeric.sum(axis=1, skipna=True)  # 每个 cell 的突变数
    col_sums = I_numeric.sum(axis=0, skipna=True)  # 每个突变被多少 cell 支持
    
    # -------------------
    # Step 5: 布局 GridSpec
    #   - 左边：行条形图
    #   - 中间：热图
    #   - 上面：列条形图
    #   - 右边：cell type 注释条
    #   - 底部：legend
    # -------------------
    fig = plt.figure(figsize=(12, 10))
    gs = fig.add_gridspec(5, 6,
                          width_ratios=[0.3, 0.3, 3, 0.05, 0.05, 0.05],
                          height_ratios=[0.5, 0.05, 3, 0.3, 0.3],
                          wspace=0.05, hspace=0.05)
    
    ax_row_bar = fig.add_subplot(gs[2, 0])   # 左边行条形图
    ax_heatmap = fig.add_subplot(gs[2, 2])   # 中间热图
    ax_celltype_bar = fig.add_subplot(gs[2, 3])  # 右边 cell type 注释条
    ax_col_bar = fig.add_subplot(gs[0, 2])   # 上方列条形图
    ax_dummy = fig.add_subplot(gs[0, 0]); ax_dummy.axis("off")  # 占位
    
    # -------------------
    # Step 6: 绘制热图
    #   - 颜色映射：0=浅蓝, 1=深红, NA=白色
    # -------------------
    cmap = ListedColormap(["#D4E8F0", "#7D2224", "white"])
    bounds = [0, 0.5, 1.5, 2.5]
    norm = BoundaryNorm(bounds, cmap.N)
    
    im = ax_heatmap.imshow(I_plot, aspect="auto", cmap=cmap,
                           interpolation="nearest", norm=norm)
    ax_heatmap.set_xlim(-0.5, I_plot.shape[1]-0.5)
    ax_heatmap.set_ylim(I_plot.shape[0]-0.5, -0.5)
    ax_heatmap.set_yticks([])
    
    # 设置横坐标 mutation 名称
    ax_heatmap.set_xticks(range(len(I_plot.columns)))
    ax_heatmap.set_xticklabels(I_plot.columns, rotation=90, fontsize=6, ha='center')
    
    # -------------------
    # Step 7: 给横坐标 mutation label 上色（根据分组）
    # -------------------
    unique_groups = sorted(set(mutation_group.values()))
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
    # Step 8: 列条形图（每个突变被多少 cell 支持）
    # -------------------
    ax_col_bar.bar(range(len(col_sums)), col_sums.values,
                   color="#7D2224", alpha=0.7, align="center")
    ax_col_bar.set_xlim(ax_heatmap.get_xlim())
    ax_col_bar.set_xticks([])
    ax_col_bar.tick_params(axis="y", labelsize=8)
    ax_col_bar.set_ylabel("#Mutations", fontsize=10)
    
    # -------------------
    # Step 9: 行条形图（每个 cell 有多少突变）
    # -------------------
    ax_row_bar.barh(range(len(row_sums)), row_sums.values,
                    color="#7D2224", alpha=0.7, align="center")
    ax_row_bar.set_ylim(ax_heatmap.get_ylim())
    ax_row_bar.set_yticks([])
    ax_row_bar.set_xlabel("#Cells", fontsize=10)
    ax_row_bar.invert_xaxis()
    
    # -------------------
    # Step 10: 绘制 cell type 注释条（右边）
    # -------------------
    celltypes = df_celltype["cell_type"].astype("category")
    type_codes = celltypes.cat.codes
    unique_types = celltypes.cat.categories
    cmap_types = plt.cm.get_cmap("tab20", len(unique_types))
    
    ax_celltype_bar.imshow(np.array(type_codes)[:, None], aspect="auto",
                           cmap=cmap_types, origin="upper")
    ax_celltype_bar.set_xticks([])
    ax_celltype_bar.set_yticks([])
    
    # -------------------
    # Step 11: 添加 Legend
    #   - 突变值 (0,1,NA)
    #   - cell types
    #   - mutation groups
    # -------------------
    # 突变值 legend
    heatmap_handles = [Patch(facecolor=c, label=l) 
                       for c, l in zip(["#D4E8F0", "#7D2224", "white"],
                                       ["0 (No Mutation)", "1 (Mutation)", "NA (Missing)"])]
    fig.legend(handles=heatmap_handles, loc="lower center", ncol=3,
               bbox_to_anchor=(0.5, -0.03), frameon=False, fontsize=9,
               title="Mutation Values", title_fontsize=10)
    
    # cell type legend
    celltype_handles = [Patch(facecolor=cmap_types(i), label=label) 
                        for i, label in enumerate(unique_types)]
    fig.legend(handles=celltype_handles, loc="lower center",
               ncol=min(len(unique_types), 5),
               bbox_to_anchor=(0.5, -0.12), frameon=False, fontsize=9,
               title="Cell Types", title_fontsize=10)
    
    # mutation group legend
    group_handles = [Patch(facecolor=color, label=f'Group {group_id}')
                     for group_id, color in group_colors.items()]
    fig.legend(handles=group_handles, loc="lower center",
               ncol=min(len(group_handles), 5),
               bbox_to_anchor=(0.5, -0.20), frameon=False, fontsize=9,
               title="Mutation Groups", title_fontsize=10)
    
    # -------------------
    # Step 12: 保存图像
    # -------------------
    plt.suptitle("Heatmap of Mutations with Cell Type Bar", fontsize=14, y=0.95)
    plt.tight_layout()
    plt.subplots_adjust(bottom=0)
    plt.savefig(pdf_file, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"保存完成: {pdf_file}")

def select_founder_mutations(I_selected_and_sorted: pd.DataFrame, mutation_group: dict):
    """
    根据 mutation_group 分组，从 I_selected_and_sorted 里挑选每组 founder mutation。
    Founder mutation 定义为该组中 Mutant Cell Fraction (MCF) 最大的突变。
    
    Parameters
    ----------
    I_selected_and_sorted : pd.DataFrame
        突变矩阵 (cells x mutations)，元素为 {1, 0, NA}。
    mutation_group : dict
        {mutation: group_id}，每个突变所属的分组。
    
    Returns
    -------
    dict
        {group_id: founder_mutation}
    """
    # 结果 dict
    founder_dict = {}
    
    # 转换成 DataFrame 方便操作
    mut2group = pd.Series(mutation_group)
    
    for group_id, muts in mut2group.groupby(mut2group):
        muts = list(muts.index)  # 当前组的突变
        mcf_values = {}
        mcn_values = {}
        for mut in muts:
            col = I_selected_and_sorted[mut]
            num_mutant = (col == 1).sum()     # 突变细胞数
            num_covered = col.notna().sum()   # 有测序覆盖的细胞数
            mcf = num_mutant / num_covered if num_covered > 0 else 0
            mcf_values[mut] = mcf
            mcn_values[mut] = num_mutant
        
        # 找到 MCF 最大的突变
        founder_mut = max(mcn_values, key=mcn_values.get)
        founder_dict[group_id] = founder_mut
    
    return founder_dict


##### 基于当前分组检查每一个突变是否分对了 group

def calculate_intersection_counts_under_backbone_nodes_scaffold(mutation_list_under_backbone_nodes, M_current, I_attached, new_mut):
    """
    计算new_mut与每个backbone node下mutation list的共现数量
    
    参数:
    mutation_list_under_backbone_nodes: dict, backbone node到mutation list的映射
    M_current: DataFrame, 细胞在backbone nodes的状态 (0/1)
    I_attached: DataFrame, 细胞在其他突变上的状态 (0/1/NaN)
    new_mut: str, 要检查的新突变
    
    返回:
    dict: 每个backbone node对应的共现数量
    """
    
    # 检查new_mut是否在I_attached中
    if new_mut not in I_attached.columns:
        raise ValueError(f"突变 {new_mut} 不在 I_attached 数据框中")
    
    intersection_counts = {}
    
    for backbone_node, mutation_list in mutation_list_under_backbone_nodes.items():
        # 1. 找到当前backbone node为1的细胞
        backbone_cells = M_current[M_current[backbone_node] == 1].index
        
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
                # 获取当前mutation在new_mut阳性细胞中的状态
                mut_status = I_attached.loc[new_mut_positive_cells, mutation]
                # 计算同时为1的数量（忽略NaN值）
                intersection_count = (mut_status == 1).sum()
                total_intersection += intersection_count
        
        intersection_counts[backbone_node] = total_intersection
    
    return intersection_counts

def find_best_backbone_node_scaffold(mutation_list_under_backbone_nodes, M_current, intersection_counts_under_backbone_node):
    """
    简化版本：只返回最佳backbone node
    """
    max_count = max(intersection_counts_under_backbone_node.values())
    max_nodes = [node for node, count in intersection_counts_under_backbone_node.items() 
                if count == max_count]
    
    if len(max_nodes) == 1:
        return max_nodes[0]
    
    # 平局处理
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
            # 如果没有细胞，归一化值为0
            if best_normalized < 0:  # 还没有找到有效节点
                best_node = node
                best_normalized = 0
    
    return best_node

def find_best_backbone_for_new_mutation_scaffold(mutation_list_under_backbone_nodes, M_current, I_attached, new_mut):
    """
    完整函数：计算共现数量并找到最佳backbone node
    """
    # 第一步：计算共现数量
    intersection_counts = calculate_intersection_counts_under_backbone_nodes_scaffold(
        mutation_list_under_backbone_nodes, M_current, I_attached, new_mut
    )
    
    # 第二步：找到最佳backbone node
    best_backbone = find_best_backbone_node_scaffold(
        mutation_list_under_backbone_nodes, M_current, intersection_counts
    )
    
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
# TreeNode: 最小化树对象
# -------------------------
class TreeNode:
    def __init__(self, name: str):
        self.name = name            # 通常是 mutation id 或 "ROOT"
        self.parent: Optional['TreeNode'] = None
        self.children: List['TreeNode'] = []
    
    def add_child(self, child: 'TreeNode'):
        child.parent = self
        self.children.append(child)
    
    def remove_child(self, child: 'TreeNode'):
        # child 存在性检查
        if child in self.children:
            self.children.remove(child)
            child.parent = None
    
    def is_root(self):
        return self.parent is None
    
    def traverse(self):
        # yield nodes in pre-order
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
        # 在 parent->child 边上插入 new node，返回新节点
        # 假定 child 是 parent 的 child
        if child not in parent.children:
            raise ValueError("child not child of parent")
        new_node = TreeNode(new_node_name)
        # replace child by new_node under parent
        parent.remove_child(child)
        parent.add_child(new_node)
        # attach child under new_node
        new_node.add_child(child)
        return new_node
    
    def add_new_parent_for_children(self, children_list: List['TreeNode'], new_node_name: str) -> 'TreeNode':
        # 把 children_list 从它们原父节点移出，放入 new_parent；new_parent 接到原先的共同父节点或 root（取第一个孩子的parent）
        parents = set(c.parent for c in children_list)
        if len(parents) != 1:
            # 如果孩子们来自不同父节点，语义上允许（你可以定义策略），这里我们只要求同父
            pass
        old_parent = children_list[0].parent
        new_node = TreeNode(new_node_name)
        if old_parent is None:
            # attach to root (if allowed)
            pass
        else:
            # replace these children in old_parent with new_node
            for c in children_list:
                old_parent.remove_child(c)
            old_parent.add_child(new_node)
        # attach children to new_node
        for c in children_list:
            new_node.add_child(c)
        return new_node
    
    def add_leaf(self, parent: 'TreeNode', new_leaf_name: str) -> 'TreeNode':
        new_leaf = TreeNode(new_leaf_name)
        parent.add_child(new_leaf)
        return new_leaf
    
    def copy(self) -> 'TreeNode':
        # 深拷贝整棵树（保留 names）
        mapping = {}
        def _copy(node):
            n2 = TreeNode(node.name)
            mapping[node] = n2
            for c in node.children:
                n2.add_child(_copy(c))
            return n2
        return _copy(self)
    
    def all_nodes(self) -> List['TreeNode']:
        """返回整棵树的所有节点对象"""
        return list(self.traverse())
    
    def all_names(self) -> List[str]:
        """返回整棵树所有节点名字"""
        return [n.name for n in self.traverse()]
    
    def all_names_no_root(self) -> List[str]:
        return [n.name for n in self.traverse() if n.name != "ROOT"]
    
    def all_edges(self) -> List[tuple]:
        """返回整棵树所有边 (parent_name, child_name)"""
        return [(n.parent.name, n.name) for n in self.traverse() if n.parent]
    
    def to_string(self, level=0):
        """将树结构转换为字符串"""
        indent = "  " * level
        result = f"{indent}└─ {self.name}\n"
        for child in self.children:
            result += child.to_string(level + 1)
        return result
    
    def save_to_file(self, filename):
        """将树结构保存到文件"""
        with open(filename, 'w', encoding='utf-8') as f:
            f.write(self.to_string())
    
    def __str__(self):
        return self.to_string()

def print_tree(node, level=0):
    """递归打印树结构"""
    indent = "  " * level
    print(f"{indent}└─ {node.name}")
    for child in node.children:
        print_tree(child, level + 1)

# # 打印整个树
# print_tree(T)


##### 树 TreeNode 格式保存和读取

import json

# 假设 T_full 是你的 TreeNode 对象
def tree_to_dict(node):
    # 假设你已经有一个 TreeNode 类或类似的结构
    tree_dict = {}
    tree_dict['name'] = node.name
    if node.children:
        tree_dict['children'] = [tree_to_dict(child) for child in node.children]
    return tree_dict

# # 将 TreeNode 转换为字典格式
# tree_dict = tree_to_dict(T_full)

# # 保存为 JSON 文件
# with open(os.path.join(phylo_dir, 'tree_node.json'), 'w') as f:
#     json.dump(tree_dict, f, indent=4)


# import json

# # 从 JSON 文件加载树的字典
# with open(os.path.join(phylo_dir, 'tree_node.json'), 'r') as f:
#     tree_dict_loaded = json.load(f)

# # 假设你有一个 TreeNode 类，可以从字典重建树
# def dict_to_tree(tree_dict_loaded):
#     node = TreeNode(tree_dict_loaded['name'])
#     if 'children' in tree_dict_loaded:
#         node.children = [dict_to_tree(child) for child in tree_dict_loaded['children']]
#     return node

# # 重建 TreeNode 对象
# T_full_loaded = dict_to_tree(tree_dict_loaded)

def print_tree_dict(tree, prefix=""):
    print(prefix + "└─ " + tree["name"])
    if "children" in tree:
        for child in tree["children"]:
            print_tree_dict(child, prefix + "  ")

# tree_json_file = "/storage/douyanmeiLab/yangqing/tools/PhyloMosaicGenie/pmg/src/phylosolid/samples/phylo_wrapped_up/sample_P6_merged/mutation_integrator.germline_no/phylo_v11.okVersion_tripleFPs_but_backbone.send_to_YD/final_cleaned_tree_node.json"

# with open(tree_json_file, 'r') as f:
#     tree_json = json.load(f)

# print_tree_dict(tree_json)




def find_paths_to_leaves_tree_node(start_node):
    """
    找到从指定节点到所有叶子节点的路径（适用于TreeNode对象）
    
    Args:
        start_node: TreeNode起始节点
    
    Returns:
        list: 包含所有从起始节点到叶子节点路径的列表
    """
    def dfs(current_node, current_path, all_paths):
        """深度优先搜索找到所有到叶子节点的路径"""
        current_path = current_path + [current_node.name]
        
        if not current_node.children:  # 如果是叶子节点
            all_paths.append(current_path)
        else:
            for child in current_node.children:
                dfs(child, current_path, all_paths)
    
    all_paths = []
    dfs(start_node, [], all_paths)
    return all_paths

def get_all_descendants_tree_node(start_node):
    """
    获取指定节点的所有后代节点（适用于TreeNode对象）
    
    Args:
        start_node: TreeNode起始节点
    
    Returns:
        set: 所有后代节点名称的集合
    """
    def get_descendants(node):
        """递归获取所有后代节点"""
        descendants = set()
        for child in node.children:
            descendants.add(child.name)
            descendants |= get_descendants(child)
        return descendants
    
    return get_descendants(start_node)

def get_subtree_nodes_tree_node(start_node):
    """获取从指定节点开始的子树所有节点名称（适用于TreeNode对象）"""
    def collect_nodes(node):
        nodes = {node.name}
        for child in node.children:
            nodes |= collect_nodes(child)
        return nodes
    return collect_nodes(start_node)

# # 使用示例
# # 1. 首先找到目标节点
# target_node = tree_data.find('chr16_265963_C_T')  # 使用TreeNode的find方法

# if target_node:
#     # 2. 找到从目标节点到所有叶子节点的路径
#     paths = find_paths_to_leaves_tree_node(target_node)
#     print("从 chr16_265963_C_T 到叶子节点的路径:")
#     for i, path in enumerate(paths):
#         print(f"路径 {i+1}: {' -> '.join(path)}")
    
#     # 3. 获取所有后代节点
#     descendants = get_all_descendants_tree_node(target_node)
#     print(f"\nchr16_265963_C_T 的所有后代节点: {descendants}")
    
#     # 4. 获取子树所有节点
#     subtree_nodes = get_subtree_nodes_tree_node(target_node)
#     print(f"\n从 chr16_265963_C_T 开始的子树所有节点: {subtree_nodes}")
# else:
#     print("未找到节点 chr16_265963_C_T")




# -------------------------
# MutationMatrix: 包装 pandas DataFrame
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
# 构建一个初始的 backbone tree (T_B, M_B)
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

def impute_backbone_clones(selected_df, backbone_mutations, mutation_group):
    """
    根据分组信息计算backbone clones
    
    参数:
    selected_df: 包含突变数据的DataFrame (cells × mutations)
    group_mutations: 所有突变的列表，按顺序排列
    mutation_group: 突变到分组的映射字典
    
    返回:
    backbone_matrix: 每个backbone clone对应的细胞二元向量
    cell_assignments: 每个细胞最终分配的backbone clone
    """
        
    # 将NA值替换为0
    selected_df = selected_df.fillna(0)
    
    # 找出每个group对应的backbone mutation
    group_to_backbone = {}
    for mut in backbone_mutations:
        group = mutation_group[mut]
        group_to_backbone[group] = mut
    
    # 为每个backbone mutation计算imputed clone vector
    backbone_vectors = {}
    sum_vectors = {}
    for group, backbone_mutation in group_to_backbone.items():
        # 找到该group中的所有mutation
        group_muts = [mut for mut, g in mutation_group.items() if g == group]
        
        # 对这些mutation的vector进行加和
        group_df = selected_df[group_muts]
        sum_vector = group_df.sum(axis=1)
        sum_vectors[backbone_mutation] = sum_vector
        # 将>0的值设为1，得到二元vector
        imputed_vector = (sum_vector > 0).astype(int)
        backbone_vectors[backbone_mutation] = imputed_vector
    
    # 创建backbone matrix (cells × backbone_clones)
    backbone_matrix = pd.DataFrame(backbone_vectors)
    sum_matrix = pd.DataFrame(sum_vectors)
    
    # 处理细胞分配：找到每个细胞属于哪个backbone clone
    cell_assignments = {}
    conflict_resolution_stats = defaultdict(int)
    
    for cell in selected_df.index:
        # 获取该细胞在所有backbone clone中的值
        cell_backbone_values = backbone_matrix.loc[cell]
        
        # 找出该细胞属于哪些backbone clone (值为1的)
        belonging_clones = cell_backbone_values[cell_backbone_values == 1].index.tolist()
        
        if len(belonging_clones) == 0:
            # 不属于任何backbone clone
            cell_assignments[cell] = None
            conflict_resolution_stats['unassigned'] += 1
        
        elif len(belonging_clones) == 1:
            # 只属于一个backbone clone，直接分配
            cell_assignments[cell] = belonging_clones[0]
            conflict_resolution_stats['unique_assignment'] += 1
        
        else:
            # 属于多个backbone clone，需要解决冲突
            conflict_resolution_stats['conflicts'] += 1
            
            # 计算该细胞在每个conflicting clone中的mutation number
            clone_scores = {}
            for clone in belonging_clones:
                mutation_count = sum_vectors[clone].loc[cell]
                clone_scores[clone] = mutation_count
            
            # 找到最高分
            max_score = max(clone_scores.values())
            
            # 找出所有得分最高的clone
            top_clones = [clone for clone, score in clone_scores.items() if score == max_score]
            
            if len(top_clones) == 1:
                # 只有一个最高分，正常分配
                assigned_clone = top_clones[0]
                conflict_resolution_stats['resolved_unique_max'] += 1
            else:
                # 多个clone得分相同，需要二级判断
                conflict_resolution_stats['ties'] += 1
                conflict_resolution_stats[f'tie_score_{max_score}'] += 1
                
                # 二级判断策略：选择backbone mutation本身表达最强的
                backbone_expr = {clone: selected_df.loc[cell, clone] for clone in top_clones}
                assigned_clone = max(backbone_expr.items(), key=lambda x: x[1])[0]
                conflict_resolution_stats['resolved_by_expression'] += 1
            
            cell_assignments[cell] = assigned_clone
            conflict_resolution_stats[f'resolved_to_{assigned_clone}'] += 1
    
    # 创建最终的backbone assignment vector
    final_assignment = pd.Series(cell_assignments, name='backbone_assignment')
    
    # 创建互斥的backbone matrix (每个细胞只属于一个clone)
    exclusive_backbone_matrix = pd.DataFrame(0, 
                                           index=selected_df.index, 
                                           columns=backbone_matrix.columns)
    
    for cell, assignment in cell_assignments.items():
        if assignment is not None:
            exclusive_backbone_matrix.loc[cell, assignment] = 1
    
    print("冲突解决统计:")
    for stat, count in conflict_resolution_stats.items():
        print(f"  {stat}: {count} cells")
    
    return exclusive_backbone_matrix, final_assignment


def impute_backbone_clones_weighted(selected_df, backbone_mutations, mutation_group, backbone_weight=2.0):
    """
    根据分组信息计算backbone clones，允许细胞跨backbone分配
    
    参数:
    selected_df: 包含突变数据的DataFrame (cells × mutations)
    backbone_mutations: 所有backbone突变的列表
    mutation_group: 突变到分组的映射字典
    backbone_weight: backbone mutation的权重系数 (默认2.0)
                    weight=1: 所有mutation平等
                    weight>1: 更重视backbone mutation
                    weight<1: 更重视其他mutation
    
    返回:
    exclusive_backbone_matrix: 每个细胞只属于一个backbone clone的矩阵
    final_assignment: 每个细胞分配的backbone clone
    """
    
    # 将NA值替换为0
    selected_df = selected_df.fillna(0)
    
    # 找出每个group对应的backbone mutation
    group_to_backbone = {}
    for mut in backbone_mutations:
        group = mutation_group[mut]
        group_to_backbone[group] = mut
    
    # 为每个backbone mutation计算加权分数
    weighted_scores = {}
    sum_vectors = {}  # 保存原始求和向量用于统计
    
    for group, backbone_mutation in group_to_backbone.items():
        # 找到该group中的所有mutation
        group_muts = [mut for mut, g in mutation_group.items() if g == group]
        group_df = selected_df[group_muts]
        
        # 计算原始求和向量（用于统计）
        sum_vector = group_df.sum(axis=1)
        sum_vectors[backbone_mutation] = sum_vector
        
        # 计算加权分数：backbone mutation有更高权重
        weighted_vector = group_df.copy()
        if backbone_mutation in weighted_vector.columns:
            weighted_vector[backbone_mutation] = weighted_vector[backbone_mutation] * backbone_weight
        
        weighted_score = weighted_vector.sum(axis=1)
        weighted_scores[backbone_mutation] = weighted_score
    
    # 创建分数矩阵
    weighted_matrix = pd.DataFrame(weighted_scores)
    sum_matrix = pd.DataFrame(sum_vectors)
    
    # 分配细胞到backbone clones
    cell_assignments = {}
    assignment_stats = defaultdict(int)
    
    for cell in selected_df.index:
        # 获取该细胞在所有backbone clone中的加权分数
        cell_scores = weighted_matrix.loc[cell]
        max_score = cell_scores.max()
        
        # 如果最高分为0，不分配任何clone
        if max_score == 0:
            cell_assignments[cell] = None
            assignment_stats['unassigned'] += 1
            continue
        
        # 找出所有得分最高的clone
        top_clones = cell_scores[cell_scores == max_score].index.tolist()
        
        if len(top_clones) == 1:
            # 只有一个最高分，直接分配
            assigned_clone = top_clones[0]
            assignment_stats['unique_assignment'] += 1
        else:
            # 多个clone得分相同，需要二级判断
            assignment_stats['ties'] += 1
            # 选择backbone mutation本身表达最强的
            backbone_expr = {clone: selected_df.loc[cell, clone] for clone in top_clones}
            assigned_clone = max(backbone_expr.items(), key=lambda x: x[1])[0]
            assignment_stats['resolved_by_expression'] += 1
        
        # 记录分配类型统计
        original_backbones = [bm for bm in backbone_mutations if selected_df.loc[cell, bm] > 0]
        
        if not original_backbones:
            assignment_stats['assigned_without_original_backbone'] += 1
        elif assigned_clone in original_backbones:
            assignment_stats['assigned_to_original_backbone'] += 1
        else:
            assignment_stats['switched_to_other_backbone'] += 1
            assignment_stats[f'switched_from_{original_backbones[0]}_to_{assigned_clone}'] += 1
        
        cell_assignments[cell] = assigned_clone
    
    # 创建最终的backbone assignment vector
    final_assignment = pd.Series(cell_assignments, name='backbone_assignment')
    
    # 创建互斥的backbone matrix (每个细胞只属于一个clone)
    exclusive_backbone_matrix = pd.DataFrame(0, 
                                           index=selected_df.index, 
                                           columns=backbone_mutations)
    
    for cell, assignment in cell_assignments.items():
        if assignment is not None:
            exclusive_backbone_matrix.loc[cell, assignment] = 1
    
    # 打印统计信息
    print("=== Backbone Clone分配统计 ===")
    print(f"总细胞数: {len(selected_df)}")
    for stat, count in assignment_stats.items():
        if not stat.startswith('switched_from_'):
            print(f"  {stat}: {count}")
    
    # 打印backbone切换统计
    switch_stats = {k: v for k, v in assignment_stats.items() if k.startswith('switched_from_')}
    if switch_stats:
        print("\nBackbone切换详情:")
        for switch, count in switch_stats.items():
            print(f"  {switch}: {count}")
    
    return exclusive_backbone_matrix, final_assignment




import pandas as pd
import numpy as np
from collections import defaultdict, Counter
from typing import Dict, List, Tuple

def find_intersecting_mutations(I: pd.DataFrame, backbone_mutation: str, min_N11: int = 1) -> List[str]:
    """找到与backbone mutation有共现的所有突变"""
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
                                       remove_zeros: bool = False) -> Tuple[pd.DataFrame, pd.Series, dict]:
    """
    基于整个突变矩阵和backbone mutations进行全面的imputation
    
    参数:
    I_somatic: 包含所有突变数据的DataFrame (cells × mutations)，可以是原始数据
    backbone_mutations: backbone mutations列表
    min_N11: 最小共现次数阈值
    remove_zeros: 是否自动移除全0的行和列
    
    返回:
    exclusive_backbone_matrix: 每个细胞在backbone clones中的分配
    final_assignment: 每个细胞的最终backbone assignment
    intersection_sets: 每个backbone mutation对应的intersection突变集合
    """
    
    # 创建数据副本以避免修改原始数据
    I = I_somatic.copy()
    
    # 将NA值替换为0
    I = I.fillna(0)
    
    print(f"原始数据形状: {I.shape}")
    print(f"原始细胞数: {len(I.index)}")
    print(f"原始突变数: {len(I.columns)}")
    
    # 可选：自动移除全0的行和列
    if remove_zeros:
        # 移除全0列
        I = I.loc[:, (I != 0).any(axis=0)]
        # 移除全0行
        I = I.loc[(I != 0).any(axis=1)]
        print(f"清理后数据形状: {I.shape}")
        print(f"清理后细胞数: {len(I.index)}")
        print(f"清理后突变数: {len(I.columns)}")
    
    # 识别有突变的细胞（至少有一个突变值为1）
    has_mutation = (I == 1).any(axis=1)
    mutated_cells = I.index[has_mutation]
    non_mutated_cells = I.index[~has_mutation]
    
    print(f"\n细胞统计:")
    print(f"  有突变的细胞: {len(mutated_cells)}")
    print(f"  无突变的细胞: {len(non_mutated_cells)}")
    
    # 只对有突变的细胞进行backbone分配处理
    I_mutated = I.loc[mutated_cells]
    
    # 为每个backbone mutation找到intersecting mutations
    intersection_sets = {}
    backbone_vectors = {}
    sum_vectors = {}
    
    print("\n为每个backbone mutation寻找intersecting mutations...")
    for backbone_mut in backbone_mutations:
        if backbone_mut not in I.columns:
            print(f"警告: backbone mutation {backbone_mut} 不在数据列中，跳过")
            continue
            
        # 找到与该backbone mutation有共现的所有突变
        intersecting_muts = find_intersecting_mutations(I_mutated, backbone_mut, min_N11)
        intersection_sets[backbone_mut] = intersecting_muts
        
        print(f"Backbone {backbone_mut}: 找到 {len(intersecting_muts)} 个intersecting mutations")
        
        # 总是包含backbone mutation本身 + intersecting mutations
        all_relevant_muts = [backbone_mut] + intersecting_muts
        
        # 确保所有突变都在数据中
        available_muts = [mut for mut in all_relevant_muts if mut in I_mutated.columns]
        
        # 计算这些突变的加和向量
        group_df = I_mutated[available_muts]
        sum_vector = group_df.sum(axis=1)
        sum_vectors[backbone_mut] = sum_vector
        
        # 二值化：只要任一相关突变表达，就认为该backbone clone活跃
        imputed_vector = (sum_vector > 0).astype(int)
        backbone_vectors[backbone_mut] = imputed_vector
    
    # 创建backbone matrix（只在有突变的细胞上）
    backbone_matrix = pd.DataFrame(backbone_vectors, index=mutated_cells)
    
    # 计算每个突变的mutant cell number（用于平局时使用）
    mutation_prevalence = I_mutated.sum(axis=0)
    
    # 处理细胞分配
    cell_assignments = {}
    conflict_resolution_stats = defaultdict(int)
    unassigned_cells_mutations = {}  # 记录未分配细胞的突变信息
    
    print(f"\n处理 {len(mutated_cells)} 个有突变细胞的分配...")
    for i, cell in enumerate(mutated_cells):
        if i > 0 and i % 1000 == 0:
            print(f"  已处理 {i}/{len(mutated_cells)} 个细胞...")
        
        # 获取该细胞在所有backbone clone中的值
        cell_backbone_values = backbone_matrix.loc[cell]
        
        # 找出该细胞属于哪些backbone clone (值为1的)
        belonging_clones = cell_backbone_values[cell_backbone_values == 1].index.tolist()
        
        if len(belonging_clones) == 0:
            # 细胞有突变，但不属于任何backbone clone
            cell_assignments[cell] = None
            conflict_resolution_stats['mutated_but_unassigned'] += 1
            
            # 记录这个细胞的突变信息用于调试
            cell_mutations = I_mutated.columns[I_mutated.loc[cell] == 1].tolist()
            unassigned_cells_mutations[cell] = cell_mutations
        
        elif len(belonging_clones) == 1:
            # 只属于一个backbone clone，直接分配
            cell_assignments[cell] = belonging_clones[0]
            conflict_resolution_stats['unique_assignment'] += 1
        
        else:
            # 属于多个backbone clone，需要解决冲突
            conflict_resolution_stats['conflicts'] += 1
            
            # 第一级判断：计算在每个conflicting clone中的支持突变数
            clone_scores = {}
            for clone in belonging_clones:
                support_count = sum_vectors[clone].loc[cell]
                clone_scores[clone] = support_count
            
            # 找到最高分
            max_score = max(clone_scores.values())
            
            # 找出所有得分最高的clone
            top_clones = [clone for clone, score in clone_scores.items() if score == max_score]
            
            if len(top_clones) == 1:
                # 只有一个最高分，正常分配
                assigned_clone = top_clones[0]
                conflict_resolution_stats['resolved_unique_max'] += 1
            else:
                # 多个clone得分相同，需要二级判断
                conflict_resolution_stats['ties'] += 1
                
                # 二级判断：比较支持突变的prevalence
                prevalence_scores = {}
                for clone in top_clones:
                    # 获取该backbone clone的所有相关突变（包含backbone本身）
                    supporting_muts = [clone] + intersection_sets[clone]
                    available_muts = [mut for mut in supporting_muts if mut in mutation_prevalence.index]
                    # 计算支持突变的平均prevalence
                    avg_prevalence = mutation_prevalence[available_muts].mean()
                    prevalence_scores[clone] = avg_prevalence
                
                # 找到prevalence最高的clone
                max_prevalence = max(prevalence_scores.values())
                top_prevalence_clones = [clone for clone, score in prevalence_scores.items() if score == max_prevalence]
                
                if len(top_prevalence_clones) == 1:
                    assigned_clone = top_prevalence_clones[0]
                    conflict_resolution_stats['resolved_by_prevalence'] += 1
                else:
                    # 三级判断：选择当前imputed clone size最大的
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
    
    # 无突变的细胞保持未分配
    for cell in non_mutated_cells:
        cell_assignments[cell] = None
        conflict_resolution_stats['no_mutation'] += 1
    
    # 关键修复：确保所有原始细胞都在分配字典中
    all_original_cells = set(I_somatic.index)
    assigned_cells = set(cell_assignments.keys())
    missing_cells = all_original_cells - assigned_cells
    
    if missing_cells:
        print(f"警告: 发现 {len(missing_cells)} 个细胞未被分配，将它们设为未分配")
        for cell in missing_cells:
            cell_assignments[cell] = None
            conflict_resolution_stats['missing_cells_set_to_unassigned'] += 1
    
    # 创建最终的backbone assignment vector（所有原始细胞）
    final_assignment = pd.Series(cell_assignments, name='backbone_assignment')
    
    # 创建互斥的backbone matrix（所有原始细胞）
    exclusive_backbone_matrix = pd.DataFrame(0, 
                                           index=I_somatic.index,  # 使用原始索引
                                           columns=backbone_mutations)
    
    # 只为有分配的有突变细胞设置值
    for cell, assignment in cell_assignments.items():
        if assignment is not None:
            exclusive_backbone_matrix.loc[cell, assignment] = 1
    
    # 分析与backbone mutations没有交集的突变
    print(f"\n=== 与backbone mutations无交集的突变分析 ===")
    
    # 找出所有与任何backbone mutation都没有交集的突变
    all_mutations = set(I_mutated.columns)
    mutations_with_intersection = set()
    
    for backbone_mut, intersecting_muts in intersection_sets.items():
        mutations_with_intersection.add(backbone_mut)
        mutations_with_intersection.update(intersecting_muts)
    
    mutations_without_intersection = all_mutations - mutations_with_intersection
    
    print(f"总突变数: {len(all_mutations)}")
    print(f"与backbone mutations有交集的突变数: {len(mutations_with_intersection)}")
    print(f"与backbone mutations无交集的突变数: {len(mutations_without_intersection)}")
    
    if mutations_without_intersection:
        print("\n与backbone mutations无交集的突变列表:")
        total_cells_in_isolated_muts = 0
        for mut in sorted(mutations_without_intersection):
            mut_cell_count = (I_mutated[mut] == 1).sum()
            total_cells_in_isolated_muts += mut_cell_count
            print(f"  {mut}: {mut_cell_count}个细胞")
        print(f"这些孤立突变总共影响 {total_cells_in_isolated_muts} 个细胞")
    
    # 分析未分配细胞的突变组成
    if unassigned_cells_mutations:
        print(f"\n=== 未分配细胞详细分析 ===")
        print(f"未分配细胞数: {len(unassigned_cells_mutations)}")
        
        # 统计这些细胞的突变分布
        all_unassigned_mutations = []
        for cell, muts in unassigned_cells_mutations.items():
            all_unassigned_mutations.extend(muts)
        
        mutation_counts = Counter(all_unassigned_mutations)
        
        print(f"未分配细胞涉及 {len(mutation_counts)} 个不同的突变")
        
        # 分类统计
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
        
        print(f"\n未分配细胞中的突变类型:")
        print(f"  - Backbone mutations: {len(backbone_muts_in_unassigned)}个")
        if backbone_muts_in_unassigned:
            for mut, count in sorted(backbone_muts_in_unassigned, key=lambda x: x[1], reverse=True):
                print(f"    {mut}: {count}个细胞")
        
        print(f"  - 与backbone有交集的突变: {len(intersecting_muts_in_unassigned)}个")
        if intersecting_muts_in_unassigned:
            for mut, count in sorted(intersecting_muts_in_unassigned, key=lambda x: x[1], reverse=True)[:10]:
                print(f"    {mut}: {count}个细胞")
        
        print(f"  - 与backbone无交集的孤立突变: {len(isolated_muts_in_unassigned)}个")
        if isolated_muts_in_unassigned:
            print(f"    (显示前20个)")
            for mut, count in sorted(isolated_muts_in_unassigned, key=lambda x: x[1], reverse=True)[:20]:
                print(f"    {mut}: {count}个细胞")
        
        # 验证：检查有多少未分配细胞只有孤立突变
        cells_with_only_isolated_muts = 0
        cells_with_backbone_or_intersecting = 0
        
        for cell, muts in unassigned_cells_mutations.items():
            has_backbone_or_intersecting = any(mut in mutations_with_intersection for mut in muts)
            if has_backbone_or_intersecting:
                cells_with_backbone_or_intersecting += 1
            else:
                cells_with_only_isolated_muts += 1
        
        print(f"\n未分配细胞验证:")
        print(f"  - 只有孤立突变的细胞: {cells_with_only_isolated_muts}")
        print(f"  - 有backbone或交集突变的细胞: {cells_with_backbone_or_intersecting}")
        
        # 如果有细胞有backbone或交集突变但仍未分配，需要进一步调查
        if cells_with_backbone_or_intersecting > 0:
            print(f"\n警告: 有 {cells_with_backbone_or_intersecting} 个细胞含有backbone或交集突变但仍未分配")
            print("这可能是由于imputation逻辑问题，需要检查:")
            
            problematic_cells = []
            for cell, muts in unassigned_cells_mutations.items():
                if any(mut in mutations_with_intersection for mut in muts):
                    problematic_cells.append((cell, muts))
            
            # 显示前几个问题细胞作为示例
            for cell, muts in problematic_cells[:3]:
                backbone_muts_in_cell = [mut for mut in muts if mut in backbone_mutations]
                intersecting_muts_in_cell = [mut for mut in muts if mut in mutations_with_intersection and mut not in backbone_mutations]
                print(f"  细胞 {cell}:")
                print(f"    Backbone mutations: {backbone_muts_in_cell}")
                print(f"    Intersecting mutations: {intersecting_muts_in_cell[:5]}...")  # 只显示前5个
    
    # 详细验证
    total_assigned_cells = exclusive_backbone_matrix.sum().sum()
    total_input_cells = len(I_somatic.index)
    total_mutated_cells = len(mutated_cells)
    total_final_assignment_cells = len(final_assignment)
    
    print(f"\n=== 最终验证 ===")
    print(f"原始输入细胞数: {total_input_cells}")
    print(f"最终分配向量细胞数: {total_final_assignment_cells}")
    print(f"有突变的细胞数: {total_mutated_cells}")
    print(f"无突变的细胞数: {len(non_mutated_cells)}")
    print(f"分配到backbone的细胞数: {total_assigned_cells}")
    print(f"有突变但未分配的细胞数: {conflict_resolution_stats['mutated_but_unassigned']}")
    print(f"分配率: {total_assigned_cells/total_mutated_cells*100:.2f}%")
    
    # 修改断言为更友好的错误信息
    if total_input_cells != total_final_assignment_cells:
        raise ValueError(f"细胞数量不匹配! 输入{total_input_cells} != 最终分配{total_final_assignment_cells}")
    
    if total_mutated_cells != total_assigned_cells + conflict_resolution_stats['mutated_but_unassigned']:
        raise ValueError(f"有突变细胞数不匹配! {total_mutated_cells} != {total_assigned_cells} + {conflict_resolution_stats['mutated_but_unassigned']}")
    
    print("\n冲突解决统计:")
    for stat, count in sorted(conflict_resolution_stats.items()):
        if count > 0:
            print(f"  {stat}: {count} cells")
    
    # 输出每个backbone clone的最终大小
    print("\n最终backbone clone大小:")
    total_backbone_cells = 0
    for backbone_mut in backbone_mutations:
        if backbone_mut in exclusive_backbone_matrix.columns:
            clone_size = exclusive_backbone_matrix[backbone_mut].sum()
            intersecting_count = len(intersection_sets.get(backbone_mut, []))
            print(f"  {backbone_mut}: {clone_size} cells, {intersecting_count} intersecting mutations")
            total_backbone_cells += clone_size
    
    print(f"所有backbone clone总细胞数: {total_backbone_cells}")
    print(f"有突变但未分配到任何backbone的细胞: {conflict_resolution_stats['mutated_but_unassigned']}")
    
    return exclusive_backbone_matrix, final_assignment, intersection_sets




# -------------------------
# merge_identical_columns(matrix)
# 输入：原始 mutation-by-cell 矩阵（每列一个 mutation）
# 输出：合并了完全相同列的矩阵，列名变成 "mut1|mut2|mut3"
# 同时返回一个映射字典，用于逆向拆分
# split_merged_columns(matrix, mapping)
# 输入：合并后的矩阵 + 映射字典
# 输出：恢复成原始 mutation-by-cell 矩阵（每列一个 mutation）
# -------------------------
def merge_duplicate_columns(df: pd.DataFrame) -> Tuple[pd.DataFrame, Dict[str,str]]:
    """
    将原始 df（cells x mutations）中完全相同的列合并成 'a|b'。
    返回 (merged_df, col_map)：
      - merged_df: cells x merged_cols
      - col_map: 原始_col -> merged_col_name
    注意：if values contain pd.NA/np.nan we treat them as None in the signature.
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
        # choose first column's data as representative
        merged_cols[merged_name] = df[cols[0]].copy()
        for c in cols:
            col_map[c] = merged_name
    merged_df = pd.DataFrame(merged_cols, index=df.index)
    return merged_df

def split_merged_columns(merged_matrix: pd.DataFrame, mut_list: list):
    """
    根据mut_list拆分合并的列
    """
    out = {}
    for mut in mut_list:
        # 找到包含该mutation的合并列
        found = False
        for merged_col in merged_matrix.columns:
            if mut in merged_col.split("|"):
                out[mut] = merged_matrix[merged_col].copy()
                found = True
                break
        if not found:
            # 如果没找到，创建全NA列
            out[mut] = pd.Series([pd.NA] * len(merged_matrix), index=merged_matrix.index)
    
    return pd.DataFrame(out, index=merged_matrix.index)[mut_list]

# # 测试数据
# df = pd.DataFrame({
#     'A': [1, 2, np.nan],
#     'B': [1, 2, np.nan],  # 与A相同
#     'C': [1, 2, None],    # 与A相同
#     'D': [1, 2, 3]        # 不同
# })
# mut_list = ['A', 'B', 'C', 'D']
# # 执行
# merged = merge_duplicate_columns(df)
# recovered = split_merged_columns(merged, mut_list)




# -------------------------
# find potiential positions that is intersected with new_mut based on backbone tree for scaffold tree
#   - intersection_nodes
#   - parent_dict
# -------------------------

# 配置日志
logger = logging.getLogger(__name__)

def find_intersection_positions_within_group_directly(T_current: TreeNode, new_mut: str, matrix, mutation_group, min_overlap=1):
    """
    基于交集分析的优化版本，直接找到相关位置，只考虑与当前 mutation 同一组的 clone 下的节点
    """
    # 1. 找到所有与目标突变有交集的节点
    intersection_nodes = find_all_intersect_muts_from_tree_by_matrix_scaffold(
        T_current, matrix, new_mut, min_overlap
    )
    
    logger.info(f"Found {len(intersection_nodes)} intersection nodes for {new_mut}: {intersection_nodes}")
    
    if len(intersection_nodes) == 0:
        logger.debug(f"No intersection nodes found for {new_mut}")
        return []  # 没有交集，返回空列表
    
    # 2. 构建树的父子关系字典
    tree_parent_dict = build_tree_parent_dict_scaffold(T_current)
    
    # 3. 获取新突变所在的组
    target_group = mutation_group[new_mut]
    
    # 4. 找到所有相关路径上的节点，仅限于同一组的 clone 路径
    all_path_nodes = get_all_path_nodes_with_group_filter(intersection_nodes, tree_parent_dict, mutation_group, target_group)
    
    logger.info(f"Found {len(all_path_nodes)} path nodes for {new_mut} in the same mutation group {target_group}")
    
    # 5. 预先创建基础树的深拷贝
    base_tree_copy = deepcopy(T_current)
    
    # 6. 只在这些相关节点上生成候选位置
    candidate_positions = []
    
    for node_name in all_path_nodes:
        if node_name == "ROOT":
            continue
            
        node = base_tree_copy.find(node_name)
        if node is None:
            logger.warning(f"Node {node_name} not found in tree")
            continue
        
        # --- 1) 放在 node 上 ---
        candidate_positions.append(_create_on_node_candidate_fast_scaffold(base_tree_copy, node, new_mut))
        
        # --- 2) 新 leaf ---
        candidate_positions.append(_create_new_leaf_candidate_fast_scaffold(base_tree_copy, node, new_mut))
        
        # --- 3) 放在每条 edge 上 ---
        for child in node.children:
            if child.name in all_path_nodes:  # 只考虑路径上的子节点
                candidate_positions.append(_create_on_edge_candidate_fast_scaffold(base_tree_copy, node, child, new_mut))
        
        # --- 4) 新 parent merge 多个子节点 ---
        if len(node.children) >= 2:
            # 只考虑路径上的子节点组合
            path_children = [child for child in node.children if child.name in all_path_nodes]
            if len(path_children) >= 2:
                # 限制组合数量避免爆炸
                for r in range(2, min(4, len(path_children) + 1)):
                    for combo in combinations(path_children, r):
                        candidate_positions.append(_create_merge_candidate_fast_scaffold(base_tree_copy, node, combo, new_mut))
    
    logger.info(f"Generated {len(candidate_positions)} candidate positions for {new_mut}")
    return candidate_positions


def find_new_leaf_positions_for_target_node(T_current: TreeNode, new_mut: str, matrix, target_node, min_overlap=1):
    """
    基于交集分析的优化版本，但只返回与目标节点相关的 new_leaf position
    """
    # 1. 找到所有与目标突变有交集的节点
    intersection_nodes = find_all_intersect_muts_from_tree_by_matrix_scaffold(
        T_current, matrix, new_mut, min_overlap
    )
    
    logger.info(f"Found {len(intersection_nodes)} intersection nodes for {new_mut}: {intersection_nodes}")
    
    if len(intersection_nodes) == 0:
        logger.debug(f"No intersection nodes found for {new_mut}")
        return []  # 没有交集，返回空列表
    
    # 2. 构建树的父子关系字典
    tree_parent_dict = build_tree_parent_dict_scaffold(T_current)
    
    # 3. 找到所有相关路径上的节点
    all_path_nodes = find_all_path_nodes_scaffold(intersection_nodes, tree_parent_dict)
    
    logger.info(f"Found {len(all_path_nodes)} path nodes for {new_mut}")
    
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
            candidate_positions.append(_create_on_node_candidate_fast_scaffold(base_tree_copy, node, new_mut))
            
            continue  # ROOT 的其他类型位置仍然跳过
        
        node = base_tree_copy.find(node_name)
        if node is None:
            logger.warning(f"Node {node_name} not found in tree")
            continue
        
        # --- 1) 放在 node 上 ---
        candidate_positions.append(_create_on_node_candidate_fast_scaffold(base_tree_copy, node, new_mut))
        
        # --- 2) 新 leaf ---
        candidate_positions.append(_create_new_leaf_candidate_fast_scaffold(base_tree_copy, node, new_mut))
        
        # --- 3) 放在每条 edge 上 ---
        for child in node.children:
            if child.name in all_path_nodes:  # 只考虑路径上的子节点
                candidate_positions.append(_create_on_edge_candidate_fast_scaffold(base_tree_copy, node, child, new_mut))
        
        # --- 4) 新 parent merge 多个子节点 ---
        if len(node.children) >= 2:
            # 只考虑路径上的子节点组合
            path_children = [child for child in node.children if child.name in all_path_nodes]
            if len(path_children) >= 2:
                # 限制组合数量避免爆炸
                for r in range(2, min(4, len(path_children) + 1)):
                    for combo in combinations(path_children, r):
                        candidate_positions.append(_create_merge_candidate_fast_scaffold(base_tree_copy, node, combo, new_mut))
    
    logger.info(f"Generated {len(candidate_positions)} candidate positions for {new_mut}")
    
    # 6. 筛选出 anchor 是 target_node 且 placement_type 是 new_leaf 的 position
    target_positions = [
        pos for pos in candidate_positions 
        if pos['placement_type'] == 'new_leaf' and pos['anchor'] == target_node.name
    ]
    
    logger.info(f"Filtered to {len(target_positions)} target positions for {new_mut} under node {target_node.name}")
    
    return target_positions


def find_intersection_positions_within_tree_directly_scaffold(T_current: TreeNode, new_mut: str, matrix, min_overlap=1):
    """
    基于交集分析的优化版本，直接找到相关位置
    """
    # 1. 找到所有与目标突变有交集的节点
    intersection_nodes = find_all_intersect_muts_from_tree_by_matrix_scaffold(
        T_current, matrix, new_mut, min_overlap
    )
    
    logger.info(f"Found {len(intersection_nodes)} intersection nodes for {new_mut}: {intersection_nodes}")
    
    if len(intersection_nodes) == 0:
        logger.debug(f"No intersection nodes found for {new_mut}")
        return []  # 没有交集，返回空列表
    
    # 2. 构建树的父子关系字典
    tree_parent_dict = build_tree_parent_dict_scaffold(T_current)
    
    # 3. 找到所有相关路径上的节点
    all_path_nodes = find_all_path_nodes_scaffold(intersection_nodes, tree_parent_dict)
    
    logger.info(f"Found {len(all_path_nodes)} path nodes for {new_mut}")
    
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
            candidate_positions.append(_create_on_node_candidate_fast_scaffold(base_tree_copy, node, new_mut))
            
            continue  # ROOT 的其他类型位置仍然跳过
        
        node = base_tree_copy.find(node_name)
        if node is None:
            logger.warning(f"Node {node_name} not found in tree")
            continue
        
        # --- 1) 放在 node 上 ---
        candidate_positions.append(_create_on_node_candidate_fast_scaffold(base_tree_copy, node, new_mut))
        
        # --- 2) 新 leaf ---
        candidate_positions.append(_create_new_leaf_candidate_fast_scaffold(base_tree_copy, node, new_mut))
        
        # --- 3) 放在每条 edge 上 ---
        for child in node.children:
            if child.name in all_path_nodes:  # 只考虑路径上的子节点
                candidate_positions.append(_create_on_edge_candidate_fast_scaffold(base_tree_copy, node, child, new_mut))
        
        # --- 4) 新 parent merge 多个子节点 ---
        if len(node.children) >= 2:
            # 只考虑路径上的子节点组合
            path_children = [child for child in node.children if child.name in all_path_nodes]
            if len(path_children) >= 2:
                # 限制组合数量避免爆炸
                for r in range(2, min(4, len(path_children) + 1)):
                    for combo in combinations(path_children, r):
                        candidate_positions.append(_create_merge_candidate_fast_scaffold(base_tree_copy, node, combo, new_mut))
    
    logger.info(f"Generated {len(candidate_positions)} candidate positions for {new_mut}")
    return candidate_positions

def filter_candidate_positions_from_target_node(candidate_positions, subtree_nodes):
    """
    筛选出anchor在指定子树节点集合中的候选位置
    
    Args:
        candidate_positions: 候选位置列表
        subtree_nodes: 子树节点名称集合
    
    Returns:
        list: 筛选后的候选位置列表
    """
    filtered_positions = []
    
    for pos in candidate_positions:
        # 检查这个候选位置的anchor是否在子树节点集合中
        if pos.get('anchor') in subtree_nodes:
            filtered_positions.append(pos)
    
    print(f"从 {len(candidate_positions)} 个候选位置中筛选出 {len(filtered_positions)} 个在目标子树内的位置")
    return filtered_positions




def get_all_path_nodes_with_group_filter(intersection_nodes, tree_parent_dict, mutation_group, target_group):
    """
    获取与交集节点相关的路径节点，但仅限于 mutation_group 中属于 target_group 的 clone 下的路径
    """
    all_path_nodes = set()
    all_path_nodes.add('ROOT')  # 总是包含 ROOT
    
    for node in intersection_nodes:
        path = get_path_to_root_scaffold(node, tree_parent_dict)
        all_path_nodes.update(path)
    
    # 交集节点之间的路径
    intersection_list = list(intersection_nodes)
    for i in range(len(intersection_list)):
        for j in range(i + 1, len(intersection_list)):
            path_between = get_path_between_nodes_scaffold(intersection_list[i], intersection_list[j], tree_parent_dict)
            all_path_nodes.update(path_between)
    
    # 只保留属于 target_group 的节点路径
    all_path_nodes = {node for node in all_path_nodes if valid_node_with_group(node, mutation_group, target_group)}
    
    return all_path_nodes


def valid_node_with_group(node, mutation_group, target_group):
    """
    判断一个节点是否属于目标 mutation_group 的 group
    """
    node_muts = node.split("|")
    for mut in node_muts:
        if mut in mutation_group and mutation_group[mut] == target_group:
            return True
    return False


def build_tree_parent_dict_scaffold(tree):
    """直接从树构建父子关系字典"""
    parent_dict = {}
    for node in tree.traverse():
        for child in node.children:
            parent_dict[child.name] = node.name
    return parent_dict


def _create_on_node_candidate_fast_scaffold(base_tree_copy, node, new_mut):
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
        "nodes": _extract_nodes_info_sacffold(new_tree),
        "edges": _extract_edges_info_scaffold(new_tree)
    }


def _create_new_leaf_candidate_fast_scaffold(base_tree_copy, node, new_mut):
    """快速创建新叶子的候选"""
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
        "nodes": _extract_nodes_info_sacffold(new_tree),
        "edges": _extract_edges_info_scaffold(new_tree)
    }


def _create_merge_candidate_fast_scaffold(base_tree_copy, parent_node, children_combo, new_mut):
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
        "nodes": _extract_nodes_info_sacffold(new_tree),
        "edges": _extract_edges_info_scaffold(new_tree)
    }


def _extract_nodes_info_sacffold(tree):
    """提取节点信息"""
    return [{"name": n.name,
             "parent": n.parent.name if n.parent else None,
             "children": [c.name for c in n.children]} 
            for n in tree.traverse()]


def _extract_edges_info_scaffold(tree):
    """提取边信息"""
    return [(n.parent.name, n.name) for n in tree.traverse() if n.parent]


def build_parent_dict_from_candidates_scaffold(candidate_positions):
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


def build_lineage_parent_dict_from_tree(tree, anchor):
    """从当前树中构建从anchor到ROOT的lineage父子关系
    
    返回格式:
    {'chr17_7389869_G_A': 'chr19_18929325_C_T', 'chr7_142870944_G_A': 'chr15_74920412_G_A', ...}
    """
    """从当前树中构建从anchor到ROOT的lineage父子关系，但去掉ROOT"""
    lineage_parent_dict = {}
    
    def find_path(node, target, current_path):
        if node.name == target:
            # 构建路径上的父子关系
            for i in range(len(current_path)-1, 0, -1):
                lineage_parent_dict[current_path[i]] = current_path[i-1]
            if current_path:  # 确保anchor本身也在关系中
                lineage_parent_dict[target] = current_path[-1] if current_path else None
            return True
        
        for child in node.children:
            if find_path(child, target, current_path + [node.name]):
                return True
        return False
    
    # 从ROOT开始查找
    find_path(tree, anchor, [])
    
    # 过滤掉包含ROOT的键值对
    filtered_dict = {child: parent for child, parent in lineage_parent_dict.items() 
                    if child != "ROOT" and parent != "ROOT"}
    
    return filtered_dict


# def find_all_intersect_muts_from_tree_by_matrix_scaffold(tree, matrix, target_mut, min_overlap=1):
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

def find_all_intersect_muts_from_tree_by_matrix_scaffold(tree, matrix, target_mut, min_overlap=1):
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

def find_all_path_nodes_scaffold(intersection_nodes, tree_parent_dict):
    """
    找到所有在连接交集节点的路径上的节点
    """
    all_path_nodes = set()
    all_path_nodes.add('ROOT')  # 总是包含 ROOT
    
    # 对于每个交集节点，找到从 ROOT 到该节点的路径
    for node in intersection_nodes:
        path_to_root = get_path_to_root_scaffold(node, tree_parent_dict)
        all_path_nodes.update(path_to_root)
    
    # 找到连接不同交集节点的路径
    intersection_list = list(intersection_nodes)
    for i in range(len(intersection_list)):
        for j in range(i + 1, len(intersection_list)):
            path_between = get_path_between_nodes_scaffold(intersection_list[i], intersection_list[j], tree_parent_dict)
            all_path_nodes.update(path_between)
    
    return all_path_nodes


def get_path_to_root_scaffold(node, tree_parent_dict):
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


def get_path_between_nodes_scaffold(node1, node2, tree_parent_dict):
    """找到两个节点之间的路径"""
    # 找到从 node1 到 ROOT 的路径
    path1 = get_path_to_root_scaffold(node1, tree_parent_dict)
    # 找到从 node2 到 ROOT 的路径  
    path2 = get_path_to_root_scaffold(node2, tree_parent_dict)
    
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
# penalty function
# 期望基因型 G 的生成：给 placement type 返回 G (cells -> 0/1)
# placement_type: 'on_node' | 'as_child' | 'on_edge' | 'new_parent_merge' | 'new_leaf'
# 简化策略（可替换为更严格的贝叶斯推断）：
#   - on_node: 期望为 node 子树内任一 mutation 存在的细胞（用 posterior>0.5 判断）
#   - as_child: 期望为 node 子树内的子集（这里我们用 node 子树的细胞作为 superset，再由后验权重调整）
#   - on_edge: 同 on_node（edge 放置后视为分割 child subset）
#   - new_parent_merge: union of children's cell sets
#   - new_leaf: subset of parent subtree (approx: empty set initially)
# -------------------------

def select_overlapping_with_minimum(results):
    """
    选择所有与最小罚分候选重叠的位置
    """
    if not results:
        return []
    
    # 找到最小penalty_min
    min_penalty = min(r['penalty_min'] for r in results)
    
    # 找到所有penalty_min等于最小值的候选
    min_candidates = [r for r in results if r['penalty_min'] == min_penalty]
    
    # 计算这些最小候选的最大penalty_max
    min_group_max = max(r['penalty_max'] for r in min_candidates)
    
    # 选择所有与最小候选组重叠的位置
    selected = []
    for result in results:
        # 如果候选的penalty_min <= 最小候选组的最大penalty_max，则选中
        if result['penalty_min'] <= min_group_max:
            selected.append(result['position'])
    
    return selected

def compute_binary_penalty_for_positions(new_mut, refined_positions, M_current, I_selected):
    """
    计算 new_mut 在所有 candidate positions 的 binary penalty 上下限。
    
    Parameters
    ----------
    new_mut : str
        待加入树的 mutation
    refined_positions : list of dict
        候选位置列表，每个 dict 含有 placement_type, anchor, nodes, edges 等信息
    M_current : pd.DataFrame
        backbone clone matrix (0/1)
    I_selected : pd.DataFrame
        observed genotype matrix (含NA)
    
    Returns
    -------
    penalty_df : pd.DataFrame
        每个位置的 penalty min/max
    """
    def compute_N_vectors(imputed_vec, observed_vec):
        """计算 N_10 和 N_01"""
        imputed_vec = imputed_vec.astype(int)
        obs_vec = observed_vec.copy()
        obs_vec = obs_vec.fillna(-1).astype(int)
        
        N_10 = np.sum((obs_vec == 1) & (imputed_vec == 0))
        N_01 = np.sum((obs_vec == 0) & (imputed_vec == 1))
        return N_10, N_01
    
    results = []
    new_mut_bin_vector = I_selected[new_mut].replace({pd.NA: np.nan}).fillna(0).astype(int)
    
    for idx, pos in enumerate(refined_positions):
        placement_type = pos['placement_type']
        anchor = pos['anchor']
        
        # 生成 imputed vector
        imputed_vec = pd.Series(0, index=M_current.index)  # 默认全 0
        
        if placement_type == 'on_node':
            # new_mut 跟 anchor node 的 cells 一致
            if anchor not in M_current.columns:
                raise ValueError(f"Anchor {anchor} not in M_current")
            imputed_vec = M_current[anchor]
        
        elif placement_type == 'new_leaf':
            # 获取父节点和子节点的信息
            parent = anchor
            if parent == 'ROOT':
                vec_parent = pd.Series(1, index=M_current.index)  # 根节点是 1
            else:
                if parent not in M_current.columns:
                    raise ValueError(f"Parent {parent} not in M_current")
                vec_parent = M_current[parent]
            
            # 新突变的交集，确保 new_mut 只在父节点下为 1
            imputed_vec = new_mut_bin_vector & vec_parent
            
            # 获取子节点信息
            child_nodes = [node['name'] for node in pos['nodes'] if node['parent'] == parent]
            child_nodes = [i for i in child_nodes if i != new_mut]
            
            # 确保 new_mut 与子节点互斥
            for child in child_nodes:
                if child not in M_current.columns:
                    raise ValueError(f"Child {child} not in M_current")
                vec_child = M_current[child]
                imputed_vec &= ~vec_child  # 与子节点互斥
            
            # 确保只有在 new_mut 为 1 的地方，父节点链上才为 1
            imputed_vec &= new_mut_bin_vector
            
        elif placement_type == 'on_edge':
            parent = anchor
            child = pos['meta']['child']
            
            # 处理 parent 向量
            if parent == 'ROOT':
                vec_parent = pd.Series(1, index=M_current.index)
            else:
                if parent not in M_current.columns:
                    raise ValueError(f"Parent {anchor} not in M_current")
                vec_parent = M_current[anchor]
            
            # 处理目标 child 向量
            if child not in M_current.columns:
                raise ValueError(f"Child {child} not in M_current")
            
            vec_child = M_current[child]
            
            # 获取 parent 下的所有其他 child (siblings)
            sibling_nodes = [node['name'] for node in pos['nodes'] 
                            if node['parent'] == parent and node['name'] != child]
            sibling_nodes = [i for i in sibling_nodes if i != new_mut]
            
            # 初始化 imputed_vec：包含目标 child，且在 parent 范围内
            imputed_vec = vec_child.copy()
            
            # 添加新突变在 parent 范围内的部分
            imputed_vec = (imputed_vec | (new_mut_bin_vector & vec_parent)).astype(int)
            
            # 关键：排除所有 sibling 节点的细胞（互斥）
            for sibling in sibling_nodes:
                if sibling not in M_current.columns:
                    raise ValueError(f"Sibling {sibling} not in M_current")
                vec_sibling = M_current[sibling]
                imputed_vec = (imputed_vec & ~vec_sibling).astype(int)
            
            # 最终确保在 parent 范围内
            imputed_vec = (imputed_vec & vec_parent).astype(int)
        
        elif placement_type == 'new_parent_merge':
            parent = anchor
            merge_children = pos['meta']['merge_children']
            
            # 处理 parent 向量
            if parent == 'ROOT':
                vec_parent = pd.Series(1, index=M_current.index)
            else:
                if parent not in M_current.columns:
                    raise ValueError(f"Parent {parent} not in M_current for edge placement")
                vec_parent = M_current[parent]
            
            # 处理要合并的 children 向量
            vec_children = pd.Series(0, index=M_current.index)
            for c in merge_children:
                if c not in M_current.columns:
                    raise ValueError(f"Child {c} not in M_current for edge placement")
                vec_children = ((vec_children == 1) | (M_current[c] == 1)).astype(int)
            
            # 初始化 imputed_vec：包含要合并的 children，且在 parent 范围内
            imputed_vec = vec_children.copy()
            
            # 添加新突变在 parent 范围内的部分
            imputed_vec = (imputed_vec | (new_mut_bin_vector & vec_parent)).astype(int)
            
            # 关键：排除所有其他 sibling 节点的细胞（互斥）
            # 获取 parent 下的所有其他 child (siblings)，不包括要合并的 children
            sibling_nodes = [node['name'] for node in pos['nodes'] 
                            if node['parent'] == parent and node['name'] not in merge_children]
            sibling_nodes = [i for i in sibling_nodes if i != new_mut]
            
            for sibling in sibling_nodes:
                if sibling not in M_current.columns:
                    raise ValueError(f"Sibling {sibling} not in M_current")
                vec_sibling = M_current[sibling]
                imputed_vec = (imputed_vec & ~vec_sibling).astype(int)
            
            # 最终确保在 parent 范围内
            imputed_vec = (imputed_vec & vec_parent).astype(int)
        
        else:
            raise ValueError(f"Unknown placement_type: {placement_type}")
        
        # 对应 observed vector
        if new_mut not in I_selected.columns:
            observed_vec = pd.Series(np.nan, index=M_current.index)  # 如果没有观测数据，全 NA
        else:
            observed_vec = I_selected[new_mut]
        
        N_10, N_01 = compute_N_vectors(imputed_vec, observed_vec)
        
        # 计算 penalty
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
    
    # 选择所有可能最优的位置。 筛选逻辑：选择下限较低且上限不超过其它位置下限的候选位置
    selected_positions = select_overlapping_with_minimum(results)
    
    return selected_positions, df_penalty

# selected_positions, df_penalty = compute_binary_penalty_for_positions('chr17_7389869_G_A', refined_positions, M_current, I_selected)
# df_penalty[['position_index', 'placement_type', 'anchor', 'N_10', 'N_01', 'penalty_min', 'penalty_max']]
# #    position_index placement_type             anchor  N_10  N_01  penalty_min  penalty_max
# # 0               0        on_node  chr1_39034563_T_A     0    25         1.25          2.5
# # 1               1       new_leaf  chr1_39034563_T_A     0     0         0.00          0.0


##### Bayesian penalty
def compute_bayesian_penalty_each_pos(input_binary_vec, posterior_vec, imputed_vec, fnfp_ratio=0.1, omega_NA=0.001):
    """
    基于数据先验分布的贝叶斯罚分计算
    """
    # 1. 计算先验概率
    non_nan_values = input_binary_vec[~np.isnan(input_binary_vec)]
    count_1 = np.sum(non_nan_values == 1)
    count_0 = np.sum(non_nan_values == 0)
    total_count = count_1 + count_0
    
    if total_count > 0:
        p_mutation = count_1 / total_count      # 突变先验概率
        p_wildtype = count_0 / total_count      # 野生型先验概率
    else:
        p_mutation = p_wildtype = 0.5
    
    # 2. 基于先验概率设置权重
    # NA→1惩罚：与野生型概率成正比（野生型越多，NA→1越不合理）
    weight_na_to_1 = omega_NA * p_wildtype  
    # NA→0惩罚：与突变概率成正比（突变越多，NA→0越不合理）
    weight_na_to_0 = omega_NA * p_mutation
    
    # 3. 计算各种情况的计数
    fp_count = np.sum((input_binary_vec == 1) & (imputed_vec == 0))
    fn_count = np.sum((input_binary_vec == 0) & (imputed_vec == 1))
    na_to_1_count = np.sum(np.isnan(input_binary_vec) & (imputed_vec == 1))
    na_to_0_count = np.sum(np.isnan(input_binary_vec) & (imputed_vec == 0))
    
    # 4. 计算罚分
    fp_penalty = fp_count
    fn_penalty = fn_count * fnfp_ratio
    na_to_1_penalty = na_to_1_count * weight_na_to_1
    na_to_0_penalty = na_to_0_count * weight_na_to_0
    
    # 5. 增强计算 FP 罚分
    
    total_penalty = fp_penalty + fn_penalty + na_to_1_penalty + na_to_0_penalty
    
    return total_penalty, weight_na_to_1, weight_na_to_0


def compute_bayesian_penalty_each_chain_mut_by_pos(input_binary_vec, posterior_vec, imputed_vec, fnfp_ratio, weight_na_to_1, weight_na_to_0):
    """
    基于数据先验分布的贝叶斯罚分计算
    """
    # 1-2. 直接导入 weight_na_to_1 和 weight_na_to_0
    
    # 3. 计算各种情况的计数
    fp_count = np.sum((input_binary_vec == 1) & (imputed_vec == 0))
    fn_count = np.sum((input_binary_vec == 0) & (imputed_vec == 1))
    na_to_1_count = np.sum(np.isnan(input_binary_vec) & (imputed_vec == 1))
    na_to_0_count = np.sum(np.isnan(input_binary_vec) & (imputed_vec == 0))
    
    # 4. 计算罚分
    fp_penalty = fp_count
    fn_penalty = fn_count * fnfp_ratio
    na_to_1_penalty = na_to_1_count * weight_na_to_1
    na_to_0_penalty = na_to_0_count * weight_na_to_0
    
    # 5. 增强计算 FP 罚分
    
    total_penalty = fp_penalty + fn_penalty + na_to_1_penalty + na_to_0_penalty
    
    return total_penalty


def compute_bayesian_penalty_for_positions_scaffold(
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
            lineage_conflict_nodes = get_all_conflict_nodes_outside_lineage_scaffold(parent, build_lineage_parent_dict_from_tree(T_current, anchor), M_current.columns)
            
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
            lineage_conflict_nodes = get_all_conflict_nodes_outside_lineage_scaffold(parent, build_lineage_parent_dict_from_tree(T_current, anchor), M_current.columns)
            
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
            lineage_conflict_nodes = get_all_conflict_nodes_outside_lineage_scaffold(parent, build_lineage_parent_dict_from_tree(T_current, anchor), M_current.columns)
            
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
            lineage_conflict_nodes = get_all_conflict_nodes_outside_lineage_scaffold(parent, build_lineage_parent_dict_from_tree(T_current, anchor), M_current.columns, exclude_nodes=merge_children)
            
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
        full_mutnode_chain = get_full_mutnode_chain_with_anchor_scaffold(anchor, parent_dict)
        
        # 计算new_mut本身的罚分
        posterior_vec = P_selected[new_mut]
        input_binary_vec = I_selected[new_mut]
        
        # 基于实际NA翻转比例计算惩罚
        new_mut_penalty, actual_na_flip_ratio, refined_ω_NA, φ_adjusted, weight_na_to_1, weight_na_to_0 = compute_dynamic_penalty_scaffold(
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
        BIC_penalty = φ_adjusted * np.log(N_nodes)  # 使用动态调整后的φ
        root_penalty = 0
        if anchor == 'ROOT':
            root_penalty = np.log(N_nodes) * 0.5
        
        # 基础总罚分
        base_total_penalty = total_chain_penalty + log_N_nodes_penalty + BIC_penalty + merge_penalty + root_penalty
        
        # -------------------------
        # 添加基于intersection模式和层级的精细调整
        # -------------------------
        intersection_penalty = compute_intersection_based_penalty_scaffold(
            new_mut, pos, intersection_nodes, M_current, I_selected, na_ratio, mut_ratio, actual_na_flip_ratio
        )
        
        hierarchy_penalty = compute_hierarchy_penalty_scaffold(
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
                logger.info(f"Skipping position with all imputed_vec values as 0: {row['position']}")
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
        # 修正：改进的M_current更新逻辑
        # -------------------------
        anchor = final_position['anchor']
        
        # 获取从锚点到ROOT的完整突变链
        full_mutnode_chain = get_full_mutnode_chain_with_anchor_scaffold(anchor, parent_dict)
        
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


def get_all_conflict_nodes_outside_lineage_scaffold(anchor, parent_dict, all_columns, exclude_nodes=None):
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


def get_full_mutnode_chain_with_anchor_scaffold(anchor, parent_dict):
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


def compute_dynamic_penalty_scaffold(input_binary_vec, posterior_vec, imputed_vec, fnfp_ratio, base_ω_NA, base_φ, na_ratio, mut_ratio, placement_type, N_nodes):
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
    φ_adjusted = compute_dynamic_bic_penalty_scaffold(
        base_φ, na_ratio, mut_ratio, actual_na_flip_ratio, placement_type, N_nodes
    )
    
    # 使用调整后的权重计算惩罚
    penalty, weight_na_to_1, weight_na_to_0 = compute_bayesian_penalty_each_pos(
        input_binary_vec, posterior_vec, imputed_vec, fnfp_ratio, refined_ω_NA
    )
    
    return penalty, actual_na_flip_ratio, refined_ω_NA, φ_adjusted, weight_na_to_1, weight_na_to_0


def compute_dynamic_bic_penalty_scaffold(base_φ, na_ratio, mut_ratio, actual_na_flip_ratio, placement_type, N_nodes):
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


def compute_intersection_based_penalty_scaffold(new_mut, position, intersection_nodes, M_current, I_selected, na_ratio, mut_ratio, actual_na_flip_ratio):
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


def compute_hierarchy_penalty_scaffold(new_mut, position, M_current, I_selected, parent_dict, 
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
        anchor_children = find_children_of_node_scaffold(anchor, M_current.columns, parent_dict)
        
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


def find_children_of_node_scaffold(node, all_columns, parent_dict):
    """找到节点的直接children"""
    children = []
    for col in all_columns:
        if col == node or col == 'ROOT':
            continue
        if parent_dict.get(col) == node:
            children.append(col)
    return children




def add_new_mutation_to_tree_independent(new_mut, T_current, final_position):
    placement_type = final_position['placement_type']
    anchor = final_position['anchor']
    meta = final_position['meta']
    
    # 获取锚点节点
    anchor_node = T_current.find(anchor)
    if not anchor_node:
        raise ValueError(f"Anchor node {anchor} not found in the tree.")
    
    # 根据 placement_type 选择不同的插入策略
    if placement_type == 'on_node':
        # 如果是 on_node，只需在 anchor 的节点名称后加上 |new_mut
        new_mut_name = new_mut  # 直接用 new_mut 作为突变名
        anchor_node.name = f"{anchor_node.name}|{new_mut_name}"
        # 返回更新后的树结构
        return T_current
    
    elif placement_type == 'new_leaf':
        # 如果是 new_leaf，在 anchor 的节点下添加一个新叶子节点
        new_leaf_name = new_mut  # 直接用 new_mut 作为新叶子节点名称
        new_leaf = TreeNode(new_leaf_name)
        anchor_node.add_child(new_leaf)
        # 返回更新后的树结构
        return T_current
    
    elif placement_type == 'on_edge':
        # 如果是 on_edge，找到 meta 中的 child，然后插入一个新节点
        child_name = meta.get('child')
        child_node = T_current.find(child_name)
        if not child_node:
            raise ValueError(f"Child node {child_name} not found in the tree.")
        
        # 创建新节点，并将其插入锚点和子节点之间
        new_edge_name = new_mut  # 直接用 new_mut 作为新节点名称
        new_node = TreeNode(new_edge_name)
        T_current.insert_on_edge(anchor_node, child_node, new_edge_name)
        # 返回更新后的树结构
        return T_current
    
    elif placement_type == 'new_parent_merge':
        # 如果是 new_parent_merge，将多个子节点合并到一个新父节点下
        children_to_merge = meta.get('merge_children', [])
        children_nodes = [T_current.find(child) for child in children_to_merge]
        if None in children_nodes:
            raise ValueError(f"One or more merge children not found in the tree.")
        
        new_parent_name = new_mut  # 直接用 new_mut 作为新父节点的名称
        new_parent = anchor_node.add_new_parent_for_children(children_nodes, new_parent_name)
        # 返回更新后的树结构
        return T_current
    
    else:
        raise ValueError(f"Unknown placement type: {placement_type}")


def add_new_mutation_to_tree_conflict_free(new_mut, T_current, pos, M_current, parent_dict, I_selected):
    """
    将新突变添加到树中，确保与现有树结构没有冲突
    基于 compute_bayesian_penalty_for_positions_scaffold 中的逻辑
    """
    placement_type = pos['placement_type']
    anchor = pos['anchor']
    meta = pos.get('meta', {})
    
    # 获取锚点节点
    anchor_node = T_current.find(anchor)
    if not anchor_node:
        raise ValueError(f"Anchor node {anchor} not found in the tree.")
    
    # 获取新突变的二进制向量
    new_mut_bin_vector = I_selected[new_mut].replace({pd.NA: np.nan}).fillna(0).astype(int)
    
    # 默认 imputed vector
    imputed_vec = pd.Series(0, index=M_current.index)
    
    # 根据 placement_type 计算无冲突的 imputed_vec
    if placement_type == 'on_node':
        # 在节点上放置，直接使用锚点节点的向量
        imputed_vec = M_current[anchor].astype(int)
        
        # 更新节点名称
        anchor_node.name = f"{anchor_node.name}|{new_mut}"
        
    elif placement_type == 'new_leaf':
        parent = anchor
        vec_parent = M_current[parent] if parent != 'ROOT' else pd.Series(1, index=M_current.index)
        
        # 获取直系sibling冲突节点
        sibling_nodes = [n['name'] for n in pos['nodes'] if n['parent'] == parent and n['name'] != new_mut]
        
        # 获取lineage之外的所有冲突节点
        lineage_conflict_nodes = get_all_conflict_nodes_outside_lineage_scaffold(parent, build_lineage_parent_dict_from_tree(T_current, anchor), M_current.columns)
        
        # 合并所有冲突节点（去重）
        all_conflict_nodes = list(set(sibling_nodes + lineage_conflict_nodes))
        
        # 构建冲突向量
        vec_conflicts = pd.Series(0, index=M_current.index)
        for conflict in all_conflict_nodes:
            vec_conflicts |= M_current[conflict]
        
        # 正确的逻辑：先排除所有冲突，再与parent交集
        new_mut_cleaned = new_mut_bin_vector & ~vec_conflicts
        imputed_vec = new_mut_cleaned.astype(int)
        
        # 添加新叶子节点
        new_leaf = TreeNode(new_mut)
        anchor_node.add_child(new_leaf)
        
    elif placement_type == 'on_edge':
        parent = anchor
        child = pos['meta']['child']
        vec_parent = M_current[parent] if parent != 'ROOT' else pd.Series(1, index=M_current.index)
        vec_child = M_current[child]
        
        # 获取直系sibling冲突节点
        sibling_nodes = [n['name'] for n in pos['nodes'] if n['parent']==parent and n['name'] not in [child,new_mut]]
        
        # 获取lineage之外的所有冲突节点
        lineage_conflict_nodes = get_all_conflict_nodes_outside_lineage_scaffold(parent, build_lineage_parent_dict_from_tree(T_current, anchor), M_current.columns)
        
        # 合并所有冲突节点（去重）
        all_conflict_nodes = list(set(sibling_nodes + lineage_conflict_nodes))
        
        # 构建冲突向量
        vec_conflicts = pd.Series(0, index=M_current.index)
        for conflict in all_conflict_nodes:
            vec_conflicts |= M_current[conflict]
        
        # 正确的逻辑：child ∪ (清理后的new_mut ∩ parent)
        new_mut_cleaned = new_mut_bin_vector & ~vec_conflicts
        imputed_vec = (vec_child | new_mut_cleaned).astype(int)
        
        # 在边上插入新节点
        new_node = TreeNode(new_mut)
        T_current.insert_on_edge(anchor_node, child_node, new_mut)
        
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
        lineage_conflict_nodes = get_all_conflict_nodes_outside_lineage_scaffold(parent, build_lineage_parent_dict_from_tree(T_current, anchor), M_current.columns, exclude_nodes=merge_children)
        
        # 合并所有冲突节点（去重）
        all_conflict_nodes = list(set(sibling_nodes + lineage_conflict_nodes))
        
        # 构建冲突向量
        vec_conflicts = pd.Series(0, index=M_current.index)
        for conflict in all_conflict_nodes:
            vec_conflicts |= M_current[conflict]
        
        # 正确的逻辑：children ∪ (清理后的new_mut ∩ parent)
        new_mut_cleaned = new_mut_bin_vector & ~vec_conflicts
        imputed_vec = (vec_children | new_mut_cleaned).astype(int)
        
        # 创建新父节点并合并子节点
        new_parent = anchor_node.add_new_parent_for_children(children_nodes, new_mut)
        
    else:
        raise ValueError(f"Unknown placement type: {placement_type}")
    
    # 更新突变矩阵 M_current
    # 首先添加新列
    M_current[new_mut] = imputed_vec
    
    # 对于在imputed_vec中为1的细胞，确保其父节点链也为1
    if placement_type != 'on_node':  # on_node 不需要这个，因为已经继承了父节点的状态
        full_mutnode_chain = get_full_mutnode_chain_with_anchor_scaffold(anchor, parent_dict)
        cells_with_final_one = imputed_vec[imputed_vec == 1].index.tolist()
        
        if len(cells_with_final_one) > 0:
            for cell in cells_with_final_one:
                for mutation in full_mutnode_chain:
                    if M_current.loc[cell, mutation] == 0:
                        M_current.loc[cell, mutation] = 1
    
    return T_current, M_current, imputed_vec


def get_all_conflict_nodes_outside_lineage_scaffold(anchor, parent_dict, all_columns, exclude_nodes=None):
    """
    获取lineage之外的所有可能冲突节点
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


def WriteTfile(out_prefix, matrix, rownames, colnames, judge): # writes matrix output as an integer matrix, uses the format given in the documentary. 
    matrix_output = matrix.astype(int)
    df_output = pd.DataFrame(matrix_output)
    df_output.index = rownames
    df_output.columns = colnames
    df_output.index.name = "cellIDxmutID"
    is_cf = scp.ul.is_conflict_free_gusfield(df_output)
    if judge == "yes":
        if is_cf:
            print("Current tree is conflict-free and output binary matrix and plot phylogenetic tree !")
            df_output.to_csv(out_prefix+".CFMatrix", sep="\t")
            tree = scp.ul.to_tree(df_output)
            scp.pl.clonal_tree(tree, output_file=out_prefix+".tree_scphylo.pdf")
        else:
            print("Current tree is not conflict-free !")
    else:
        print("Only output binary matrix !")
        df_output.to_csv(out_prefix+".CFMatrix", sep="\t")

def find_flipping_spots(series_in_bin, series_phylogeny, condition_in_bin, condition_phylogeny):
    """
    Find a list of eligible line names (Spots)
    """
    return series_in_bin[(series_in_bin == condition_in_bin) & (series_phylogeny == condition_phylogeny)].index.tolist()


def integrate_mutations_to_scaffold_within_group(sorted_attached_mutations, T_current, M_current, I_attached, P_attached, 
                                   mutation_group, ω_NA, fnfp_ratio, φ, logger):
    """
    处理外部突变并将其整合到scaffold进化树中
    
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
    mutation_group : list
        突变分组信息
    ω_NA : float
        NA值的权重参数
    fnfp_ratio : float
        假阴性假阳性比率
    φ : float
        贝叶斯罚分参数
    logger : logging.Logger
        日志记录器
    
    Returns:
    --------
    tuple : (external_mutations, T_current, M_current)
        处理后的外部突变列表、更新后的树、更新后的矩阵
    """
    
    external_mutations = []
    
    for new_mut in tqdm(sorted_attached_mutations, desc="Processing mutations", unit="mutation"):
        logger.info(f"Processing mutation: {new_mut}")
        
        # 找到交集节点
        intersection_nodes = find_all_intersect_muts_from_tree_by_matrix_scaffold(T_current, I_attached, new_mut)
        if len(intersection_nodes) == 0:
            external_mutations.append(new_mut)
            logger.info(f"Mutation {new_mut} added to external_mutations (no intersection found)")
            continue
        
        # 使用优化方法获取候选位置
        refined_positions = find_intersection_positions_within_group_directly(T_current, new_mut, I_attached, mutation_group, min_overlap=1)
        parent_dict = build_parent_dict_from_candidates_scaffold(refined_positions)
        
        # 检查是否找到候选位置
        if len(refined_positions) == 0:
            external_mutations.append(new_mut)
            logger.warning(f"Mutation {new_mut} added to external_mutations (no candidate positions found despite having intersection nodes)")
            continue
                
        # 计算贝叶斯罚分并更新 M_current
        final_position, final_imputed_vec, df_penalty_score, M_current = compute_bayesian_penalty_for_positions_scaffold(
            new_mut, refined_positions, T_current, M_current, I_attached, P_attached, parent_dict, intersection_nodes, 
            ω_NA=ω_NA, fnfp_ratio=fnfp_ratio, φ=φ
        )
        logger.info(f"The new_mut should be placed on position: {final_position['placement_type']}.")
        
        # 更新 M_current
        if final_position['placement_type'] == 'on_node':
            mut_in_mtx = final_position['anchor']
            M_current = M_current.rename(columns={mut_in_mtx: mut_in_mtx + '|' + new_mut})
        else:
            M_current[new_mut] = final_imputed_vec
        
        # 更新 T_current
        T_current = add_new_mutation_to_tree_independent(new_mut, T_current, final_position)
        
        # 打印当前树的结构
        logger.info(f"Updated tree after mutation {new_mut}:")
        print_tree(T_current)
        
        # 检查冲突
        if scp.ul.is_conflict_free_gusfield(M_current):
            logger.info(f"Current M_current is conflict-free and shaped as: {M_current.shape}")
        else:
            raise ValueError(f"Current M_current is conflict !!! Break!!!.")
    
    return external_mutations, T_current, M_current




def process_subtree_mutations_to_specific_node(
    subtree_groups, target_node_names, T_current, M_current, I_selected, P_selected, 
    ω_NA, fnfp_ratio, φ, logger, root_mutations=None
):
    """
    通过子树组处理外部突变，支持多突变组和单突变组的分别处理，挂到指定节点下
    """
    
    if root_mutations is None:
        root_mutations = []
    
    external_mutations = []
    
    # 查找目标节点
    if target_node_names is None:
        raise ValueError(f"Target node '{target_node_names}' not found in the tree")
    
    logger.info(f"Found target node: {target_node_names}")
    
    # 分离子树组和单元素组
    multi_mut_subtree_groups = [g for g in subtree_groups if len(g) > 1]
    singleton_subtree_groups = [g for g in subtree_groups if len(g) == 1]
    
    logger.info(f"Found {len(multi_mut_subtree_groups)} multi-mutation groups and {len(singleton_subtree_groups)} singleton groups to attach under {target_node_names}")
    
    ##### 处理长度 >1 的子树组
    for group_idx, group in enumerate(tqdm(multi_mut_subtree_groups, desc="Processing multiple subtrees")):
        logger.info(f"Building subtree {group_idx+1}/{len(multi_mut_subtree_groups)} for mutations: {group}")
        
        # 根据 I_selected 中每个 mutation 的 1 的个数排序（降序）
        sorted_group = sorted(group, key=lambda subtree_mut: I_selected[subtree_mut].sum(), reverse=True)
        # 重排序，把 sorted_group 中跟 T_current 上 node 有交集的突变先放到 list 最前面
        # 简洁版本
        viable_mutations = [mut for mut in sorted_group 
                           if len(find_all_intersect_muts_from_tree_by_matrix_scaffold(T_current, I_selected, mut)) > 0]
        non_viable_mutations = [mut for mut in sorted_group if mut not in viable_mutations]
        
        if not viable_mutations:
            external_mutations.extend(sorted_group)
            logger.info(f"All mutations in group {group_idx+1} added to external_mutations")
            continue
        
        sorted_group = viable_mutations + non_viable_mutations
        
        # 记录当前组的起始节点（第一个突变创建的新节点）
        current_group_start_node = None
        reattached_mutations = []
        for idx, subtree_mut in enumerate(tqdm(sorted_group, desc="Processing mutations in group")):
            T_rollback = copy.deepcopy(T_current)
            M_rollback = M_current.copy()
            
            logger.info(f"Processing mutation {idx+1}/{len(sorted_group)}: {subtree_mut}")
            
            if idx == 0:
                # 找到 intersection nodes
                intersection_nodes = find_all_intersect_muts_from_tree_by_matrix_scaffold(T_current, I_selected, subtree_mut)
                
                # 第一个 mutation 要在当前最晚 node 上找到罚分更小的 node 挂到目标节点的 new_leaf 上
                candidate_positions = find_intersection_positions_within_tree_directly_scaffold(
                    T_current, subtree_mut, I_selected, min_overlap=1
                )
                parent_dict = build_parent_dict_from_candidates_scaffold(candidate_positions)
                potential_positions = [p for i,p in enumerate(candidate_positions) if p['placement_type'] == 'new_leaf']
                # T_current = add_new_mutation_to_tree_independent(subtree_mut, T_current, final_position_list[0])
                final_position, final_imputed_vec, df_penalty_score, M_current = compute_bayesian_penalty_for_positions_scaffold(
                    subtree_mut, potential_positions, T_current, M_current, I_selected, P_selected, parent_dict, set(list(intersection_nodes)+target_node_names), 
                    ω_NA=ω_NA, fnfp_ratio=fnfp_ratio, φ=φ
                )
                M_current[subtree_mut] = final_imputed_vec
                T_current = add_new_mutation_to_tree_independent(subtree_mut, T_current, final_position)
                
                # 记录这个新创建的节点作为后续突变的搜索起点
                current_group_start_node = T_current.find(subtree_mut)
                logger.info(f"Created new leaf node: {subtree_mut}, will use as start node for remaining mutations in group")
            
            else:
                
                # 后续突变从第一个突变创建的节点开始向下搜索
                if current_group_start_node is None:
                    logger.error(f"No start node found for group {group}, skipping mutation {subtree_mut}")
                    reattached_mutations.append(subtree_mut)
                    continue
                
                # 使用优化方法获取候选位置（从当前组的起始节点开始）
                candidate_positions = find_intersection_positions_within_tree_directly_scaffold(
                    T_current, subtree_mut, I_selected, min_overlap=1
                )
                parent_dict = build_parent_dict_from_candidates_scaffold(candidate_positions)
                subtree_nodes = get_subtree_nodes_tree_node(current_group_start_node)
                potential_positions = filter_candidate_positions_from_target_node(candidate_positions, subtree_nodes)
                
                # 检查是否找到候选位置
                if len(potential_positions) == 0:
                    reattached_mutations.append(subtree_mut)
                    logger.info(f"Mutation {subtree_mut} added to reattached_mutations (no candidate positions found)")
                    continue
                
                # 计算贝叶斯罚分并更新 M_current
                final_position, final_imputed_vec, df_penalty_score, M_current = compute_bayesian_penalty_for_positions_scaffold(
                    subtree_mut, potential_positions, T_current, M_current, I_selected, P_selected, parent_dict, subtree_nodes, 
                    ω_NA=ω_NA, fnfp_ratio=fnfp_ratio, φ=φ
                )
                logger.info(f"Mutation {subtree_mut} should be placed on position: {final_position['placement_type']}")
                
                # 更新 M_current
                if final_position['placement_type'] == 'on_node':
                    mut_in_mtx = final_position['anchor']
                    if mut_in_mtx == 'ROOT':
                        root_mutations.append(subtree_mut)
                    else:
                        # 重命名矩阵列，合并突变
                        M_current = M_current.rename(columns={mut_in_mtx: mut_in_mtx + '|' + subtree_mut})
                        T_current = add_new_mutation_to_tree_independent(subtree_mut, T_current, final_position)
                else:
                    # 创建新节点或新边
                    M_current[subtree_mut] = final_imputed_vec
                    T_current = add_new_mutation_to_tree_independent(subtree_mut, T_current, final_position)
            
            # 打印当前树结构
            logger.info(f"Tree after adding {subtree_mut}:")
            print_tree(T_current)
            
            # 检查冲突
            if not scp.ul.is_conflict_free_gusfield(M_current):
                logger.warning(f"Conflict detected after adding {subtree_mut}, rolling back")
                
                # 回滚操作：从矩阵中移除这个突变
                T_current = copy.deepcopy(T_rollback)
                M_current = M_rollback.copy()
                
                # 把这个突变放到 reattached_mutations
                reattached_mutations.append(subtree_mut)
                logger.info(f"Mutation {subtree_mut} added to reattached_mutations due to conflict")
                continue  # 跳过这个突变，继续处理下一个
        
        # 处理重新挂载的突变（第一轮未处理的）
        second_reattached_mutations = []
        if len(reattached_mutations) > 0:
            logger.info(f"Processing {len(reattached_mutations)} reattached mutations for group {group_idx+1}")
            
            sorted_reattached_mutations = [i for i in I_selected.columns if i in reattached_mutations]
            for subtree_mut in tqdm(sorted_reattached_mutations, desc="Processing re-attached mutations"):
                T_rollback = copy.deepcopy(T_current)
                M_rollback = M_current.copy()
                
                logger.info(f"Processing re-attached mutation: {subtree_mut}")
                
                # 仍然从当前组的起始节点开始搜索
                if current_group_start_node is None:
                    second_reattached_mutations.append(subtree_mut)
                    logger.info(f"Mutation {subtree_mut} added to second_reattached_mutations (no start node available)")
                    continue
                
                # 使用优化方法获取候选位置（从当前组的起始节点开始）
                candidate_positions = find_intersection_positions_within_tree_directly_scaffold(
                    T_current, subtree_mut, I_selected, min_overlap=1
                )
                parent_dict = build_parent_dict_from_candidates_scaffold(candidate_positions)
                subtree_nodes = get_subtree_nodes_tree_node(current_group_start_node)
                potential_positions = filter_candidate_positions_from_target_node(candidate_positions, subtree_nodes)
                
                if len(potential_positions) == 0:
                    second_reattached_mutations.append(subtree_mut)
                    logger.info(f"Mutation {subtree_mut} added to second_reattached_mutations (no candidate positions found)")
                    continue
                                
                # 计算贝叶斯罚分并更新 M_current
                final_position, final_imputed_vec, df_penalty_score, M_current = compute_bayesian_penalty_for_positions_scaffold(
                    subtree_mut, potential_positions, T_current, M_current, I_selected, P_selected, parent_dict, subtree_nodes, 
                    ω_NA=ω_NA, fnfp_ratio=fnfp_ratio, φ=φ
                )
                logger.info(f"The subtree_mut should be placed on position: {final_position['placement_type']}.")
                
                if final_position['placement_type'] == 'on_node':
                    mut_in_mtx = final_position['anchor']
                    if mut_in_mtx == 'ROOT':
                        root_mutations.append(subtree_mut)
                    else:
                        M_current = M_current.rename(columns={mut_in_mtx: mut_in_mtx + '|' + subtree_mut})
                        T_current = add_new_mutation_to_tree_independent(subtree_mut, T_current, final_position)
                else:
                    M_current[subtree_mut] = final_imputed_vec
                    T_current = add_new_mutation_to_tree_independent(subtree_mut, T_current, final_position)
                
                logger.info(f"Updated tree after re-attaching mutation {subtree_mut}:")
                print_tree(T_current)
                
                if not scp.ul.is_conflict_free_gusfield(M_current):
                    logger.warning(f"Conflict detected after adding {subtree_mut}, rolling back")
                    
                    # 回滚操作：从矩阵中移除这个突变
                    T_current = copy.deepcopy(T_rollback)
                    M_current = M_rollback.copy()
                    
                    # 把这个突变放到external_mutations
                    second_reattached_mutations.append(subtree_mut)
                    logger.info(f"Mutation {subtree_mut} added to second_reattached_mutations due to conflict")
                    continue  # 跳过这个突变，继续处理下一个
        
        # 记录仍然未处理的突变
        if second_reattached_mutations:
            logger.warning(f"Group {group_idx+1} has {len(second_reattached_mutations)} mutations still remaining: {second_reattached_mutations}")
            external_mutations.extend(second_reattached_mutations)
    
    ##### 处理长度 =1 的单元素组
    logger.info(f"Processing {len(singleton_subtree_groups)} singleton groups")
    # for group_idx, group in enumerate(tqdm(singleton_subtree_groups, desc="Processing singleton subtrees")):
    for group in tqdm(singleton_subtree_groups, desc="Processing singleton subtrees"):
        T_rollback = copy.deepcopy(T_current)
        M_rollback = M_current.copy()
        
        subtree_mut = group[0]
        logger.info(f"Attaching singleton mutation directly to target node {target_node_names}: {subtree_mut}")
                
        # 找到 intersection nodes
        intersection_nodes = find_all_intersect_muts_from_tree_by_matrix_scaffold(T_current, I_selected, subtree_mut)
        # 检查是否有intersection
        if len(intersection_nodes) == 0:
            logger.warning(f"No intersection found for singleton mutation {subtree_mut}, skipping")
            external_mutations.append(subtree_mut)
            continue
        
        # 第一个 mutation 要在当前最晚 node 上找到罚分更小的 node 挂到目标节点的 new_leaf 上
        candidate_positions = find_intersection_positions_within_tree_directly_scaffold(
            T_current, subtree_mut, I_selected, min_overlap=1
        )
        parent_dict = build_parent_dict_from_candidates_scaffold(candidate_positions)
        potential_positions = [p for i,p in enumerate(candidate_positions) if p['placement_type'] == 'new_leaf']
        # T_current = add_new_mutation_to_tree_independent(subtree_mut, T_current, final_position_list[0])
        final_position, final_imputed_vec, df_penalty_score, M_current = compute_bayesian_penalty_for_positions_scaffold(
            subtree_mut, potential_positions, T_current, M_current, I_selected, P_selected, parent_dict, set(list(intersection_nodes)+target_node_names), 
            ω_NA=ω_NA, fnfp_ratio=fnfp_ratio, φ=φ
        )
        M_current[subtree_mut] = final_imputed_vec
        T_current = add_new_mutation_to_tree_independent(subtree_mut, T_current, final_position)
        
        logger.info(f"Tree after adding singleton {subtree_mut}:")
        print_tree(T_current)
        
        if not scp.ul.is_conflict_free_gusfield(M_current):
            logger.warning(f"Conflict detected after adding {subtree_mut}, rolling back")
            
            # 回滚操作：从矩阵中移除这个突变
            T_current = copy.deepcopy(T_rollback)
            M_current = M_rollback.copy()
            
            # 把这个突变放到external_mutations
            external_mutations.append(subtree_mut)
            logger.info(f"Mutation {subtree_mut} added to external_mutations due to conflict")
            continue  # 跳过这个突变，继续处理下一个
    
    logger.info("All external mutations have been processed successfully.")
    logger.info(f"Remained mutations count: {len(external_mutations)}")
    
    return external_mutations, T_current, M_current, root_mutations




# -------------------------
# Sub-grouping within Backbone Groups
# -------------------------

def find_latest_hub_node(root, hub_mutations):
    """
    找到包含 hub_mutations 的所有最晚节点，并返回这些节点及其中的 hub_mutations
    
    Args:
        root: TreeNode 根节点
        hub_mutations: List[str] 要查找的hub突变列表
    
    Returns:
        Tuple[List[TreeNode], List[List[str]]]: (最晚节点列表, 对应节点中的hub_mutations列表)
    """
    hub_nodes = []
    
    # 遍历所有非ROOT节点
    for node in root.traverse():
        if node.name == "ROOT":
            continue
            
        # 解析节点中的突变
        if "|" in node.name:
            node_mutations = node.name.split("|")
        else:
            node_mutations = [node.name]
        
        # 找出该节点中包含的hub_mutations
        hub_in_node = [mut for mut in node_mutations if mut in hub_mutations]
        
        if hub_in_node:
            hub_nodes.append((node, hub_in_node))
    
    if not hub_nodes:
        return [], []
    
    # 计算所有节点的深度
    def calculate_depth(node):
        """计算节点深度（ROOT深度为0）"""
        depth = 0
        current = node
        while current.parent:
            depth += 1
            current = current.parent
        return depth
    
    # 计算所有hub节点的深度
    hub_nodes_with_depth = [(node, hub_muts, calculate_depth(node)) 
                           for node, hub_muts in hub_nodes]
    
    # 找到最大深度
    max_depth = max(depth for _, _, depth in hub_nodes_with_depth)
    
    # 找到所有深度等于最大深度的节点
    latest_nodes = [(node, hub_muts) for node, hub_muts, depth in hub_nodes_with_depth 
                   if depth == max_depth]
    
    # 分离节点列表和突变列表
    if latest_nodes:
        nodes, mutations = zip(*latest_nodes)
        return list(nodes), list(mutations)
    else:
        return [], []

# # 使用修改后的函数
# latest_nodes, found_mutations_list = find_latest_hub_node(T_current, hub_mutations)

# if latest_nodes:
#     print(f"找到 {len(latest_nodes)} 个最晚包含hub突变的节点:")
#     for i, (node, mutations) in enumerate(zip(latest_nodes, found_mutations_list)):
#         print(f"节点 {i+1}: {node.name}")
#         print(f"  包含的hub突变: {mutations}")
# else:
#     print("没有找到包含hub突变的节点")




# 更简洁的版本（使用向量化操作）
def calculate_group_cooccurrence_fraction(sorted_parallel_groups, backbone_mut, I_selected):
    """
    简化版本，使用向量化操作提高效率
    """
    # 找到backbone突变为1的cells
    backbone_pos_mask = I_selected[backbone_mut] == 1
    backbone_cells = I_selected[backbone_pos_mask].index
    
    if len(backbone_cells) == 0:
        return {group_id: 0.0 for group_id in sorted_parallel_groups.keys()}
    
    group_fractions = {}
    
    for group_id, mutation_list in sorted_parallel_groups.items():
        # 检查突变是否存在
        valid_mutations = [mut for mut in mutation_list if mut in I_selected.columns]
        if not valid_mutations:
            group_fractions[group_id] = 0.0
            continue
        
        # 提取数据并计算
        group_data = I_selected.loc[backbone_cells, valid_mutations]
        
        # 使用max(axis=1)来找到每个cell是否有至少一个突变为1（忽略NA）
        # 先填充NA为0，然后找最大值
        group_data_filled = group_data.fillna(0)
        has_mutation = (group_data_filled == 1).any(axis=1)
        
        fraction = has_mutation.mean() if len(has_mutation) > 0 else 0.0
        group_fractions[group_id] = fraction
    
    return group_fractions

# # 使用示例：
# group_fractions = calculate_group_cooccurrence_fraction(
#     sorted_parallel_groups, backbone_mut, I_selected
# )
# # 输出示例:
# # {0: 0.85, 1: 0.72, 2: 0.63}

def find_max_fraction_group(group_fractions):
    """
    找到比例最大的group key
    
    参数:
    group_fractions: dict, {group_id: fraction}
    
    返回:
    tuple: (max_key, max_value)
    """
    if not group_fractions:
        return None, 0.0
    
    max_key = max(group_fractions, key=group_fractions.get)
    max_value = group_fractions[max_key]
    
    return max_key, max_value

# # 使用示例
# max_key, max_value = find_max_fraction_group(group_fractions)
# print(f"最大比例的group: {max_key}, 比例: {max_value:.3f}")
# # 输出: 最大比例的group: 1, 比例: 0.069


def calculate_cv_for_single_mutation(df_reads_group, mutation, cv_thresh, logger):
    """
    计算单个突变的CV统计信息
    
    Parameters:
    -----------
    df_reads_group : pd.DataFrame
        包含覆盖度数据的DataFrame
    mutation : str
        突变名称
    cv_thresh : float
        CV阈值
    logger : logging.Logger
        日志记录器
        
    Returns:
    --------
    dict : 包含CV统计信息的字典
    """
    try:
        cov_data = df_reads_group[mutation].dropna()
        if len(cov_data) == 0:
            return {'cov_median': 0, 'cov_CV': float('inf'), 'cov_mean': 0, 'cov_std': 0, 'pass_CV': False}
        
        cov_median = np.median(cov_data)
        cov_mean = np.mean(cov_data)
        cov_std = np.std(cov_data)
        cov_CV = cov_std / cov_mean if cov_mean > 0 else float('inf')
        
        return {
            'cov_median': cov_median,
            'cov_CV': cov_CV,
            'cov_mean': cov_mean,
            'cov_std': cov_std,
            'pass_CV': cov_CV <= cv_thresh
        }
    except Exception as e:
        logger.warning(f"Error calculating CV for mutation {mutation}: {str(e)}")
        return {'cov_median': 0, 'cov_CV': float('inf'), 'cov_mean': 0, 'cov_std': 0, 'pass_CV': False}


def generate_subgrouping_report(complete_mutation_hierarchy, backbones_of_group, subgroup_backbone_mutations, 
                              cv_stats_all, outputpath, sampleid, logger, high_cv_mutations=None):
    """Generate detailed report for subgrouping results including high CV mutations."""
    
    # 按主分组统计
    main_group_summary = {}
    for group_id in backbones_of_group.keys():
        group_mutations = [mut for mut, info in complete_mutation_hierarchy.items() 
                          if info['main_group'] == group_id]
        backbone_mut = backbones_of_group[group_id]
        sub_groups = set([info['sub_group'] for info in complete_mutation_hierarchy.values() 
                         if info['main_group'] == group_id])
        
        main_group_summary[group_id] = {
            'backbone_mutation': backbone_mut,
            'total_mutations': len(group_mutations),
            'subgroups_count': len(sub_groups),
            'mutations': group_mutations
        }
    
    # 生成报告
    report_lines = []
    report_lines.append("=" * 80)
    report_lines.append("MUTATION SUB-GROUPING FINAL REPORT")
    report_lines.append("=" * 80)
    report_lines.append(f"Sample: {sampleid}")
    report_lines.append(f"Total main groups: {len(backbones_of_group)}")
    report_lines.append(f"Total subgroups: {len(subgroup_backbone_mutations)}")
    report_lines.append(f"Total mutations with hierarchy: {len(complete_mutation_hierarchy)}")
    
    # 添加高CV突变信息
    if high_cv_mutations is not None:
        report_lines.append(f"High CV mutations (excluded from analysis): {len(high_cv_mutations)}")
    else:
        report_lines.append(f"High CV mutations (excluded from analysis): 0")
    
    # 添加CV统计信息（如果可用）
    if cv_stats_all is not None and not cv_stats_all.empty:
        passed_cv = len(cv_stats_all[cv_stats_all['pass_CV'] == True])
        total_cv = len(cv_stats_all)
        report_lines.append(f"CV filtering: {passed_cv}/{total_cv} mutations passed CV threshold")
    
    report_lines.append("")
    
    # 详细的主分组信息
    for group_id, summary in main_group_summary.items():
        report_lines.append(f"Main Group {group_id}:")
        report_lines.append(f"  Backbone mutation: {summary['backbone_mutation']}")
        report_lines.append(f"  Total mutations: {summary['total_mutations']}")
        report_lines.append(f"  Number of subgroups: {summary['subgroups_count']}")
        
        # 按子分组显示
        subgroup_mutations = {}
        for mut in summary['mutations']:
            sub_group = complete_mutation_hierarchy[mut]['sub_group']
            if sub_group not in subgroup_mutations:
                subgroup_mutations[sub_group] = []
            subgroup_mutations[sub_group].append(mut)
        
        for sub_group, muts in subgroup_mutations.items():
            sub_backbone = subgroup_backbone_mutations.get(sub_group, "None")
            backbone_indicator = " (SUB_BACKBONE)" if sub_backbone in muts else ""
            report_lines.append(f"    Subgroup {sub_group}{backbone_indicator}: {len(muts)} mutations")
            for mut in muts:
                flags = []
                if mut == summary['backbone_mutation']:
                    flags.append("MAIN_BACKBONE")
                if mut == sub_backbone:
                    flags.append("SUB_BACKBONE")
                flag_str = " [" + ", ".join(flags) + "]" if flags else ""
                report_lines.append(f"      - {mut}{flag_str}")
        report_lines.append("")
    
    # 添加高CV突变详细列表
    if high_cv_mutations and len(high_cv_mutations) > 0:
        report_lines.append("=" * 80)
        report_lines.append("HIGH CV MUTATIONS (EXCLUDED FROM ANALYSIS)")
        report_lines.append("=" * 80)
        report_lines.append(f"Total high CV mutations: {len(high_cv_mutations)}")
        report_lines.append("")
        
        # 按主分组组织高CV突变
        high_cv_by_group = {}
        for mut in high_cv_mutations:
            # 尝试从mutation_group中获取分组信息，或者根据其他逻辑确定分组
            group_id = "Unknown"
            # 这里可以根据您的数据结构调整如何确定高CV突变的分组
            # 例如：如果mutation_group包含所有突变的信息
            # group_id = mutation_group.get(mut, "Unknown")
            high_cv_by_group.setdefault(group_id, []).append(mut)
        
        for group_id, muts in high_cv_by_group.items():
            if group_id == "Unknown":
                report_lines.append(f"Mutations with unknown group assignment:")
            else:
                report_lines.append(f"Group {group_id}:")
            
            for mut in sorted(muts):
                # 尝试获取CV统计信息（如果可用）
                cv_info = ""
                if cv_stats_all is not None and not cv_stats_all.empty:
                    mut_cv_stats = cv_stats_all[cv_stats_all.index == mut]
                    if not mut_cv_stats.empty:
                        cv_value = mut_cv_stats.iloc[0]['cov_CV']
                        cv_threshold = mut_cv_stats.iloc[0].get('cv_threshold', 'Unknown')
                        cv_info = f" (CV: {cv_value:.3f}, threshold: {cv_threshold})"
                
                report_lines.append(f"  - {mut}{cv_info}")
            report_lines.append("")
        
        report_lines.append("Note: High CV mutations were excluded from sub-grouping and tree building")
        report_lines.append("due to high coverage variability (CV above threshold).")
        report_lines.append("")
    
    # 添加总结统计
    report_lines.append("=" * 80)
    report_lines.append("SUMMARY STATISTICS")
    report_lines.append("=" * 80)
    
    total_analyzed = len(complete_mutation_hierarchy)
    total_high_cv = len(high_cv_mutations) if high_cv_mutations else 0
    total_considered = total_analyzed + total_high_cv
    
    report_lines.append(f"Total mutations considered: {total_considered}")
    report_lines.append(f"  - Successfully analyzed and placed in tree: {total_analyzed}")
    report_lines.append(f"  - Excluded due to high CV: {total_high_cv}")
    
    if total_considered > 0:
        analyzed_percent = (total_analyzed / total_considered) * 100
        high_cv_percent = (total_high_cv / total_considered) * 100
        report_lines.append(f"Success rate: {analyzed_percent:.1f}% ({total_analyzed}/{total_considered})")
        report_lines.append(f"High CV exclusion rate: {high_cv_percent:.1f}% ({total_high_cv}/{total_considered})")
    
    # 保存报告
    report_path = os.path.join(outputpath, f"{sampleid}.subgrouping_report.txt")
    with open(report_path, 'w') as f:
        f.write("\n".join(report_lines))
    
    logger.info(f"Subgrouping report saved to: {report_path}")
    
    # 同时保存高CV突变的单独列表文件
    if high_cv_mutations and len(high_cv_mutations) > 0:
        high_cv_path = os.path.join(outputpath, f"{sampleid}.high_cv_mutations.txt")
        with open(high_cv_path, 'w') as f:
            f.write("High CV Mutations (excluded from analysis)\n")
            f.write("=" * 50 + "\n\n")
            for mut in sorted(high_cv_mutations):
                f.write(f"{mut}\n")
        logger.info(f"High CV mutations list saved to: {high_cv_path}")


def perform_subgrouping_within_backbone_groups_and_build_initial_scaffold_tree(sorted_I_resolved, sorted_P_resolved, T_current, M_current, mutation_group, backbones_of_group, df_reads_resolved, df_features_new, outputpath, sampleid, logger, params, M_B, root_mutations, cutoff_mcf_for_graph, cutoff_mcn_for_graph):
    """
    Perform sub-grouping within each backbone group to identify sub-clones.
    Includes CV-based filtering before subgrouping.
    
    Parameters:
    -----------
    sorted_I_resolved : pd.DataFrame
        Resolved binary mutation matrix
    sorted_P_resolved : pd.DataFrame
        Resolved posterior mutation matrix
    T_current : object
        Current tree structure
    M_current : object  
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
    logger : logging.Logger
        Logger instance
    params : dict
        Including CV threshold for filtering mutations before subgrouping
    M_B : pd.DataFrame
        Mutation backbone matrix
    root_mutations : list
        List of root mutations
        
    Returns:
    --------
    Tuple containing:
        - complete_mutation_hierarchy: Complete hierarchy information for all mutations
        - subgroup_backbone_mutations: Backbone mutations for each subgroup
        - subgroup_details: Detailed information about each subgroup
        - cv_stats_all: CV statistics for all mutations considered
        - T_current: Updated tree structure
        - M_current: Updated mutation matrix
        - root_mutations: Updated root mutations
        - external_mutations: Mutations that failed to be added to the tree
        - high_cv_mutations: Mutations with CV above threshold
    """
    logger.info("Performing sub-grouping within backbone groups...")
    logger.info(f"Using CV threshold: {params['cv_thresh']} for pre-subgrouping filtration")
    
    # 存储最终的子分组结果
    mutation_subgroups = {}
    subgroup_backbone_mutations = {}
    subgroup_details = {}
    all_cv_stats = []
    external_mutations = []
    high_cv_mutations = []  # 新增：存储高CV突变
    
    df_cv_stats_within_subgroup = pd.DataFrame(columns=['cov_median', 'cov_CV', 'cov_mean', 'cov_std', 'pass_CV', 'group_id', 'backbone_mut'])
    
    # 定义全局参数（需要根据实际情况设置）
    ω_NA = params['general_weight_NA'] if params['general_weight_NA'] else 0.001
    fnfp_ratio = params['fnfp_ratio']
    φ = params['phi']
    
    for group_id, backbone_mut in backbones_of_group.items():
        logger.info(f"Processing subgroup for backbone group {group_id} with backbone mutation {backbone_mut}")
        
        dir_subgroup = outputpath+"/backbone_clone_"+str(group_id)+"."+str(backbone_mut)
        os.makedirs(dir_subgroup, exist_ok=True)
        
        # 获取当前group的所有突变（不包括backbone mutation）
        group_muts = [mut for mut, gid in mutation_group.items() 
                     if gid == group_id and mut != backbone_mut]
        
        # 生成对应的 subgroup 的数据矩阵
        I_group = sorted_I_resolved.loc[list(M_B[M_B[backbone_mut]==1].index), group_muts]
        df_reads_group = df_reads_resolved.loc[['bulk']+list(M_B[M_B[backbone_mut]==1].index), group_muts]
        
        if len(group_muts) < 1:
            ### condition_0: 这个 backbone group 下只有一个 backbone mutation 没有其他的可挂突变
            logger.info(f"Group {group_id} has no non-backbone mutations, skipping sub-grouping")
            continue
        
        elif len(group_muts) == 1:
            ### condition_1: 可能这下面只有一个突变就直接往上挂
            logger.info(f"Group {group_id} has only 1 non-backbone mutations, directly hang mutation to tree")
            
            # 首先检查这个突变的CV值
            single_mut = group_muts[0]
            cv_stats_single = calculate_cv_for_single_mutation(df_reads_group, single_mut, params['cv_thresh'], logger)
            
            # 如果CV值高于阈值，添加到high_cv_mutations
            if cv_stats_single['cov_CV'] > params['cv_thresh']:
                logger.info(f"Single mutation {single_mut} has high CV ({cv_stats_single['cov_CV']:.3f}), adding to high_cv_mutations")
                high_cv_mutations.append(single_mut)
                continue  # 跳过挂树操作
            
            # 如果CV值正常，继续挂树
            mutation_subgroups[single_mut] = f"{group_id}_0"
            subgroup_backbone_mutations[f"{group_id}_0"] = single_mut
            subgroup_details[f"{group_id}_0"] = {
                'main_group': group_id,
                'backbone_mutation': backbone_mut,
                'mutations': group_muts,
                'is_trivial': True,
                'cv_filter_applied': False
            }
            # 一个突变直接挂树
            subtree_groups = [group_muts]
            latest_nodes, found_mutations_list = find_latest_hub_node(T_current, [backbone_mut])
            target_node_names = [node.name for node in latest_nodes]  # 获取所有目标节点名称
            external_mutations_temp, T_current, M_current, root_mutations = process_subtree_mutations_to_specific_node(
                subtree_groups, target_node_names, T_current, M_current, sorted_I_resolved, sorted_P_resolved, 
                ω_NA, fnfp_ratio, φ, logger, root_mutations
            )            
            external_mutations.extend(external_mutations_temp)
            continue
        
        logger.info(f"Group {group_id} has {len(group_muts)} non-backbone mutations for sub-grouping")
        
        # Step 1: 计算CV值并过滤
        logger.info(f"Calculating CV values for {len(group_muts)} mutations in group {group_id}")
        low_cv_mutations, cv_stats_group = calculate_cv_for_subgrouping(df_reads_group, group_muts, params['cv_thresh'], logger)
        
        # 新增：将高CV突变添加到high_cv_mutations列表
        high_cv_in_group = [mut for mut in group_muts if mut not in low_cv_mutations]
        high_cv_mutations.extend(high_cv_in_group)
        logger.info(f"Group {group_id}: {len(high_cv_in_group)} mutations with high CV added to high_cv_mutations: {high_cv_in_group}")
        
        # 添加组信息到CV统计中
        cv_stats_group['group_id'] = group_id
        cv_stats_group['backbone_mut'] = backbone_mut
        df_cv_stats_within_subgroup = pd.concat([df_cv_stats_within_subgroup, cv_stats_group], ignore_index=True)
        
        # 如果低CV突变数量不足，跳过子分组
        if len(low_cv_mutations) < 2:
            ### condition_2: 可能有多个突变但是不足以执行 graph，就按照 maf 顺序往上挂
            logger.info(f"Group {group_id}: only {len(low_cv_mutations)} mutations passed CV filter, skipping sub-grouping")
            
            # 只处理通过CV过滤的突变，高CV突变已经添加到high_cv_mutations
            for i, mut in enumerate(low_cv_mutations):  # 修改：只遍历low_cv_mutations
                subgroup_id = f"{group_id}_{i}"
                mutation_subgroups[mut] = subgroup_id
                subgroup_backbone_mutations[subgroup_id] = mut
                subgroup_details[subgroup_id] = {
                    'main_group': group_id,
                    'backbone_mutation': backbone_mut,
                    'mutations': [mut],
                    'is_trivial': True,
                    'cv_filter_applied': True,
                    'passed_cv_filter': True  # 这些突变都通过了CV过滤
                }
            
            # 按MAF顺序挂树（只挂低CV突变）
            if low_cv_mutations:  # 确保有突变需要挂树
                subtree_groups = [[mut] for mut in low_cv_mutations]  # 每个突变单独成组
                latest_nodes, found_mutations_list = find_latest_hub_node(T_current, [backbone_mut])
                target_node_names = [node.name for node in latest_nodes]  # 获取所有目标节点名称
                external_mutations_temp, T_current, M_current, root_mutations = process_subtree_mutations_to_specific_node(
                    subtree_groups, target_node_names, T_current, M_current, sorted_I_resolved, sorted_P_resolved, 
                    ω_NA, fnfp_ratio, φ, logger, root_mutations
                )
                external_mutations.extend(external_mutations_temp)
            continue
        
        logger.info(f"Group {group_id}: {len(low_cv_mutations)}/{len(group_muts)} mutations passed CV filter, proceeding with sub-grouping")
        
        # 从原始矩阵中提取通过CV过滤的突变数据
        I_group_low_cv = I_group[low_cv_mutations].copy()
        
        # 进一步过滤掉全为0或全为NA的突变
        valid_mutations = []
        for mut in low_cv_mutations:
            mut_data = I_group_low_cv[mut]
            if (mut_data == 1).any() and (mut_data == 0).any():
                valid_mutations.append(mut)
        
        # 新增：将无效的突变添加到external_mutations
        invalid_mutations = [mut for mut in low_cv_mutations if mut not in valid_mutations]
        external_mutations.extend(invalid_mutations)
        logger.info(f"Group {group_id}: {len(invalid_mutations)} mutations failed quality filter, added to external_mutations")
        
        I_group_final = I_group_low_cv[valid_mutations]
        
        if I_group_final.shape[1] < 2:
            ### condition_2: 可能有多个突变但是不足以执行 graph，就按照 maf 顺序往上挂
            logger.info(f"After quality filtration, group {group_id} has {I_group_final.shape[1]} mutations, skipping sub-grouping")
            
            # 只处理有效的低CV突变
            for i, mut in enumerate(valid_mutations):
                subgroup_id = f"{group_id}_{i}"
                mutation_subgroups[mut] = subgroup_id
                subgroup_backbone_mutations[subgroup_id] = mut
                subgroup_details[subgroup_id] = {
                    'main_group': group_id,
                    'backbone_mutation': backbone_mut,
                    'mutations': [mut],
                    'is_trivial': True,
                    'cv_filter_applied': True,
                    'passed_cv_filter': True,
                    'passed_quality_filter': True
                }
            
            # 按MAF顺序挂树（只挂有效的低CV突变）
            if valid_mutations:
                subtree_groups = [[mut] for mut in valid_mutations]
                latest_nodes, found_mutations_list = find_latest_hub_node(T_current, [backbone_mut])
                target_node_names = [node.name for node in latest_nodes]  # 获取所有目标节点名称
                external_mutations_temp, T_current, M_current, root_mutations = process_subtree_mutations_to_specific_node(
                    subtree_groups, target_node_names, T_current, M_current, sorted_I_resolved, sorted_P_resolved, 
                    ω_NA, fnfp_ratio, φ, logger, root_mutations
                )
                external_mutations.extend(external_mutations_temp)
            continue
        
        logger.info(f"Group {group_id}: {I_group_final.shape[1]} mutations available for Leiden sub-grouping")
        
        try:
            # Step 2: 计算突变间的相关性权重
            logger.info(f"Calculating correlation weights for group {group_id}...")
            clone_weights_sub, pair_weights_sub = get_correlation_graph_elements(I_group_final, 100, 42, cutoff_mcf_for_graph, cutoff_mcn_for_graph)
            
            if not pair_weights_sub or len(pair_weights_sub) < 2:
                ### condition_2: 可能有多个突变但是不足以执行 graph，就按照 maf 顺序往上挂
                # 如果图结构不足以进行社区检测，直接按MAF顺序挂树
                if not pair_weights_sub:
                    logger.warning(f"Group {group_id}: pair_weights is empty, all clones are singleton mutations")
                else:
                    logger.warning(f"Group {group_id}: insufficient edges ({len(pair_weights_sub)}) for community detection")
                
                # 但仍然需要处理这些突变 - 按MAF顺序挂树
                logger.info(f"Group {group_id}: hanging {len(valid_mutations)} mutations to tree in MAF order")
                
                # 创建子分组信息
                for i, mut in enumerate(valid_mutations):
                    subgroup_id = f"{group_id}_{i}"
                    mutation_subgroups[mut] = subgroup_id
                    subgroup_backbone_mutations[subgroup_id] = mut
                    subgroup_details[subgroup_id] = {
                        'main_group': group_id,
                        'backbone_mutation': backbone_mut,
                        'mutations': [mut],
                        'is_trivial': True,
                        'cv_filter_applied': True,
                        'passed_cv_filter': True,
                        'passed_quality_filter': True,
                        'community_detection_skipped': True,
                        'skip_reason': 'insufficient_graph_structure'
                    }
                
                # 关键：执行挂树操作, 将每个突变作为单独的子组挂树
                subtree_groups = [[mut] for mut in valid_mutations]
                latest_nodes, found_mutations_list = find_latest_hub_node(T_current, [backbone_mut])
                target_node_names = [node.name for node in latest_nodes]  # 获取所有目标节点名称
                external_mutations_temp, T_current, M_current, root_mutations = process_subtree_mutations_to_specific_node(
                    subtree_groups, target_node_names, T_current, M_current, sorted_I_resolved, sorted_P_resolved, 
                    ω_NA, fnfp_ratio, φ, logger, root_mutations
                )
                external_mutations.extend(external_mutations_temp)
                continue
            
            else:
                # Step 3: 使用Leiden算法进行子分组
                logger.info(f"Performing Leiden sub-grouping for group {group_id}...")
                mutation_subgroup, partition_sub, G_ig_sub = leiden_mutation_groups(
                    clone_weights_sub, pair_weights_sub, 
                    dir_subgroup + "/" + sampleid + f".group_{group_id}_subgraph.pdf", 
                    params['resolution_of_graph']
                )
                external_mutations_outgroup = [i for i in valid_mutations if i not in list(mutation_subgroup.keys())]
                external_mutations.extend(external_mutations_outgroup)
                
                # Step 4: 检测当前 graph 中的 hub clusters
                hub_clusters, cluster_degrees = detect_hub_clusters(G_ig_sub, mutation_subgroup)
                
                if len(hub_clusters) == 0:
                    ### condition_3: 可能像 epi non-tumor 能分 group 但是分不出来 hub 也是直接按照 maf 往上挂
                    if len(set(mutation_subgroup.values())) <= 2:
                        # 直接这里面的所有的突变就按照 maf 排序 within subclone 内往树上挂
                        sorted_subgroup_mutations_but_backbone = [i for i in sorted_I_resolved.columns if i in list(mutation_subgroup.keys())]
                        external_mutations_temp, T_current, M_current = integrate_mutations_to_scaffold_within_group(
                            sorted_attached_mutations=sorted_subgroup_mutations_but_backbone,
                            T_current=T_current,
                            M_current=M_current,
                            I_attached=sorted_I_resolved,
                            P_attached=sorted_P_resolved,
                            mutation_group=mutation_group,
                            ω_NA=ω_NA,
                            fnfp_ratio=fnfp_ratio,
                            φ=φ,
                            logger=logger
                        )
                        external_mutations.extend(external_mutations_temp)
                    else:
                        # 多个cluster并列或者只有一个 cluster 了，按cluster分组挂树
                        sorted_parallel_groups = {v: sorted([k for k in mutation_subgroup if mutation_subgroup[k] == v],
                                                           key=lambda x: list(sorted_I_resolved.columns).index(x))
                                                for v in set(mutation_subgroup.values())}
                        subtree_groups = list(sorted_parallel_groups.values())
                        latest_nodes, found_mutations_list = find_latest_hub_node(T_current, [backbone_mut])
                        target_node_names = [node.name for node in latest_nodes]  # 获取所有目标节点名称
                        external_mutations_temp, T_current, M_current, root_mutations = process_subtree_mutations_to_specific_node(
                            subtree_groups, target_node_names, T_current, M_current, sorted_I_resolved, sorted_P_resolved, 
                            ω_NA, fnfp_ratio, φ, logger, root_mutations
                        )
                        external_mutations.extend(external_mutations_temp)
                
                elif len(hub_clusters) == 1:
                    ### condition_4: 可能像 epi tumor 第一层先分出来一个 hub 直接挂，再往下再分出来相互独立的 cluster 这些挂到树上会平行的分别成 subtrees
                    hub_mutations = [i for i,m in mutation_subgroup.items() if m==hub_clusters[0]]
                    # 要把这个为一个 hub cluster 拿出来先按照 maf 挂到树上去，之后再做一次 graph 分群
                    sorted_subgroup_mutations_which_hub = [i for i in sorted_I_resolved.columns if i in hub_mutations]
                    external_mutations_temp, T_current, M_current = integrate_mutations_to_scaffold_within_group(
                        sorted_attached_mutations=sorted_subgroup_mutations_which_hub,
                        T_current=T_current,
                        M_current=M_current,
                        I_attached=sorted_I_resolved,
                        P_attached=sorted_P_resolved,
                        mutation_group=mutation_group,
                        ω_NA=ω_NA,
                        fnfp_ratio=fnfp_ratio,
                        φ=φ,
                        logger=logger
                    )
                    external_mutations.extend(external_mutations_temp)
                    
                    # 去掉上面的 hub clusters 中的 mutations 之后在剩余的突变中再次利用 graph 划分 group 直至分成并列的 clone 或者直至不能再分
                    I_group_dehub = I_group_final.drop(columns=hub_mutations, errors='ignore').copy()
                    
                    if I_group_dehub.shape[1] >= 2:
                        clone_weights_dehub, pair_weights_dehub = get_correlation_graph_elements(I_group_dehub, 100, 42, cutoff_mcf_for_graph, cutoff_mcn_for_graph)
                        
                        if pair_weights_dehub and len(pair_weights_dehub) >= 2:
                            mutation_dehubgroup, partition_dehub, G_ig_dehub = leiden_mutation_groups(
                                clone_weights_dehub, pair_weights_dehub, 
                                dir_subgroup + "/" + sampleid + f".group_{group_id}_for_group_exclude_hub_cluster.pdf", 
                                params['resolution_of_graph']
                            )
                            hub_clusters_dehub, cluster_degrees_dehub = detect_hub_clusters(G_ig_dehub, mutation_dehubgroup)
                            
                            if all(np.array(list(cluster_degrees_dehub.values())) == 0):
                                # parallel_groups = {v: [k for k in mutation_dehubgroup if mutation_dehubgroup[k] == v] for v in set(mutation_dehubgroup.values())}
                                sorted_parallel_groups = {v: sorted([k for k in mutation_dehubgroup if mutation_dehubgroup[k] == v],
                                                                   key=lambda x: list(sorted_I_resolved.columns).index(x))
                                                        for v in set(mutation_dehubgroup.values())}
                                subtree_groups = list(sorted_parallel_groups.values())
                                latest_nodes, found_mutations_list = find_latest_hub_node(T_current, hub_mutations)
                                target_node_names = [node.name for node in latest_nodes]  # 获取所有目标节点名称
                                external_mutations_temp, T_current, M_current, root_mutations = process_subtree_mutations_to_specific_node(
                                    subtree_groups, target_node_names, T_current, M_current, sorted_I_resolved, sorted_P_resolved, 
                                    ω_NA, fnfp_ratio, φ, logger, root_mutations
                                )
                                external_mutations.extend(external_mutations_temp)
                            
                            else:
                                external_mutations.extend(list(I_group_dehub.columns))
                
                else:
                    ### condition_5: 可能直接没有 hub cluster 直接就是相互独立的 cluster 就直接平行的挂树
                    # 这个应就是直接很多个 cluster 并列了（目前看起来只有 backbone cluster 的格式）
                    if all(np.array(list(cluster_degrees.values())) == 0):
                        # 创建sorted_parallel_groups
                        sorted_parallel_groups = {v: sorted([k for k in mutation_subgroup if mutation_subgroup[k] == v],
                                                           key=lambda x: list(sorted_I_resolved.columns).index(x))
                                                for v in set(mutation_subgroup.values())}
                        
                        # 计算每一个分组中 cells 占比
                        group_fractions = calculate_group_cooccurrence_fraction(sorted_parallel_groups, backbone_mut, I_group_final)
                        max_key, max_value = find_max_fraction_group(group_fractions)
                        
                        # 分情况处理
                        if max_value > 0.9:
                            # 要再分一次
                            I_group_maxsub = I_group_final[sorted_parallel_groups[max_key]].copy()
                            
                            if I_group_maxsub.shape[1] >= 2:
                                clone_weights_maxsub, pair_weights_maxsub = get_correlation_graph_elements(I_group_maxsub, 100, 42, cutoff_mcf_for_graph, cutoff_mcn_for_graph)
                                
                                if pair_weights_maxsub and len(pair_weights_maxsub) >= 2:
                                    mutation_maxsubgroup, partition_maxsub, G_ig_maxsub = leiden_mutation_groups(
                                        clone_weights_maxsub, pair_weights_maxsub, 
                                        dir_subgroup + "/" + sampleid + f".group_{group_id}_for_more90percent_group.pdf", 
                                        params['resolution_of_graph']
                                    )
                                    hub_clusters_maxsub, cluster_degrees_maxsub = detect_hub_clusters(G_ig_maxsub, mutation_maxsubgroup)
                                    
                                    if all(np.array(list(cluster_degrees_maxsub.values())) == 0) or len(cluster_degrees_maxsub) == 1:
                                        sorted_parallel_groups_inner = {v: sorted([k for k in mutation_maxsubgroup if mutation_maxsubgroup[k] == v],
                                                                                  key=lambda x: list(sorted_I_resolved.columns).index(x))
                                                                      for v in set(mutation_maxsubgroup.values())}
                                        subtree_groups = list(sorted_parallel_groups_inner.values())
                                        latest_nodes, found_mutations_list = find_latest_hub_node(T_current, [backbone_mut])
                                        target_node_names = [node.name for node in latest_nodes]  # 获取所有目标节点名称
                                        external_mutations_temp, T_current, M_current, root_mutations = process_subtree_mutations_to_specific_node(
                                            subtree_groups, target_node_names, T_current, M_current, sorted_I_resolved, sorted_P_resolved, 
                                            ω_NA, fnfp_ratio, φ, logger, root_mutations
                                        )
                                        external_mutations.extend(external_mutations_temp)
                        
                        else:
                            # 直接一个一个 subtree 往上挂
                            subtree_groups = list(sorted_parallel_groups.values())
                            latest_nodes, found_mutations_list = find_latest_hub_node(T_current, [backbone_mut])
                            target_node_names = [node.name for node in latest_nodes]  # 获取所有目标节点名称
                            external_mutations_temp, T_current, M_current, root_mutations = process_subtree_mutations_to_specific_node(
                                subtree_groups, target_node_names, T_current, M_current, sorted_I_resolved, sorted_P_resolved, 
                                ω_NA, fnfp_ratio, φ, logger, root_mutations
                            )
                            external_mutations.extend(external_mutations_temp)
            
            # Step 4: 对子分组矩阵进行排序
            I_group_sorted, mut_df_sorted_sub, subgroup_to_muts, final_order_sub = sort_I_hierarchical_freeze_ones_fixed(
                I_group_final, mutation_subgroup
            )
            
            # Step 5: 选择每个子分组的early mutation
            subgroup_backbones = select_founder_mutations(I_group_sorted, mutation_subgroup)
            
            # 存储子分组结果
            for mut, subgroup_id in mutation_subgroup.items():
                full_subgroup_id = f"{group_id}_{subgroup_id}"
                mutation_subgroups[mut] = full_subgroup_id
            
            # 存储子分组的backbone mutations
            for subgroup_id, sub_backbone in subgroup_backbones.items():
                full_subgroup_id = f"{group_id}_{subgroup_id}"
                subgroup_backbone_mutations[full_subgroup_id] = sub_backbone
            
            # 存储子分组详细信息
            for subgroup_id in set(mutation_subgroup.values()):
                full_subgroup_id = f"{group_id}_{subgroup_id}"
                subgroup_muts = [mut for mut, sid in mutation_subgroup.items() if sid == subgroup_id]
                subgroup_details[full_subgroup_id] = {
                    'main_group': group_id,
                    'backbone_mutation': backbone_mut,
                    'subgroup_backbone': subgroup_backbones.get(subgroup_id),
                    'mutations': subgroup_muts,
                    'is_trivial': False,
                    'cv_filter_applied': True,
                    'passed_cv_filter': all(mut in low_cv_mutations for mut in subgroup_muts),
                    'passed_quality_filter': all(mut in valid_mutations for mut in subgroup_muts)
                }
                
            logger.info(f"Group {group_id}: found {len(subgroup_backbones)} subgroups with backbones: {list(subgroup_backbones.values())}")
            
            # 可视化子分组结果（如果有相应的数据）
            try:
                # 注意：这里需要根据实际情况调整，因为df_celltype_sub可能未定义
                plot_heatmap_with_celltype_by_your_sorting(
                    I_group_sorted, 
                    None,  # 替换为实际的celltype数据
                    mutation_subgroup,
                    list(mut_df_sorted_sub['mutation']),
                    os.path.join(dir_subgroup, f"{sampleid}.group_{group_id}_subgroup_heatmap.pdf")
                )
            except Exception as e:
                logger.warning(f"Could not generate heatmap for group {group_id}: {str(e)}")
            
        except Exception as e:
            logger.error(f"Error in sub-grouping for group {group_id}: {str(e)}")
            # 如果子分组失败，为每个有效的低CV突变创建独立的子分组
            for i, mut in enumerate(valid_mutations):
                subgroup_id = f"{group_id}_{i}"
                mutation_subgroups[mut] = subgroup_id
                subgroup_backbone_mutations[subgroup_id] = mut
                subgroup_details[subgroup_id] = {
                    'main_group': group_id,
                    'backbone_mutation': backbone_mut,
                    'mutations': [mut],
                    'is_trivial': True,
                    'cv_filter_applied': True,
                    'passed_cv_filter': True,
                    'passed_quality_filter': True
                }
            
            # 按MAF顺序挂树（只挂有效的低CV突变）
            if valid_mutations:
                subtree_groups = [[mut] for mut in valid_mutations]
                latest_nodes, found_mutations_list = find_latest_hub_node(T_current, [backbone_mut])
                target_node_names = [node.name for node in latest_nodes]  # 获取所有目标节点名称
                external_mutations_temp, T_current, M_current, root_mutations = process_subtree_mutations_to_specific_node(
                    subtree_groups, target_node_names, T_current, M_current, sorted_I_resolved, sorted_P_resolved, 
                    ω_NA, fnfp_ratio, φ, logger, root_mutations
                )
                external_mutations.extend(external_mutations_temp)
    
    # 合并所有CV统计数据
    cv_stats_all = pd.concat(all_cv_stats, axis=0) if all_cv_stats else pd.DataFrame()
    
    # 创建完整的分组层次结构
    complete_mutation_hierarchy = {}
    
    # 添加backbone mutations到层次结构中
    for group_id, backbone_mut in backbones_of_group.items():
        complete_mutation_hierarchy[backbone_mut] = {
            'main_group': group_id,
            'sub_group': f"{group_id}_backbone",
            'is_backbone': True,
            'is_sub_backbone': False
        }
    
    # 添加子分组的backbone mutations
    for full_subgroup_id, sub_backbone in subgroup_backbone_mutations.items():
        main_group_id = int(full_subgroup_id.split('_')[0])
        complete_mutation_hierarchy[sub_backbone] = {
            'main_group': main_group_id,
            'sub_group': full_subgroup_id,
            'is_backbone': False,
            'is_sub_backbone': True
        }
    
    # 添加其他突变
    for mut, full_subgroup_id in mutation_subgroups.items():
        if mut not in complete_mutation_hierarchy:  # 避免覆盖backbone mutations
            main_group_id = int(full_subgroup_id.split('_')[0])
            complete_mutation_hierarchy[mut] = {
                'main_group': main_group_id,
                'sub_group': full_subgroup_id,
                'is_backbone': False,
                'is_sub_backbone': False
            }
    
    # 保存分组结果和CV统计
    df_mutation_hierarchy = pd.DataFrame.from_dict(complete_mutation_hierarchy, orient='index')
    df_mutation_hierarchy.reset_index(inplace=True)
    df_mutation_hierarchy.rename(columns={'index': 'mutation'}, inplace=True)
    df_mutation_hierarchy.to_csv(os.path.join(outputpath, f"{sampleid}.mutation_hierarchy.csv"), index=False)
    
    if not df_cv_stats_within_subgroup.empty:
        df_cv_stats_within_subgroup.to_csv(os.path.join(outputpath, f"{sampleid}.subgroup_cv_statistics.csv"), index=False)
    
    # 生成分组报告（包含高CV突变信息）
    generate_subgrouping_report(complete_mutation_hierarchy, backbones_of_group, subgroup_backbone_mutations, 
                              df_cv_stats_within_subgroup, outputpath, sampleid, logger, high_cv_mutations)
    
    logger.info(f"Sub-grouping completed: {len(complete_mutation_hierarchy)} mutations in hierarchy")
    logger.info(f"External mutations (failed to add to tree): {len(external_mutations)}")
    logger.info(f"High CV mutations (excluded from analysis): {len(high_cv_mutations)}")
    if high_cv_mutations:
        logger.info(f"High CV mutations list: {high_cv_mutations}")
    
    # 修改返回语句，添加high_cv_mutations
    return complete_mutation_hierarchy, subgroup_backbone_mutations, subgroup_details, df_cv_stats_within_subgroup, T_current, M_current, root_mutations, external_mutations, high_cv_mutations




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
    outputpath: str,
    sampleid: str,
    immune_mutations: list,
    df_celltype: Optional[str] = None
) -> Tuple[Any, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Main function for building scaffold tree from mutation data.
    
    Parameters:
    -----------
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
    df_reads : pd.DataFrame
        Reads dataframe
    params : Dict[str, Any]
        Parameters dictionary
    outputpath : str
        Output directory path
    sampleid : str
        Sample ID
    df_celltype : pd.DataFrame
        Celltype dataframe
        
    Returns:
    --------
    Tuple containing:
        - T_scaffold: Final scaffold tree
        - M_scaffold: Final imputed mutation matrix
        - df_flipping_spots: Dataframe of flipping spots
        - df_total_flipping_count: Dataframe of total flipping counts
        - final_cleaned_bin_withNA3: Final cleaned binary matrix with NA=3
    """
    
    # Set up logger
    logger = logging.getLogger(__name__)
    
    # ------------------------------
    # Step 1: Check Celltype Data
    # ------------------------------
    logger.info("Running scaffold builder ...")
    
    logger.info(f"Celltype data loaded: {df_celltype.shape[0]} cells")
    
    if is_filter_quality=="yes":
        # ------------------------------
        # Step 2: Initial Filtration
        # ------------------------------
        logger.info("Starting initial filtration (Step 3.1) ...")
        kept_cells, kept_mutations, P_sub, V_sub, A_sub, C_sub, I_sub = initial_filter(P_somatic, V_somatic, A_somatic, C_somatic, I_somatic, params)
        
        # ------------------------------
        # Step 3: Coverage-Based Filtration
        # ------------------------------
        logger.info("Running coverage-based filtration (Step 3.2) ...")
        scaffold_mutations, df_summary_filtration = coverage_filters(kept_mutations, df_reads_somatic, df_celltype, params, outputpath)
        
        if not scaffold_mutations:
            raise ValueError("No mutations passed coverage filtration after removing germline variants")
        
        logger.info(f"Final scaffold set has {len(scaffold_mutations)} mutations")
    else:
        P_sub, V_sub, A_sub, C_sub, I_sub = P_somatic.copy(), V_somatic.copy(), A_somatic.copy(), C_somatic.copy(), I_somatic.copy()
        scaffold_mutations = list(I_somatic.columns)
    
    # ------------------------------
    # Step 4: Create Final Scaffold Matrices
    # ------------------------------
    logger.info("Creating final scaffold matrices ...")
    I_scaffold = I_sub[scaffold_mutations].copy()
    P_scaffold = P_sub[scaffold_mutations].copy()
    V_scaffold = V_sub[scaffold_mutations].copy()
    A_scaffold = A_sub[scaffold_mutations].copy()
    C_scaffold = C_sub[scaffold_mutations].copy()
    df_reads_scaffold = df_reads_somatic.loc[I_scaffold.index.insert(0, 'bulk'), scaffold_mutations].copy()
    
    # ------------------------------
    # Step 5: Resolve Immune Mutations
    # ------------------------------
    logger.info("Resolving immune mutations ...")
    I_resolved, P_resolved, V_resolved, A_resolved, C_resolved, df_reads_resolved, spots_to_split = resolved_spots_by_immune_mutations(
        I_scaffold, immune_mutations, P_scaffold, V_scaffold, A_scaffold, C_scaffold, df_reads_scaffold, p_threshold=0.5
    )
    I_somatic_resolved, P_somatic_resolved = split_spots_by_immune_mutations_scaffold(spots_to_split, [i for i in immune_mutations if i in I_somatic.columns], I_somatic, P_somatic)
    I_somatic_resolved_withNA3 = I_somatic_resolved.replace({np.nan: 3}).astype(int)
    I_somatic_resolved_withNA3.to_csv(os.path.join(outputpath, "I_somatic_resolved_withNA3.txt"), sep="\t")
    final_cleaned_I_somatic_resolved = I_somatic_resolved.loc[:, (I_somatic_resolved == 1).any(axis=0)]
    final_cleaned_I_somatic_resolved = final_cleaned_I_somatic_resolved.loc[(final_cleaned_I_somatic_resolved == 1).any(axis=1)]
    final_cleaned_I_somatic_resolved_withNA3 = final_cleaned_I_somatic_resolved.replace({np.nan: 3}).astype(int)
    final_cleaned_I_somatic_resolved_withNA3.to_csv(os.path.join(outputpath, "final_cleaned_I_somatic_resolved_withNA3.txt"), sep="\t")
    
    # ------------------------------
    # Step 6: Mutation Grouping with Leiden Algorithm
    # ------------------------------
    logger.info("Performing mutation grouping using Leiden algorithm ...")
    if is_filter_quality=="yes":
        cutoff_mcf_for_graph = 0.05
        cutoff_mcn_for_graph = 5
    else:
        cutoff_mcf_for_graph = 0
        cutoff_mcn_for_graph = 0
    
    clone_weights, pair_weights = get_correlation_graph_elements(I_resolved, 100, 42, cutoff_mcf_for_graph, cutoff_mcn_for_graph)
    mutation_group, partition, G_ig = leiden_mutation_groups(clone_weights, pair_weights, outputpath + "/" + sampleid + ".graph_for_mut_grouping.pdf", params['resolution_of_graph'])
    group_mutations = list(mutation_group.keys())
    I_selected_and_sorted, mut_df_sorted, group_to_muts, final_order = sort_I_hierarchical_freeze_ones_fixed(I_resolved, mutation_group)
    df_celltype_sub = df_celltype[df_celltype['barcode'].isin(I_selected_and_sorted.index)].copy()
    plot_heatmap_with_celltype_by_your_sorting(I_selected_and_sorted, df_celltype_sub, mutation_group, list(mut_df_sorted['mutation']), os.path.join(outputpath, sampleid+".heatmap_with_celltype_right_in_I_selected_and_sorted_after_graph_grouping.pdf"))
    logger.info(f"Total mutation groups: {len(mutation_group)}")
    
    # ------------------------------
    # Step 7: Select Backbone Mutations
    # ------------------------------
    logger.info("Selecting backbone mutations ...")
    backbones_of_group = select_founder_mutations(I_selected_and_sorted, mutation_group)
    backbone_mutations = list(backbones_of_group.values())
    # Check if the group assigned for mutation is correct
    mutation_list_under_backbone_mutations = {backbone: [mutation for mutation, value in mutation_group.items() if value == mutation_group[backbone]] for backbone in backbone_mutations}
    misgroup_mutations = []
    for check_mut in [i for i in list(mutation_group.keys()) if i not in backbone_mutations]:
        graph_key = next(key for key, mutations in mutation_list_under_backbone_mutations.items() if check_mut in mutations)
        best_backbone, intersection_counts = find_best_backbone_for_new_mutation_scaffold(mutation_list_under_backbone_mutations, I_somatic_resolved, I_somatic_resolved, check_mut)
        if graph_key != best_backbone:
            misgroup_mutations.append(check_mut)
    
    checked_mutation_group = {k: v for k, v in mutation_group.items() if k not in misgroup_mutations}
    mutation_list_under_backbone_mutations = {backbone: [mutation for mutation, value in checked_mutation_group.items() if value == checked_mutation_group[backbone]] for backbone in backbone_mutations}
    
    logger.info(f"Found backbone mutations: {list(backbones_of_group.values())}")
    
    # ------------------------------
    # Step 8: Build Backbone Tree
    # ------------------------------
    dir_backbone = os.path.join(outputpath, "phylo_backbone_tree")
    if not os.path.exists(dir_backbone):
        os.makedirs(dir_backbone)
    
    logger.info("Building backbone tree ...")
    T_B = build_backbone_tree(backbone_mutations)
    I_B = I_selected_and_sorted[backbone_mutations]
    I_B_withNA3 = I_B.replace({np.nan: 3}).astype(int)
    I_B_withNA3.to_csv(os.path.join(dir_backbone, "I_B_withNA3.txt"), sep="\t")
        
    M_B, cell_assignments = impute_backbone_clones(I_selected_and_sorted, backbone_mutations, checked_mutation_group)
    
    scp.ul.is_conflict_free_gusfield(M_B)
    WriteTfile(os.path.join(dir_backbone, "M_backbone_basedPivots.filtered_sites_inferred"), 
               M_B, M_B.index.tolist(), M_B.columns.tolist(), judge="yes")
    
    # Clean M_B: remove all zeros columns(muts) or rows(cells)
    final_cleaned_M_B = M_B.loc[:, (M_B != 0).any(axis=0)]  # 移除全0列
    final_cleaned_M_B = final_cleaned_M_B.loc[(final_cleaned_M_B != 0).any(axis=1)]  # 移除全0行
    WriteTfile(os.path.join(dir_backbone, "final_cleaned_M_B_basedPivots.filtered_sites_inferred"), 
               final_cleaned_M_B, final_cleaned_M_B.index.tolist(), final_cleaned_M_B.columns.tolist(), judge="yes")
    kept_rows = final_cleaned_M_B.index
    kept_cols = final_cleaned_M_B.columns
    final_cleaned_I_B_withNA3 = I_B_withNA3.loc[kept_rows, kept_cols]
    final_cleaned_I_B_withNA3.to_csv(os.path.join(dir_backbone, "final_cleaned_I_B_withNA3_for_circosPlot.txt"), sep="\t")
    
    # ------------------------------
    # Step 8: Sub-grouping within Backbone Groups
    # Step 9: Calculate Penalty and Refine Placement
    # ------------------------------
    
    # 为后续的操作准备数据
    I_selected = I_selected_and_sorted[list(checked_mutation_group.keys())]
    P_selected = P_resolved.loc[I_selected.index, I_selected.columns]
    V_selected = V_resolved.loc[I_selected.index, I_selected.columns]
    A_selected = A_resolved.loc[I_selected.index, I_selected.columns]
    C_selected = C_resolved.loc[I_selected.index, I_selected.columns]
    df_reads_selected = df_reads_resolved.loc[I_selected.index.insert(0, 'bulk'), I_selected.columns]
    
    logger.info("Calculating penalties and refining placements ...")
    T_current = copy.deepcopy(T_B)
    M_current = M_B.copy()
    M_current.insert(0, 'ROOT', 1)
    
    # 计算罚分时的 NA 权重设置
    ω_NA = params['general_weight_NA'] if params['general_weight_NA'] else 0.001
    fnfp_ratio = params['fnfp_ratio']
    φ = params['phi']
    
    # Refining positions for each mutation
    group_but_backbone_mutations = [i for i in list(checked_mutation_group.keys()) if i not in backbone_mutations]
    all_mutations = I_somatic.columns.tolist()
    
    # resorted_I_selected_and_sorted, sorting_stats_of_resorted_I_selected_and_sorted = reorder_columns_by_mutant_stats(I_selected_and_sorted, df_features_new)
    # sorted_selected_but_backbone_mutations = [i for i in resorted_I_selected_and_sorted.columns if i in group_but_backbone_mutations]
    
    ##### 鉴定 hub clusters 并且同时利用 penalty function 实现突变挂树
    root_mutations = []
    sorted_I_resolved, sorting_stats_of_I_resolved = reorder_columns_by_mutant_stats(I_resolved, df_features_new)
    sorted_P_resolved, sorting_stats_of_P_resolved = reorder_columns_by_mutant_stats(P_resolved, df_features_new)
    results_of_subgrouping = perform_subgrouping_within_backbone_groups_and_build_initial_scaffold_tree(
        sorted_I_resolved, 
        sorted_P_resolved, 
        T_current, 
        M_current, 
        checked_mutation_group, 
        backbones_of_group, 
        df_reads_resolved, 
        df_features_new, 
        outputpath, 
        sampleid, 
        logger, 
        params, 
        M_B, 
        root_mutations, 
        cutoff_mcf_for_graph, 
        cutoff_mcn_for_graph
    )
    complete_mutation_hierarchy, subgroup_backbone_mutations, subgroup_details, df_cv_stats_within_subgroup, T_current, M_current, root_mutations, external_mutations, high_cv_mutations = results_of_subgrouping
        
    # ------------------------------
    # Step 10: Calculate Penalty and Refine Placement
    # ------------------------------
    # 能放到树上的就是属于 scaffold mutations 的了
    # 如果还有未处理的突变，再次处理
    remained_mutations = []
    if len(external_mutations) > 0:
        logger.info(f"Reprocessing {len(external_mutations)} external mutations...")
        remained_mutations, T_current, M_current = integrate_mutations_to_scaffold_within_group(
            sorted_attached_mutations=external_mutations,
            T_current=T_current,
            M_current=M_current,
            I_attached=I_selected,
            P_attached=P_selected,
            mutation_group=checked_mutation_group,
            ω_NA=ω_NA,
            fnfp_ratio=fnfp_ratio,
            φ=φ,
            logger=logger
        )
    else:
        logger.info("There is no external scaffold mutations after first building.")
    
    # Final scaffold tree and matrix    
    T_scaffold = copy.deepcopy(T_current)
    M_scaffold_initial = M_current.copy()
    M_scaffold_initial = M_scaffold_initial.drop(columns=['ROOT'], errors='ignore')
    mutations_on_T_scaffold = M_scaffold_initial.columns.to_series().apply(lambda x: x.split("|")).explode().unique().tolist()
    M_scaffold = split_merged_columns(M_scaffold_initial, mutations_on_T_scaffold)
    
    logger.info("Final scaffold tree:")
    print_tree(T_scaffold)
    logger.info(f"Final scaffold matrix shape: {M_scaffold.shape}")
    
    # ------------------------------
    # Step 11: Output Results
    # ------------------------------
    logger.info("Outputting results ...")
    
    # Create output directory
    dir_scaffold = os.path.join(outputpath, "phylo_scaffold_tree")
    if not os.path.exists(dir_scaffold):
        os.makedirs(dir_scaffold)
    
    # Prepare binary matrix with NA=3
    I_selected_withNA3 = I_selected.replace({np.nan: 3}).astype(int)
    I_selected_withNA3.to_csv(os.path.join(dir_scaffold, "I_scaffold_withNA3.txt"), sep="\t")
    WriteTfile(os.path.join(dir_scaffold, "M_scaffold_basedPivots.filtered_sites_inferred"), 
               M_scaffold, M_scaffold.index.tolist(), M_scaffold.columns.tolist(), judge="yes")
    
    # Clean M_scaffold: remove all zeros columns(muts) or rows(cells)
    final_cleaned_M_scaffold = M_scaffold.loc[:, (M_scaffold != 0).any(axis=0)]  # 移除全0列
    final_cleaned_M_scaffold = final_cleaned_M_scaffold.loc[(final_cleaned_M_scaffold != 0).any(axis=1)]  # 移除全0行
    
    # 获取保留的行列名
    kept_rows = final_cleaned_M_scaffold.index
    kept_cols = final_cleaned_M_scaffold.columns
    
    # 从 I_selected_withNA3 提取
    final_cleaned_I_selected_withNA3 = I_selected_withNA3.loc[kept_rows, kept_cols]
    
    WriteTfile(os.path.join(dir_scaffold, "final_cleaned_M_scaffold_basedPivots.filtered_sites_inferred"), 
               final_cleaned_M_scaffold, final_cleaned_M_scaffold.index.tolist(), final_cleaned_M_scaffold.columns.tolist(), judge="yes")
    final_cleaned_I_selected_withNA3.to_csv(os.path.join(dir_scaffold, "final_cleaned_I_scaffold_withNA3_for_circosPlot.txt"), sep="\t")
    
    # ------------------------------
    # Step 12: Identify Flipping Spots
    # ------------------------------
    logger.info("Identifying flipping spots ...")
    df_I_selected_withNA3_for_flipping = final_cleaned_I_selected_withNA3.copy()
    df_phylogeny = final_cleaned_M_scaffold.copy()
    
    # get false_negative_flipping spots
    false_negative_flipping_spots = df_I_selected_withNA3_for_flipping.apply(
        lambda col: find_flipping_spots(col, df_phylogeny[col.name], condition_in_bin=0, condition_phylogeny=1)
    )
    
    # get NAto1 spots
    NAto1_flipping_spots = df_I_selected_withNA3_for_flipping.apply(
        lambda col: find_flipping_spots(col, df_phylogeny[col.name], condition_in_bin=3, condition_phylogeny=1)
    )
    
    # get false_positive_flipping spots
    false_positive_flipping_spots = df_I_selected_withNA3_for_flipping.apply(
        lambda col: find_flipping_spots(col, df_phylogeny[col.name], condition_in_bin=1, condition_phylogeny=0)
    )
    
    # get NAto0 spots
    NAto0_flipping_spots = df_I_selected_withNA3_for_flipping.apply(
        lambda col: find_flipping_spots(col, df_phylogeny[col.name], condition_in_bin=3, condition_phylogeny=0)
    )
    
    # process na list
    if false_negative_flipping_spots.empty:
        false_negative_flipping_spots = {col: [] for col in df_I_selected_withNA3_for_flipping.columns}
    
    if false_positive_flipping_spots.empty:
        false_positive_flipping_spots = {col: [] for col in df_I_selected_withNA3_for_flipping.columns}
    
    if NAto1_flipping_spots.empty:
        NAto1_flipping_spots = {col: [] for col in df_I_selected_withNA3_for_flipping.columns}
    
    if NAto0_flipping_spots.empty:
        NAto0_flipping_spots = {col: [] for col in df_I_selected_withNA3_for_flipping.columns}
    
    # get flipping_spots dataframe
    df_flipping_spots = pd.DataFrame({
        'Mutation': df_I_selected_withNA3_for_flipping.columns,
        'false_negative_flipping_spots': [', '.join(false_negative_flipping_spots.get(col, [])) for col in df_I_selected_withNA3_for_flipping.columns],
        'false_positive_flipping_spots': [', '.join(false_positive_flipping_spots.get(col, [])) for col in df_I_selected_withNA3_for_flipping.columns],
        'NAto1_flipping_spots': [', '.join(NAto1_flipping_spots.get(col, [])) for col in df_I_selected_withNA3_for_flipping.columns],
        'NAto0_flipping_spots': [', '.join(NAto0_flipping_spots.get(col, [])) for col in df_I_selected_withNA3_for_flipping.columns]
    })
    df_flipping_spots.to_csv(os.path.join(dir_scaffold, "df_flipping_spots.txt"), sep="\t", index=False)
    
    # ------------------------------
    # Step 12: Calculate Total Flipping Counts
    # ------------------------------
    logger.info("Calculating total flipping counts ...")
    total_FN_flipping = ((df_I_selected_withNA3_for_flipping == 0) & (df_phylogeny == 1)).sum().sum()
    total_FP_flipping = ((df_I_selected_withNA3_for_flipping == 1) & (df_phylogeny == 0)).sum().sum()
    total_NAto0 = ((df_I_selected_withNA3_for_flipping == 3) & (df_phylogeny == 0)).sum().sum()
    total_NAto1 = ((df_I_selected_withNA3_for_flipping == 3) & (df_phylogeny == 1)).sum().sum()
    
    logger.info(f"Total False Negative flipping count: {total_FN_flipping}")
    logger.info(f"Total False Positive flipping count: {total_FP_flipping}")
    logger.info(f"Total NA to 0 flipping count: {total_NAto0}")
    logger.info(f"Total NA to 1 flipping count: {total_NAto1}")
    
    df_total_flipping_count = pd.DataFrame({
        'total_flipping_False_Negative': [total_FN_flipping],
        'total_flipping_False_Positive': [total_FP_flipping],
        'total_flipping_NA_to_0': [total_NAto0],
        'total_flipping_NA_to_1': [total_NAto1]
    })
    df_total_flipping_count.to_csv(os.path.join(dir_scaffold, "df_total_flipping_count.txt"), sep="\t", index=False)
    
    logger.info("Scaffold building completed successfully!")
    
    # 将 TreeNode 转换为字典格式并且保存为 JSON 文件
    tree_dict = tree_to_dict(T_scaffold)
    with open(os.path.join(dir_scaffold, 'T_scaffold.json'), 'w') as f:
        json.dump(tree_dict, f, indent=4)
    
    # 直接返回字符串文本格式并且保存
    T_scaffold.save_to_file(os.path.join(dir_scaffold, 'T_scaffold.txt'))
    
    return T_scaffold, M_scaffold, df_flipping_spots, df_total_flipping_count, final_cleaned_I_selected_withNA3, final_cleaned_M_scaffold, backbone_mutations, checked_mutation_group, spots_to_split, list(checked_mutation_group.keys()), remained_mutations, high_cv_mutations






























#!/usr/bin/env python3

"""
germline_filter.py

Identify putative heterozygous germline variants
based on Methods Section 2 (MCF, J_r, U_r, S_r^FP, gap detection).

Inputs (from data_loader.load_all):
  - P : posterior matrix (cells x muts)
  - V : mutant allele frequency matrix (cells x muts)
  - C : coverage matrix (cells x muts)
  - A : mutant allele count matrix (cells x muts)

Main function:
  identify_germline_variants(P, V, C, A, mcf_cutoff=0.05)

Returns:
  A Python set of mutation IDs (columns in P/V/C/A) identified as germline.
"""
import numpy as np
import pandas as pd
from typing import Set, Dict, List, Tuple, Optional
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.colors import ListedColormap
import logging

logger = logging.getLogger(__name__)

# -------------------------
# Helper functions
# -------------------------

def calculate_mutant_fraction(I_detected, df_reads_detected):
    # 计算每列中 mutant_cell_number（1的个数）
    mutant_cell_number = (I_detected == 1).sum(axis=0)
    
    # 计算每列中非 NaN 元素的个数
    coverage_cell_number = df_reads_detected.notna().sum(axis=0)
    
    # 计算 mutant_cell_fraction (mutant_cell_number / coverage_cell_number)
    mutant_cell_fraction = mutant_cell_number / coverage_cell_number
    
    return mutant_cell_number, mutant_cell_fraction

def build_binary_I(P: pd.DataFrame, V: pd.DataFrame, C: pd.DataFrame, p_thresh: float = 0.5) -> pd.DataFrame:
    """二值化矩阵 I，根据公式(1)"""
    I = pd.DataFrame(np.nan, index=P.index, columns=P.columns)
    covered = (C >= 1)
    mask_mut = (V > 0) & (P > p_thresh) & covered
    I[mask_mut] = 1
    mask_ref = (V == 0) & (P <= p_thresh) & covered
    I[mask_ref] = 0
    # 对 coverage>=1 但未赋值的 NA 置为 0
    I[covered & I.isna()] = 0
    return I

def pairwise_counts(I: pd.DataFrame, j1: str, j2: str) -> Dict[str,int]:
    """计算两个突变的 N11, N10, N01, N00（排除 NA）"""
    a = I[j1]; b = I[j2]
    valid = (~a.isna()) & (~b.isna())
    a = a[valid]; b = b[valid]
    N11 = int(((a==1)&(b==1)).sum())
    N10 = int(((a==1)&(b==0)).sum())
    N01 = int(((a==0)&(b==1)).sum())
    N00 = int(((a==0)&(b==0)).sum())
    return dict(N11=N11, N10=N10, N01=N01, N00=N00)

def pairwise_counts_for_two_columns(col1: pd.Series, col2: pd.Series) -> Dict[str, int]:
    """计算两个列的 N11, N10, N01, N00（排除 NA）"""
    # 确保两个列的索引一致
    if not col1.index.equals(col2.index):
        raise ValueError("Columns do not have the same index")
    # 过滤掉 NA 值
    valid = (~col1.isna()) & (~col2.isna())
    col1 = col1[valid]
    col2 = col2[valid]
    # 计算 N11, N10, N01, N00
    N11 = int(((col1 == 1) & (col2 == 1)).sum())
    N10 = int(((col1 == 1) & (col2 == 0)).sum())
    N01 = int(((col1 == 0) & (col2 == 1)).sum())
    N00 = int(((col1 == 0) & (col2 == 0)).sum())
    return dict(N11=N11, N10=N10, N01=N01, N00=N00)

def jaccard_index(I, j1, j2):
    """对称的Jaccard """
    counts = pairwise_counts(I, j1, j2)
    N11, N10, N01 = counts['N11'], counts['N10'], counts['N01']
    denominator = N11 + N10 + N01
    return N11 / denominator if denominator > 0 else 0.0

def f_fraction(I: pd.DataFrame, j1: str, j2: str) -> float:
    """f(j1,j2) = N11 / |S(j1)|"""
    counts = pairwise_counts(I,j1,j2)
    N11 = counts['N11']
    S1 = int((I[j1]==1).sum())
    return N11/S1 if S1>0 else 0.0

def are_mutations_correlated(I: pd.DataFrame, j1: str, j2: str) -> bool:
    """
    判断两个突变是否相关，根据标准：
    - N11(j1,j2) ≥ 3 ∧ J(j1,j2) ≥ 0.2
    - 或 N11(j1,j2) ≥ 3 ∧ 0 < J(j1,j2) < 0.2 ∧ max(f(j1,j2), f(j2,j1)) ≥ 0.9
    """
    N11_threshold = 1
    J_val_threshold = 0.08
    f_fraction_threshold = 0.5
    counts = pairwise_counts(I, j1, j2)
    N11 = counts['N11']
    # 首先检查N11是否满足最小细胞数要求
    if N11 < N11_threshold:
        return False
    # 计算Jaccard指数
    J_val = jaccard_index(I, j1, j2)
    # 第一个条件：Jaccard指数 ≥ 0.2
    if J_val >= J_val_threshold:
        return True
    # 第二个条件：0 < Jaccard指数 < 0.2 且 max(f(j1,j2), f(j2,j1)) ≥ 0.9
    if 0 < J_val < J_val_threshold:
        f_j1j2 = f_fraction(I, j1, j2)
        f_j2j1 = f_fraction(I, j2, j1)
        if max(f_j1j2, f_j2j1) >= f_fraction_threshold:
            return True
    return False

def build_J_r(I: pd.DataFrame, r: str) -> Set[str]:
    """
    构建与参考突变r相关的突变集合 J_r
    J_r = { j≠r | j与r相关 }
    """
    J = set()
    for j in I.columns:
        if j == r:
            continue
        if are_mutations_correlated(I, r, j):
            J.add(j)
    return J

def infer_U_r(I: pd.DataFrame, r: str, J_r: Set[str]) -> Set[str]:
    """U_r = { cells with I[r]=1 or q_i ≥ q_threshold }"""
    if len(J_r) < 3:
        return set(I.index[I[r] == 1]), 'unknown'
    # 计算每个细胞在 J_r 上的 mutant count
    row_sums = I[list(J_r)].apply(lambda x: x[x == 1].count(), axis=1)
    # 如果整体上限太低（q_max <= 2）或者低值占比较高，就不扩展
    value_counts = row_sums.value_counts().sort_index()
    low_counts_ratio = (value_counts.get(1,0) + value_counts.get(2,0)) / len(row_sums)
    if row_sums.max() <= 2 or low_counts_ratio > 0.5:
        # 不扩展，只取 r=1 的细胞
        return set(I.index[I[r] == 1]), low_counts_ratio
    # 否则使用累积百分比阈值方法
    value_counts_desc = row_sums.value_counts().sort_index(ascending=False)
    cumulative_percent = value_counts_desc.cumsum() / len(row_sums) * 100
    for i, (value, percent) in enumerate(cumulative_percent.items()):
        if percent > 10:
            if i == 0:
                q_threshold = value + 1
            else:
                q_threshold = list(cumulative_percent.index)[i-1]
            break
    else:
        q_threshold = cumulative_percent.index[-1]
    # 构建 mask
    mask = (I[r] == 1) | (row_sums >= q_threshold)
    return set(I.index[mask]), low_counts_ratio

def compute_S_r_FP(I: pd.DataFrame, r: str) -> float:
    """计算 S_r^FP = 平均 S_j^FP(r) over J_r_plus"""
    J_r = build_J_r(I,r)
    J_plus = set(J_r)|{r}
    if not J_plus:
        return 0.0
    U_r, low_counts_ratio = infer_U_r(I,r,J_r)
    scores=[]
    for j in I.columns:
        if j==r: continue
        S_j = set(I.index[I[j]==1])
        if len(S_j)==0: 
            continue
        N_in = len([c for c in S_j if c in U_r])
        N_out = len(S_j)-N_in
        if j in J_plus:
            N_FP = N_out
        else:
            N_FP = N_in
        scores.append(N_FP/len(S_j))
    mean_score = float(np.mean(scores)) if scores else 0.0
    std_score = float(np.std(scores, ddof=0)) if scores else 0.0
    cv_score = std_score / mean_score if mean_score > 0 else 0.0
    return mean_score, std_score, cv_score, U_r, low_counts_ratio


def plot_heatmap_with_germline_mutations(I_raw, germline_mutations, pdf_file):
    """
    绘制带有突变矩阵的热图，germline_mutations 中的突变放在最左边，横坐标刻度标红，带有行/列突变数统计的柱状图，并在下方显示 legend。

    Parameters
    ----------
    I_raw : pd.DataFrame
        原始突变矩阵 (cell x mutation)，元素为 {0,1,NA}
    germline_mutations : set
        包含 germline 突变的集合，将其放在热图的最左边
    pdf_file : str
        保存图像的 PDF 文件路径
    """
    
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap, BoundaryNorm
    from matplotlib.patches import Patch
    
    # -------------------
    # Step 1: 移动 germline_mutations 到最左边
    # -------------------
    germline_mutations_list = list(germline_mutations)
    # 确保 germline_mutations 只包含在 I_raw 的列中
    germline_mutations_in_data = [mut for mut in germline_mutations_list if mut in I_raw.columns]
    
    # 将 germline_mutations 放到数据框的最左边
    I_sorted = I_raw[germline_mutations_in_data + [col for col in I_raw.columns if col not in germline_mutations_in_data]]
    
    # -------------------
    # Step 2: 转换矩阵数值（NA → 2，用于单独上色）
    # -------------------
    I_numeric = I_sorted.fillna(np.nan).apply(pd.to_numeric, errors="coerce")
    I_plot = I_numeric.copy()
    I_plot = I_plot.where(~I_plot.isna(), 2)  # 把 NA 填为 2
    
    # -------------------
    # Step 3: 计算行/列突变数量统计（用于条形图）
    # -------------------
    row_sums = I_numeric.sum(axis=1, skipna=True)  # 每个 cell 的突变数
    col_sums = I_numeric.sum(axis=0, skipna=True)  # 每个突变被多少 cell 支持
    
    # -------------------
    # Step 4: 布局 GridSpec
    #   - 左边：行条形图
    #   - 中间：热图
    #   - 上面：列条形图
    # -------------------
    fig = plt.figure(figsize=(12, 10))
    gs = fig.add_gridspec(5, 6,
                          width_ratios=[0.3, 0.3, 3, 0.05, 0.05, 0.05],
                          height_ratios=[0.5, 0.05, 3, 0.3, 0.3],
                          wspace=0.05, hspace=0.05)
    
    ax_row_bar = fig.add_subplot(gs[2, 0])   # 左边行条形图
    ax_heatmap = fig.add_subplot(gs[2, 2])   # 中间热图
    ax_col_bar = fig.add_subplot(gs[0, 2])   # 上方列条形图
    ax_dummy = fig.add_subplot(gs[0, 0]); ax_dummy.axis("off")  # 占位
    
    # -------------------
    # Step 5: 绘制热图
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
    # Step 6: 给横坐标 mutation label 上色（germline_mutations 的 label 使用红色）
    # -------------------
    for label in ax_heatmap.get_xticklabels():
        mut_name = label.get_text()
        if mut_name in germline_mutations_in_data:
            label.set_color('red')
        else:
            label.set_color('black')
    
    # -------------------
    # Step 7: 列条形图（每个突变被多少 cell 支持）
    # -------------------
    ax_col_bar.bar(range(len(col_sums)), col_sums.values,
                   color="#7D2224", alpha=0.7, align="center")
    ax_col_bar.set_xlim(ax_heatmap.get_xlim())
    ax_col_bar.set_xticks([])
    ax_col_bar.tick_params(axis="y", labelsize=8)
    ax_col_bar.set_ylabel("#Mutations", fontsize=10)
    
    # -------------------
    # Step 8: 行条形图（每个 cell 有多少突变）
    # -------------------
    ax_row_bar.barh(range(len(row_sums)), row_sums.values,
                    color="#7D2224", alpha=0.7, align="center")
    ax_row_bar.set_ylim(ax_heatmap.get_ylim())
    ax_row_bar.set_yticks([])
    ax_row_bar.set_xlabel("#Cells", fontsize=10)
    ax_row_bar.invert_xaxis()
    
    # -------------------
    # Step 9: 添加 Legend
    #   - 突变值 (0,1,NA)
    # -------------------
    # 突变值 legend
    heatmap_handles = [Patch(facecolor=c, label=l) 
                       for c, l in zip(["#D4E8F0", "#7D2224", "white"],
                                       ["0 (No Mutation)", "1 (Mutation)", "NA (Missing)"])]
    fig.legend(handles=heatmap_handles, loc="lower center", ncol=3,
               bbox_to_anchor=(0.5, -0.03), frameon=False, fontsize=9,
               title="Mutation Values", title_fontsize=10)
    
    # -------------------
    # Step 10: 保存图像
    # -------------------
    plt.suptitle("Heatmap of Mutations with Germline Mutations Highlighted", fontsize=14, y=0.95)
    plt.tight_layout()
    plt.subplots_adjust(bottom=0)
    plt.savefig(pdf_file, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"保存完成: {pdf_file}")

def update_germline_status(final_df, I, mcn_cutoff=5):
    # 筛选出 'germline_determined' 列为 'germline' 的行
    germline_mutations = final_df[final_df['germline_determined'] == 'germline']
    
    for idx, row in germline_mutations.iterrows():
        # 获取 mutation id
        mutation_id = idx
        
        # 获取对应的 mutation 列
        mutation_column = I[mutation_id]
        
        # 计算该突变列中值为 1 的个数
        mutant_cell_count = np.sum(mutation_column == 1)
        
        # 如果小于 mcn_cutoff，更新为 'non-germline'
        if mutant_cell_count <= mcn_cutoff:
            final_df.loc[idx, 'germline_determined'] = 'non-germline'
    
    return final_df


def calculate_prob_threshold(probs_trimmed):
    mean = probs_trimmed.mean()
    std = probs_trimmed.std()
    
    # 计算初步的阈值
    threshold = mean + 2.6 * std
    
    # 如果计算值大于 0.9，则改为 mean + 2 * std
    if threshold > 0.9:
        print(f"阈值大于 0.9，改为 {mean + 2 * std} (mean + 2 * std)")
        threshold = mean + 2 * std
    
    # 如果计算值大于 0.9，则改为 mean + 1.5 * std
    if threshold > 0.9:
        print(f"阈值大于 0.9，改为 {mean + 1.5 * std} (mean + 1.5 * std)")
        threshold = mean + 1.5 * std
    
    # 如果计算值大于 0.9，则改为 mean + std
    if threshold > 0.9:
        print(f"阈值大于 0.9，改为 {mean + std} (mean + std)")
        threshold = mean + std
    
    # 如果还是大于 0.9，则设置为 0.9
    if threshold > 0.9:
        print("阈值依然大于 0.9，改为 0.9")
        threshold = 0.9
    
    # 如果小于 0.1，则设置为 0.1
    if threshold < 0.1:
        print("阈值小于 0.1，改为 0.1")
        threshold = 0.1
    
    return threshold


def update_features_matrix(I, df_reads, df_features, mcf_cutoff):
    """
    Build features matrix and identify candidate founders
    
    Parameters:
    -----------
    I : pandas.DataFrame
        Binary matrix
    df_reads : pandas.DataFrame
        Reads data
    df_features : pandas.DataFrame
        Existing features matrix
    mcf_cutoff : float
        Minimum mutant cell fraction cutoff for candidates
        
    Returns:
    --------
    df_features_new : pandas.DataFrame
        Updated features matrix with mutant cell fraction data
    """
    
    # Step 1. Filter detected mutations and cells
    rows_with_ones = (I == 1).any(axis=1)
    cols_with_ones = (I == 1).any(axis=0)
    I_detected = I.loc[rows_with_ones, cols_with_ones]
    df_reads_detected = df_reads.loc[I_detected.index, I_detected.columns]
    empty_mutations = I.columns[~cols_with_ones]
    if not empty_mutations.empty:
        logger.warning(f"The following mutations are not detected in any cell: {', '.join(empty_mutations)}")
    
    # Step 2. Calculate mutant fractions
    mutant_cell_number_detected, mutant_cell_fraction_detected = calculate_mutant_fraction(I_detected, df_reads_detected)
    mutant_cell_number_input, mutant_cell_fraction_input = calculate_mutant_fraction(I, df_reads)
    
    # Step 3. Build new features matrix
    df_features_new = pd.concat([df_features, 
                             pd.DataFrame([mutant_cell_fraction_detected], index=['mutant_cell_fraction_detected']),
                             pd.DataFrame([mutant_cell_fraction_input], index=['mutant_cell_fraction_input'])])
    
    return df_features_new, empty_mutations


def add_mutation_proportions_to_features(df_features, df_cells):
    """
    计算每个突变在细胞数据框中的0、1、NA比例，并添加到特征数据框中
    
    Parameters:
    -----------
    df_features : DataFrame
        突变特征数据框（df_features_new），用于存储结果
    df_cells : DataFrame
        细胞基因型数据框（I_attached），用于计算
    """
    zero_props = []
    one_props = []
    na_props = []
    
    for mutation in df_features.columns:
        if mutation in df_cells.columns:
            # 获取该突变列的值计数（包括NA）
            value_counts = df_cells[mutation].value_counts(dropna=False)
            total_cells = len(df_cells[mutation])
            
            # 计算各种值的比例
            zero_count = value_counts.get(0, 0)
            one_count = value_counts.get(1, 0)
            na_count = value_counts.get(np.nan, 0) if np.nan in value_counts.index else 0
            
            zero_prop = zero_count / total_cells if total_cells > 0 else 0
            one_prop = one_count / total_cells if total_cells > 0 else 0
            na_prop = na_count / total_cells if total_cells > 0 else 0
        else:
            # 如果突变不在df_cells中，设为0或NaN
            zero_prop = 0
            one_prop = 0
            na_prop = 1  # 或者设为1，表示完全缺失
    
        zero_props.append(zero_prop)
        one_props.append(one_prop)
        na_props.append(na_prop)
    
    # 将结果添加到df_features中
    df_features.loc['zero_prop_detected'] = zero_props
    df_features.loc['one_prop_detected'] = one_props
    df_features.loc['na_prop_detected'] = na_props
    
    return df_features


def reorder_columns_by_mutant_stats(df_values, df_features_new, 
                                    min_cell_threshold=30, bin_size=5, 
                                    descending=True, return_stats=True):
    """
    最优化的列重排序函数：按mutant cell number分组，组内按mutant cell fraction排序
    
    Parameters:
    -----------
    df_values : DataFrame
        包含0,1,NA的原始数据框 (rows: cells, columns: mutations)
    df_features_new : DataFrame
        包含突变统计信息的数据框
    min_cell_threshold : int
        最小细胞数阈值，大于等于此值的突变单独作为高优先级组
    bin_size : int
        阈值以下的分组间隔大小
    descending : bool
        True: 从大到小排序 (高mutant cell number在前)  
        False: 从小到大排序
    return_stats : bool
        是否返回排序统计信息
    
    Returns:
    --------
    df_reordered : DataFrame
        重新排序列后的数据框
    sorting_stats : DataFrame (可选)
        列的排序统计信息
    """
    
    # 1. 获取两个数据框列的交集
    common_columns = list(set(df_values.columns) & set(df_features_new.columns))
    print(f"原始df_values列数: {len(df_values.columns)}")
    print(f"原始df_features_new列数: {len(df_features_new.columns)}")
    print(f"共同列数: {len(common_columns)}")
    
    if len(common_columns) == 0:
        raise ValueError("两个数据框没有共同的列！")
    
    # 2. 筛选共同列
    df_values_common = df_values[common_columns]
    
    # 3. 提取关键统计信息（只针对共同列）
    mutant_cell_num = df_features_new[common_columns].loc['mutant_cellnum'].astype(int)
    mutant_cell_frac = df_features_new[common_columns].loc['mutant_cell_fraction'].astype(float)
    
    # 4. 创建排序统计DataFrame
    stats_df = pd.DataFrame({
        'column_name': mutant_cell_num.index,
        'mutant_cell_num': mutant_cell_num.values,
        'mutant_cell_frac': mutant_cell_frac.values
    })
    
    # 5. 定义分组逻辑
    def create_mutant_group(num):
        """创建mutant cell number分组标签"""
        if num >= min_cell_threshold:
            return f'≥{min_cell_threshold}'
        else:
            lower = (num // bin_size) * bin_size
            upper = lower + bin_size - 1
            return f'{lower:02d}-{upper:02d}'
    
    stats_df['mutant_group'] = stats_df['mutant_cell_num'].apply(create_mutant_group)
    
    # 6. 定义分组排序顺序
    # 高mutant cell number的组在前
    high_priority_groups = [f'≥{min_cell_threshold}']
    
    # 低mutant cell number的组，从大到小
    low_priority_groups = []
    for i in range(min_cell_threshold - bin_size, -1, -bin_size):
        lower = i
        upper = i + bin_size - 1
        if lower >= 0:
            low_priority_groups.append(f'{lower:02d}-{upper:02d}')
    
    group_order = high_priority_groups + low_priority_groups
    
    # 7. 转换为有序分类变量
    stats_df['mutant_group'] = pd.Categorical(
        stats_df['mutant_group'], 
        categories=group_order, 
        ordered=True
    )
    
    # 8. 排序：先按分组，组内按mutant cell fraction
    if descending:
        # 从大到小：高mutant number + 高fraction在前
        stats_df_sorted = stats_df.sort_values(
            ['mutant_group', 'mutant_cell_frac'], 
            ascending=[True, False]  # 分组用分类顺序，分数降序
        )
    else:
        # 从小到大：低mutant number + 低fraction在前  
        stats_df_sorted = stats_df.sort_values(
            ['mutant_group', 'mutant_cell_frac'], 
            ascending=[True, True]   # 分组用分类顺序，分数升序
        )
    
    # 9. 获取排序后的列名
    sorted_columns = stats_df_sorted['column_name'].tolist()
    
    # 10. 重新排列数据框列（只针对共同列）
    df_reordered = df_values_common[sorted_columns]
    
    # 11. 重置索引以便查看
    stats_df_sorted = stats_df_sorted.reset_index(drop=True)
    stats_df_sorted['final_order'] = stats_df_sorted.index + 1
    
    print(f"最终重排序列数: {len(sorted_columns)}")
    
    if return_stats:
        return df_reordered, stats_df_sorted
    else:
        return df_reordered




# -------------------------
# Main function
# -------------------------

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
from typing import Set, Optional, Tuple
from sklearn.linear_model import LogisticRegression
from matplotlib.backends.backend_pdf import PdfPages

def identify_germline_variants(
    P: pd.DataFrame, V: pd.DataFrame, C: pd.DataFrame, df_reads: pd.DataFrame, df_features_new: pd.DataFrame, 
    p_thresh: float = 0.5, mcf_cutoff: float = 0.05, mcn_cutoff: int = 5, 
    outputpath: Optional[str] = None,
    sampleid: Optional[str] = None,
    df_labeled: Optional[pd.DataFrame] = None
) -> Tuple[pd.DataFrame, Set[str]]:
    """
    Identify germline mutations via logistic regression on mean/std/cv,
    and visualize scatter plots.
    
    df_labeled: optional, labeled training dataset for logistic regression.
    """
    # Step 1. Build binary matrix
    I = build_binary_I(P, V, C, p_thresh)
    
    # Step 2. Candidate founder set    
    candidates = df_features_new.loc['mutant_cell_fraction_detected'][df_features_new.loc['mutant_cell_fraction_detected'] > mcf_cutoff].index.tolist()
    if not candidates:
        return pd.DataFrame(), set()
    
    # Step 3. Compute stats
    S_r_scores, S_r_std, S_r_cv, S_r_fn, S_r_lcr = {}, {}, {}, {}, {}
    for r in I.columns:
    # for r in candidates:
        mean_score, std_score, cv_score, U_r, low_counts_ratio = compute_S_r_FP(I, r)
        S_r_scores[r] = mean_score
        S_r_std[r] = std_score
        S_r_cv[r] = cv_score
        S_r_fn[r] = len(U_r) - len(I[I[r] == 1])
        S_r_lcr[r] = low_counts_ratio
    
    stats_df = pd.DataFrame({
        "FP_mean": pd.Series(S_r_scores),
        "FP_std": pd.Series(S_r_std),
        "FP_cv": pd.Series(S_r_cv),
        "FN_num": pd.Series(S_r_fn),
        "low_counts_ratio": pd.Series(S_r_lcr)
    })
    
    # Step 4. 拆分数据框
    # 过滤掉 FP_mean = FP_std = FP_cv = 0 的突变，放入 non_germline 数据框
    stats_df_non_germline = stats_df[(stats_df['FP_mean'] == 0) &
                                     (stats_df['FP_std'] == 0) &
                                     (stats_df['FP_cv'] == 0)].copy()
    stats_df_candidates = stats_df.drop(stats_df_non_germline.index)  # 剩下的突变
    
    if len(stats_df_candidates)==0:
        return pd.DataFrame(), set()
    
    # Step 5. 标记 non-germline 数据框
    stats_df_non_germline['germline_prob'] = np.nan  # 添加空列
    stats_df_non_germline['germline_pred'] = np.nan  # 添加空列
    stats_df_non_germline['germline_determined'] = 'non-germline'
    
    # Step 6. 对 stats_df_candidates 进行 Logistic Regression 处理
    pred_germline_mutations = set()
    prob_threshold = None
    
    if df_labeled is None:
        # 获取项目根目录（假设当前文件在 src/phylosolid/germline_filter/ 下）
        current_file_dir = os.path.dirname(os.path.abspath(__file__))
        # 向上到项目根目录：src/phylosolid/germline_filter -> src/phylosolid -> src -> 根目录
        project_root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(current_file_dir))))
        
        # 构建 resource 文件夹下的路径
        criteria_path = os.path.join(
            project_root,
            "resource",
            "stats_data_by_merged_3samples_withLabels.csv"
        )
        
        if os.path.exists(criteria_path):
            df_labeled = pd.read_csv(criteria_path)
            print(f"Loaded criteria file from: {criteria_path}")  # 可选的调试信息
        else:
            print(f"Warning: Criteria file not found at {criteria_path}")  # 可选的警告信息
    
    if df_labeled is not None and not df_labeled.empty:
        if sampleid is not None and "sampleid" in df_labeled.columns:
            df_labeled = df_labeled[df_labeled["sampleid"] != sampleid]
        
        if not df_labeled.empty:
            # 准备训练数据
            X_train = df_labeled[['FP_mean', 'FP_std', 'FP_cv']].values
            y_train = (df_labeled['label'] == 'germline').astype(int).values
            
            from sklearn.linear_model import LogisticRegression
            clf = LogisticRegression(class_weight="balanced", random_state=42)
            clf.fit(X_train, y_train)
            
            # 预测当前数据集
            X_test = stats_df_candidates[['FP_mean', 'FP_std', 'FP_cv']].values
            probs = clf.predict_proba(X_test)[:, 1]
            stats_df_candidates['germline_prob'] = probs
            
            # 自动计算 cutoff
            probs_sorted = np.sort(probs)
            n_remove = int(len(probs_sorted) * 0.05)
            probs_trimmed = probs_sorted[n_remove:-n_remove] if n_remove > 0 else probs_sorted
            # prob_threshold = calculate_prob_threshold(probs_trimmed)
            prob_threshold = max(probs_trimmed.mean() + 1.5 * probs_trimmed.std(), 0.1)
            
            stats_df_candidates['germline_pred'] = (probs > prob_threshold).astype(int)
            pred_germline_mutations = set(stats_df_candidates[stats_df_candidates['germline_pred'] == 1].index.tolist())
            
            print("Identified prob_threshold: ", str(prob_threshold))
            
            # 可视化 PDF
            if outputpath is not None:
                os.makedirs(outputpath, exist_ok=True)
                pdf_file = os.path.join(outputpath, "logreg_scatter_plots.pdf")
                from matplotlib.backends.backend_pdf import PdfPages
                import matplotlib.pyplot as plt
                
                with PdfPages(pdf_file) as pdf:  # pdf 只在这个块内部可用
                    pairs = [('FP_mean', 'FP_std'), ('FP_mean', 'FP_cv'), ('FP_std', 'FP_cv')]
                    for xcol, ycol in pairs:
                        plt.figure(figsize=(7, 6))
                        scatter = plt.scatter(stats_df_candidates[xcol], stats_df_candidates[ycol],
                                              c=stats_df_candidates['germline_prob'],
                                              cmap='coolwarm', s=40, edgecolors='k')
                        plt.colorbar(scatter, label="Predicted germline probability")
                        plt.xlabel(xcol)
                        plt.ylabel(ycol)
                        plt.title(f"Logistic regression: {xcol} vs {ycol}\nCutoff={prob_threshold:.3f}")
                        plt.tight_layout()
                        pdf.savefig()  # 一定要在 with 块内部
                        plt.close()
    
    # 新增列 'germline_determined'，根据行名是否在 pred_germline_mutations 中决定值
    stats_df_candidates['germline_determined'] = stats_df_candidates.index.to_series().apply(
        lambda x: 'germline' if x in pred_germline_mutations and candidates else 'non-germline'
    )
    
    # Step 7. 合并两个数据框    
    merged_df = pd.concat([stats_df_candidates, stats_df_non_germline], axis=0)
    # Step 7. 更新 'germline_determined' 列
    print("Updating germline status...")
    final_df = update_germline_status(merged_df, I, mcn_cutoff)
    print("Germline status updated.")
    
    # Step 8. 保存输出文件
    if outputpath is not None:
        os.makedirs(outputpath, exist_ok=True)
        final_df.to_csv(os.path.join(outputpath, "S_r_FP_stats_df.csv"))
    
    final_germline_mutations = set(final_df[final_df['germline_determined'] == 'germline'].index)
    
    print("Identified germline mutations:\n", final_germline_mutations)
    return final_df, final_germline_mutations




# -------------------------
# 3.2 Coverage-based filtration
# -------------------------

def filter_scaffold_muts_by_na_proportion_germline(filtered_sites, df_reads, df_celltype, na_prop_thresh=0.9):
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


def get_total_reads_withoutNAcells_germline(x):
    if pd.isna(x):
        return np.nan
    try:
        _, total = x.split('/')
        return int(total)
    except:
        return np.nan

def get_total_reads_withNAcells_germline(x):
    if pd.isna(x):
        return 0  # Treat NA as 0 reads
    try:
        _, total = x.split('/')
        return int(total)
    except:
        return 0  # In case of any parsing issues, treat as 0


def coverage_filters_germline(kept_mutations, df_reads, df_celltype, params, outputpath):
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
    step1_mutations, df_NA_prop = filter_scaffold_muts_by_na_proportion_germline(
        kept_mutations, df_reads, df_celltype, na_prop_thresh
    )
    logger.info("Section 3.2.1) Selection of ubiquitously expressed regions across cell (types)")
    print("=====> Step1 (coverage-based) mutations:", len(step1_mutations))
    
    # --- Step2: CV filter ---
    # Convert the reads to total read counts, treating NA as 0
    reads_matrix_withoutNAcells = df_reads.drop(index='bulk', errors='ignore').applymap(get_total_reads_withoutNAcells_germline)
    reads_matrix_withNAcells = df_reads.drop(index='bulk', errors='ignore').applymap(get_total_reads_withNAcells_germline)
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
    # final_scaffold_mutations = list(set(step1_mutations) | set(step2_mutations))
    final_scaffold_mutations = list(set(step2_mutations))
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



# -------------------------
# 按照 两两算 jaccard index 然后做 leiden graph
# -------------------------

def compute_clone_and_pair_weights_germline(muts, corr_cache, n_shuffle=100):
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


import matplotlib.pyplot as plt
import networkx as nx
def plot_mutation_graph_germline(G_ig, mutation_group, pdf_file, figsize=(8,8), edge_scale=0.2, seed=42):
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


from typing import List, Tuple, Dict
from collections import defaultdict
import itertools
import numpy as np
import pandas as pd

def cal_jaccard_index_by_pairs_for_graph_elements(I_S: pd.DataFrame):
    """
    Compute clone_weights and pair_jacidx, filtering out low-support singleton mutations.
    
    Args:
        I_S: binary matrix (cells x mutations), values 0/1/NA
        n_shuffle: number of random permutations per reference mutation
        seed: random seed
        min_frac: minimum mutant cell fraction for singleton mutations
        min_cells: minimum mutant cell number for singleton mutations
    
    Returns:
        clone_weights: dict mapping clone tuples to weight
        pair_jacidx: dict mapping mutation pairs to weight
    """
    muts = list(I_S.columns)
    n_mut = len(muts)
    
    # Step1: precompute pairwise jaccard index cache
    jacidx_cache = {}
    for u, v in itertools.combinations(muts, 2):
        jacidx = jaccard_index(I_S, u, v)
        jacidx_cache[(u, v)] = jacidx
        # jacidx_cache[(v, u)] = jacidx
    
    for m in muts:
        jacidx_cache[(m, m)] = 1.0
    
    # Step2: 将结果变成 leiden graph 输入格式
    pair_jacidx = defaultdict(float)
    for (var1, var2), weight in jacidx_cache.items():
        if var1 != var2:  # 跳过对角线
            sorted_pair = tuple(sorted((var1, var2)))
            pair_jacidx[sorted_pair] = weight
                
    return pair_jacidx


import igraph as ig
import leidenalg
def leiden_mutation_groups_using_jaccard_index(pair_jacidx, pdf_file, resolution=1.0, seed=42):
    """
    根据 clone_weights 和 pair_jacidx 构建加权共现图，并使用 Leiden 算法划分 mutation group。
    
    Parameters
    ----------
    clone_weights : dict
        {tuple(mutations): weight}  每个 clone 的全局权重，包括单节点 clone
    pair_jacidx : dict
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
    for clone in pair_jacidx.keys():
        all_mutations.update(clone)
    
    # 2. 构建 igraph 图
    G_ig = ig.Graph()
    G_ig.add_vertices(list(all_mutations))  # 所有 mutation 作为节点
    
    # 添加边（只考虑长度>=2的 clone）
    for (m1, m2), w in pair_jacidx.items():
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
    plot_mutation_graph_germline(G_ig, mutation_group, pdf_file)
    
    return mutation_group, partition, G_ig


def get_correlation_graph_elements_germline(I_S: pd.DataFrame, n_shuffle: int = 100, seed: int = 42, cutoff_mcf_for_graph: float = 0.05, cutoff_mcn_for_graph: int = 5) -> Tuple[Dict[Tuple[str], float], Dict[Tuple[str,str], float]]:
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
    clone_weights, pair_weights = compute_clone_and_pair_weights_germline(muts, corr_cache, n_shuffle=n_shuffle)
    
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


import igraph as ig
import leidenalg
def leiden_mutation_groups_germline(clone_weights, pair_weights, pdf_file, resolution=1.0, seed=42):
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
    plot_mutation_graph_germline(G_ig, mutation_group, pdf_file)
    
    return mutation_group, partition, G_ig




##### 在每一个 graph 划分的 groups 中找到 hub group
def detect_hub_clusters_germline(G_ig, mutation_group):
    """
    基于加权度中心性检测hub群体
    """
    # 1. 构建群体级别的图（就是你例子中的方法）
    cluster_graph = build_cluster_graph_germline(G_ig, mutation_group)
    
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


def build_cluster_graph_germline(G_ig, mutation_group):
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




# def identify_germline_variants(
#     P: pd.DataFrame, V: pd.DataFrame, C: pd.DataFrame, df_reads: pd.DataFrame, df_features_new: pd.DataFrame, 
#     p_thresh: float = 0.5, mcf_cutoff: float = 0.05, mcn_cutoff: int = 5, 
#     outputpath: Optional[str] = None,
#     sampleid: Optional[str] = None,
#     df_labeled: Optional[pd.DataFrame] = None
# ):



# pair_jacidx = cal_jaccard_index_by_pairs_for_graph_elements(I_add_germline)
# mutation_group, partition, G_ig = leiden_mutation_groups_using_jaccard_index(pair_jacidx, outputpath + "/" + sampleid + ".graph_for_mut_grouping.pdf")

# group_mutations = list(mutation_group.keys())

# I_selected_and_sorted, mut_df_sorted, group_to_muts, final_order = sort_I_hierarchical_freeze_ones_fixed(I_add_germline, mutation_group)
# df_celltype_sub = df_celltype[df_celltype['barcode'].isin(I_selected_and_sorted.index)].copy()
# plot_heatmap_with_celltype_by_your_sorting(I_selected_and_sorted, df_celltype_sub, mutation_group, list(mut_df_sorted['mutation']), os.path.join(outputpath, sampleid+".heatmap_with_celltype_right_in_I_selected_and_sorted_after_graph_grouping.pdf"))
# logger.info(f"Total mutation groups: {len(mutation_group)}")

# hub_clusters, cluster_degrees = detect_hub_clusters(G_ig, mutation_group)




# logger.info("Performing mutation grouping using Leiden algorithm ...")
# clone_weights, pair_weights = get_correlation_graph_elements(I_add_germline, 100, 42)
# mutation_group, partition, G_ig = leiden_mutation_groups(clone_weights, pair_weights, outputpath + "/" + sampleid + ".graph_for_mut_grouping.pdf")
# group_mutations = list(mutation_group.keys())
# I_selected_and_sorted, mut_df_sorted, group_to_muts, final_order = sort_I_hierarchical_freeze_ones_fixed(I_add_germline, mutation_group)
# df_celltype_sub = df_celltype[df_celltype['barcode'].isin(I_selected_and_sorted.index)].copy()
# plot_heatmap_with_celltype_by_your_sorting(I_selected_and_sorted, df_celltype_sub, mutation_group, list(mut_df_sorted['mutation']), os.path.join(outputpath, sampleid+".heatmap_with_celltype_right_in_I_selected_and_sorted_after_graph_grouping.pdf"))
# logger.info(f"Total mutation groups: {len(mutation_group)}")

# hub_clusters, cluster_degrees = detect_hub_clusters(G_ig, mutation_group)




# -------------------------
# Demo
# -------------------------
if __name__ == "__main__":
    import pandas as pd
    
    # 构造一个简单的 toy 矩阵（2个细胞 × 3个位点）
    P = pd.DataFrame([[0.9, 0.2, 0.8], [0.1, 0.8, 0.3]], 
                     index=["cell1", "cell2"], 
                     columns=["mut1", "mut2", "mut3"])
    M = pd.DataFrame([[1, 0, 1], [0, 1, 0]], index=P.index, columns=P.columns)
    C = pd.DataFrame([[10, 10, 10], [10, 10, 10]], index=P.index, columns=P.columns)
    A = (M * C).astype(int)
    
    # 构造一个简单的 labeled dataset，用于训练 logistic regression
    df_labeled = pd.DataFrame({
        'FP_mean': [0.8, 0.1, 0.5],
        'FP_std': [0.05, 0.02, 0.1],
        'FP_cv': [0.0625, 0.2, 0.2],
        'label': ['germline', 'mosaic', 'germline']
    }, index=['mut1','mut2','mut3'])
    
    # 调用函数
    stats_df, final_germline_mutations = identify_germline_variants(P, M, C, 
                                                       p_thresh=0.5, 
                                                       mcf_cutoff=0.05, 
                                                       outputpath=None, 
                                                       plot=False,
                                                       df_labeled=df_labeled)
    
    print("Stats dataframe:")
    print(stats_df)
    print("Putative germline variants:")
    print(final_germline_mutations)




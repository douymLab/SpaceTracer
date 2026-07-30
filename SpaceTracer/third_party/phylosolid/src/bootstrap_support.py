#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Bootstrap Branch Support Calculator for PhyloSOLID

用法:
    python bootstrap_support.py --input data.csv --branches branches.json --n_bootstrap 500 --output results.json

或者直接在Python中:
    from bootstrap_support import run_bootstrap
    results = run_bootstrap(df, branch_defs, n_bootstrap=500)
"""

import numpy as np
import pandas as pd
import random
import json
import pickle
import argparse
import time
from collections import defaultdict
from typing import Dict, List, Callable, Optional, Tuple
import warnings
warnings.filterwarnings('ignore')


# ============================================================
# 核心类
# ============================================================

class BootstrapSupportCalculator:
    """
    Bootstrap分支支持度计算器
    
    参数:
        mutation_matrix_df: pandas DataFrame, 行=细胞, 列=突变, 值=0/1/NaN
        branch_definitions: dict, 分支名称 -> 突变列表
        n_bootstrap: int, Bootstrap重复次数 (默认: 1000)
        random_seed: int, 随机种子 (默认: 42)
        verbose: bool, 是否显示进度 (默认: True)
    """
    
    def __init__(
        self,
        mutation_matrix_df: pd.DataFrame,
        branch_definitions: Dict[str, List[str]],
        n_bootstrap: int = 1000,
        random_seed: int = 42,
        verbose: bool = True
    ):
        # 输入验证
        if not isinstance(mutation_matrix_df, pd.DataFrame):
            raise ValueError("mutation_matrix_df 必须是 pandas DataFrame")
        if not branch_definitions:
            raise ValueError("branch_definitions 不能为空")
        if n_bootstrap < 10:
            raise ValueError("n_bootstrap 至少为 10")
        
        # 存储数据
        self.original_df = mutation_matrix_df.copy()
        self.branch_definitions = branch_definitions
        self.n_bootstrap = n_bootstrap
        self.random_seed = random_seed
        self.verbose = verbose
        
        # 提取元数据
        self.mutation_names = list(mutation_matrix_df.columns)
        self.cell_names = list(mutation_matrix_df.index)
        self.n_mutations = len(self.mutation_names)
        self.n_cells = len(self.cell_names)
        
        # 存储结果
        self.branch_support = {}
        self.branch_counts = {}
        self.branch_presence_records = {branch: [] for branch in branch_definitions}
        self.total_valid_bootstraps = 0
        self.failed_bootstraps = 0
        
        # 设置随机种子
        random.seed(random_seed)
        np.random.seed(random_seed)
        
        if self.verbose:
            self._print_init_info()
    
    def _print_init_info(self):
        """打印初始化信息"""
        print("\n" + "="*70)
        print("Bootstrap 分支支持度计算器")
        print("="*70)
        print(f"  突变数量: {self.n_mutations}")
        print(f"  细胞数量: {self.n_cells}")
        print(f"  分支数量: {len(self.branch_definitions)}")
        print(f"  Bootstrap次数: {self.n_bootstrap}")
        print(f"  随机种子: {self.random_seed}")
        print("\n分支定义:")
        for branch_name, mutations in self.branch_definitions.items():
            print(f"  {branch_name}: {len(mutations)} 个突变")
        print("="*70)
    
    def bootstrap_resample(self) -> List[str]:
        """
        从所有突变中有放回地抽取 n_mutations 个突变
        
        返回:
            抽取的突变名称列表 (可能有重复)
        """
        return random.choices(self.mutation_names, k=self.n_mutations)
    
    def create_bootstrap_matrix(self, sampled_mutations: List[str]) -> pd.DataFrame:
        """
        根据抽到的突变列表创建bootstrap矩阵
        
        参数:
            sampled_mutations: 抽到的突变名称列表 (可能有重复)
        
        返回:
            bootstrap DataFrame
        """
        bootstrap_df = pd.DataFrame(index=self.cell_names)
        
        for mut_name in sampled_mutations:
            if mut_name in self.original_df.columns:
                bootstrap_df[mut_name] = self.original_df[mut_name].values
            else:
                bootstrap_df[mut_name] = np.nan
        
        return bootstrap_df
    
    def check_branch_in_tree(
        self,
        tree_result,
        branch_name: str,
        branch_mutations: List[str],
        tolerance: float = 0.3
    ) -> bool:
        """
        检查分支是否在bootstrap树中存在
        
        参数:
            tree_result: PhyloSOLID建树结果
            branch_name: 分支名称
            branch_mutations: 分支包含的突变列表
            tolerance: 允许的外来突变比例 (默认: 0.3)
        
        返回:
            True/False
        """
        if not tree_result:
            return False
        
        branch_set = set(branch_mutations)
        
        # ============================================================
        # 格式1: 树是列表的列表 [[mut1, mut2], [mut3, mut4, mut5]]
        # ============================================================
        if isinstance(tree_result, (list, tuple)):
            for cluster in tree_result:
                if isinstance(cluster, (list, tuple, set)):
                    cluster_set = set(cluster)
                    if branch_set.issubset(cluster_set):
                        extra = cluster_set - branch_set
                        if len(extra) / max(len(cluster_set), 1) < tolerance:
                            return True
        
        # ============================================================
        # 格式2: 树是Newick字符串
        # ============================================================
        elif isinstance(tree_result, str):
            # 检查所有分支突变是否在字符串中
            for mut in branch_mutations:
                if mut not in tree_result:
                    return False
            
            # 查找包含这些突变的括号对
            import re
            # 匹配括号内的内容
            matches = re.findall(r'\(([^()]+)\)', tree_result)
            for match in matches:
                muts_in_cluster = set(re.findall(r'[A-Za-z0-9_\-\.]+', match))
                if branch_set.issubset(muts_in_cluster):
                    extra = muts_in_cluster - branch_set
                    if len(extra) / max(len(muts_in_cluster), 1) < tolerance:
                        return True
        
        # ============================================================
        # 格式3: 树是字典
        # ============================================================
        elif isinstance(tree_result, dict):
            for key in ['clusters', 'branches', 'groups', 'clades']:
                if key in tree_result:
                    for cluster in tree_result[key]:
                        cluster_set = set(cluster)
                        if branch_set.issubset(cluster_set):
                            extra = cluster_set - branch_set
                            if len(extra) / max(len(cluster_set), 1) < tolerance:
                                return True
        
        return False
    
    def run_bootstrap(
        self,
        phylosolid_func: Callable,
        branch_tolerance: float = 0.3
    ) -> Dict[str, float]:
        """
        运行Bootstrap分析
        
        参数:
            phylosolid_func: PhyloSOLID建树函数，接受DataFrame，返回树结构
            branch_tolerance: 允许的外来突变比例
        
        返回:
            Dict[分支名称, 支持度百分比]
        """
        # 初始化计数
        support_counts = {branch: 0 for branch in self.branch_definitions}
        self.branch_presence_records = {branch: [] for branch in self.branch_definitions}
        
        start_time = time.time()
        
        if self.verbose:
            print(f"\n开始 Bootstrap 分析 (共 {self.n_bootstrap} 次)...")
            print("-"*70)
        
        for i in range(self.n_bootstrap):
            # 1. 有放回抽取突变
            sampled_mutations = self.bootstrap_resample()
            
            # 2. 创建bootstrap矩阵
            bootstrap_df = self.create_bootstrap_matrix(sampled_mutations)
            
            # 3. 运行 PhyloSOLID 建树
            try:
                tree_result = phylosolid_func(bootstrap_df)
                self.total_valid_bootstraps += 1
            except Exception as e:
                self.failed_bootstraps += 1
                if self.verbose and self.failed_bootstraps <= 5:
                    print(f"  ⚠️ 第 {i+1} 次建树失败: {e}")
                continue
            
            # 4. 检查每个分支
            for branch_name, branch_muts in self.branch_definitions.items():
                # 检查该分支的突变是否被抽到了 (至少50%被抽到)
                sampled_branch_muts = [m for m in branch_muts if m in sampled_mutations]
                
                if len(sampled_branch_muts) < len(branch_muts) * 0.5:
                    self.branch_presence_records[branch_name].append(False)
                    continue
                
                # 检查分支是否在树中存在
                exists = self.check_branch_in_tree(
                    tree_result, branch_name, branch_muts, branch_tolerance
                )
                self.branch_presence_records[branch_name].append(exists)
                
                if exists:
                    support_counts[branch_name] += 1
            
            # 5. 进度显示
            if self.verbose and (i + 1) % max(1, self.n_bootstrap // 10) == 0:
                elapsed = time.time() - start_time
                remaining = (elapsed / (i + 1)) * (self.n_bootstrap - i - 1)
                print(f"  完成 {i+1}/{self.n_bootstrap} 次 ({100*(i+1)/self.n_bootstrap:.0f}%)")
                print(f"    耗时: {elapsed:.1f}s, 预计剩余: {remaining:.1f}s")
                # 显示当前支持度
                for branch_name in self.branch_definitions:
                    valid_count = len(self.branch_presence_records[branch_name])
                    if valid_count > 0:
                        current = (support_counts[branch_name] / valid_count) * 100
                        print(f"    {branch_name}: {current:.1f}%")
                print("-"*70)
        
        # 6. 计算最终支持度
        for branch_name in self.branch_definitions:
            valid_count = len(self.branch_presence_records[branch_name])
            if valid_count > 0:
                support = (support_counts[branch_name] / valid_count) * 100
            else:
                support = 0.0
            self.branch_support[branch_name] = support
            self.branch_counts[branch_name] = support_counts[branch_name]
        
        # 7. 打印最终统计
        if self.verbose:
            elapsed = time.time() - start_time
            print(f"\nBootstrap 完成!")
            print(f"  总耗时: {elapsed:.1f}s ({elapsed/60:.1f}min)")
            print(f"  成功建树: {self.total_valid_bootstraps}/{self.n_bootstrap}")
            if self.failed_bootstraps > 0:
                print(f"  建树失败: {self.failed_bootstraps} 次")
        
        return self.branch_support
    
    def print_results(self):
        """打印Bootstrap支持度结果"""
        if not self.branch_support:
            print("⚠️ 请先运行 run_bootstrap()")
            return
        
        print("\n" + "="*70)
        print("Bootstrap 分支支持度结果")
        print("="*70)
        
        # 按支持度排序
        sorted_branches = sorted(
            self.branch_support.items(),
            key=lambda x: x[1],
            reverse=True
        )
        
        for branch_name, support in sorted_branches:
            # 星号可视化 (每5%一个星号)
            stars = '*' * int(support / 5)
            
            # 支持度评级
            if support >= 95:
                level = "⭐⭐⭐ 非常稳健"
            elif support >= 80:
                level = "⭐⭐ 稳健"
            elif support >= 70:
                level = "⭐ 可接受"
            elif support >= 50:
                level = "⚠️ 中等支持度"
            else:
                level = "❌ 低支持度"
            
            print(f"\n{branch_name}:")
            print(f"  支持度: {support:.1f}%  {stars}")
            print(f"  评级: {level}")
            print(f"  包含突变数: {len(self.branch_definitions[branch_name])}")
            if len(self.branch_definitions[branch_name]) <= 10:
                print(f"  突变: {self.branch_definitions[branch_name]}")
            else:
                preview = self.branch_definitions[branch_name][:5]
                print(f"  突变 (前5个): {preview} ... (共{len(self.branch_definitions[branch_name])}个)")
        
        # 统计摘要
        print("\n" + "-"*70)
        print("摘要统计:")
        total = len(self.branch_support)
        high = sum(1 for v in self.branch_support.values() if v >= 70)
        medium = sum(1 for v in self.branch_support.values() if 50 <= v < 70)
        low = sum(1 for v in self.branch_support.values() if v < 50)
        
        print(f"  总分支数: {total}")
        print(f"  ✅ 高支持度 (>=70%): {high} 个 ({100*high/total:.1f}%)")
        print(f"  📊 中等支持度 (50-69%): {medium} 个 ({100*medium/total:.1f}%)")
        print(f"  ❌ 低支持度 (<50%): {low} 个 ({100*low/total:.1f}%)")
        print("="*70)
    
    def save_results(self, filename: str = "bootstrap_results.pkl"):
        """保存结果到文件"""
        results = {
            'branch_definitions': self.branch_definitions,
            'branch_support': self.branch_support,
            'branch_counts': self.branch_counts,
            'branch_presence_records': self.branch_presence_records,
            'n_bootstrap': self.n_bootstrap,
            'total_valid': self.total_valid_bootstraps,
            'failed': self.failed_bootstraps,
            'random_seed': self.random_seed,
            'n_mutations': self.n_mutations,
            'n_cells': self.n_cells,
            'mutation_names': self.mutation_names,
            'cell_names': self.cell_names
        }
        
        with open(filename, 'wb') as f:
            pickle.dump(results, f)
        
        print(f"\n结果已保存至: {filename}")
        
        # 同时导出CSV
        csv_file = filename.replace('.pkl', '.csv')
        df = pd.DataFrame({
            'Branch': list(self.branch_support.keys()),
            'Support_%': list(self.branch_support.values()),
            'Count': [self.branch_counts[b] for b in self.branch_support],
            'Num_Mutations': [len(self.branch_definitions[b]) for b in self.branch_support]
        })
        df.to_csv(csv_file, index=False)
        print(f"支持度数据已导出至: {csv_file}")
        
        return results
    
    def load_results(self, filename: str):
        """从文件加载结果"""
        with open(filename, 'rb') as f:
            results = pickle.load(f)
        
        self.branch_support = results['branch_support']
        self.branch_counts = results['branch_counts']
        self.branch_presence_records = results['branch_presence_records']
        self.total_valid_bootstraps = results['total_valid']
        self.failed_bootstraps = results['failed']
        
        print(f"结果已加载: {filename}")
        return results


# ============================================================
# 辅助函数
# ============================================================

def load_data(csv_file: str, index_col: int = 0) -> pd.DataFrame:
    """加载数据"""
    df = pd.read_csv(csv_file, index_col=index_col)
    print(f"加载数据: {df.shape[0]} 行 × {df.shape[1]} 列")
    return df


def load_branch_definitions(json_file: str) -> Dict[str, List[str]]:
    """从JSON文件加载分支定义"""
    with open(json_file, 'r') as f:
        return json.load(f)


def save_branch_definitions(branch_defs: Dict[str, List[str]], filename: str = "branches.json"):
    """保存分支定义到JSON"""
    with open(filename, 'w') as f:
        json.dump(branch_defs, f, indent=2)
    print(f"分支定义已保存至: {filename}")


# ============================================================
# 主运行函数
# ============================================================

def run_bootstrap(
    mutation_matrix_df: pd.DataFrame,
    branch_definitions: Dict[str, List[str]],
    phylosolid_func: Callable,
    n_bootstrap: int = 1000,
    random_seed: int = 42,
    branch_tolerance: float = 0.3,
    save_path: str = "bootstrap_results.pkl",
    verbose: bool = True
) -> Dict[str, float]:
    """
    运行Bootstrap分析的便捷函数
    
    参数:
        mutation_matrix_df: pandas DataFrame
        branch_definitions: 分支定义字典
        phylosolid_func: PhyloSOLID建树函数
        n_bootstrap: Bootstrap次数
        random_seed: 随机种子
        branch_tolerance: 分支检查容差
        save_path: 结果保存路径
        verbose: 是否显示进度
    
    返回:
        Dict[分支名称, 支持度百分比]
    """
    # 创建计算器
    calculator = BootstrapSupportCalculator(
        mutation_matrix_df=mutation_matrix_df,
        branch_definitions=branch_definitions,
        n_bootstrap=n_bootstrap,
        random_seed=random_seed,
        verbose=verbose
    )
    
    # 运行Bootstrap
    support_values = calculator.run_bootstrap(
        phylosolid_func=phylosolid_func,
        branch_tolerance=branch_tolerance
    )
    
    # 打印结果
    calculator.print_results()
    
    # 保存结果
    calculator.save_results(save_path)
    
    return support_values


# ============================================================
# 命令行接口
# ============================================================

def main():
    """命令行入口"""
    parser = argparse.ArgumentParser(
        description="Bootstrap Branch Support Calculator for PhyloSOLID"
    )
    parser.add_argument(
        '--input', '-i',
        type=str,
        required=True,
        help='输入数据文件 (CSV格式，行=细胞，列=突变)'
    )
    parser.add_argument(
        '--branches', '-b',
        type=str,
        required=True,
        help='分支定义文件 (JSON格式)'
    )
    parser.add_argument(
        '--output', '-o',
        type=str,
        default='bootstrap_results.pkl',
        help='输出文件路径 (默认: bootstrap_results.pkl)'
    )
    parser.add_argument(
        '--n_bootstrap', '-n',
        type=int,
        default=1000,
        help='Bootstrap重复次数 (默认: 1000)'
    )
    parser.add_argument(
        '--seed',
        type=int,
        default=42,
        help='随机种子 (默认: 42)'
    )
    parser.add_argument(
        '--tolerance',
        type=float,
        default=0.3,
        help='分支检查容差，允许的外来突变比例 (默认: 0.3)'
    )
    parser.add_argument(
        '--quiet',
        action='store_true',
        help='安静模式，不显示进度'
    )
    
    args = parser.parse_args()
    
    # 加载数据
    print(f"加载数据: {args.input}")
    df = load_data(args.input)
    
    # 加载分支定义
    print(f"加载分支定义: {args.branches}")
    branch_defs = load_branch_definitions(args.branches)
    
    # 这里需要用户定义自己的PhyloSOLID函数
    # 示例: 使用简单的聚类作为占位
    def dummy_phylosolid(matrix_df):
        """示例建树函数 - 请替换为你的PhyloSOLID"""
        from scipy.cluster.hierarchy import linkage, fcluster
        from scipy.spatial.distance import pdist
        
        data = matrix_df.fillna(0).values
        if data.shape[1] < 2:
            return []
        
        try:
            dist = pdist(data.T)
            linkage_matrix = linkage(dist, method='average')
            clusters = fcluster(linkage_matrix, t=0.6, criterion='distance')
            
            cluster_dict = defaultdict(list)
            mut_names = matrix_df.columns.tolist()
            for i, label in enumerate(clusters):
                cluster_dict[label].append(mut_names[i])
            
            return list(cluster_dict.values())
        except:
            return []
    
    print("\n⚠️ 注意: 当前使用的是 dummy_phylosolid 函数!")
    print("   请替换为你的 PhyloSOLID 建树函数")
    
    # 运行Bootstrap
    results = run_bootstrap(
        mutation_matrix_df=df,
        branch_definitions=branch_defs,
        phylosolid_func=dummy_phylosolid,
        n_bootstrap=args.n_bootstrap,
        random_seed=args.seed,
        branch_tolerance=args.tolerance,
        save_path=args.output,
        verbose=not args.quiet
    )
    
    # 打印最终结果
    print("\n" + "="*70)
    print("最终结果:")
    for branch, support in results.items():
        print(f"  {branch}: {support:.1f}%")
    print("="*70)


# ============================================================
# 使用示例
# ============================================================

def example_usage():
    """使用示例"""
    print("\n" + "="*70)
    print("Bootstrap Support Calculator - 使用示例")
    print("="*70)
    
    # 1. 准备数据
    np.random.seed(42)
    n_cells = 50
    n_mutations = 30
    
    data = np.random.binomial(1, 0.3, (n_cells, n_mutations))
    df = pd.DataFrame(
        data,
        index=[f'Cell_{i}' for i in range(n_cells)],
        columns=[f'M{i+1}' for i in range(n_mutations)]
    )
    
    # 添加NaN
    mask = np.random.random((n_cells, n_mutations)) < 0.1
    df[mask] = np.nan
    
    # 2. 定义分支
    branch_defs = {
        'Branch_A': [f'M{i+1}' for i in range(10)],
        'Branch_B': [f'M{i+1}' for i in range(10, 20)],
        'Branch_C': [f'M{i+1}' for i in range(20, 30)]
    }
    
    # 3. 定义建树函数 (示例)
    def my_phylosolid(matrix_df):
        from scipy.cluster.hierarchy import linkage, fcluster
        from scipy.spatial.distance import pdist
        
        data = matrix_df.fillna(0).values
        if data.shape[1] < 2:
            return []
        
        try:
            dist = pdist(data.T)
            linkage_matrix = linkage(dist, method='average')
            clusters = fcluster(linkage_matrix, t=0.6, criterion='distance')
            
            cluster_dict = defaultdict(list)
            mut_names = matrix_df.columns.tolist()
            for i, label in enumerate(clusters):
                cluster_dict[label].append(mut_names[i])
            
            return list(cluster_dict.values())
        except:
            return []
    
    # 4. 运行Bootstrap (使用100次作为演示)
    # 实际使用时将 n_bootstrap 改为 500 或 1000
    results = run_bootstrap(
        mutation_matrix_df=df,
        branch_definitions=branch_defs,
        phylosolid_func=my_phylosolid,
        n_bootstrap=100,  # 示例用100次
        random_seed=42,
        save_path="example_results.pkl"
    )
    
    return results


if __name__ == "__main__":
    # 如果直接运行，执行示例
    example_usage()

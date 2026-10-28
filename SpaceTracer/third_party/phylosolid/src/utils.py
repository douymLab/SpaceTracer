import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import pdist, squareform
from collections import defaultdict
import random

class BootstrapSupportCalculator:
    """
    计算系统发育树分支的 Bootstrap 支持度
    """
    
    def __init__(self, mutation_matrix, mutation_names=None, n_bootstrap=1000, random_seed=42):
        """
        参数:
        - mutation_matrix: 突变矩阵，行=突变，列=样本/细胞，值=0/1 (1表示该突变存在)
        - mutation_names: 突变名称列表，如 ['M1', 'M2', ...]
        - n_bootstrap: Bootstrap 重复次数
        - random_seed: 随机种子，确保结果可重复
        """
        self.original_matrix = np.array(mutation_matrix)
        self.n_mutations = self.original_matrix.shape[0]
        self.n_samples = self.original_matrix.shape[1]
        self.mutation_names = mutation_names if mutation_names else [f'M{i+1}' for i in range(self.n_mutations)]
        self.n_bootstrap = n_bootstrap
        self.random_seed = random_seed
        
        # 构建原始树，识别第一层级分支
        self.original_tree = None
        self.main_branches = None  # 存储第一层级分支的突变集合
        
    def build_tree_and_get_branches(self, distance_matrix=None, threshold=0.5):
        """
        基于距离矩阵构建树，并识别第一层级分支
        
        参数:
        - distance_matrix: 突变间的距离矩阵 (如果为None，则基于共现模式计算)
        - threshold: 聚类阈值，用于确定第一层级分支
        """
        if distance_matrix is None:
            # 基于Jaccard距离计算突变间距离
            distance_matrix = self._compute_jaccard_distance()
        
        # 使用UPGMA进行层次聚类
        condensed_dist = squareform(distance_matrix)
        linkage_matrix = linkage(condensed_dist, method='average')
        
        # 根据阈值切割树，得到第一层级分支
        # 这里我们取距离矩阵中位数的某个倍数作为阈值
        if threshold is None:
            threshold = np.median(distance_matrix[distance_matrix > 0]) * 0.8
        
        cluster_labels = fcluster(linkage_matrix, t=threshold, criterion='distance')
        
        # 提取每个分支包含的突变
        branches = defaultdict(list)
        for i, label in enumerate(cluster_labels):
            branches[label].append(i)
        
        # 只保留包含至少2个突变的分支
        self.main_branches = {f'Branch_{k}': sorted(indices) 
                             for k, indices in branches.items() if len(indices) >= 2}
        
        return self.main_branches
    
    def _compute_jaccard_distance(self):
        """计算突变间的Jaccard距离"""
        n = self.n_mutations
        distance_matrix = np.zeros((n, n))
        
        for i in range(n):
            for j in range(i+1, n):
                # 计算两个突变的共现情况
                intersection = np.sum(self.original_matrix[i] & self.original_matrix[j])
                union = np.sum(self.original_matrix[i] | self.original_matrix[j])
                
                if union == 0:
                    distance = 1.0  # 两个突变都不存在，距离最大
                else:
                    jaccard_sim = intersection / union
                    distance = 1 - jaccard_sim  # Jaccard距离
                
                distance_matrix[i, j] = distance
                distance_matrix[j, i] = distance
        
        return distance_matrix
    
    def bootstrap_resample(self):
        """
        有放回地抽取突变，生成bootstrap样本
        """
        # 有放回抽取，数量等于原始突变数
        sampled_indices = random.choices(range(self.n_mutations), k=self.n_mutations)
        bootstrap_matrix = self.original_matrix[sampled_indices]
        return bootstrap_matrix, sampled_indices
    
    def compute_branch_support(self):
        """
        计算每个第一层级分支的Bootstrap支持度
        """
        # 首先构建原始树，识别第一层级分支
        if self.main_branches is None:
            self.build_tree_and_get_branches()
        
        print(f"识别到 {len(self.main_branches)} 个第一层级分支:")
        for branch_name, mutations in self.main_branches.items():
            mutation_labels = [self.mutation_names[i] for i in mutations]
            print(f"  {branch_name}: {mutation_labels}")
        
        # 存储每个分支的支持度计数
        support_counts = {branch_name: 0 for branch_name in self.main_branches}
        
        # 设置随机种子
        random.seed(self.random_seed)
        np.random.seed(self.random_seed)
        
        print(f"\n开始 Bootstrap 分析 (重复 {self.n_bootstrap} 次)...")
        
        for bootstrap_iter in range(self.n_bootstrap):
            if (bootstrap_iter + 1) % 100 == 0:
                print(f"  完成 {bootstrap_iter + 1}/{self.n_bootstrap} 次")
            
            # 生成bootstrap样本
            bootstrap_matrix, sampled_indices = self.bootstrap_resample()
            
            # 检查每个主分支是否在bootstrap树中出现
            for branch_name, original_mutations in self.main_branches.items():
                if self._check_branch_exists(bootstrap_matrix, original_mutations, sampled_indices):
                    support_counts[branch_name] += 1
        
        # 计算支持度百分比
        support_values = {}
        for branch_name, count in support_counts.items():
            support_percent = (count / self.n_bootstrap) * 100
            support_values[branch_name] = support_percent
        
        return support_values
    
    def _check_branch_exists(self, bootstrap_matrix, original_mutations, sampled_indices):
        """
        检查原始分支的突变在bootstrap树中是否仍然聚在一起
        
        参数:
        - bootstrap_matrix: bootstrap重采样后的突变矩阵
        - original_mutations: 原始分支中的突变索引列表
        - sampled_indices: 重采样时抽取的原始索引
        """
        # 获取bootstrap矩阵中对应原始突变的行
        # 注意：bootstrap矩阵中的行对应于sampled_indices中的索引
        
        # 首先找出哪些原始突变被抽到了bootstrap样本中
        # 以及它们在bootstrap矩阵中的位置
        mutation_positions = []
        for orig_idx in original_mutations:
            # 在sampled_indices中找到所有等于orig_idx的位置
            positions = [i for i, idx in enumerate(sampled_indices) if idx == orig_idx]
            if positions:
                mutation_positions.extend(positions)
        
        # 如果分支中的突变在bootstrap样本中丢失太多，认为分支不存在
        if len(mutation_positions) < len(original_mutations) * 0.5:
            return False
        
        # 基于bootstrap矩阵计算这些突变之间的距离
        if len(mutation_positions) < 2:
            return False
        
        # 提取这些突变的数据
        sub_matrix = bootstrap_matrix[mutation_positions]
        
        # 计算它们之间的平均距离（Jaccard距离）
        distances = []
        for i in range(len(sub_matrix)):
            for j in range(i+1, len(sub_matrix)):
                intersection = np.sum(sub_matrix[i] & sub_matrix[j])
                union = np.sum(sub_matrix[i] | sub_matrix[j])
                if union == 0:
                    dist = 1.0
                else:
                    dist = 1 - (intersection / union)
                distances.append(dist)
        
        avg_distance = np.mean(distances) if distances else 1.0
        
        # 计算这些突变与其他突变之间的平均距离
        other_distances = []
        all_positions = list(range(len(bootstrap_matrix)))
        other_positions = [p for p in all_positions if p not in mutation_positions]
        
        if other_positions:
            for pos in mutation_positions:
                for other in other_positions:
                    intersection = np.sum(bootstrap_matrix[pos] & bootstrap_matrix[other])
                    union = np.sum(bootstrap_matrix[pos] | bootstrap_matrix[other])
                    if union == 0:
                        dist = 1.0
                    else:
                        dist = 1 - (intersection / union)
                    other_distances.append(dist)
            
            avg_other_distance = np.mean(other_distances) if other_distances else 1.0
            
            # 如果分支内平均距离小于分支间平均距离，认为分支存在
            return avg_distance < avg_other_distance
        else:
            # 如果没有其他突变可以比较，认为分支存在
            return True
    
    def print_results(self, support_values):
        """
        打印Bootstrap支持度结果
        """
        print("\n" + "="*60)
        print("Bootstrap 分支支持度结果")
        print("="*60)
        
        for branch_name, support in sorted(support_values.items(), key=lambda x: x[1], reverse=True):
            mutations = [self.mutation_names[i] for i in self.main_branches[branch_name]]
            stars = '*' * int(support / 5)  # 每5%显示一个星号
            print(f"{branch_name}: {support:.1f}%  {stars}")
            print(f"  包含突变: {mutations}")
            print()
        
        # 显示支持度摘要
        high_support = sum(1 for v in support_values.values() if v >= 70)
        medium_support = sum(1 for v in support_values.values() if 50 <= v < 70)
        low_support = sum(1 for v in support_values.values() if v < 50)
        
        print("-"*60)
        print(f"摘要:")
        print(f"  高支持度 (>=70%): {high_support} 个分支")
        print(f"  中等支持度 (50-69%): {medium_support} 个分支")
        print(f"  低支持度 (<50%): {low_support} 个分支")
        print("="*60)


# ============================================================
# 示例使用代码
# ============================================================

def example_usage():
    """
    示例：如何使用BootstrapSupportCalculator
    """
    
    # 创建示例数据
    # 假设我们有30个突变 (M1-M30)
    # 3个主分支: Branch1 (M1-M15), Branch2 (M16-M22), Branch3 (M23-M30)
    # 10个样本/细胞
    
    np.random.seed(42)
    
    n_mutations = 30
    n_samples = 10
    
    # 创建突变矩阵
    mutation_matrix = np.zeros((n_mutations, n_samples), dtype=int)
    
    # 为每个分支生成特征模式
    # Branch1: M1-M15 在样本0-4中存在
    for i in range(0, 15):
        mutation_matrix[i, :5] = np.random.binomial(1, 0.8, 5)
        mutation_matrix[i, 5:] = np.random.binomial(1, 0.1, 5)
    
    # Branch2: M16-M22 在样本3-7中存在
    for i in range(15, 22):
        mutation_matrix[i, 3:8] = np.random.binomial(1, 0.8, 5)
        mutation_matrix[i, :3] = np.random.binomial(1, 0.1, 3)
        mutation_matrix[i, 8:] = np.random.binomial(1, 0.1, 2)
    
    # Branch3: M23-M30 在样本6-9中存在
    for i in range(22, 30):
        mutation_matrix[i, 6:10] = np.random.binomial(1, 0.8, 4)
        mutation_matrix[i, :6] = np.random.binomial(1, 0.1, 6)
    
    # 添加一些随机噪声
    noise_mask = np.random.random((n_mutations, n_samples)) < 0.05
    mutation_matrix[noise_mask] = 1 - mutation_matrix[noise_mask]
    
    # 突变名称
    mutation_names = [f'M{i+1}' for i in range(n_mutations)]
    
    # 创建计算器实例
    calculator = BootstrapSupportCalculator(
        mutation_matrix=mutation_matrix,
        mutation_names=mutation_names,
        n_bootstrap=1000,  # 可以改为100进行快速测试
        random_seed=42
    )
    
    # 计算Bootstrap支持度
    support_values = calculator.compute_branch_support()
    
    # 打印结果
    calculator.print_results(support_values)
    
    return calculator, support_values


if __name__ == "__main__":
    # 运行示例
    calculator, support_values = example_usage()



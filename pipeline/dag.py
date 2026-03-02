#!/usr/bin/env python3
"""
PipelineDAG - 你的原始代码 + 以下修复：
  1. 修复 'pileup' → 'mpileup' 的key不一致问题
  2. 修复 start_from/stop_at 在有分支时不正确（不能用list切片）
  3. 修复 prior step 的 name 字段与 key 不一致
  4. 完善 only_steps 的依赖自动补全逻辑
"""

from typing import List, Dict, Set, Optional
from dataclasses import dataclass
from enum import Enum


class StepStatus(Enum):
    PENDING   = "pending"
    RUNNING   = "running"
    COMPLETED = "completed"
    FAILED    = "failed"
    SKIPPED   = "skipped"


@dataclass
class Step:
    """描述一个步骤的元数据（不含执行逻辑）"""
    name: str
    description: str
    dependencies: List[str]   # 上游依赖的step name列表
    produces: List[str]       # 产出的数据/文件标识
    optional: bool = False
    parallel: bool = False


class PipelineDAG:
    """
    基于DAG的流程调度器

    职责：
    - 定义步骤间的依赖关系（纯元数据）
    - 生成正确的拓扑排序执行计划
    - 计算可并行执行的步骤组
    - 查询上下游关系
    """

    STEPS = {
        # 'cluster': Step(
        #     name='cluster',
        #     description='聚类分析',
        #     dependencies=[],
        #     produces=['clustered_bam']
        # ),
        'bam_processing': Step(
            name='bam_processing',
            description='BAM过滤',
            dependencies=[], #'cluster'
            produces=['filtered_bam']
        ),
        'cell_num': Step(
            name='cell_num',
            description='cell number per spot calculation',
            dependencies=['bam_processing'],
            produces=['cell_num_results']
        ),
        'mpileup': Step(
            name='mpileup',
            description='mpileup生成',
            dependencies=['bam_processing'],
            produces=['pileup_results']
        ),
        # ────────────────────────────────────────────
        # BUG FIX: 原来写的 dependencies=['pileup']
        # 但字典的key是 'mpileup'，导致依赖永远找不到
        # ────────────────────────────────────────────
        'umi_combine': Step(
            name='umi_combine',
            description='UMI合并',
            dependencies=['mpileup'],          # ← 修正: 'pileup' → 'mpileup'
            produces=['umi_combined']
        ),
        'prior': Step(
            name='prior',                      # ← 修正: 原来name='prior_error'但key='prior'
            description='先验错误率计算',
            dependencies=['mpileup'],          # ← 修正: 'pileup' → 'mpileup'
            produces=['prior_error_results']
        ),
        # 'error': Step(
        #     name='error',
        #     description='错误率计算',
        #     dependencies=['umi_combine'],
        #     produces=['error_results']
        # ),
        'genotyping': Step(
            name='genotyping',
            description='基因分型',
            dependencies=['umi_combine', 'prior'],  # 汇合两个分支
            produces=['genotype_results']
        ),
        # 'phasing': Step(
        #     name='phasing',
        #     description='单倍型phasing',
        #     dependencies=['genotyping'],
        #     produces=['phasing_results']
        # ),
        'spatial_feature': Step(
            name='spatial_feature',
            description='提取空间特征',
            dependencies=['genotyping'],
            produces=['spatial_features']
        ),
        # 'feature': Step(
        #     name='feature',
        #     description='特征整合',
        #     dependencies=['spatial_feature', 'prior'],
        #     produces=['feature_results']
        # ),
        # 'hfdr': Step(
        #     name='hfdr',
        #     description='hFDR校正',
        #     dependencies=['error', 'genotyping', 'phasing', 'feature'],
        #     produces=['hfdr_results']
        # ),
        # 'final_filter': Step(
        #     name='final_filter',
        #     description='最终过滤',
        #     dependencies=['hfdr'],
        #     produces=['final_variants']
        # ),
    }

    def __init__(self):
        self.steps = self.STEPS
        self._validate_dag()

    # ─────────────────────────────────────────────
    # 原有方法：保持你的接口不变
    # ─────────────────────────────────────────────

    def _validate_dag(self):
        """验证DAG是否有环"""
        visited = set()
        recursion_stack = set()

        def dfs(step_name):
            if step_name in recursion_stack:
                raise ValueError(f"检测到循环依赖: {step_name}")
            if step_name in visited:
                return

            visited.add(step_name)
            recursion_stack.add(step_name)

            step = self.steps.get(step_name)
            if step:
                for dep in step.dependencies:
                    if dep in self.steps:
                        dfs(dep)

            recursion_stack.remove(step_name)

        for step_name in self.steps:
            dfs(step_name)

    def get_execution_plan(self,
                           start_from: Optional[str] = None,
                           stop_at: Optional[str] = None,
                           only_steps: Optional[List[str]] = None) -> List[str]:
        """
        生成执行计划（拓扑排序）

        BUG FIX: 原来用 all_steps[start_idx:end_idx] 线性切片，
        在有分支的DAG中会漏掉或错误包含步骤。
        现在改为基于依赖关系的子图提取。
        """
        if only_steps:
            # 补全 only_steps 的依赖，然后拓扑排序
            return self._get_subgraph_with_deps(only_steps)

        if start_from or stop_at:
            # 提取 [start_from, stop_at] 之间的子图
            return self._get_range_subgraph(start_from, stop_at)

        # 全量执行：对所有step做拓扑排序
        return self._topological_sort(list(self.steps.keys()))

    def _get_range_subgraph(self,
                            start_from: Optional[str],
                            stop_at: Optional[str]) -> List[str]:
        """
        提取从 start_from 到 stop_at 的子图

        策略：
        - 先确定 stop_at 及其所有上游 → 这是终止集合
        - 再从终止集合中去掉 start_from 的上游（这些已完成，跳过）
        - 对剩余部分做拓扑排序
        """
        all_steps = set(self.steps.keys())

        # 确定终止范围
        if stop_at:
            if stop_at not in self.steps:
                raise ValueError(f"Unknown step: '{stop_at}'")
            # stop_at 及其所有上游
            included = self.get_upstream(stop_at) | {stop_at}
        else:
            included = all_steps

        # 去掉 start_from 的上游（已完成，不需要跑）
        if start_from:
            if start_from not in self.steps:
                raise ValueError(f"Unknown step: '{start_from}'")
            excluded = self.get_upstream(start_from)  # start_from本身要跑，不排除
            included -= excluded

        return self._topological_sort(list(included))

    def _get_subgraph_with_deps(self, step_names: List[str]) -> List[str]:
        """
        获取指定步骤 + 它们的所有依赖，做拓扑排序
        （与原来的 _get_subgraph 相同，保留你的原始实现）
        """
        from collections import deque

        needed = set(step_names)
        queue = deque(step_names)

        while queue:
            step_name = queue.popleft()
            if step_name in self.steps:
                for dep in self.steps[step_name].dependencies:
                    if dep not in needed:
                        needed.add(dep)
                        queue.append(dep)

        return self._topological_sort(list(needed))

    def _topological_sort(self, step_names: List[str]) -> List[str]:
        """
        Kahn算法拓扑排序（你的原始实现，保持不变）
        """
        from collections import deque, defaultdict

        graph = defaultdict(list)
        in_degree = defaultdict(int)
        step_set = set(step_names)

        for step_name in step_names:
            step = self.steps[step_name]
            for dep in step.dependencies:
                if dep in step_set:
                    graph[dep].append(step_name)
                    in_degree[step_name] += 1

        queue = deque([s for s in step_names if in_degree[s] == 0])
        result = []

        while queue:
            step = queue.popleft()
            result.append(step)
            for neighbor in graph[step]:
                in_degree[neighbor] -= 1
                if in_degree[neighbor] == 0:
                    queue.append(neighbor)

        if len(result) != len(step_names):
            cycle = set(step_names) - set(result)
            raise ValueError(f"检测到循环依赖: {cycle}")

        return result

    def get_parallel_groups(self, step_names: List[str]) -> List[List[str]]:
        """
        获取可以并行执行的步骤组（你的原始实现，保持不变）
        """
        from collections import defaultdict

        levels = {}

        def get_level(step_name):
            if step_name in levels:
                return levels[step_name]
            step = self.steps[step_name]
            if not step.dependencies:
                levels[step_name] = 0
                return 0
            # 只考虑在 step_names 里的依赖
            relevant_deps = [d for d in step.dependencies if d in step_names]
            if not relevant_deps:
                levels[step_name] = 0
                return 0
            max_dep_level = max(get_level(dep) for dep in relevant_deps)
            levels[step_name] = max_dep_level + 1
            return max_dep_level + 1

        for step_name in step_names:
            get_level(step_name)

        groups = defaultdict(list)
        for step_name, level in levels.items():
            if step_name in step_names:
                groups[level].append(step_name)

        return [groups[i] for i in sorted(groups.keys())]

    def get_upstream(self, step_name: str) -> Set[str]:
        """获取所有上游依赖（递归，你的原始实现）"""
        upstream = set()

        def collect(step):
            if step in self.steps:
                for dep in self.steps[step].dependencies:
                    if dep not in upstream:
                        upstream.add(dep)
                        collect(dep)

        collect(step_name)
        return upstream

    def get_downstream(self, step_name: str) -> Set[str]:
        """获取所有下游步骤（递归，你的原始实现）"""
        downstream = set()

        for s, step in self.steps.items():
            if step_name in step.dependencies:
                downstream.add(s)
                downstream.update(self.get_downstream(s))

        return downstream
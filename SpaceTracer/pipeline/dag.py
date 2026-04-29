#!/usr/bin/env python3
"""
Pipeline DAG
"""

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Set


@dataclass
class DAGStep:
    name: str
    dependencies: List[str] = field(default_factory=list)
    produces: List[str] = field(default_factory=list)
    run_level: Optional[str] = None
    output_level: Optional[str] = None


class PipelineDAG:
 
    def __init__(self):
        self.steps: Dict[str, DAGStep] = {
            "cluster": DAGStep(
                name="cluster",
                dependencies=[],
                produces=["cluster_file"],
                run_level="sample",
                output_level="sample"
            ),
            "bam_processing": DAGStep(
                name="bam_processing",
                dependencies=[],
                produces=["filtered_bam"],
                run_level="sample",
                output_level="sample"
            ),
            "mpileup": DAGStep(
                name="mpileup",
                dependencies=["bam_processing"],
                produces=["pileup_results"],
                run_level="sample",
                output_level="chunk"
            ),
            "umi_combine": DAGStep(
                name="umi_combine",
                dependencies=["mpileup"],
                produces=["umi_combined"],
                run_level="chunk",
                output_level="chunk"
            ),
            "cell_num": DAGStep(
                name="cell_num",
                dependencies=["cluster","umi_combine"],
                produces=["spot_count_file", "error_count_file"],
                run_level="sample",
                output_level="sample"
            ),
            "prior": DAGStep(
                name="prior",
                dependencies=["umi_combine"],
                produces=["prior_file"],
                run_level="chrom",
                output_level="chrom"
            ),
            "genotyping": DAGStep(
                name="genotyping",
                dependencies=["cluster","prior", "umi_combine","cell_num"],
                produces=["genotype_results"],
                run_level="chunk",
                output_level="chunk"
            ),
            "spatial_feature": DAGStep(
                name="spatial_feature",
                dependencies=["genotyping"],
                produces=["spatial_feature_results"],
                run_level="chunk",
                output_level="chunk"
            ),
            "mappability_feature": DAGStep(
                name="mappability_feature",
                dependencies=["genotyping"],
                produces=["mappability_feature_results"],
                run_level="chrom",
                output_level="chrom"
            ),
            "read_feature": DAGStep(
                name="read_feature",
                dependencies=["genotyping"],
                produces=["read_feature_results"],
                run_level="chunk",
                output_level="chunk"
            ),
            "RNA_feature": DAGStep(
                name="RNA_feature",
                dependencies=["genotyping"],
                produces=["RNA_feature_results"],
                run_level="sample",
                output_level="sample"
            ),
            "phasing":DAGStep(
                name="phasing",
                dependencies=["RNA_feature"],
                produces=["phasing_results","cluster_events"],
                run_level="sample",
                output_level="sample"
            ),
            "merge_feature": DAGStep(
                name="merge_feature",
                dependencies=[
                    "spatial_feature",
                    "mappability_feature",
                    "read_feature",
                    "RNA_feature",
                    "phasing"
                ],
                produces=["merged_features"],
                run_level="sample",
                output_level="sample"
            ),
            "mutation_prediction": DAGStep(
                name="mutation_predict",
                dependencies=["merge_feature"],
                produces=["final_vcf"],
                run_level="sample",
                output_level="sample"
            ),
            # "phylogeny":DAGStep(
            #     name="phylogeny",
            #     dependencies=["mutation_prediction"],
            #     produces=["phylogeny_results"],
            #     run_level="sample",
            #     output_level="sample"
            # ),
        }

    # ─────────────────────────────────────────────────────────────
    # query
    # ─────────────────────────────────────────────────────────────

    def has_step(self, step_name: str) -> bool:
        return step_name in self.steps

    def get_dependencies(self, step_name: str) -> List[str]:
        self._require_step(step_name)
        return list(self.steps[step_name].dependencies)

    def get_upstream(self, step_name: str) -> Set[str]:
        self._require_step(step_name)

        visited = set()

        def dfs(s: str):
            for dep in self.steps[s].dependencies:
                if dep not in visited:
                    visited.add(dep)
                    dfs(dep)

        dfs(step_name)
        return visited

    def _require_step(self, step_name: str) -> None:
        if step_name not in self.steps:
            raise ValueError(f"Unknown step: {step_name}")

    # ─────────────────────────────────────────────────────────────
    # topo sort
    # ─────────────────────────────────────────────────────────────

    def _topological_sort(self, step_names: List[str]) -> List[str]:
        subset = set(step_names)
        in_degree = {s: 0 for s in subset}
        graph = {s: [] for s in subset}

        for step_name in subset:
            for dep in self.steps[step_name].dependencies:
                if dep in subset:
                    graph[dep].append(step_name)
                    in_degree[step_name] += 1

        queue = [s for s in step_names if in_degree[s] == 0]
        result = []

        while queue:
            current = queue.pop(0)
            result.append(current)

            for nxt in graph[current]:
                in_degree[nxt] -= 1
                if in_degree[nxt] == 0:
                    queue.append(nxt)

        if len(result) != len(step_names):
            raise ValueError(f"Cyclic dependency detected in steps: {step_names}")

        return result

    # ─────────────────────────────────────────────────────────────
    # run
    # ─────────────────────────────────────────────────────────────

    def get_execution_plan(
        self,
        start_from: Optional[str] = None,
        stop_at: Optional[str] = None,
        only_steps: Optional[List[str]] = None,
    ) -> List[str]:
      
        all_steps = self._topological_sort(list(self.steps.keys()))

        if only_steps:
            unique_steps = list(dict.fromkeys(only_steps))
            for step_name in unique_steps:
                self._require_step(step_name)
            return self._topological_sort_subset(unique_steps)

        if start_from is not None:
            self._require_step(start_from)
        if stop_at is not None:
            self._require_step(stop_at)

        start_idx = 0
        stop_idx = len(all_steps) - 1

        if start_from is not None:
            start_idx = all_steps.index(start_from)

        if stop_at is not None:
            stop_idx = all_steps.index(stop_at)

        if start_idx > stop_idx:
            raise ValueError(
                f"Invalid range: start_from='{start_from}' is after stop_at='{stop_at}'"
            )

        return all_steps[start_idx: stop_idx + 1]

    def _topological_sort_subset(self, subset_steps: List[str]) -> List[str]:
     
        subset = set(subset_steps)
        in_degree = {s: 0 for s in subset}
        graph = {s: [] for s in subset}

        for step_name in subset:
            for dep in self.steps[step_name].dependencies:
                if dep in subset:
                    graph[dep].append(step_name)
                    in_degree[step_name] += 1

        queue = [s for s in subset_steps if in_degree[s] == 0]
        result = []

        while queue:
            current = queue.pop(0)
            result.append(current)

            for nxt in graph[current]:
                in_degree[nxt] -= 1
                if in_degree[nxt] == 0:
                    queue.append(nxt)

        if len(result) != len(subset_steps):
            raise ValueError(f"Cyclic dependency detected in only_steps: {subset_steps}")

        return result

    # ─────────────────────────────────────────────────────────────
    # smart plan 
    # ─────────────────────────────────────────────────────────────

    def get_missing_dependencies(
        self,
        step_name: str,
        completed_steps: Set[str]
    ) -> List[str]:
       
        self._require_step(step_name)

        upstream = self.get_upstream(step_name)
        missing = [s for s in upstream if s not in completed_steps]

        if not missing:
            return []

        return self._topological_sort(missing)

    # ─────────────────────────────────────────────────────────────
    # parallel groups
    # ─────────────────────────────────────────────────────────────

    def get_parallel_groups(self, steps_to_run: List[str]) -> List[List[str]]:
        
        remaining = list(steps_to_run)
        completed = set()
        groups: List[List[str]] = []

        while remaining:
            current_group = []

            for step_name in remaining:
                deps_in_plan = [
                    dep for dep in self.steps[step_name].dependencies
                    if dep in steps_to_run
                ]
                if all(dep in completed for dep in deps_in_plan):
                    current_group.append(step_name)

            if not current_group:
                raise ValueError(
                    f"Cannot build parallel groups, possible cycle in plan: {steps_to_run}"
                )

            groups.append(current_group)

            for step_name in current_group:
                completed.add(step_name)
                remaining.remove(step_name)

        return groups

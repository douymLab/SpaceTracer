#!/usr/bin/env python3
"""
Pipeline Orchestrator
"""

import json
import logging
import os
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple

from SpaceTracer.pipeline.checkpoint import CheckpointManager
from SpaceTracer.pipeline.registry import get_step_class
from SpaceTracer.pipeline.dag import PipelineDAG
from SpaceTracer.utils.parallel import memory_checkpoint

logger = logging.getLogger(__name__)
os.environ["NUMBA_THREADING_LAYER"] = "workqueue"


_PATH_LIKE_EXTENSIONS = {
    ".bam", ".txt", ".tsv", ".vcf", ".bed",
    ".mpileup", ".json", ".fa", ".fasta", ".gz", ".parquet",
}


def _resolve_paths(override: Dict) -> Dict:
    """
    Walk override dict and convert string values that look like file paths
    into Path objects. Everything else is passed through unchanged.
    """
    out = {}
    for k, v in override.items():
        if isinstance(v, str):
            p = Path(v)
            if p.suffix in _PATH_LIKE_EXTENSIONS or ("/" in v or "\\" in v):
                out[k] = str(p)
                continue
        out[k] = v
    return out


class PipelineOrchestrator:

    def __init__(
        self,
        config: Dict,
        *,
        force: bool = False,
        chunk_index: Optional[int] = None,
    ):
        self.config = config
        self.force = force
        self.chunk_index = chunk_index
        self.disabled=False # 默认打开CheckpointManager

        # cell number
        cell_num_str = self.config.get("steps", {}).get("cell_number","")
        if cell_num_str and os.path.exists(cell_num_str):
            cell_num = cell_num_str
        else:
            cell_num = int(cell_num_str) if cell_num_str else 0
        self.config["cell_num"] = cell_num

        # cluster
        cluster_str = self.config.get("steps", {}).get("cluster", {}).get("cluster_file")
        cluster_path = Path(cluster_str) if cluster_str else None
        if cluster_path and cluster_path.exists():
            self.config["cluster"] = str(cluster_path)
        else:
            if self.config.get("sequence_type") == "visium":
                self.config["cluster"] = 0 # 0 means the cluster stpe will run
            else:
                pass

        output_dir = Path(self.config.get("output_dir"))
        if chunk_index is not None:
            self.work_dir = output_dir / f"chunk_{chunk_index:04d}"
        else:
            self.work_dir = output_dir
        self.work_dir.mkdir(parents=True, exist_ok=True)

        self.checkpoint = CheckpointManager(self.work_dir, disabled=self.disabled)
        # self.validator = Validator(self.config.get("run", {}).get("skip_validation", False))
        self.skip_validation = self.config.get("run", {}).get("skip_validation", False)
        self.dag = PipelineDAG()
        self._step_instances: Dict[str, Any] = {}

        self.stats = {
            "start_time": None,
            "end_time": None,
            "step_times": {},
            "step_status": {},
            "chunk_index": chunk_index,
        }

    # ─────────────────────────────────────────────────────────────
    # 主入口
    # ─────────────────────────────────────────────────────────────

    def run(
        self,
        *,
        start_from: Optional[str] = None,
        stop_at: Optional[str] = None,
        only_steps: Optional[List[str]] = None,
        enable_parallel: bool = True,
        context_override: Optional[Dict[str, Any]] = None,
    ) -> Dict:
        logger.info("=" * 70)
        label = f"chunk {self.chunk_index:04d}" if self.chunk_index is not None else "full run"
        logger.info(f"Starting Pipeline ({label})")
        logger.info("=" * 70)

        if only_steps and (start_from is not None or stop_at is not None):
            raise ValueError("only_steps cannot be used together with start_from/stop_at")

        self.stats["start_time"] = time.time()

        #build context
        context = self._build_base_context() # base
        self.base_keys=list(context.keys()) # base keys will not be validated
        context = self._load_completed_outputs(context) # load from pipeline check

        if context_override: # override from outside
            resolved = _resolve_paths(context_override)
            context.update(resolved)

            logger.info(f"Context override keys: {sorted(resolved.keys())}")

        # Run
        if only_steps: # only 1 step
            steps_to_run = self._get_only_steps_plan(only_steps)
            logger.info(
                f"Only-steps plan ({len(steps_to_run)} steps): {' → '.join(steps_to_run)}"
            )
            context = self._run_only_steps(steps_to_run, context)
        else: # build run steps 
            requested_steps = self.dag.get_execution_plan(
                start_from=start_from,
                stop_at=stop_at,
                only_steps=None,
            )

            if self.force:
                steps_to_run = requested_steps
                logger.info("Force mode: running all requested steps")
            else:
                steps_to_run = self._smart_plan(requested_steps, context_override)

            logger.info(
                f"Execution plan ({len(steps_to_run)} steps): "
                f"{' → '.join(steps_to_run)}"
            )

            if enable_parallel:
                context = self._run_with_parallel_groups(steps_to_run, context)
            else:
                context = self._run_sequential(steps_to_run, context)

        self.stats["end_time"] = time.time()
        elapsed = self.stats["end_time"] - self.stats["start_time"]

        logger.info(f"Pipeline finished in {elapsed:.2f}s ({label})")

        return {
            "final_vcf": context.get("final_variants"),
            "elapsed_time": elapsed,
            "chunk_index": self.chunk_index,
            "stats": self.stats,
        }

    # ─────────────────────────────────────────────────────────────
    # context
    # ─────────────────────────────────────────────────────────────

    def _build_base_context(self) -> Dict[str, Any]:
        return {
            "bam_file": self.config["bam_file"],
            "config": self.config,
            "chunk_index": self.chunk_index,
        }

    def _load_completed_outputs(self, context: Dict[str, Any]) -> Dict[str, Any]:
        context = self.checkpoint.load_all_completed_outputs(context)
        return context

    # ─────────────────────────────────────────────────────────────
    # plan
    # ─────────────────────────────────────────────────────────────

    def _smart_plan(
        self,
        requested_steps: List[str],
        context_override: Optional[Dict[str, Any]],
    ) -> List[str]:
        logger.info("Analyzing checkpoint and outputs...")

        completed = set()

        for step_name in self.dag.steps.keys():
            is_complete = self.checkpoint.is_complete(step_name)
            logger.info(f"  Checking {step_name}: complete={is_complete}")
            if is_complete:
                completed.add(step_name)

        logger.info(f"Completed steps: {completed}")
        logger.info(f"Requested steps: {requested_steps}")

        if context_override:
            for step_name in self.dag.steps.keys():
                step = self.dag.steps[step_name]
                produces = getattr(step, "produces", [])
                if produces and all(key in context_override for key in produces):
                    completed.add(step_name)
                    logger.info(f"  ✓ {step_name} (provided by context_override)")

        needed = set()

        for step_name in requested_steps:
            logger.info(f"  Processing {step_name}: in_completed={step_name in completed}")

            if step_name in completed:
                logger.info(f"  ⊘ {step_name} (skipping, already complete)")
                continue

            needed.add(step_name)

            missing_deps = self.dag.get_missing_dependencies(step_name, completed)
            if missing_deps:
                logger.info(f"  + {step_name} requires: {missing_deps}")
                needed.update(missing_deps)

        logger.info(f"Steps to run: {needed}")

        if not needed:
            logger.info("All requested steps are complete!")
            return []

        return self.dag._topological_sort(list(needed))

    # def _get_only_steps_plan(self, only_steps: List[str]) -> List[str]:
    #     if not only_steps:
    #         return []

    #     unknown = [s for s in only_steps if s not in self.dag.steps]
    #     if unknown:
    #         raise ValueError(f"Unknown steps in only_steps: {unknown}")

    #     unique_steps = list(dict.fromkeys(only_steps))
    #     return self._topological_sort_subset(unique_steps)

    def _get_only_steps_plan(self, only_steps: List[str]) -> List[str]:
        return self.dag.get_execution_plan(only_steps=only_steps)

    def _topological_sort_subset(self, subset_steps: List[str]) -> List[str]:
        """
        只在 subset 内部排序，不补集合外依赖
        """
        subset = set(subset_steps)
        in_degree = {s: 0 for s in subset}
        graph = {s: [] for s in subset}

        for step_name in subset:
            node = self.dag.steps[step_name]
            deps = getattr(node, "dependencies", [])
            for dep in deps:
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
            raise ValueError(f"Cyclic dependency detected within only_steps: {subset_steps}")

        return result

    # ─────────────────────────────────────────────────────────────
    # merge
    # ─────────────────────────────────────────────────────────────

    def _merge_new_outputs(self, context: Dict[str, Any], step_name: str, new_outputs: Dict[str, Any]) -> None:
        if new_outputs:
            for key, value in new_outputs.items():
                if key not in context:
                    context[key] = value
                elif context[key] == value:
                    continue
                else:
                    raise ValueError(
                        f"Output key conflict for step '{step_name}': "
                        f"key='{key}', existing='{context[key]}', new='{value}'"
                    )

    # ─────────────────────────────────────────────────────────────
    # step execution
    # ─────────────────────────────────────────────────────────────

    def _execute_step(self, step_name: str, step, context: Dict[str, Any]) -> Dict[str, Any]:
        logger.info(f"[{step_name}] starting ...")

        t0 = time.time()
        # logger.info(f"@@@@@@@@ base_keys: {self.base_keys}")
        
        # context_dict={k: v for k, v in context.items() if k not in self.base_keys}
        # logger.info(f"====== raw context: {context_dict}")

        try:
            
            context = step.execute(context,self.skip_validation)
            # context_dict={k: v for k, v in context.items() if k not in self.base_keys}
            # logger.info(f"++++++++ run context: {context_dict}")

            new_outputs = step.get_outputs(context)

            memory_checkpoint(step_name,logger=logger, top_n=8, do_gc=True)
            
            # logger.info(f"-------- only outputs: {context_dict}")

            if not isinstance(new_outputs, dict):
                raise TypeError(
                    f"Step '{step_name}' execute() must return dict, got {type(new_outputs)}"
                )

            filtered_outputs = {k: v for k, v in new_outputs.items() if k not in self.base_keys}
            # logger.info(f"******** only filtered outputs: {new_outputs}")

            # 简单输出存在性校验
            # self._validate_step_outputs_exist(step_name, filtered_outputs) # no need to validate here

            elapsed = time.time() - t0
            self.stats["step_times"][step_name] = elapsed   
            self.stats["step_status"][step_name] = "success"

            logger.info(f"[{step_name}] ✓ {elapsed:.2f}s")
            return filtered_outputs

        except Exception as exc:
            self.stats["step_status"][step_name] = "failed"
            self.checkpoint.mark_failed(step_name, str(exc))
            logger.error(f"[{step_name}] ✗ {exc}", exc_info=True)

            return {}
            # raise

    def _validate_step_outputs_exist(self, step_name: str, new_outputs: Dict[str, Any]) -> None:
        for key, value in new_outputs.items():
            if key not in self.base_keys:
                if value is None:
                    raise ValueError(f"Step '{step_name}' output '{key}' is None")

                if not isinstance(value,int):
                    path = Path(value)
                    if not path.exists():
                        raise FileNotFoundError(
                            f"Step '{step_name}' output '{key}' does not exist: {value}"
                        )

    # ─────────────────────────────────────────────────────────────
    # sequential
    # ─────────────────────────────────────────────────────────────

    def _run_sequential(self, steps_to_run: List[str], context: Dict[str, Any]) -> Dict[str, Any]:
        for step_name in steps_to_run:
            logger.info("=" * 70)
            logger.info(f"Running step: {step_name}")
            logger.info("=" * 70)

            if not self.force and self.checkpoint.is_complete(step_name):
                logger.info(f"Step {step_name} already complete, skipping")
                context = self.checkpoint.load_outputs_to_context(step_name, context)
                self.stats["step_status"][step_name] = "skipped"
                continue

            step = self._get_step_instance(step_name, context)
            new_outputs = self._execute_step(step_name, step, context)
            self._merge_new_outputs(context, step_name, new_outputs)
            self.checkpoint.mark_complete(step_name, new_outputs)

            if not self.skip_validation:
                if not new_outputs:
                    raise ValueError(f"Output validation failed: {step_name}")

        return context

    # ─────────────────────────────────────────────────────────────
    # parallel
    # ─────────────────────────────────────────────────────────────

    def _run_with_parallel_groups(
        self,
        steps_to_run: List[str],
        context: Dict[str, Any],
    ) -> Dict[str, Any]:
        groups = self.dag.get_parallel_groups(steps_to_run)
        # logger.info(f'{groups}')
        max_workers = self.config.get("runtime", {}).get("max_parallel", 4)

        for group in groups:
            if len(group) == 1:
                step_name = group[0]
                check=self.checkpoint.is_complete(step_name)
                # logger.info(f"{self.force}  {check}")
                if not self.force and self.checkpoint.is_complete(step_name):
                    logger.info(f"Step {step_name} already complete, skipping")
                    # context = self.checkpoint.load_outputs_to_context(step_name, context)
                    self.stats["step_status"][step_name] = "skipped"
                    continue

                step = self._get_step_instance(step_name, context)

                new_outputs = self._execute_step(step_name, step, context)

                self._merge_new_outputs(context, step_name, new_outputs)
                if new_outputs:
                    self.checkpoint.mark_complete(step_name, new_outputs)
                else:
                    self.checkpoint.mark_failed(step_name, "No outputs.")
                    raise RuntimeError(f"Step '{step_name}' has no outputs!")


                if not self.skip_validation:
                    if not new_outputs:
                        raise ValueError(f"Output validation failed: {step_name}")

            else:
                logger.info(f"Running in parallel: {group}")

                futures = {}
                base_context = context.copy()

                with ThreadPoolExecutor(max_workers=min(len(group), max_workers)) as pool:
                    for step_name in group:
                        if not self.force and self.checkpoint.is_complete(step_name):
                            logger.info(f"Step {step_name} already complete, skipping")
                            # context = self.checkpoint.load_outputs_to_context(step_name, context)

                            self.stats["step_status"][step_name] = "skipped"
                            continue

                        step = self._get_step_instance(step_name, base_context)
                        fut = pool.submit(
                            self._execute_step,
                            step_name,
                            step,
                            base_context.copy()
                        )
                        futures[fut] = step_name

                    not_finish_step=[]
                    for fut in as_completed(futures):
                        step_name = futures[fut]
                        new_outputs = fut.result()

                        if new_outputs:
                            self._merge_new_outputs(context, step_name, new_outputs)
                            self.checkpoint.mark_complete(step_name, new_outputs)
                        else:
                            not_finish_step.append(step_name)
                    
                    if not self.skip_validation and not_finish_step:
                        raise ValueError(f'The following steps have not finished: {not_finish_step}!')
                        # if not self.config.get("skip_validation", False):
                        #     if not self.validator.validate_step_output(step_name, context):
                        #         raise ValueError(f"Output validation failed: {step_name}")

        return context

    # ─────────────────────────────────────────────────────────────
    # only_steps
    # ─────────────────────────────────────────────────────────────

    def _check_step_inputs_available(self, step, context: Dict[str, Any]) -> Tuple[List[str], List[Tuple[str, str]]]:
        missing_keys = []
        missing_files = []

        try:
            inputs = step.get_inputs(context)
        except KeyError as e:
            missing_keys.append(str(e).strip("'"))
            return missing_keys, missing_files

        for key, path in inputs.items():
            if path is None:
                missing_files.append((key, str(path)))
                continue
            if not Path(path).exists():
                missing_files.append((key, str(path)))

        return missing_keys, missing_files

    def _run_only_steps(self, only_steps_plan: List[str], context: Dict[str, Any]) -> Dict[str, Any]:
        for step_name in only_steps_plan:
            logger.info("=" * 70)
            logger.info(f"Running only-step: {step_name}")
            logger.info("=" * 70)

            if not self.force and self.checkpoint.is_complete(step_name):
                logger.info(f"Step {step_name} already complete, skipping")
                # context = self.checkpoint.load_outputs_to_context(step_name, context)
                self.stats["step_status"][step_name] = "skipped"
                continue

            step = self._get_step_instance(step_name, context)

            missing_keys, missing_files = self._check_step_inputs_available(step, context)
            if missing_keys or missing_files:
                messages = []
                if missing_keys:
                    messages.append(f"missing context keys: {missing_keys}")
                if missing_files:
                    messages.append(f"missing files: {missing_files}")

                raise ValueError(
                    f"Step '{step_name}' cannot run in only_steps mode; "
                    + "; ".join(messages)
                )

            new_outputs = self._execute_step(step_name, step, context)
            self._merge_new_outputs(context, step_name, new_outputs)

            if new_outputs:
                self.checkpoint.mark_complete(step_name, new_outputs)
            else:
                self.checkpoint.mark_failed(step_name,"The output is empty, please check the log info.")
                    
            if not self.skip_validation:
                if not new_outputs:
                    raise ValueError(f"Output validation failed: {step_name}")

        return context

    # ─────────────────────────────────────────────────────────────
    # step instance
    # ─────────────────────────────────────────────────────────────

    def _get_step_instance(self, step_name: str, context: Dict[str, Any]):
        cls = get_step_class(step_name)
        if cls is None:
            raise ValueError(
                f"No implementation class registered for '{step_name}'. "
            )
        return cls(step_name, context, self.checkpoint)

    # ─────────────────────────────────────────────────────────────
    # diagnostics
    # ─────────────────────────────────────────────────────────────

    def show_plan(
        self,
        start_from: Optional[str] = None,
        stop_at: Optional[str] = None,
        only_steps: Optional[List[str]] = None,
    ):
        if only_steps:
            plan = self._get_only_steps_plan(only_steps)
        else:
            plan = self.dag.get_execution_plan(start_from, stop_at, None)

        groups = self.dag.get_parallel_groups(plan)
        label = (
            f"chunk {self.chunk_index:04d}"
            if self.chunk_index is not None else "full run"
        )

        print(f"\n📋  Execution plan  [{label}]")
        print("─" * 52)
        for i, group in enumerate(groups):
            if len(group) == 1:
                print(f"  layer {i}:  {group[0]}")
            else:
                print(f"  layer {i}:  {' | '.join(group)}   ← parallel")
        print("─" * 52)
        print(f"  {len(plan)} steps  ·  {len(groups)} layers")

    # ─────────────────────────────────────────────────────────────
    # CLI helper
    # ─────────────────────────────────────────────────────────────

    @classmethod
    def from_cli_args(cls, args) -> "PipelineOrchestrator":
        override = None
        raw = getattr(args, "context_override", None)
        if raw:
            try:
                override = json.loads(raw)
            except json.JSONDecodeError as exc:
                raise ValueError(f"--context-override is not valid JSON: {exc}") from exc

        instance = cls(
            config=args.config,
            force=getattr(args, "force", False),
            chunk_index=getattr(args, "chunk_index", None),
        )
        instance._pending_context_override = override
        return instance

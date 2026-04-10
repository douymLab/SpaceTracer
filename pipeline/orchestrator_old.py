#!/usr/bin/env python3
"""
PipelineOrchestrator  —  two external-input interfaces

Interface 1 │ context_override  (in run())
            │   Inject or overwrite any context key before the first step
            │   runs.  Lets CLI / Snakemake skip upstream computation by
            │   supplying pre-existing files directly.
            │
Interface 2 │ chunk_index  (in __init__)
            │   Mark this Orchestrator instance as handling one split chunk.
            │   Automatically isolates work_dir + checkpoint so N parallel
            │   Snakemake/SLURM jobs never collide.
            │   Also stored in context['chunk_index'] for steps to use.

──────────────────────────────────────────────────────────────────────────
Snakemake integration sketch
──────────────────────────────────────────────────────────────────────────

  # Snakefile
  CHUNKS = list(range(100))   # produced by a split step

  rule run_chunk:
      input:
          chunk = "splits/chunk_{index}.mpileup"
      output:
          result = "results/chunk_{index}/hfdr_results.tsv"
      params:
          index = lambda wc: int(wc.index)
      shell:
          \"\"\"
          mutation-caller run-chunk \\
              --bam        {input.bam} \\
              --output     results/chunk_{wildcards.index}/ \\
              --chunk-index {params.index} \\
              --context-override '{{"filtered_candidates": "{input.chunk}"}}' \\
              --start-from umi_combine \\
              --stop-at    hfdr
          \"\"\"

  # The key insight:
  #   --chunk-index   → isolates paths (work_dir / checkpoint)
  #   --context-override → injects the pre-split input file into context
  #   --start-from    → skips the steps already done by the split step
"""

import json
import logging
import os
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any, Dict, List, Optional, Set

from SpaceTracer.pipeline.checkpoint import CheckpointManager
from SpaceTracer.pipeline.validator import Validator
from SpaceTracer.pipeline.dag import PipelineDAG
from SpaceTracer.steps.step1_bam_processing import BamProcessingStep
from SpaceTracer.steps.step2_mpileup import MpileupStep
from SpaceTracer.steps.step3_UMI_combine import UMICombineStep
from SpaceTracer.steps.step3_cell_number import CellNumStep
from SpaceTracer.steps.step3_get_prior import PriorCalculator
from SpaceTracer.steps.step4_genotyping import GenotypingStep
from SpaceTracer.steps.step5_RNA_level_feature import RNAFeatureStep
from SpaceTracer.steps.step5_read_feature import ReadFeatureStep
from SpaceTracer.steps.step5_spatial_feature import SpatialFeatureStep
from SpaceTracer.steps.step5_mappability_feature import MappabilityFeatureStep
from SpaceTracer.steps.step6_merge_all_features import MergeFeatureStep
from SpaceTracer.utils.get_genome_info import GenomeDetails

logger = logging.getLogger(__name__)


# ─────────────────────────────────────────────────────────────────────────────
# Orchestrator
# ─────────────────────────────────────────────────────────────────────────────

class PipelineOrchestrator:

    # Subclass or register concrete Step classes here


    # ── Construction ──────────────────────────────────────────────────────────

    def __init__(
        self,
        config: Dict,
        *,
        force: bool = False,
        # ═══════════════════════════════════════════════════════════════
        # INTERFACE 2 — chunk_index
        #
        # Identifies which split file this job is processing.
        # None  →  single full run  (normal mode)
        # int   →  one chunk of a scatter-gather run  (Snakemake mode)
        #
        # What it does automatically:
        #   • work_dir  = <output_dir>/work/           (None)
        #              = <output_dir>/work/chunk_0003/ (chunk_index=3)
        #   • checkpoint is scoped to that sub-dir, so jobs are isolated
        #   • stored in context['chunk_index'] for steps to read
        # ═══════════════════════════════════════════════════════════════
        chunk_index: Optional[int] = None,
    ):
        self.config=config
        STEP_CLASSES = {
            'bam_processing':BamProcessingStep,
            'mpileup': MpileupStep,
            'umi_combine': UMICombineStep,
            'cell_num': CellNumStep,
            'prior': PriorCalculator,
            'genotyping': GenotypingStep,
            'spatial_feature': SpatialFeatureStep,
            'mappability_feature': MappabilityFeatureStep,
            'read_feature': ReadFeatureStep,
            'RNA_feature': RNAFeatureStep,
            'merge_feature':MergeFeatureStep
        }

        # cluster?
        cluster_str=self.config.get('steps').get('cluster').get('cluster_file')
        cluster_path=Path(cluster_str) if cluster_str else ""
        print(cluster_path)
        if cluster_path and cluster_path.exists():
            self.config['cluster']=cluster_str
        else:
            if self.config["sequence_type"]=="visium":
                self.config['cluster']=0 # 0 means we'll run cluster function
            else:
                raise ValueError(f'Unsupported sequence_type for cluster (Only "visium" supported). Please provide cluster file in steps-cluster-cluster_file.')
        
        # cell number?
        cell_num_str=self.config.get('steps').get('cell_number')
        if os.path.exists(cell_num_str):
            cell_num=cell_num_str
        else:
            cell_num = int(cell_num_str) if cell_num_str else 0 # we wil run cell number function, if cell_num_str is empty
        self.config['cell_num']=cell_num


        # keep_intermediates
        self.keep_intermediates = self.config.get('run').get('keep_intermediates',False)
        
        self.config       = config
        self.force        = force
        self.chunk_index  = chunk_index
        self.STEP_CLASSES =STEP_CLASSES
        # work_dir is chunk-aware — parallel jobs never share a directory
        output_dir=Path(self.config.get('output_dir'))
        if chunk_index is not None:
            self.work_dir =  output_dir / f'chunk_{chunk_index:04d}'
        else:
            self.work_dir = output_dir
        self.work_dir.mkdir(parents=True, exist_ok=True)


        # Checkpoint scoped to work_dir so each chunk has its own state
        self.checkpoint = CheckpointManager(self.work_dir, disabled=force)
        self.validator   = Validator(self.config.get('run').get('skip_validation', False))
        self.dag         = PipelineDAG()
        self._step_instances: Dict = {}
        self.stats = {
            'start_time':  None,
            'end_time':    None,
            'step_times':  {},
            'step_status': {},
            'chunk_index': chunk_index,
        }

    # ── Primary run interface ─────────────────────────────────────────────────
        
    def run(
        self,
        *,
        start_from: Optional[str] = None,
        stop_at: Optional[str] = None,
        only_steps: Optional[List[str]] = None,
        enable_parallel: bool = False,
        context_override: Optional[Dict[str, Any]] = None,
    ) -> Dict:
        
        logger.info("=" * 70)
        label = f"chunk {self.chunk_index:04d}" if self.chunk_index is not None else "full run"
        logger.info(f"Starting Pipeline ({label})")
        logger.info("=" * 70)
        
        self.stats['start_time'] = time.time()
        
        # ── Step 1: 确定用户想跑哪些步骤 ──
        requested_steps = self.dag.get_execution_plan(
            start_from=start_from,
            stop_at=stop_at,
            only_steps=only_steps,
        )
        
        # ── Step 2: 智能计划（如果 force=False）──
        if self.force:
            steps_to_run = requested_steps
            logger.info("Force mode: running all requested steps")
        else:
            steps_to_run = self._smart_plan(requested_steps, context_override)
        
        logger.info(
            f"Execution plan ({len(steps_to_run)} steps): "
            f"{' → '.join(steps_to_run)}"
        )
        
        # ── Step 3: 构建 context ──
        context = self._build_base_context()
        
        # 加载已完成步骤的输出
        if not self.force:
            context = self._load_completed_outputs(context)
        
        # 应用 context_override
        if context_override:
            resolved = _resolve_paths(context_override)
            context.update(resolved)
            logger.info(f"Context override keys: {sorted(resolved.keys())}")
        
        # ── Step 4: 执行 ──
        if enable_parallel:
            context = self._run_with_parallel_groups(steps_to_run, context)
        else:
            context = self._run_sequential(steps_to_run, context)
        
        self.stats['end_time'] = time.time()
        elapsed = self.stats['end_time'] - self.stats['start_time']
        
        logger.info(f"Pipeline finished in {elapsed:.2f}s ({label})")
        
        return {
            'final_vcf': context.get('final_variants'),
            'elapsed_time': elapsed,
            'chunk_index': self.chunk_index,
            'stats': self.stats,
        }

    def _smart_plan(self, 
                    requested_steps: List[str],
                    context_override: Optional[Dict]) -> List[str]:
        logger.info("Analyzing checkpoint and outputs...")
        
        # ── 1. 找出已完成的步骤 ──
        completed = set()
        
        for step_name in self.dag.steps.keys():
            is_complete = self.checkpoint.check_outputs_exist(step_name)
            logger.info(f"  Checking {step_name}: complete={is_complete}")  # ← 添加
            
            if is_complete:
                completed.add(step_name)
                logger.debug(f"  ✓ {step_name} (outputs exist)")
        
        logger.info(f"Completed steps: {completed}")  # ← 添加
        logger.info(f"Requested steps: {requested_steps}")  # ← 添加
        
        # ── 2. 检查 context_override 提供了哪些输出 ──
        if context_override:
            for step_name in self.dag.steps.keys():
                step = self.dag.steps[step_name]
                if all(key in context_override for key in step.produces):
                    completed.add(step_name)
                    logger.info(f"  ✓ {step_name} (provided by context_override)")
        
        # ── 3. 计算需要运行的步骤 ──
        needed = set()
        
        for step_name in requested_steps:
            logger.info(f"  Processing {step_name}: in_completed={step_name in completed}")  # ← 添加
            
            if step_name in completed:
                logger.info(f"  ⊘ {step_name} (skipping, already complete)")
                continue
            
            needed.add(step_name)
            
            missing_deps = self.dag.get_missing_dependencies(step_name, completed)
            if missing_deps:
                logger.info(f"  + {step_name} requires: {missing_deps}")
                needed.update(missing_deps)
        
        logger.info(f"Steps to run: {needed}")  # ← 添加
        
        # ── 4. 拓扑排序 ──
        if not needed:
            logger.info("All requested steps are complete!")
            return []
        
        return self.dag._topological_sort(list(needed))

    # ── Context construction ──────────────────────────────────────────────────

    def _build_base_context(self) -> Dict:
        """
        Default context.  chunk_index is always injected here so every
        step can read context['chunk_index'] without special handling.
        """

        context = {
            'bam_file': self.config['bam_file'],
            'regions_file': self.config['regions_file'],
            'config': self.config,
            'chunk_index': self.chunk_index   # None for non-chunked runs

        }
        return context

    def _load_completed_outputs(self, context: Dict) -> Dict:
        """
        加载所有已完成步骤的输出到 context
        
        新增方法
        """
        for step_name in self.checkpoint.get_completed_steps():
            context = self.checkpoint.load_outputs_to_context(step_name, context)
            logger.debug(f"Loaded outputs from {step_name}")
        
        return context

    # ── Class-level factory: build from parsed CLI args ───────────────────────

    @classmethod
    def from_cli_args(cls, args) -> 'PipelineOrchestrator':
        """
        Build an Orchestrator from an argparse Namespace (or similar).
        Decodes --context-override JSON string automatically.

        Attaches decoded override as ._pending_context_override so the
        caller can pass it directly into run():

            orchestrator = PipelineOrchestrator.from_cli_args(args)
            orchestrator.run(
                start_from=args.start_from,
                context_override=orchestrator._pending_context_override,
            )
        """
        override = None
        raw = getattr(args, 'context_override', None)
        if raw:
            try:
                override = json.loads(raw)
            except json.JSONDecodeError as exc:
                raise ValueError(
                    f"--context-override is not valid JSON: {exc}"
                ) from exc

        instance = cls(
            bam_file=args.bam,
            output_dir=args.output,
            config=args.config,
            regions_file=getattr(args, 'regions', None),
            force=getattr(args, 'force', False),
            chunk_index=getattr(args, 'chunk_index', None),  # ← Interface 2
        )
        instance._pending_context_override = override        # ← Interface 1
        return instance

    # ── Sequential execution ──────────────────────────────────────────────────

    def _run_sequential(self, steps_to_run: List[str], context: Dict) -> Dict:
        """顺序执行步骤"""
        for step_name in steps_to_run:
            logger.info(f"{'='*70}")
            logger.info(f"Running step: {step_name}")
            logger.info(f"{'='*70}")
            
            try:
                # 创建步骤实例
                step = self._get_step_instance(step_name, context)
                context = self._execute_step(step_name, step, context)
                
                # 执行步骤
                # context = step.run(context)
                
                # ← 修改这里：传入 step_instance
                self.checkpoint.mark_complete(step_name, context, step_instance=step)
                
            except Exception as e:
                logger.error(f"Step {step_name} failed: {e}", exc_info=True)
                self.checkpoint.mark_failed(step_name, str(e))
                raise
        
        return context


    # ── Parallel-group execution ──────────────────────────────────────────────

    def _run_with_parallel_groups(
        self, steps_to_run: List[str], context: Dict
    ) -> Dict:
        groups = self.dag.get_parallel_groups(steps_to_run)
        max_workers = self.config.get('runtime', {}).get('max_parallel', 4)

        for group in groups:
            if len(group) == 1:
                step_name = group[0]
                step = self._get_step_instance(step_name, context)
                if self.checkpoint.is_complete(step_name):
                    context = step.load_outputs(context)
                    self.stats['step_status'][step_name] = 'skipped'
                else:
                    context = self._execute_step(step_name, step, context)
            else:
                logger.info(f"Running in parallel: {group}")
                futures: Dict = {}
                with ThreadPoolExecutor(
                    max_workers=min(len(group), max_workers)
                ) as pool:
                    for step_name in group:
                        if self.checkpoint.is_complete(step_name):
                            self.stats['step_status'][step_name] = 'skipped'
                            continue
                        step = self._get_step_instance(step_name, context)
                        futures[pool.submit(
                            self._execute_step, step_name, step, context.copy()
                        )] = step_name

                for fut in as_completed(futures):
                    sn = futures[fut]
                    try:
                        context.update(fut.result())
                    except Exception as exc:
                        logger.error(f"[{sn}] ✗ {exc}")
                        raise

        return context

    # ── Single step execution ─────────────────────────────────────────────────

    def _execute_step(self, step_name: str, step, context: Dict) -> Dict:
        logger.info(f"[{step_name}] starting …")
        t0 = time.time()
        try:
            context = step.execute(context)
            elapsed = time.time() - t0
            self.stats['step_times'][step_name]  = elapsed
            self.stats['step_status'][step_name] = 'success'
            logger.info(f"[{step_name}] ✓  {elapsed:.2f}s")
            self.checkpoint.mark_complete(step_name, context)
            if not self.config.get('skip_validation', False):
                if not self.validator.validate_step_output(step_name, context):
                    raise ValueError(f"Output validation failed: {step_name}")
        except Exception as exc:
            self.stats['step_status'][step_name] = 'failed'
            self.checkpoint.mark_failed(step_name, str(exc))
            logger.error(f"[{step_name}] ✗  {exc}")
            raise
        return context

    # ── Step instance cache ───────────────────────────────────────────────────
    def _get_step_instance(self, step_name: str, context: Dict):
        if step_name not in self._step_instances:
            print(self.STEP_CLASSES)
            cls = self.STEP_CLASSES.get(step_name)
            if cls is None:
                raise ValueError(
                    f"No implementation class registered for '{step_name}'. "
                    "Add it to STEP_CLASSES."
                )
            self._step_instances[step_name] = cls(
                step_name, context, self.checkpoint
            )
        return self._step_instances[step_name]

    # ── Diagnostics ───────────────────────────────────────────────────────────

    def show_plan(
        self,
        start_from=None, stop_at=None, only_steps=None
    ):
        plan   = self.dag.get_execution_plan(start_from, stop_at, only_steps)
        groups = self.dag.get_parallel_groups(plan)
        label  = (f"chunk {self.chunk_index:04d}"
                    if self.chunk_index is not None else "full run")
        print(f"\n📋  Execution plan  [{label}]")
        print("─" * 52)
        for i, group in enumerate(groups):
            if len(group) == 1:
                print(f"  layer {i}:  {group[0]}")
            else:
                print(f"  layer {i}:  {' | '.join(group)}   ← parallel")
        print("─" * 52)
        print(f"  {len(plan)} steps  ·  {len(groups)} layers\n")

    def get_execution_plan(self,
                        start_from: Optional[str] = None,
                        stop_at: Optional[str] = None,
                        only_steps: Optional[List[str]] = None,
                        # ← 新增参数
                        completed_steps: Optional[Set[str]] = None) -> List[str]:
        """
        生成执行计划（拓扑排序）
        
        新增功能：
        - completed_steps: 已完成的步骤集合
        - 自动过滤掉已完成的步骤
        """
        if only_steps:
            plan = self._get_subgraph_with_deps(only_steps)
        elif start_from or stop_at:
            plan = self._get_range_subgraph(start_from, stop_at)
        else:
            plan = self._topological_sort(list(self.steps.keys()))
        
        # ← 新增：过滤已完成的步骤
        if completed_steps:
            plan = [s for s in plan if s not in completed_steps]
        
        return plan
    
    def get_missing_dependencies(self, 
                                step_name: str,
                                completed_steps: Set[str]) -> List[str]:
        """
        获取缺失的依赖步骤
        
        新增方法
        
        Args:
            step_name: 目标步骤
            completed_steps: 已完成的步骤集合
        
        Returns:
            需要运行的依赖步骤列表（拓扑排序）
        """
        # 获取所有上游依赖
        all_deps = self.get_upstream(step_name)
        
        # 找出未完成的依赖
        missing = all_deps - completed_steps
        
        # 拓扑排序
        if missing:
            return self._topological_sort(list(missing))
        return []

# ─────────────────────────────────────────────────────────────────────────────
# Helper
# ─────────────────────────────────────────────────────────────────────────────

_PATH_LIKE_EXTENSIONS = {
    '.bam', '.txt', '.tsv', '.vcf', '.bed',
    '.mpileup', '.json', '.fa', '.fasta', '.gz',
}

def _resolve_paths(override: Dict) -> Dict:
    """
    Walk override dict and convert string values that look like file paths
    into Path objects.  Everything else is passed through unchanged.
    """
    out = {}
    for k, v in override.items():
        if isinstance(v, str):
            p = Path(v)
            if p.suffix in _PATH_LIKE_EXTENSIONS or ('/' in v or '\\' in v):
                out[k] = p
                continue
        out[k] = v
    return out


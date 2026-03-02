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
from typing import Any, Dict, List, Optional

from SpaceTracer.pipeline.checkpoint import CheckpointManager
from SpaceTracer.pipeline.validator import Validator
from SpaceTracer.pipeline.dag import PipelineDAG
from SpaceTracer.steps.step1_bam_processing import BamProcessingStep
from SpaceTracer.steps.step2_mpileup import MpileupStep
from SpaceTracer.steps.step3_UMI_combine import UMICombineStep
from SpaceTracer.steps.step3_cell_number import CellNumStep
from SpaceTracer.steps.step3_get_prior import PriorCalculator
from SpaceTracer.steps.step4_genotyping import GenotypingStep
from SpaceTracer.steps.step5_spatial_feature import SpatialFeatureStep
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
        resume: bool = False,
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
            'spatial_feature': SpatialFeatureStep
        }

        # cell number?
        cell_num_str=self.config.get('steps').get('cell_number')
        if os.path.exists(cell_num_str):
            cell_num=cell_num_str
        else:
            cell_num = int(cell_num_str) if cell_num_str else 0
        self.config['cell_num']=cell_num

        # cluster?
        cluster_str=self.config.get('steps').get('cluster')
        cluster_path=Path(cluster_str) if cluster_str else ""
        if cluster_path and cluster_path.exists:
            self.config['cluster']=cluster_str
        
        # keep_intermediates
        self.keep_intermediates = self.config.get('run').get('keep_intermediates',False)
        
        self.config       = config
        self.resume       = resume
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

        Genome=GenomeDetails(self.config.get('genome'),self.config.get('genome_fasta'))
        self.genome_details=Genome._get_genome_details()

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
        stop_at:    Optional[str] = None,
        only_steps: Optional[List[str]] = None,
        enable_parallel: bool = False,
        # ═══════════════════════════════════════════════════════════════
        # INTERFACE 1 — context_override
        #
        # A dict of key→value pairs injected into context *before* the
        # first step executes.  Values here win over every default.
        #
        # Use cases
        # ─────────
        # A) Re-enter pipeline at a later step using a pre-built file:
        #      context_override = {
        #          "filtered_candidates": "/data/chunk_003.mpileup"
        #      }
        #    Combine with start_from="umi_combine" to skip upstream steps.
        #
        # B) Override a step's input in an automated test:
        #      context_override = {
        #          "umi_combined": "tests/fixtures/tiny_umi.tsv"
        #      }
        #
        # C) Snakemake passes the chunk file path for each scatter job:
        #      context_override = json.loads(args.context_override)
        #
        # Path strings are automatically resolved to Path objects.
        # Other types (int, bool, list …) are passed through unchanged.
        # ═══════════════════════════════════════════════════════════════
        context_override: Optional[Dict[str, Any]] = None,
    ) -> Dict:

        logger.info("=" * 70)
        label = f"chunk {self.chunk_index:04d}" \
                if self.chunk_index is not None else "full run"
        logger.info(f"Starting Pipeline ({label})")
        logger.info("=" * 70)

        self.stats['start_time'] = time.time()

        steps_to_run = self.dag.get_execution_plan(
            start_from=start_from,
            stop_at=stop_at,
            only_steps=only_steps,
        )
        logger.info(
            f"Execution plan ({len(steps_to_run)} steps): "
            f"{' → '.join(steps_to_run)}"
        )

        # ── Build context ─────────────────────────────────────────────────────
        context = self._build_base_context()

        # ── Apply Interface 1 overrides ───────────────────────────────────────
        if context_override:
            resolved = _resolve_paths(context_override)
            context.update(resolved)
            logger.info(f"Context override keys: {sorted(resolved.keys())}")
        # ──────────────────────────────────────────────────────────────────────

        if enable_parallel:
            context = self._run_with_parallel_groups(steps_to_run, context)
        else:
            context = self._run_sequential(steps_to_run, context)

        self.stats['end_time'] = time.time()
        elapsed = self.stats['end_time'] - self.stats['start_time']

        logger.info(f"Pipeline finished in {elapsed:.2f}s ({label})")

        return {
            'final_vcf':   context.get('final_variants'),
            'elapsed_time': elapsed,
            'chunk_index': self.chunk_index,
            'stats':       self.stats,
        }

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
            resume=getattr(args, 'resume', False),
            force=getattr(args, 'force', False),
            chunk_index=getattr(args, 'chunk_index', None),  # ← Interface 2
        )
        instance._pending_context_override = override        # ← Interface 1
        return instance

    # ── Sequential execution ──────────────────────────────────────────────────

    def _run_sequential(self, steps_to_run: List[str], context: Dict) -> Dict:
        for step_name in steps_to_run:
            step = self._get_step_instance(step_name, context)
            if self.resume and not self.force:
                if self.checkpoint.is_complete(step_name):
                    logger.info(f"[{step_name}] skipped (checkpoint)")
                    context = step.load_outputs(context)
                    self.stats['step_status'][step_name] = 'skipped'
                    continue
            context = self._execute_step(step_name, step, context)
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
                if self.resume and self.checkpoint.is_complete(step_name):
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
                        if self.resume and self.checkpoint.is_complete(step_name):
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
#!/usr/bin/env python3
"""
CLI command: mutation-caller run-chunk
For Snakemake/SLURM scatter-gather integration.

Usage examples
──────────────
# Normal full run (no chunking)
mutation-caller run \
    --bam sample.bam \
    --output results/ \
    --genome hg38

# Single chunk run (called by Snakemake rule)
mutation-caller run-chunk \
    --bam      sample.bam \
    --output   results/ \
    --chunk-index 3 \
    --context-override '{"filtered_candidates": "splits/chunk_0003.mpileup"}' \
    --start-from umi_combine \
    --stop-at    hfdr
"""

import click
import json
import logging
from pathlib import Path

from SpaceTracer.pipeline.orchestrator_old import PipelineOrchestrator
from SpaceTracer.config.config_loader import load_config

logger = logging.getLogger(__name__)


# ─────────────────────────────────────────────────────────────────────────────
# run-chunk  command
# ─────────────────────────────────────────────────────────────────────────────

@click.command('run-chunk')
@click.option('--bam', required=True, type=click.Path(exists=True))
@click.option('--output', '-o', required=True, type=click.Path())
@click.option('--genome', default='hg38')
@click.option('--genome-fasta', type=click.Path(exists=True))
@click.option('--config', 'config_file', type=click.Path(exists=True))
@click.option('--start-from', default=None)
@click.option('--stop-at',    default=None)
@click.option('--only-steps', default=None,
                help='Comma-separated step names')
@click.option('--resume', is_flag=True)
@click.option('--force',  is_flag=True)
# ── Interface 2 ──────────────────────────────────────────────────────────────
@click.option(
    '--chunk-index', type=int, default=None,
    help=(
        'Index of the split chunk this job is processing. '
        'Isolates work_dir to work/chunk_XXXX/ and scopes the checkpoint. '
        'Required when running inside a Snakemake scatter rule.'
    ),
)
# ── Interface 1 ──────────────────────────────────────────────────────────────
@click.option(
    '--context-override', default=None,
    help=(
        'JSON string of key→value pairs injected into the pipeline context '
        'before execution. Lets you skip upstream steps by supplying '
        'pre-existing output files directly. '
        'Example: \'{"filtered_candidates": "splits/chunk_003.mpileup"}\''
    ),
)
# ─────────────────────────────────────────────────────────────────────────────
def run_chunk(
    bam, output, genome, genome_fasta, config_file,
    start_from, stop_at, only_steps,
    resume, force,
    chunk_index,        # Interface 2
    context_override,   # Interface 1
):
    """
    Run one chunk of the pipeline (scatter-gather mode).

    Designed to be called by a Snakemake rule for each split input file.
    Combines --chunk-index (path isolation) with --context-override
    (inject the split file as a step input) so the pipeline resumes
    from the correct step with the correct data.
    """
    output_dir = Path(output)
    output_dir.mkdir(parents=True, exist_ok=True)

    config = load_config(
        genome=genome,
        genome_fasta=genome_fasta,
        custom_config=config_file,
    )

    # ── Decode Interface 1 from CLI string ────────────────────────────────────
    override_dict = None
    if context_override:
        try:
            override_dict = json.loads(context_override)
        except json.JSONDecodeError as exc:
            raise click.BadParameter(
                f"--context-override is not valid JSON: {exc}",
                param_hint='--context-override',
            ) from exc

    only_steps_list = only_steps.split(',') if only_steps else None

    # ── Build orchestrator (Interface 2 via chunk_index) ──────────────────────
    pipeline = PipelineOrchestrator(
        bam_file=bam,
        output_dir=output_dir,
        config=config,
        resume=resume,
        force=force,
        chunk_index=chunk_index,    # ← Interface 2
    )

    pipeline.show_plan(
        start_from=start_from,
        stop_at=stop_at,
        only_steps=only_steps_list,
    )

    # ── Run (Interface 1 via context_override) ────────────────────────────────
    results = pipeline.run(
        start_from=start_from,
        stop_at=stop_at,
        only_steps=only_steps_list,
        context_override=override_dict,   # ← Interface 1
    )

    label = f"chunk {chunk_index:04d}" if chunk_index is not None else "run"
    click.echo(f"✓ {label} completed in {results['elapsed_time']:.1f}s")


# ─────────────────────────────────────────────────────────────────────────────
# Snakemake integration example  (not Python, shown as a docstring)
# ─────────────────────────────────────────────────────────────────────────────

SNAKEFILE_EXAMPLE = '''
# Snakefile
# ─────────────────────────────────────────────────────────────────────────
# Pattern:
#   1. One "split" rule  → produces N chunk files
#   2. One "run_chunk" rule  → one SLURM job per chunk
#   3. One "merge" rule  → combines all chunk results
# ─────────────────────────────────────────────────────────────────────────

configfile: "config.yaml"

BAM    = config["bam"]
GENOME = config["genome"]
N_CHUNKS = config.get("n_chunks", 50)
CHUNKS = [f"{i:04d}" for i in range(N_CHUNKS)]

# ── Rule 0: split mpileup into chunks ─────────────────────────────────────
rule split:
    input:
        bam = BAM
    output:
        chunks = expand("splits/chunk_{index}.mpileup", index=CHUNKS),
        manifest = "splits/manifest.json"
    shell:
        """
        mutation-caller run \\
            --bam {input.bam} \\
            --genome {GENOME} \\
            --output splits/ \\
            --stop-at filter_candidates
        mutation-caller split-mpileup \\
            --input  splits/work/filter_candidates/filtered_candidates.txt \\
            --n-chunks {N_CHUNKS} \\
            --output splits/
        """

# ── Rule 1: run pipeline on each chunk ────────────────────────────────────
rule run_chunk:
    input:
        chunk = "splits/chunk_{index}.mpileup",
        bam   = BAM,
    output:
        result = "results/chunk_{index}/hfdr_results.tsv"
    params:
        index = lambda wc: int(wc.index)
    resources:
        mem_mb  = 8000,
        runtime = 120    # minutes
    shell:
        """
        mutation-caller run-chunk \\
            --bam          {input.bam} \\
            --output       results/chunk_{wildcards.index}/ \\
            --genome       {GENOME} \\
            --chunk-index  {params.index} \\
            --context-override \'{{"filtered_candidates": "{input.chunk}"}}\' \\
            --start-from   umi_combine \\
            --stop-at      hfdr
        """

# ── Rule 2: merge all chunks ───────────────────────────────────────────────
rule merge:
    input:
        results = expand("results/chunk_{index}/hfdr_results.tsv",
                         index=CHUNKS)
    output:
        final = "results/final/variants.vcf"
    shell:
        """
        mutation-caller merge \\
            --input  results/ \\
            --output {output.final}
        """
'''
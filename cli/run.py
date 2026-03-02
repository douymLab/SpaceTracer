#!/usr/bin/env python3
"""
CLI command: mutation-caller run-chunk

Context override via --set key=value (repeatable).

Examples
────────
# Snakemake scatter job: inject split file, resume from umi_combine
mutation-caller run-chunk \
    --bam         sample.bam \
    --output      results/ \
    --chunk-index 3 \
    --start-from  umi_combine \
    --stop-at     hfdr \
    --set filtered_candidates=splits/chunk_0003.mpileup

# Override multiple keys
mutation-caller run-chunk \
    --bam    sample.bam \
    --output results/ \
    --set    filtered_candidates=splits/chunk_003.mpileup \
    --set    chunk_index=3 \
    --set    min_depth=10
"""

import click
import logging
from pathlib import Path

from SpaceTracer.pipeline.orchestrator import PipelineOrchestrator
from SpaceTracer.config.config_loader import load_config

logger = logging.getLogger(__name__)


# ─────────────────────────────────────────────────────────────────────────────
# Value type inference
# ─────────────────────────────────────────────────────────────────────────────

_PATH_SUFFIXES = {
    '.bam', '.txt', '.tsv', '.vcf', '.bed',
    '.mpileup', '.json', '.fa', '.fasta', '.gz',
    '.csv', '.h5', '.pickle', '.pkl',
}

def _infer_value(raw: str):
    """
    Convert a raw string value from --set into the most appropriate type.

    Rules (checked in order):
      1. "true" / "false"  → bool
      2. Pure integer       → int
      3. Pure float         → float
      4. Looks like a path  → Path  (has a known suffix, or contains / or \\)
      5. Anything else      → str   (unchanged)

    Why Path detection matters: context keys like 'filtered_candidates'
    or 'umi_combined' are Path objects internally.  Steps call .exists()
    on them, so we need real Path objects, not strings.
    """
    # bool
    if raw.lower() == 'true':
        return True
    if raw.lower() == 'false':
        return False

    # int
    try:
        return int(raw)
    except ValueError:
        pass

    # float
    try:
        return float(raw)
    except ValueError:
        pass

    # Path  (has a file extension we recognise, or contains a path separator)
    p = Path(raw)
    if p.suffix in _PATH_SUFFIXES or '/' in raw or '\\' in raw:
        return p

    # plain string
    return raw


def _parse_set_args(set_args: tuple) -> dict:
    """
    Parse a sequence of "key=value" strings into a dict.
    Raises click.BadParameter on malformed entries.
    """
    result = {}
    for entry in set_args:
        if '=' not in entry:
            raise click.BadParameter(
                f"'{entry}' — expected format: key=value",
                param_hint='--set',
            )
        key, _, raw_value = entry.partition('=')
        key = key.strip()
        if not key:
            raise click.BadParameter(
                f"Empty key in '{entry}'",
                param_hint='--set',
            )
        result[key] = _infer_value(raw_value)
    return result


# ─────────────────────────────────────────────────────────────────────────────
# run-chunk command
# ─────────────────────────────────────────────────────────────────────────────

@click.command('run-chunk')
@click.option('--bam',          required=True, type=click.Path(exists=True))
@click.option('--output', '-o', required=True, type=click.Path())
@click.option('--genome',       default='hg38')
@click.option('--genome-fasta', type=click.Path(exists=True))
@click.option('--config', 'config_file', type=click.Path(exists=True))
@click.option('--start-from',   default=None)
@click.option('--stop-at',      default=None)
@click.option('--only-steps',   default=None, help='Comma-separated step names')
@click.option('--resume',       is_flag=True)
@click.option('--force',        is_flag=True)
# ── chunk index (path isolation) ──────────────────────────────────────────────
@click.option(
    '--chunk-index', type=int, default=None,
    help='Index of the split chunk. Isolates work_dir and checkpoint.',
)
# ── context override via --set ─────────────────────────────────────────────────
@click.option(
    '--set', 'set_args',
    multiple=True,          # can be specified more than once
    metavar='KEY=VALUE',
    help=(
        'Inject or overwrite a context key before execution. '
        'Repeatable. Values are auto-typed (int/float/bool/Path/str). '
        'Example: --set filtered_candidates=splits/chunk_003.mpileup'
    ),
)
def run_chunk(
    bam, output, genome, genome_fasta, config_file,
    start_from, stop_at, only_steps,
    resume, force,
    chunk_index,
    set_args,           # tuple of "key=value" strings
):
    """
    Run one chunk of the pipeline (scatter-gather / Snakemake mode).

    Use --set to inject pre-existing files into the context so upstream
    steps can be skipped.  Combine with --start-from to re-enter the
    pipeline at any step.
    """
    output_dir = Path(output)
    output_dir.mkdir(parents=True, exist_ok=True)

    config = load_config(
        genome=genome,
        genome_fasta=genome_fasta,
        custom_config=config_file,
    )

    # Parse --set arguments into a plain dict
    context_override = _parse_set_args(set_args)

    if context_override:
        logger.info("Context override:")
        for k, v in context_override.items():
            logger.info(f"  {k} = {v!r}  ({type(v).__name__})")

    only_steps_list = only_steps.split(',') if only_steps else None

    # Build orchestrator  (chunk_index isolates paths)
    pipeline = PipelineOrchestrator(
        bam_file=bam,
        output_dir=output_dir,
        config=config,
        resume=resume,
        force=force,
        chunk_index=chunk_index,
    )

    pipeline.show_plan(
        start_from=start_from,
        stop_at=stop_at,
        only_steps=only_steps_list,
    )

    results = pipeline.run(
        start_from=start_from,
        stop_at=stop_at,
        only_steps=only_steps_list,
        context_override=context_override or None,
    )

    label = f"chunk {chunk_index:04d}" if chunk_index is not None else "run"
    click.echo(f"✓ {label} completed in {results['elapsed_time']:.1f}s")
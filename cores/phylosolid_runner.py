# SpaceTracer/cores/phylosolid_runner.py

from pathlib import Path

from SpaceTracer.installer.config import MANIFEST_PATH
from SpaceTracer.utils.phylosolid_wrapper import run_phylosolid


def run_phylosolid_step(
    input_file,
    output_dir,
    manifest_path=MANIFEST_PATH,
    extra_args=None,
    capture_output=False,
    dry_run=False
):
    """
    Business-level wrapper for one PhyloSOLID step.

    NOTE:
    Replace '--input' and '--output' with real PhyloSOLID CLI options.
    """
    input_path = Path(input_file)
    outdir = Path(output_dir)

    if not input_path.exists():
        raise FileNotFoundError("input file not found: {}".format(input_path))

    outdir.mkdir(parents=True, exist_ok=True)

    args = [
        "--input", str(input_path),
        "--output", str(outdir)
    ]

    if extra_args:
        args.extend(extra_args)

    return run_phylosolid(
        args=args,
        manifest_path=manifest_path,
        capture_output=capture_output,
        check=True,
        dry_run=dry_run
    )


def run_phylosolid_with_raw_args(
    raw_args,
    manifest_path=MANIFEST_PATH,
    capture_output=False,
    dry_run=False
):
    """
    Direct pass-through mode.
    """
    return run_phylosolid(
        args=raw_args,
        manifest_path=manifest_path,
        capture_output=capture_output,
        check=True,
        dry_run=dry_run
    )

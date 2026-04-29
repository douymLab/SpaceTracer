# not stable

import os
from pathlib import Path
import subprocess

from SpaceTracer.installer.config import MANIFEST_PATH
from SpaceTracer.utils.phylosolid_wrapper import run_phylosolid

def run_phylosolid(manifest_record, args):
    python_bin = manifest_record["python_bin"]
    source_dir = manifest_record["source_dir"]
    module = manifest_record["module"]
    
    env = os.environ.copy()
    pythonpath = env.get("PYTHONPATH", "")
    env["PYTHONPATH"] = f"{source_dir}:{pythonpath}" if pythonpath else source_dir
    
    cmd = [python_bin, "-m", module] + args
    result = subprocess.run(cmd, env=env, cwd=source_dir)
    return result


def run_phylosolid_cmd(
    sample,
    mutation_list,
    bam,
    barcode,
    outdir,
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
    sample = sample
    mutation_list = Path(mutation_list)

    if not mutation_list.exists():
        raise FileNotFoundError("input file not found: {}".format(mutation_list))

    outdir.mkdir(parents=True, exist_ok=True)

    args = [
        "--sample", str(sample),
        "--output", str(outdir),
        "--bam", str(bam) ,
        "--barcode", str(barcode),
        "--mutation-list"
        "--barcode",

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


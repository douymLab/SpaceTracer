# SpaceTracer/cores/phylosolid_runner.py

from pathlib import Path

from SpaceTracer.installer.config import MANIFEST_PATH
from SpaceTracer.utils.phylosolid_wrapper import run_phylosolid

def run_phylosolid(manifest_record, args):
    """运行 PhyloSOLID"""
    python_bin = manifest_record["python_bin"]
    source_dir = manifest_record["source_dir"]
    module = manifest_record["module"]
    
    # 构建环境变量，避免冲突
    env = os.environ.copy()
    pythonpath = env.get("PYTHONPATH", "")
    env["PYTHONPATH"] = f"{source_dir}:{pythonpath}" if pythonpath else source_dir
    
    # 使用 subprocess 运行
    cmd = [python_bin, "-m", module] + args
    result = subprocess.run(cmd, env=env, cwd=source_dir)
    return result


def run_phylosolid_cmd(
    sample,
    mutation_list,
    bam,
    barcode,
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


# def run_phylosolid_with_raw_args(
#     raw_args,
#     manifest_path=MANIFEST_PATH,
#     capture_output=False,
#     dry_run=False
# ):
#     """
#     Direct pass-through mode.
#     """
#     return run_phylosolid(
#         args=raw_args,
#         manifest_path=manifest_path,
#         capture_output=capture_output,
#         check=True,
#         dry_run=dry_run
#     )

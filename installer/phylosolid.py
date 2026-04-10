# installer/phylosolid.py

import sys
from pathlib import Path

from SpaceTracer.installer.config import PHYLOSOLID_DIR, MANIFEST_PATH
from SpaceTracer.installer.manifest import init_phylosolid_record


def register_phylosolid(
    source_dir: Path = PHYLOSOLID_DIR,
    python_bin: str = sys.executable
):
    """
    Register local PhyloSOLID repo into installer/manifest.json
    """
    source_dir = Path(source_dir).resolve()

    if not source_dir.exists():
        raise FileNotFoundError(f"PhyloSOLID source dir not found: {source_dir}")

    init_phylosolid_record(
        source_dir=str(source_dir),
        python_bin=python_bin,
        invoke_mode="python_module",
        module="cli.main",
        cwd=str(source_dir),
        found=True,
        healthcheck="registered",
        manifest_path=MANIFEST_PATH
    )

    print(f"PhyloSOLID registered into manifest: {MANIFEST_PATH}")

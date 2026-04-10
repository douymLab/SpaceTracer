# SpaceTracer/installer/main.py

import sys
from pathlib import Path

from SpaceTracer.installer.config import MANIFEST_PATH, PHYLOSOLID_DIR
from SpaceTracer.installer.manifest import (
    init_manifest_if_needed,
    init_phylosolid_record,
    read_manifest,
)


def log(msg):
    print("[SpaceTracer Installer] {}".format(msg))


def register_phylosolid():
    source_dir = Path(PHYLOSOLID_DIR).resolve()

    if not source_dir.exists():
        raise FileNotFoundError("PhyloSOLID source dir not found: {}".format(source_dir))

    init_phylosolid_record(
        source_dir=str(source_dir),
        python_bin=sys.executable,
        invoke_mode="python_module",
        module="cli.main",
        cwd=str(source_dir),
        found=True,
        healthcheck="registered",
        manifest_path=MANIFEST_PATH,
    )

    log("PhyloSOLID registered successfully.")
    log("Manifest path: {}".format(MANIFEST_PATH))


def show_manifest():
    data = read_manifest(MANIFEST_PATH)
    log("Current manifest content:")
    print(data)


def main():
    log("Initializing manifest if needed...")
    init_manifest_if_needed(MANIFEST_PATH)

    log("Registering PhyloSOLID...")
    register_phylosolid()

    show_manifest()


if __name__ == "__main__":
    main()

# SpaceTracer/installer/config.py

from pathlib import Path

# SpaceTracer package root: .../SpaceTracer
PACKAGE_ROOT = Path(__file__).resolve().parent.parent

# External tools directory
EXTERNAL_TOOLS_DIR = PACKAGE_ROOT / "external_tools"

# PhyloSOLID source repo
PHYLOSOLID_DIR = EXTERNAL_TOOLS_DIR / "PhyloSOLID"

# Manifest path
MANIFEST_PATH = PACKAGE_ROOT / "installer" / "manifest.json"

#!/usr/bin/env bash
set -euo pipefail
PHYLOSOLID_SRC="/storage/douyanmeiLab/yangzhirui/SpaceTracer_new/SpaceTracer/external_tools/PhyloSOLID"
PYTHON_BIN="/home/douyanmeiLab/yangzhirui/conda/envs/SpaceTracer/bin/python"
export PYTHONPATH="${PHYLOSOLID_SRC}:${PYTHONPATH:-}"
cd "${PHYLOSOLID_SRC}"
exec "${PYTHON_BIN}" "/storage/douyanmeiLab/yangzhirui/SpaceTracer_new/SpaceTracer/external_tools/PhyloSOLID/cli/main.py" "$@"

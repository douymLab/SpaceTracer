#!/usr/bin/env bash
set -euo pipefail

PROJECT_NAME="SpaceTracer"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

CONDA_ENV=""
SKIP_PIP=0

print_help() {
    cat <<EOF
Usage:
  bash install.sh [options]

Options:
  --conda-env NAME    Create/use this conda environment
  --skip-pip          Skip 'pip install -e .'
  -h, --help          Show this help message

Examples:
  bash install.sh
  bash install.sh --conda-env SpaceTracer
  bash install.sh --skip-pip
EOF
}

log() {
    echo "[INFO] $*"
}

warn() {
    echo "[WARN] $*" >&2
}

err() {
    echo "[ERROR] $*" >&2
    exit 1
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --conda-env)
            [[ $# -ge 2 ]] || err "--conda-env requires a name"
            CONDA_ENV="$2"
            shift 2
            ;;
        --skip-pip)
            SKIP_PIP=1
            shift
            ;;
        -h|--help)
            print_help
            exit 0
            ;;
        *)
            err "Unknown argument: $1"
            ;;
    esac
done

log "Installing ${PROJECT_NAME}"
log "Project directory: ${SCRIPT_DIR}"

# ---------------------------
# Detect Python
# ---------------------------
if command -v python >/dev/null 2>&1; then
    PYTHON_BIN="python"
elif command -v python3 >/dev/null 2>&1; then
    PYTHON_BIN="python3"
else
    err "Python not found in PATH"
fi

log "Detected Python: $(${PYTHON_BIN} --version)"

# ---------------------------
# Optional conda setup
# ---------------------------
if [[ -n "${CONDA_ENV}" ]]; then
    if ! command -v conda >/dev/null 2>&1; then
        err "conda not found, but --conda-env was provided"
    fi

    # shellcheck disable=SC1091
    eval "$(conda shell.bash hook)"

    if conda env list | awk '{print $1}' | grep -qx "${CONDA_ENV}"; then
        log "Conda environment '${CONDA_ENV}' already exists"
    else
        if [[ -f "${SCRIPT_DIR}/environment.yml" ]]; then
            log "Creating conda environment '${CONDA_ENV}' from environment.yml"
            conda env create -n "${CONDA_ENV}" -f "${SCRIPT_DIR}/environment.yml"
        elif [[ -f "${SCRIPT_DIR}/environment.yaml" ]]; then
            log "Creating conda environment '${CONDA_ENV}' from environment.yaml"
            conda env create -n "${CONDA_ENV}" -f "${SCRIPT_DIR}/environment.yaml"
        else
            err "No environment.yml or environment.yaml found in project root"
        fi
    fi

    log "Activating conda environment '${CONDA_ENV}'"
    conda activate "${CONDA_ENV}"

    if command -v python >/dev/null 2>&1; then
        PYTHON_BIN="python"
    fi

    log "Using Python in conda env: $(${PYTHON_BIN} --version)"
fi

# ---------------------------
# Install package
# ---------------------------
if [[ "${SKIP_PIP}" -eq 0 ]]; then
    if [[ -f "${SCRIPT_DIR}/pyproject.toml" || -f "${SCRIPT_DIR}/setup.py" ]]; then
        log "Upgrading pip/setuptools/wheel"
        "${PYTHON_BIN}" -m pip install --upgrade pip setuptools wheel

        log "Installing ${PROJECT_NAME} with pip install -e ."
        "${PYTHON_BIN}" -m pip install -e "${SCRIPT_DIR}"
    else
        warn "No pyproject.toml or setup.py found, skipping pip install"
    fi
else
    log "Skipping pip install"
fi

# ---------------------------
# Done
# ---------------------------
echo
echo "======================================"
echo "${PROJECT_NAME} installation completed"
echo "======================================"
echo
if [[ -n "${CONDA_ENV}" ]]; then
    echo "To use the environment later:"
    echo "  conda activate ${CONDA_ENV}"
    echo
fi
echo "Try:"
echo "  spacetracer run --help"
echo

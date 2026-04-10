#!/usr/bin/env bash
set -euo pipefail

########################################
# SpaceTracer install.sh v1
#
# Features:
#   - multi-genome support
#   - incremental manifest update
#   - PhyloSOLID managed as source tool
#   - ANNOVAR presence + db validation
#   - genome-specific resources download
########################################

############ defaults ############
PROJECT_NAME="SpaceTracer"
PROJECT_VERSION="0.1.0"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MANIFEST_PATH="${MANIFEST_PATH:-${SCRIPT_DIR}/manifest.json}"
RESOURCE_ROOT="${RESOURCE_ROOT:-${SCRIPT_DIR}/resources}"
EXTERNAL_ROOT="${EXTERNAL_ROOT:-${SCRIPT_DIR}/external}"
BIN_DIR="${BIN_DIR:-${SCRIPT_DIR}/bin}"

GENOMES="${GENOMES:-hg38}"

ANNOVAR_DIR="${ANNOVAR_DIR:-/storage/douyanmeiLab/yangzhirui/00.Softwares/annovar}"
PHYLOSOLID_SOURCE_DIR="${PHYLOSOLID_SOURCE_DIR:-${EXTERNAL_ROOT}/PhyloSOLID}"
PYTHON_BIN="${PYTHON_BIN:-$(command -v python || true)}"

FORCE_REDOWNLOAD="${FORCE_REDOWNLOAD:-false}"
FORCE_REINSTALL="${FORCE_REINSTALL:-false}"
SKIP_PHYLOSOLID="${SKIP_PHYLOSOLID:-false}"

mkdir -p "${RESOURCE_ROOT}" "${EXTERNAL_ROOT}" "${BIN_DIR}"

############ runtime vars ############
PHYLOSOLID_FOUND="false"
PHYLOSOLID_HEALTHCHECK="not_checked"
PHYLOSOLID_INVOKE_MODE=""
PHYLOSOLID_MODULE="cli.main"
PHYLOSOLID_WRAPPER="${BIN_DIR}/phylosolid_wrapper.sh"

ANNOVAR_FOUND="false"
ANNOVAR_TABLE=""
ANNOVAR_ANNOTATE=""
ANNOVAR_HUMANDB=""

############ helpers ############
log() {
    echo "[INFO] $*"
}

warn() {
    echo "[WARN] $*" >&2
}

die() {
    echo "[ERROR] $*" >&2
    exit 1
}

check_command() {
    command -v "$1" >/dev/null 2>&1 || die "Required command not found: $1"
}

parse_args() {
    while [[ $# -gt 0 ]]; do
        case "$1" in
            --genomes)
                GENOMES="$2"
                shift 2
                ;;
            --annovar-dir)
                ANNOVAR_DIR="$2"
                shift 2
                ;;
            --resource-root)
                RESOURCE_ROOT="$2"
                shift 2
                ;;
            --manifest)
                MANIFEST_PATH="$2"
                shift 2
                ;;
            --phylosolid-source)
                PHYLOSOLID_SOURCE_DIR="$2"
                shift 2
                ;;
            --python-bin)
                PYTHON_BIN="$2"
                shift 2
                ;;
            --force-redownload)
                FORCE_REDOWNLOAD="true"
                shift
                ;;
            --force-reinstall)
                FORCE_REINSTALL="true"
                shift
                ;;
            --skip-phylosolid)
                SKIP_PHYLOSOLID="true"
                shift
                ;;
            -h|--help)
                cat <<EOF
Usage: bash install.sh [options]

Options:
  --genomes hg38,hg19
  --annovar-dir /path/to/annovar
  --resource-root /path/to/resources
  --manifest /path/to/manifest.json
  --phylosolid-source /path/to/PhyloSOLID
  --python-bin /path/to/python
  --force-redownload
  --force-reinstall
  --skip-phylosolid
EOF
                exit 0
                ;;
            *)
                die "Unknown argument: $1"
                ;;
        esac
    done
}

############ genome-specific config ############
get_annovar_table_protocols_for_genome() {
    local genome="$1"
    case "$genome" in
        hg38)
            echo "dbnsfp30a,cosmic70"
            ;;
        hg19)
            echo "dbnsfp30a,cosmic70"
            ;;
        *)
            die "Unsupported genome build for ANNOVAR protocols: $genome"
            ;;
    esac
}

get_annovar_annotate_dbtypes_for_genome() {
    local genome="$1"
    case "$genome" in
        hg38)
            echo "wgEncodeGencodeBasicV46"
            ;;
        hg19)
            echo "wgEncodeGencodeBasicV46"
            ;;
        *)
            die "Unsupported genome build for ANNOVAR dbtypes: $genome"
            ;;
    esac
}

get_annovar_db_list_for_genome() {
    local genome="$1"
    local p d
    p="$(get_annovar_table_protocols_for_genome "$genome")"
    d="$(get_annovar_annotate_dbtypes_for_genome "$genome")"
    if [[ -n "$p" && -n "$d" ]]; then
        echo "${p},${d}"
    elif [[ -n "$p" ]]; then
        echo "${p}"
    else
        echo "${d}"
    fi
}

get_editing_bed_url_for_genome() {
    local genome="$1"
    case "$genome" in
        hg38)
            echo "https://ndownloader.figshare.com/files/12345678"
            ;;
        hg19)
            echo "https://ndownloader.figshare.com/files/22345678"
            ;;
        *)
            die "Unsupported genome build for editing bed url: $genome"
            ;;
    esac
}

get_pon_url_for_genome() {
    local genome="$1"
    case "$genome" in
        hg38)
            echo "https://ndownloader.figshare.com/files/12345679"
            ;;
        hg19)
            echo "https://ndownloader.figshare.com/files/22345679"
            ;;
        *)
            die "Unsupported genome build for PON url: $genome"
            ;;
    esac
}

get_editing_bed_name_for_genome() {
    local genome="$1"
    case "$genome" in
        hg38)
            echo "COMBINED_RADAR_REDIprotal_DARNED_hg38_all_sites.bed"
            ;;
        hg19)
            echo "COMBINED_RADAR_REDIprotal_DARNED_hg19_all_sites.bed"
            ;;
        *)
            die "Unsupported genome build for editing bed name: $genome"
            ;;
    esac
}

get_pon_name_for_genome() {
    local genome="$1"
    case "$genome" in
        hg38)
            echo "PON_new_model_hg38.txt"
            ;;
        hg19)
            echo "PON_new_model_hg19.txt"
            ;;
        *)
            die "Unsupported genome build for PON name: $genome"
            ;;
    esac
}

############ file download ############
download_file_if_needed() {
    local url="$1"
    local out="$2"

    mkdir -p "$(dirname "$out")"

    if [[ -s "$out" && "${FORCE_REDOWNLOAD}" != "true" ]]; then
        log "File exists, skip download: $out"
        return 0
    fi

    log "Downloading: $url -> $out"
    if command -v curl >/dev/null 2>&1; then
        curl -L "$url" -o "$out"
    elif command -v wget >/dev/null 2>&1; then
        wget -O "$out" "$url"
    else
        die "Neither curl nor wget is available"
    fi

    [[ -s "$out" ]] || die "Downloaded file is empty: $out"
}

############ PhyloSOLID ############
check_phylosolid_source() {
    [[ -d "${PHYLOSOLID_SOURCE_DIR}" ]] || die "PhyloSOLID source dir not found: ${PHYLOSOLID_SOURCE_DIR}"
    [[ -f "${PHYLOSOLID_SOURCE_DIR}/cli/main.py" ]] || die "PhyloSOLID cli/main.py not found: ${PHYLOSOLID_SOURCE_DIR}"
}

check_phylosolid_health() {
    check_phylosolid_source
    [[ -n "${PYTHON_BIN}" ]] || die "python not found"

    if (cd "${PHYLOSOLID_SOURCE_DIR}" && "${PYTHON_BIN}" -m cli.main -h >/dev/null 2>&1); then
        PHYLOSOLID_FOUND="true"
        PHYLOSOLID_HEALTHCHECK="ok"
        PHYLOSOLID_INVOKE_MODE="python_module"
        log "PhyloSOLID healthcheck passed via: ${PYTHON_BIN} -m cli.main"
    else
        PHYLOSOLID_FOUND="false"
        PHYLOSOLID_HEALTHCHECK="failed"
        warn "PhyloSOLID healthcheck failed in source dir: ${PHYLOSOLID_SOURCE_DIR}"
        return 1
    fi
}

create_phylosolid_wrapper() {
    cat > "${PHYLOSOLID_WRAPPER}" <<EOF
#!/usr/bin/env bash
set -euo pipefail
PHYLOSOLID_SRC="${PHYLOSOLID_SOURCE_DIR}"
PYTHON_BIN="${PYTHON_BIN}"
export PYTHONPATH="\${PHYLOSOLID_SRC}:\${PYTHONPATH:-}"
cd "\${PHYLOSOLID_SRC}"
exec "\${PYTHON_BIN}" -m cli.main "\$@"
EOF
    chmod +x "${PHYLOSOLID_WRAPPER}"
    log "Created wrapper: ${PHYLOSOLID_WRAPPER}"
}

############ ANNOVAR ############
check_annovar() {
    check_command perl

    [[ -d "${ANNOVAR_DIR}" ]] || die "ANNOVAR_DIR does not exist: ${ANNOVAR_DIR}"

    ANNOVAR_TABLE="${ANNOVAR_DIR}/table_annovar.pl"
    ANNOVAR_ANNOTATE="${ANNOVAR_DIR}/annotate_variation.pl"
    ANNOVAR_HUMANDB="${ANNOVAR_DIR}/humandb"

    [[ -f "${ANNOVAR_TABLE}" ]] || die "table_annovar.pl not found in ${ANNOVAR_DIR}"
    [[ -f "${ANNOVAR_ANNOTATE}" ]] || die "annotate_variation.pl not found in ${ANNOVAR_DIR}"
    [[ -d "${ANNOVAR_HUMANDB}" ]] || die "humandb directory not found in ${ANNOVAR_DIR}"

    ANNOVAR_FOUND="true"
    log "ANNOVAR found: ${ANNOVAR_DIR}"
}

get_annovar_db_status_json_for_genome() {
    local genome="$1"
    local db_list="$2"

    python - <<PY
import os, json, glob
humandb = "${ANNOVAR_HUMANDB}"
build = "${genome}"
dbs = [x.strip() for x in "${db_list}".split(",") if x.strip()]

status = {}
for db in dbs:
    pattern = os.path.join(humandb, f"{build}_{db}*")
    files = sorted(glob.glob(pattern))
    status[db] = {
        "present": len(files) > 0,
        "matched_files": files
    }
print(json.dumps(status))
PY
}

assert_required_annovar_databases_for_genome() {
    local genome="$1"
    local db_status_json="$2"

    python - <<PY
import json, sys
genome = "${genome}"
status = json.loads('''${db_status_json}''')
missing = [db for db, info in status.items() if not info["present"]]
if missing:
    print(f"[ERROR] Missing required ANNOVAR databases for {genome}:")
    for db in missing:
        print(f"  - {db}")
    sys.exit(1)
else:
    print(f"[INFO] All required ANNOVAR databases are present for {genome}.")
PY
}

############ manifest ############
init_manifest_if_needed() {
    if [[ ! -f "${MANIFEST_PATH}" ]]; then
        log "Initializing manifest: ${MANIFEST_PATH}"
        python - <<PY
import json
manifest = {
    "project": {
        "name": "${PROJECT_NAME}",
        "version": "${PROJECT_VERSION}"
    },
    "tools": {},
    "genomes": {},
    "downloads": []
}
with open("${MANIFEST_PATH}", "w") as f:
    json.dump(manifest, f, indent=2)
PY
    fi
}

update_manifest_tools() {
    python - <<PY
import json
from datetime import datetime

path = "${MANIFEST_PATH}"
with open(path) as f:
    manifest = json.load(f)

manifest.setdefault("project", {})
manifest["project"]["name"] = "${PROJECT_NAME}"
manifest["project"]["version"] = "${PROJECT_VERSION}"
manifest["project"]["install_time"] = datetime.now().isoformat(timespec="seconds")

manifest.setdefault("tools", {})
manifest["tools"]["phylosolid"] = {
    "found": ${PHYLOSOLID_FOUND},
    "source_dir": "${PHYLOSOLID_SOURCE_DIR}",
    "python_bin": "${PYTHON_BIN}",
    "invoke_mode": "${PHYLOSOLID_INVOKE_MODE}",
    "module": "${PHYLOSOLID_MODULE}",
    "cwd": "${PHYLOSOLID_SOURCE_DIR}",
    "wrapper": "${PHYLOSOLID_WRAPPER}",
    "healthcheck": "${PHYLOSOLID_HEALTHCHECK}"
}

manifest["tools"]["annovar"] = {
    "found": ${ANNOVAR_FOUND},
    "script_dir": "${ANNOVAR_DIR}",
    "table_annovar": "${ANNOVAR_TABLE}",
    "annotate_variation": "${ANNOVAR_ANNOTATE}",
    "humandb": "${ANNOVAR_HUMANDB}"
}

with open(path, "w") as f:
    json.dump(manifest, f, indent=2)
PY
}

update_manifest_genome() {
    local genome="$1"
    local resource_dir="$2"
    local editing_bed="$3"
    local pon_file="$4"
    local table_protocols="$5"
    local annotate_dbtypes="$6"
    local db_status_json="$7"
    local editing_url="$8"
    local pon_url="$9"

    python - <<PY
import json

path = "${MANIFEST_PATH}"
with open(path) as f:
    manifest = json.load(f)

manifest.setdefault("genomes", {})
manifest["genomes"]["${genome}"] = {
    "annovar": {
        "build": "${genome}",
        "table_annovar": {
            "protocols": [x.strip() for x in "${table_protocols}".split(",") if x.strip()]
        },
        "annotate_variation": {
            "dbtypes": [x.strip() for x in "${annotate_dbtypes}".split(",") if x.strip()]
        },
        "databases": json.loads('''${db_status_json}''')
    },
    "resources": {
        "resource_dir": "${resource_dir}",
        "editing_bed": "${editing_bed}",
        "PON_file": "${pon_file}"
    }
}

manifest.setdefault("downloads", [])

def upsert_download(name, genome, url, path):
    for item in manifest["downloads"]:
        if item.get("name") == name and item.get("genome") == genome:
            item["url"] = url
            item["path"] = path
            item["downloaded"] = True
            return
    manifest["downloads"].append({
        "name": name,
        "genome": genome,
        "url": url,
        "path": path,
        "downloaded": True,
        "checksum": ""
    })

upsert_download("editing_bed", "${genome}", "${editing_url}", "${editing_bed}")
upsert_download("PON_file", "${genome}", "${pon_url}", "${pon_file}")

with open(path, "w") as f:
    json.dump(manifest, f, indent=2)
PY
}

############ per-genome install ############
process_one_genome() {
    local genome="$1"
    local genome_resource_dir="${RESOURCE_ROOT}/${genome}"

    mkdir -p "${genome_resource_dir}"

    local table_protocols
    local annotate_dbtypes
    local annovar_db_list

    table_protocols="$(get_annovar_table_protocols_for_genome "${genome}")"
    annotate_dbtypes="$(get_annovar_annotate_dbtypes_for_genome "${genome}")"
    annovar_db_list="$(get_annovar_db_list_for_genome "${genome}")"

    local editing_url
    local pon_url
    local editing_name
    local pon_name

    editing_url="$(get_editing_bed_url_for_genome "${genome}")"
    pon_url="$(get_pon_url_for_genome "${genome}")"
    editing_name="$(get_editing_bed_name_for_genome "${genome}")"
    pon_name="$(get_pon_name_for_genome "${genome}")"

    local editing_path="${genome_resource_dir}/${editing_name}"
    local pon_path="${genome_resource_dir}/${pon_name}"

    log "Processing genome: ${genome}"

    local db_status_json
    db_status_json="$(get_annovar_db_status_json_for_genome "${genome}" "${annovar_db_list}")"

    python - <<PY
import json
status = json.loads('''${db_status_json}''')
for db, info in status.items():
    if info["present"]:
        print(f"[INFO] ANNOVAR DB found for ${genome}: {db}")
    else:
        print(f"[WARN] ANNOVAR DB missing for ${genome}: {db}")
PY

    assert_required_annovar_databases_for_genome "${genome}" "${db_status_json}"

    download_file_if_needed "${editing_url}" "${editing_path}"
    download_file_if_needed "${pon_url}" "${pon_path}"

    update_manifest_genome \
        "${genome}" \
        "${genome_resource_dir}" \
        "${editing_path}" \
        "${pon_path}" \
        "${table_protocols}" \
        "${annotate_dbtypes}" \
        "${db_status_json}" \
        "${editing_url}" \
        "${pon_url}"
}

############ main ############
main() {
    parse_args "$@"

    [[ -n "${PYTHON_BIN}" ]] || die "No python found; please use --python-bin"
    init_manifest_if_needed

    if [[ "${SKIP_PHYLOSOLID}" != "true" ]]; then
        check_phylosolid_health
        create_phylosolid_wrapper
    else
        warn "Skipping PhyloSOLID check by user request"
    fi

    check_annovar
    update_manifest_tools

    IFS=',' read -r -a genome_array <<< "${GENOMES}"
    for genome in "${genome_array[@]}"; do
        genome="$(echo "$genome" | xargs)"
        [[ -n "${genome}" ]] || continue
        process_one_genome "${genome}"
    done

    log "Installation finished successfully."
    log "Manifest written to: ${MANIFEST_PATH}"
}

main "$@"

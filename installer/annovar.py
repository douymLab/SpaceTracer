import glob
from pathlib import Path
from SpaceTracer.installer.utils import die, check_command_exists, log


def check_annovar(annovar_dir: str):
    if not check_command_exists("perl"):
        die("perl not found in PATH")

    base = Path(annovar_dir)
    if not base.exists():
        die(f"ANNOVAR dir does not exist: {annovar_dir}")

    table_annovar = base / "table_annovar.pl"
    annotate_variation = base / "annotate_variation.pl"
    humandb = base / "humandb"

    if not table_annovar.exists():
        die(f"table_annovar.pl not found: {table_annovar}")
    if not annotate_variation.exists():
        die(f"annotate_variation.pl not found: {annotate_variation}")
    if not humandb.exists():
        die(f"humandb not found: {humandb}")

    log(f"ANNOVAR found: {annovar_dir}")
    return {
        "found": True,
        "script_dir": str(base.resolve()),
        "table_annovar": str(table_annovar.resolve()),
        "annotate_variation": str(annotate_variation.resolve()),
        "humandb": str(humandb.resolve())
    }


def check_annovar_databases(humandb: str, genome: str, table_protocols: list, annotate_dbtypes: list):
    all_dbs = list(table_protocols) + list(annotate_dbtypes)
    status = {}
    for db in all_dbs:
        pattern = str(Path(humandb) / f"{genome}_{db}*")
        files = sorted(glob.glob(pattern))
        status[db] = {
            "present": len(files) > 0,
            "matched_files": files
        }
    return status


def assert_annovar_databases_present(status: dict, genome: str):
    missing = [db for db, info in status.items() if not info["present"]]
    if missing:
        details = "\n".join([f"  - {x}" for x in missing])
        die(f"Missing required ANNOVAR databases for {genome}:\n{details}")

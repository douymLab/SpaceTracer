# SpaceTracer/installer/manifest.py

import json
from pathlib import Path

from SpaceTracer.installer.config import MANIFEST_PATH


def _default_manifest():
    return {
        "tools": {}
    }


def ensure_manifest_file(manifest_path=MANIFEST_PATH):
    """
    Ensure manifest JSON exists.
    """
    manifest_path = Path(manifest_path)
    manifest_path.parent.mkdir(parents=True, exist_ok=True)

    if not manifest_path.exists():
        with manifest_path.open("w", encoding="utf-8") as f:
            json.dump(_default_manifest(), f, indent=2, ensure_ascii=False)

    return manifest_path


def init_manifest_if_needed(manifest_path=MANIFEST_PATH):
    """
    Backward-compatible alias.
    """
    return ensure_manifest_file(manifest_path)


def load_manifest(manifest_path=MANIFEST_PATH):
    """
    Load manifest JSON as dict.
    """
    manifest_path = Path(manifest_path)
    ensure_manifest_file(manifest_path)

    with manifest_path.open("r", encoding="utf-8") as f:
        return json.load(f)


def read_manifest(manifest_path=MANIFEST_PATH):
    """
    Backward-compatible alias.
    """
    return load_manifest(manifest_path)


def save_manifest(data, manifest_path=MANIFEST_PATH):
    """
    Save dict to manifest JSON.
    """
    manifest_path = Path(manifest_path)
    manifest_path.parent.mkdir(parents=True, exist_ok=True)

    with manifest_path.open("w", encoding="utf-8") as f:
        json.dump(data, f, indent=2, ensure_ascii=False)

    return manifest_path


def write_manifest(data, manifest_path=MANIFEST_PATH):
    """
    Backward-compatible alias.
    """
    return save_manifest(data, manifest_path)


def get_tool_record(tool_name, manifest_path=MANIFEST_PATH):
    """
    Get one tool record from manifest.
    """
    data = load_manifest(manifest_path)
    return data.get("tools", {}).get(tool_name)


def set_tool_record(tool_name, record, manifest_path=MANIFEST_PATH):
    """
    Set one tool record.
    """
    data = load_manifest(manifest_path)

    if "tools" not in data:
        data["tools"] = {}

    data["tools"][tool_name] = record
    return save_manifest(data, manifest_path)


def update_tool_record(tool_name, updates, manifest_path=MANIFEST_PATH):
    """
    Update one tool record partially.
    """
    data = load_manifest(manifest_path)

    if "tools" not in data:
        data["tools"] = {}

    current = data["tools"].get(tool_name, {})
    current.update(updates)
    data["tools"][tool_name] = current

    return save_manifest(data, manifest_path)


def init_phylosolid_record(
    source_dir,
    python_bin="python",
    invoke_mode="python_module",
    module="cli.main",
    cwd=None,
    found=True,
    healthcheck="unknown",
    manifest_path=MANIFEST_PATH
):
    """
    Initialize/update the phylosolid record.
    """
    if cwd is None:
        cwd = source_dir

    record = {
        "found": found,
        "source_dir": str(source_dir),
        "python_bin": str(python_bin),
        "invoke_mode": str(invoke_mode),
        "module": str(module),
        "cwd": str(cwd),
        "healthcheck": str(healthcheck)
    }

    return set_tool_record("phylosolid", record, manifest_path)


def set_phylosolid_record(
    source_dir,
    python_bin="python",
    invoke_mode="python_module",
    module="cli.main",
    cwd=None,
    found=True,
    healthcheck="unknown",
    manifest_path=MANIFEST_PATH
):
    """
    Backward-compatible alias.
    """
    return init_phylosolid_record(
        source_dir=source_dir,
        python_bin=python_bin,
        invoke_mode=invoke_mode,
        module=module,
        cwd=cwd,
        found=found,
        healthcheck=healthcheck,
        manifest_path=manifest_path
    )

# SpaceTracer/utils/phylosolid_wrapper.py

import os
import subprocess
from pathlib import Path

from SpaceTracer.installer.config import MANIFEST_PATH
from SpaceTracer.installer.manifest import get_tool_record


def get_phylosolid_tool_info(manifest_path=MANIFEST_PATH):
    """
    Read PhyloSOLID metadata from manifest.
    """
    tool = get_tool_record("phylosolid", manifest_path=manifest_path)

    if tool is None:
        raise KeyError("phylosolid record not found in installer/manifest.json")

    if not tool.get("found", False):
        raise RuntimeError("phylosolid is marked as not found in manifest")

    source_dir = tool.get("source_dir")
    if not source_dir:
        raise RuntimeError("phylosolid source_dir missing in manifest")

    if not Path(source_dir).exists():
        raise FileNotFoundError("phylosolid source_dir does not exist: {}".format(source_dir))

    return tool


def build_phylosolid_command(args, manifest_path=MANIFEST_PATH):
    """
    Build command, cwd and env for running PhyloSOLID.
    """
    tool = get_phylosolid_tool_info(manifest_path=manifest_path)

    invoke_mode = tool.get("invoke_mode", "python_module")
    source_dir = tool["source_dir"]
    cwd = tool.get("cwd", source_dir)
    python_bin = tool.get("python_bin", "python")
    module = tool.get("module", "cli.main")

    env = os.environ.copy()
    old_pythonpath = env.get("PYTHONPATH", "")
    env["PYTHONPATH"] = "{}:{}".format(source_dir, old_pythonpath) if old_pythonpath else source_dir

    if invoke_mode == "python_module":
        cmd = [python_bin, "-m", module] + list(args)

    elif invoke_mode == "console_script":
        executable = tool.get("executable")
        if not executable:
            raise RuntimeError("manifest missing 'executable' for console_script mode")
        cmd = [executable] + list(args)

    elif invoke_mode == "wrapper":
        wrapper = tool.get("wrapper")
        if not wrapper:
            raise RuntimeError("manifest missing 'wrapper' for wrapper mode")
        cmd = [wrapper] + list(args)

    else:
        raise RuntimeError("unsupported invoke_mode: {}".format(invoke_mode))

    return {
        "cmd": cmd,
        "cwd": cwd,
        "env": env
    }


def run_phylosolid(
    args,
    manifest_path=MANIFEST_PATH,
    capture_output=False,
    check=True,
    dry_run=False
):
    """
    Execute PhyloSOLID command.
    """
    payload = build_phylosolid_command(args=args, manifest_path=manifest_path)
    cmd = payload["cmd"]
    cwd = payload["cwd"]
    env = payload["env"]

    print("[PhyloSOLID CMD] {}".format(" ".join(map(str, cmd))))
    print("[PhyloSOLID CWD] {}".format(cwd))
    print("[Manifest] {}".format(manifest_path))

    if dry_run:
        return {
            "cmd": cmd,
            "cwd": cwd,
            "env": env,
            "returncode": 0,
            "stdout": "",
            "stderr": ""
        }

    result = subprocess.run(
        cmd,
        cwd=cwd,
        env=env,
        text=True,
        capture_output=capture_output,
        check=check
    )
    return result


def run_phylosolid_help(manifest_path=MANIFEST_PATH, capture_output=True):
    """
    Debug helper.
    """
    return run_phylosolid(
        args=["-h"],
        manifest_path=manifest_path,
        capture_output=capture_output,
        check=True,
        dry_run=False
    )


def print_phylosolid_debug_info(manifest_path=MANIFEST_PATH):
    """
    Print current PhyloSOLID tool info.
    """
    tool = get_phylosolid_tool_info(manifest_path=manifest_path)

    print("=== PhyloSOLID tool info ===")
    for k, v in tool.items():
        print("{}: {}".format(k, v))

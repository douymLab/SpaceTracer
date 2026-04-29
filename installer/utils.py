import json
import os
import shutil
import subprocess
import sys
from pathlib import Path
from urllib.request import urlopen


def log(msg: str):
    print(f"[INFO] {msg}")


def warn(msg: str):
    print(f"[WARN] {msg}", file=sys.stderr)


def die(msg: str, code: int = 1):
    print(f"[ERROR] {msg}", file=sys.stderr)
    raise SystemExit(code)


def ensure_dir(path):
    Path(path).mkdir(parents=True, exist_ok=True)


def check_command_exists(cmd: str) -> bool:
    return shutil.which(cmd) is not None


def run_command(cmd, cwd=None, env=None, capture_output=False, check=True, text=True):
    return subprocess.run(
        cmd,
        cwd=cwd,
        env=env,
        capture_output=capture_output,
        check=check,
        text=text
    )


def read_json(path):
    with open(path) as f:
        return json.load(f)


def write_json(path, data):
    with open(path, "w") as f:
        json.dump(data, f, indent=2)


def download_file(url: str, output_path: str, force: bool = False):
    output = Path(output_path)
    ensure_dir(output.parent)

    if output.exists() and output.stat().st_size > 0 and not force:
        log(f"File exists, skip download: {output}")
        return str(output)

    log(f"Downloading: {url} -> {output}")
    with urlopen(url) as resp, open(output, "wb") as out:
        out.write(resp.read())

    if not output.exists() or output.stat().st_size == 0:
        die(f"Downloaded file is empty: {output}")

    return str(output)

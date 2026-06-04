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
    
    # 检查是否是有效的 PhyloSOLID 目录
    cli_dir = source_dir / "cli"
    if not cli_dir.exists():
        raise FileNotFoundError(f"Invalid PhyloSOLID directory: missing 'cli/' subdirectory")
    
    # 记录入口信息，运行时会动态添加路径
    init_phylosolid_record(
        source_dir=str(source_dir),
        python_bin=python_bin,
        invoke_mode="python_module_with_path",  # 自定义模式
        module="cli.main",
        cwd=str(source_dir),
        pythonpath=str(source_dir),  # 额外记录需要添加到 PYTHONPATH 的路径
        found=True,
        healthcheck="registered",
        manifest_path=MANIFEST_PATH
    )

    print(f"PhyloSOLID registered into manifest: {MANIFEST_PATH}")
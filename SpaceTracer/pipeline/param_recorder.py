#!/usr/bin/env python3
"""
Record resolved pipeline / step parameters per run.

Written next to checkpoints as: {output_dir}/.pipeline_params.json

Why a file (not log only):
- easy to diff after re-runs with changed params
- keeps history of previous runs
- readable JSON, same style as .pipeline_checkpoints.json
"""

from __future__ import annotations

import json
import logging
import os
import shutil
import tempfile
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional

logger = logging.getLogger(__name__)

# Keys that are huge / not useful for param review.
_SKIP_TOP_LEVEL = {
    "genome_details",
    "input_details",
    "resource_details",
}

_KEEP_TOP_LEVEL = (
    "sample",
    "condition",
    "genome",
    "sequence_type",
    "bin_size",
    "output_dir",
    "spaceranger_dir",
    "resource_dir",
    "regions_file",
    "bam_file",
    "tissue_position",
    "barcode_key",
    "barcode_mapping",
    "cell_info",
    "cell_number",
    "cluster_file",
    "cluster_method",
    "model_name",
    "model_dir",
    "contain_chrM",
)


def _now() -> str:
    return datetime.now().isoformat(timespec="seconds")


def _run_id(existing_ids: Optional[set] = None) -> str:
    """Unique run id; include microseconds so same-second re-runs do not collide."""
    existing_ids = existing_ids or set()
    base = datetime.now().strftime("%Y%m%d_%H%M%S_%f")
    run_id = base
    n = 1
    while run_id in existing_ids:
        run_id = f"{base}_{n}"
        n += 1
    return run_id


def to_jsonable(obj: Any) -> Any:
    if isinstance(obj, dict):
        return {str(k): to_jsonable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple, set)):
        return [to_jsonable(v) for v in obj]
    if isinstance(obj, Path):
        return str(obj)
    if isinstance(obj, (str, int, float, bool)) or obj is None:
        return obj
    return str(obj)


def snapshot_config(config: Dict[str, Any]) -> Dict[str, Any]:
    """Compact, human-readable snapshot of the resolved config."""
    snap: Dict[str, Any] = {}
    for key in _KEEP_TOP_LEVEL:
        if key in config and config[key] is not None:
            snap[key] = to_jsonable(config[key])

    run_cfg = config.get("run")
    if isinstance(run_cfg, dict):
        snap["run"] = to_jsonable(run_cfg)

    logging_cfg = config.get("logging")
    if isinstance(logging_cfg, dict):
        snap["logging"] = to_jsonable(logging_cfg)

    steps = config.get("steps")
    if isinstance(steps, dict):
        snap["steps"] = to_jsonable(steps)

    # Keep any other simple top-level scalars the user may have added.
    for key, value in config.items():
        if key in snap or key in _SKIP_TOP_LEVEL:
            continue
        if key in ("run", "logging", "steps"):
            continue
        if isinstance(value, (str, int, float, bool)) or value is None:
            snap[key] = value

    return snap


class ParamRecorder:
    """Append-friendly parameter history for SpaceTracer runs."""

    def __init__(self, output_dir: Path, *, max_history: int = 50):
        self.output_dir = Path(output_dir)
        self.params_file = self.output_dir / ".pipeline_params.json"
        self.backup_file = self.output_dir / ".pipeline_params.json.bak"
        self.max_history = max_history
        self.state = self._load()
        self.current_run_id: Optional[str] = None

    def _default_state(self) -> Dict[str, Any]:
        return {
            "meta": {
                "format_version": 1,
                "created_at": _now(),
                "updated_at": _now(),
                "note": (
                    "current = this run; history = previous runs "
                    "(kept because params may change between re-runs)."
                ),
            },
            "current": None,
            "history": [],
        }

    def _load(self) -> Dict[str, Any]:
        if not self.params_file.exists():
            return self._default_state()
        try:
            with open(self.params_file, "r", encoding="utf-8") as f:
                raw = f.read().strip()
            if not raw:
                return self._default_state()
            data = json.loads(raw)
            if not isinstance(data, dict) or "meta" not in data:
                raise ValueError("invalid params file structure")
            data.setdefault("current", None)
            data.setdefault("history", [])
            return data
        except Exception as exc:
            logger.warning(f"Failed to load {self.params_file}: {exc}")
            if self.backup_file.exists():
                try:
                    with open(self.backup_file, "r", encoding="utf-8") as f:
                        data = json.load(f)
                    logger.info("Recovered params history from backup")
                    return data
                except Exception as exc2:
                    logger.warning(f"Failed to load params backup: {exc2}")
            return self._default_state()

    def _save(self) -> None:
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.state["meta"]["updated_at"] = _now()

        fd, tmp_path = tempfile.mkstemp(
            dir=self.output_dir,
            prefix=".params_",
            suffix=".tmp",
        )
        try:
            with os.fdopen(fd, "w", encoding="utf-8") as f:
                json.dump(self.state, f, indent=2, ensure_ascii=False)
                f.write("\n")

            if self.params_file.exists():
                shutil.copy2(self.params_file, self.backup_file)

            shutil.move(tmp_path, self.params_file)
        except Exception:
            Path(tmp_path).unlink(missing_ok=True)
            raise

    def start_run(self, config: Dict[str, Any]) -> str:
        """Start a new run record. Previous current run is pushed to history."""
        existing_ids = set()
        if self.state.get("current") and self.state["current"].get("run_id"):
            existing_ids.add(self.state["current"]["run_id"])
        for old in self.state.get("history") or []:
            if isinstance(old, dict) and old.get("run_id"):
                existing_ids.add(old["run_id"])

        run_id = _run_id(existing_ids)
        self.current_run_id = run_id

        if self.state.get("current"):
            history: List[Dict[str, Any]] = list(self.state.get("history") or [])
            history.insert(0, self.state["current"])
            self.state["history"] = history[: self.max_history]

        self.state["current"] = {
            "run_id": run_id,
            "started_at": _now(),
            "config": snapshot_config(config),
            "steps": {},
        }
        self._save()
        logger.info(f"Run params recording started -> {self.params_file} (run_id={run_id})")
        return run_id

    def record_step(
        self,
        step_name: str,
        step_params: Optional[Dict[str, Any]] = None,
        *,
        extra: Optional[Dict[str, Any]] = None,
    ) -> None:
        """Record parameters used by one step in the current run."""
        if not self.state.get("current"):
            # Safety: allow recording even if start_run was missed.
            self.start_run({})

        entry: Dict[str, Any] = {
            "recorded_at": _now(),
            "params": to_jsonable(step_params or {}),
        }
        if extra:
            entry["extra"] = to_jsonable(extra)

        self.state["current"]["steps"][step_name] = entry
        self._save()
        logger.info(f"[{step_name}] params recorded in {self.params_file.name}")

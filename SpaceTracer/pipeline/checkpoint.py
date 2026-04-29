#!/usr/bin/env python3
"""
Checkpoint Manager
"""

import json
import logging
import os
import shutil
import tempfile
from pathlib import Path
from datetime import datetime
from typing import Dict, Any, Optional, List

logger = logging.getLogger(__name__)


class CheckpointManager:
   
    def __init__(self, output_dir: Path, disabled: bool = False):
        self.output_dir = Path(output_dir)
        self.checkpoint_file = self.output_dir / ".pipeline_checkpoints.json"
        self.backup_file = self.output_dir / ".pipeline_checkpoints.json.bak"
        self.disabled = disabled

        self.state = self._load_checkpoints()

    # ─────────────────────────────────────────────────────────────
    # basic helpers
    # ─────────────────────────────────────────────────────────────

    def _now(self) -> str:
        return datetime.now().isoformat(timespec="seconds")

    def _default_state(self) -> Dict[str, Any]:
        now = self._now()
        return {
            "meta": {
                "format_version": 2,
                "created_at": now,
                "updated_at": now,
                "validation": {
                    "mode": "exists"
                }
            },
            "steps": {}
        }

    def _is_valid_state(self, state: Dict[str, Any]) -> bool:
        return (
            isinstance(state, dict)
            and "meta" in state
            and "steps" in state
            and isinstance(state["steps"], dict)
        )

    # ─────────────────────────────────────────────────────────────
    # load / save
    # ─────────────────────────────────────────────────────────────

    def _load_checkpoints(self) -> Dict[str, Any]:
        if self.disabled:
            logger.info("Checkpoint is disabled")
            return self._default_state()

        if not self.checkpoint_file.exists():
            logger.info("No checkpoint file found, starting fresh")
            return self._default_state()

        try:
            state = self._load_json_file(self.checkpoint_file)
            if not self._is_valid_state(state):
                raise ValueError("Invalid checkpoint structure")
            logger.info(
                f"Loaded {len(state.get('steps', {}))} checkpoints from {self.checkpoint_file}"
            )
            return state

        except Exception as e:
            logger.warning(f"Failed to load checkpoint file: {e}")

            if self.backup_file.exists():
                try:
                    logger.info("Trying backup checkpoint...")
                    state = self._load_json_file(self.backup_file)
                    if self._is_valid_state(state):
                        self._save_state(state)
                        logger.info("Recovered checkpoint from backup")
                        return state
                except Exception as e2:
                    logger.error(f"Failed to load backup checkpoint: {e2}")

            logger.warning("Starting with empty checkpoint state")
            return self._default_state()

    def _load_json_file(self, path: Path) -> Dict[str, Any]:
        with open(path, "r") as f:
            content = f.read().strip()

        if not content:
            return self._default_state()

        data = json.loads(content)
        return data

    def _save_state(self, state: Dict[str, Any]) -> None:
        if self.disabled:
            return

        self.checkpoint_file.parent.mkdir(parents=True, exist_ok=True)

        fd, tmp_path = tempfile.mkstemp(
            dir=self.checkpoint_file.parent,
            prefix=".checkpoint_",
            suffix=".tmp"
        )

        try:
            with os.fdopen(fd, "w") as f:
                json.dump(state, f, indent=2)

            if self.checkpoint_file.exists():
                shutil.copy2(self.checkpoint_file, self.backup_file)

            shutil.move(tmp_path, self.checkpoint_file)
            logger.debug(f"Saved checkpoint to {self.checkpoint_file}")

        except Exception:
            Path(tmp_path).unlink(missing_ok=True)
            raise

    def _save_checkpoints(self) -> None:
        if self.disabled:
            return
        self.state["meta"]["updated_at"] = self._now()
        self._save_state(self.state)

    # ─────────────────────────────────────────────────────────────
    # file info / validation
    # ─────────────────────────────────────────────────────────────

    def _build_file_info(self, file_path: str) -> Dict[str, Any]:
        path = Path(file_path)
        info = {
            "path": str(path),
            "real_path": str(path.resolve()) if path.exists() else str(path),
        }

        if path.exists():
            stat = path.stat()
            info["mtime"] = stat.st_mtime
            info["size"] = stat.st_size

        return info


    def verify_output_file(self, file_info: Dict[str, Any]) -> bool:
        try:
            path = file_info.get("path")
            if not path:
                return False
            return Path(path).exists()
        except:
            return False

    # ─────────────────────────────────────────────────────────────
    # record access
    # ─────────────────────────────────────────────────────────────

    def get_step_record(self, step_name: str) -> Dict[str, Any]:
        return self.state["steps"].get(step_name, {})

    def get_step_status(self, step_name: str) -> str:
        if self.disabled:
            return "disabled"
        return self.get_step_record(step_name).get("status", "not_started")

    def get_error(self, step_name: str) -> Optional[str]:
        return self.get_step_record(step_name).get("error")

    def get_outputs(self, step_name: str) -> Dict[str, str]:
        record = self.get_step_record(step_name)
        if record.get("status") != "complete":
            return {}

        outputs = record.get("outputs", {})
        result = {}
        for key, file_info in outputs.items():
            if self.verify_output_file(file_info):
                result[key] = file_info["path"]
            elif isinstance(file_info,int):
                result[key] = file_info
        return result

    # ─────────────────────────────────────────────────────────────
    # mark status
    # ─────────────────────────────────────────────────────────────

    def mark_complete(self, step_name: str, new_outputs: Dict[str, Any]) -> None:
        if self.disabled:
            return
        logger.info(
            f"Marking step '{step_name}' as complete, outputs={list(new_outputs)}" #.keys()
        )

        outputs = {}
        for key, value in new_outputs.items():
            if value is None:
                continue
            if isinstance(value,str):
                outputs[key] = self._build_file_info(str(value))
            elif isinstance(value,int):
                outputs[key] = {"path": value}
            else:
                raise ValueError(f'Wrong output file of {key}: {value} in step[{step_name}]!')

        self.state["steps"][step_name] = {
            "status": "complete",
            "timestamp": self._now(),
            "outputs": outputs,
        }
        self._save_checkpoints()

    def mark_failed(self, step_name: str, error_message: str) -> None:
        if self.disabled:
            return

        logger.error(f"Marking step '{step_name}' as failed: {error_message}")
        self.state["steps"][step_name] = {
            "status": "failed",
            "timestamp": self._now(),
            "error": error_message,
        }
        self._save_checkpoints()

    # ─────────────────────────────────────────────────────────────
    # state checks
    # ─────────────────────────────────────────────────────────────

    def check_outputs_exist(self, step_name: str) -> bool:
        if self.disabled:
            return False

        record = self.get_step_record(step_name)
        if record.get("status") != "complete":
            return False

        outputs = record.get("outputs", {})
        if not outputs:
            logger.debug(f"No outputs recorded for {step_name}")
            return False

        for key, file_info in outputs.items():
            if not isinstance(file_info,int) and not self.verify_output_file(file_info):
                logger.warning(
                    f"Output file missing for {step_name}:{key} -> {file_info}"
                )
                return False

        return True

    def is_complete(self, step_name: str) -> bool:
        if self.disabled:
            return False
        return self.check_outputs_exist(step_name)

    def is_failed(self, step_name: str) -> bool:
        if self.disabled:
            return False
        return self.get_step_status(step_name) == "failed"

    # ─────────────────────────────────────────────────────────────
    # restore outputs
    # ─────────────────────────────────────────────────────────────

    def load_outputs_to_context(self, step_name: str, context: Dict[str, Any]) -> Dict[str, Any]:
        record = self.get_step_record(step_name)
        if record.get("status") != "complete":
            return context

        outputs = record.get("outputs", {})
        for key, file_info in outputs.items():
            if self.verify_output_file(file_info):
                context[key] = file_info["path"]
            elif isinstance(file_info,int):
                context[key] = file_info
            else:
                logger.warning(f"Output file missing: {file_info.get('path')}")

        return context

    def load_all_completed_outputs(self, context: Dict[str, Any]) -> Dict[str, Any]:
        for step_name in self.get_completed_steps():
            context = self.load_outputs_to_context(step_name, context)
        return context

    # ─────────────────────────────────────────────────────────────
    # summary / maintenance
    # ─────────────────────────────────────────────────────────────

    def get_completed_steps(self) -> List[str]:
        return [
            name for name, cp in self.state.get("steps", {}).items()
            if cp.get("status") == "complete"
        ]

    def get_failed_steps(self) -> List[str]:
        return [
            name for name, cp in self.state.get("steps", {}).items()
            if cp.get("status") == "failed"
        ]

    def clear(self) -> None:
        self.state = self._default_state()

        if self.checkpoint_file.exists():
            if self.backup_file.exists():
                self.backup_file.unlink()
            shutil.move(self.checkpoint_file, self.backup_file)

        logger.info("Cleared all checkpoints (old file moved to backup)")

    def clear_step(self, step_name: str) -> None:
        if step_name in self.state["steps"]:
            del self.state["steps"][step_name]
            self._save_checkpoints()
            logger.info(f"Cleared checkpoint for {step_name}")

    def get_all_checkpoints(self) -> Dict[str, Any]:
        return self.state.copy()

    def print_summary(self) -> None:
        logger.info("=" * 50)
        logger.info("Checkpoint Summary")
        logger.info("=" * 50)

        completed = self.get_completed_steps()
        failed = self.get_failed_steps()

        logger.info(f"Total checkpoints: {len(self.state.get('steps', {}))}")
        logger.info(f"Completed: {len(completed)}")
        logger.info(f"Failed: {len(failed)}")

        if completed:
            logger.info("Completed steps:")
            for step in completed:
                timestamp = self.state["steps"][step].get("timestamp", "unknown")
                logger.info(f"  ✓ {step} ({timestamp})")

        if failed:
            logger.info("Failed steps:")
            for step in failed:
                error = self.state["steps"][step].get("error", "unknown error")
                logger.info(f"  ✗ {step}: {error}")

        logger.info("=" * 50)

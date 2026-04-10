#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import sys
from pathlib import Path

from SpaceTracer.pipeline.orchestrator import PipelineOrchestrator
from SpaceTracer.config.config_loader import LoadConfig
from SpaceTracer.utils.logger import setup_logger


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="spacetracer run",
        description="Run SpaceTracer pipeline",
    )

    parser.add_argument(
        "--config",
        required=True,
        help="Path to config yaml file",
    )

    parser.add_argument(
        "--start-from",
        dest="start_from",
        default=None,
        help="Start running from this step",
    )

    parser.add_argument(
        "--stop-at",
        dest="stop_at",
        default=None,
        help="Stop running at this step",
    )

    parser.add_argument(
        "--force",
        action="store_true",
        help="Force rerun even if checkpoint/output exists",
    )

    return parser


def load_pipeline(config_path: str, force: bool) -> PipelineOrchestrator:
    config = LoadConfig().load_config(custom_config=config_path)
    pipeline = PipelineOrchestrator(
        config=config,
        force=force,
    )
    return pipeline


def run_command(args: argparse.Namespace):
    config_path = Path(args.config)
    if not config_path.exists():
        raise FileNotFoundError(f"Config file not found: {config_path}")

    setup_logger("INFO")

    pipeline = load_pipeline(
        config_path=str(config_path),
        force=args.force,
    )

    results = pipeline.run(
        start_from=args.start_from,
        stop_at=args.stop_at,
    )

    print("✅ Pipeline completed!")
    return results


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    run_command(args)

    # try:
    #     run_command(args)
    # except Exception as e:
    #     print(f"[SpaceTracer run] ERROR: {e}", file=sys.stderr)
    #     sys.exit(1)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
from pathlib import Path

from SpaceTracer.pipeline.orchestrator import PipelineOrchestrator
from SpaceTracer.config.config_loader import LoadConfig
from SpaceTracer.utils.logger import setup_logger


def parse_only_steps(only_steps_str: str):
    if not only_steps_str:
        return None
    steps = [step.strip() for step in only_steps_str.split(",") if step.strip()]
    return steps if steps else None


def build_run_parser(subparsers):
    parser = subparsers.add_parser(
        "run",
        help="Run SpaceTracer pipeline",
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
        "--only-steps",
        dest="only_steps",
        default=None,
        help="Comma-separated list of step names to run",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Force rerun even if checkpoint/output exists",
    )

    parser.set_defaults(func=run_command)


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

    only_steps = parse_only_steps(args.only_steps)

    if only_steps:
        pipeline.run(only_steps=only_steps)
    else:
        pipeline.run(
            start_from=args.start_from,
            stop_at=args.stop_at,
        )

    print("✅ Pipeline completed!")


def build_parser():
    parser = argparse.ArgumentParser(
        prog="spacetracer",
        description="SpaceTracer command line interface",
    )
    subparsers = parser.add_subparsers(dest="command")

    build_run_parser(subparsers)

    return parser


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)

    if not hasattr(args, "func"):
        parser.print_help()
        return 1

    return args.func(args)


if __name__ == "__main__":
    main()

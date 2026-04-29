#!/usr/bin/env python3

import argparse
import yaml

from SpaceTracer.pipeline.snakemake_generator import SnakemakeGenerator

def main():
    parser = argparse.ArgumentParser(description="Generate Snakefile from SpaceTracer DAG")
    parser.add_argument("--config", required=True, help="Path to config.yaml")
    parser.add_argument("--output", default="Snakefile", help="Output Snakefile path")
    parser.add_argument(
        "--runner-mode",
        default="orchestrator",
        choices=["orchestrator", "step"],
        help="How Snakemake rules run step commands",
    )
    parser.add_argument(
        "--final-step",
        default=None,
        help="Optional final step for rule all",
    )
    args = parser.parse_args()

    with open(args.config, "r") as f:
        config = yaml.safe_load(f)

    generator = SnakemakeGenerator(
        config=config,
        snakefile_path=args.output,
        runner_mode=args.runner_mode,
        final_step=args.final_step,
    )
    path = generator.write()
    print(f"[OK] Snakefile generated: {path}")


if __name__ == "__main__":
    main()

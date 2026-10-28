#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import math
from pathlib import Path

from SpaceTracer.config.config_loader import LoadConfig
from SpaceTracer.utils.logger import setup_logger


def positive_float(value):
    """Parse a positive floating-point CLI value."""
    parsed = float(value)
    if not math.isfinite(parsed) or parsed <= 0:
        raise argparse.ArgumentTypeError("value must be a finite number greater than 0")
    return parsed


def build_plot_parser(subparsers):
    parser = subparsers.add_parser(
        "plot",
        help="Plot SpaceTracer results",
        description="Plot SpaceTracer results",
    )
    plot_subparsers = parser.add_subparsers(dest="plot_command")

    scatter = plot_subparsers.add_parser(
        "scatter",
        help="Posterior scatter plots for PASS mutations",
        description=(
            "Plot mutation probability and mosaic likelihood "
            "from spatial_feature barcode files."
        ),
    )
    scatter.add_argument(
        "--config",
        required=True,
        help="Path to the config yaml used for spacetracer run",
    )
    scatter.add_argument(
        "--fig-size",
        dest="fig_size",
        type=positive_float,
        default=5,
        help="Scatter figure size (default: 5)",
    )
    scatter.add_argument(
        "--point-size",
        dest="point_size",
        type=positive_float,
        default=None,
        help="Point size. Default is chosen from the tissue span.",
    )
    scatter.set_defaults(func=scatter_command)

    parser.set_defaults(func=plot_help, plot_parser=parser)
    return parser


def plot_help(args):
    args.plot_parser.print_help()
    return 0


def scatter_command(args):
    config_path = Path(args.config)
    if not config_path.exists():
        raise FileNotFoundError(f"Config file not found: {config_path}")

    setup_logger("INFO")
    config = LoadConfig().load_config(custom_config=str(config_path))

    from SpaceTracer.visualization.scatter_posterior import plot_posterior_scatter

    result = plot_posterior_scatter(
        config,
        fig_size=args.fig_size,
        point_size=args.point_size,
    )
    print(
        f"Posterior scatter plots: {result['plotted']} mutations, "
        f"{len(result['missing'])} without spatial info."
    )
    for out_dir in result["out_dirs"]:
        print(f"  {out_dir}")
    return 0

"""Posterior scatter plots for SpaceTracer PASS mutations."""

import math
import os
from pathlib import Path

import matplotlib
matplotlib.use("Agg")

import numpy as np
import pandas as pd

from SpaceTracer.utils.logger import get_logger
from SpaceTracer.visualization.scatter_plot import plot_scatter_spatial

logger = get_logger(__name__)

MUTINFO_COLUMNS = [
    "barcode_name",
    "cluster",
    "in_tissue",
    "pos_x",
    "pos_y",
    "depth",
    "vaf",
    "mut_prob",
    "l_mosaic",
]
MUTINFO_SUFFIX = ".barcode.mutinfo.txt"
TISSUE_COLUMNS = [
    "barcode",
    "in_tissue",
    "array_row",
    "array_col",
    "pxl_row_in_fullres",
    "pxl_col_in_fullres",
]


def pass_mutation_list_path(config, model_name):
    sample = config["sample"]
    return (
        Path(config["output_dir"])
        / "mutation_prediction"
        / "results"
        / f"{sample}_{model_name}_total_pred_truesites_PASS_mutation_list.txt"
    )


def spatial_feature_dir(config):
    return Path(config["output_dir"]) / "spatial_feature" / str(config["sample"])


def model_names(config):
    names = config.get("model_name")
    if isinstance(names, str):
        names = [part.strip() for part in names.replace("，", ",").split(",") if part.strip()]
    elif isinstance(names, (list, tuple)):
        names = [str(name).strip() for name in names if str(name).strip()]
    else:
        names = []
    if not names:
        raise ValueError("model_name is missing from the config, so the PASS mutation list cannot be resolved.")
    return names


def read_mutation_list(path):
    mutations = []
    with open(path, "r", encoding="utf-8") as handle:
        for line in handle:
            name = line.strip()
            if name and not name.startswith("#"):
                mutations.append(name)
    return mutations


def index_mutinfo(spatial_dir):
    """Map mutation id -> barcode mutinfo file across spatial_feature chunks."""
    found = {}
    if not spatial_dir.is_dir():
        return found
    for path in spatial_dir.glob(f"**/barcode_dir/*{MUTINFO_SUFFIX}"):
        mutation = path.name[: -len(MUTINFO_SUFFIX)]
        found[mutation] = path
    return found


def read_mutinfo(path):
    df = pd.read_csv(
        path,
        sep="\t",
        header=None,
        names=MUTINFO_COLUMNS,
        comment="#",
        usecols=range(len(MUTINFO_COLUMNS)),
    )
    if df.empty:
        return df
    df["pos_x"] = pd.to_numeric(df["pos_x"], errors="coerce")
    df["pos_y"] = pd.to_numeric(df["pos_y"], errors="coerce")
    df["vaf"] = pd.to_numeric(df["vaf"], errors="coerce").fillna(0)
    df["mut_prob"] = pd.to_numeric(df["mut_prob"], errors="coerce")
    df["l_mosaic"] = pd.to_numeric(df["l_mosaic"], errors="coerce")
    return df.dropna(subset=["pos_x", "pos_y"])


def load_tissue_positions(path):
    tissue_path = Path(path)
    if tissue_path.suffix == ".parquet":
        df = pd.read_parquet(tissue_path)
    else:
        df = pd.read_csv(tissue_path)
        if "array_row" not in df.columns or "array_col" not in df.columns:
            df = pd.read_csv(tissue_path, header=None, names=TISSUE_COLUMNS)
    missing = [col for col in ("array_row", "array_col") if col not in df.columns]
    if missing:
        raise ValueError(
            f"tissue_position file is missing columns {missing}: {tissue_path}"
        )
    df["array_row"] = pd.to_numeric(df["array_row"], errors="coerce")
    df["array_col"] = pd.to_numeric(df["array_col"], errors="coerce")
    return df.dropna(subset=["array_row", "array_col"])


def axis_limits(positions):
    """Return spatial bounds and unique spot count for plot_scatter_spatial."""
    positions = positions.drop_duplicates(subset=["array_row", "array_col"])
    if positions.empty:
        raise ValueError("No valid spatial coordinates were found for scatter plotting.")
    row_min = float(positions["array_row"].min())
    row_max = float(positions["array_row"].max())
    col_min = float(positions["array_col"].min())
    col_max = float(positions["array_col"].max())
    return row_min, row_max, col_min, col_max, len(positions)


def limits_from_mutinfo(mutinfo_paths):
    frames = []
    for path in mutinfo_paths:
        df = read_mutinfo(path)
        if not df.empty:
            frames.append(df.rename(columns={"pos_x": "array_row", "pos_y": "array_col"}))
    if not frames:
        raise FileNotFoundError(
            "No barcode coordinates were found in spatial_feature mutinfo files."
        )
    return axis_limits(pd.concat(frames, ignore_index=True))


def choose_point_size(point_size, row_min, row_max, col_min, col_max, n_spots):
    if point_size is not None:
        point_size = float(point_size)
        if not math.isfinite(point_size) or point_size <= 0:
            raise ValueError("point_size must be a finite number greater than 0.")
        return point_size
    array_area = (row_max - row_min) * (col_max - col_min)
    density = n_spots / array_area if array_area else 1
    size = round((1.5 / density) ** 2, 1) if density else 6
    return float(np.clip(size, 5, 50))


def plot_posterior_scatter(config, fig_size=5, point_size=None):
    """Plot posterior and likelihood scatters for each configured PASS list.

    Returns a dict with plotted, missing, and output directories.
    """
    if not math.isfinite(fig_size) or fig_size <= 0:
        raise ValueError("fig_size must be a finite number greater than 0.")

    spatial_dir = spatial_feature_dir(config)
    mutinfo = index_mutinfo(spatial_dir)
    if not mutinfo:
        raise FileNotFoundError(
            "spatial_feature barcode files were not found under "
            f"{spatial_dir}. Run spatial_feature before plotting."
        )

    tissue_position = config.get("tissue_position")
    if tissue_position:
        limits = axis_limits(load_tissue_positions(tissue_position))
    else:
        limits = limits_from_mutinfo(mutinfo.values())
    row_min, row_max, col_min, col_max, n_spots = limits
    point_size = choose_point_size(point_size, row_min, row_max, col_min, col_max, n_spots)

    plotted = 0
    missing = []
    out_dirs = []
    for model_name in model_names(config):
        mutation_list = pass_mutation_list_path(config, model_name)
        if not mutation_list.is_file():
            raise FileNotFoundError(
                "PASS mutation list was not found: "
                f"{mutation_list}. Run mutation_prediction before plotting."
            )
        mutations = read_mutation_list(mutation_list)
        if not mutations:
            raise FileNotFoundError(
                f"PASS mutation list is empty: {mutation_list}"
            )

        plot_dir = (
            Path(config["output_dir"])
            / "plot"
            / "scatter"
            / f"{config['sample']}_{model_name}"
        )
        os.makedirs(plot_dir, exist_ok=True)
        out_dirs.append(str(plot_dir))

        for mutation in mutations:
            mutinfo_path = mutinfo.get(mutation)
            if mutinfo_path is None:
                missing.append(mutation)
                logger.warning("%s does not have spatial info", mutation)
                continue
            df = read_mutinfo(mutinfo_path)
            if df.empty:
                missing.append(mutation)
                logger.warning("%s spatial info has no coordinates", mutation)
                continue
            plot_scatter_spatial(
                df=df,
                output_dir=str(plot_dir),
                sample_name=mutation,
                column="mut_prob",
                note="prob",
                fig_size=fig_size,
                xmin=row_min,
                xmax=row_max,
                ymin=col_min,
                ymax=col_max,
                point_size=point_size,
            )
            plot_scatter_spatial(
                df=df,
                output_dir=str(plot_dir),
                sample_name=mutation,
                column="l_mosaic",
                note="likelihood",
                fig_size=fig_size,
                xmin=row_min,
                xmax=row_max,
                ymin=col_min,
                ymax=col_max,
                point_size=point_size,
            )
            plotted += 1

    if plotted == 0:
        raise FileNotFoundError(
            "None of the PASS mutations have spatial_feature barcode files under "
            f"{spatial_dir}."
        )

    logger.info(
        "Wrote posterior scatter plots for %s mutations (%s missing spatial info).",
        plotted,
        len(missing),
    )
    return {"plotted": plotted, "missing": missing, "out_dirs": out_dirs}

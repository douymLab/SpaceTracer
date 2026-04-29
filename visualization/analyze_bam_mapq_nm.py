#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import argparse
from collections import Counter

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

import pysam
from scipy.ndimage import gaussian_filter1d


def parse_args():
    parser = argparse.ArgumentParser(
        description="Analyze MAPQ and NM distributions in BAM and plot joint scatter + marginal density."
    )
    parser.add_argument(
        "-i", "--bam",
        required=True,
        help="Input BAM file"
    )
    parser.add_argument(
        "-o", "--outdir",
        required=True,
        help="Output directory"
    )
    parser.add_argument(
        "--sample",
        default=None,
        help="Sample name; default: infer from BAM filename"
    )
    parser.add_argument(
        "--max-nm",
        type=int,
        default=20,
        help="Max NM shown in plot; larger NM will be grouped into this max value in plot/stat summary table still keeps real value in pair table"
    )
    parser.add_argument(
        "--mapq-max",
        type=int,
        default=255,
        help="Max MAPQ shown in plot"
    )
    parser.add_argument(
        "--skip-secondary",
        action="store_true",
        help="Skip secondary alignments"
    )
    parser.add_argument(
        "--skip-supplementary",
        action="store_true",
        help="Skip supplementary alignments"
    )
    parser.add_argument(
        "--skip-duplicate",
        action="store_true",
        help="Skip duplicate alignments"
    )
    parser.add_argument(
        "--skip-qcfail",
        action="store_true",
        help="Skip QC-fail alignments"
    )
    parser.add_argument(
        "--skip-unmapped",
        action="store_true",
        help="Skip unmapped reads"
    )
    return parser.parse_args()


def infer_sample_name(bam_path):
    base = os.path.basename(bam_path)
    for suffix in [".bam", ".cram", ".sam"]:
        if base.endswith(suffix):
            base = base[:-len(suffix)]
    return base


def moving_density(counts, sigma=1.2):
    arr = np.asarray(counts, dtype=float)
    if arr.sum() == 0:
        return arr
    den = gaussian_filter1d(arr, sigma=sigma, mode="nearest")
    if den.sum() > 0:
        den = den / den.sum()
    return den


def analyze_bam(bam_path, max_nm, mapq_max,
                skip_secondary=False,
                skip_supplementary=False,
                skip_duplicate=False,
                skip_qcfail=False,
                skip_unmapped=False):
    pair_counter = Counter()
    mapq_counter = Counter()
    nm_counter = Counter()

    total_reads = 0
    used_reads = 0
    missing_nm = 0

    bam = pysam.AlignmentFile(bam_path, "rb")

    for read in bam.fetch(until_eof=True):
        total_reads += 1

        if skip_unmapped and read.is_unmapped:
            continue
        if skip_secondary and read.is_secondary:
            continue
        if skip_supplementary and read.is_supplementary:
            continue
        if skip_duplicate and read.is_duplicate:
            continue
        if skip_qcfail and read.is_qcfail:
            continue

        try:
            nm = read.get_tag("nM")
        except KeyError:
            missing_nm += 1
            continue

        mapq = int(read.mapping_quality)

        used_reads += 1
        pair_counter[(mapq, nm)] += 1
        mapq_counter[mapq] += 1
        nm_counter[nm] += 1

    bam.close()

    pair_rows = []
    for (mapq, nm), count in sorted(pair_counter.items()):
        pair_rows.append({
            "MAPQ": mapq,
            "NM": nm,
            "count": count
        })
    pair_df = pd.DataFrame(pair_rows)

    mapq_rows = []
    total_mapq = sum(mapq_counter.values())
    for q in sorted(mapq_counter):
        c = mapq_counter[q]
        mapq_rows.append({
            "MAPQ": q,
            "count": c,
            "fraction": c / total_mapq if total_mapq else 0
        })
    mapq_df = pd.DataFrame(mapq_rows)

    nm_rows = []
    total_nm = sum(nm_counter.values())
    for n in sorted(nm_counter):
        c = nm_counter[n]
        nm_rows.append({
            "NM": n,
            "count": c,
            "fraction": c / total_nm if total_nm else 0
        })
    nm_df = pd.DataFrame(nm_rows)

    stats = {
        "total_reads_scanned": total_reads,
        "reads_used_with_nm": used_reads,
        "reads_missing_nm": missing_nm,
        "unique_mapq_values": len(mapq_counter),
        "unique_nm_values": len(nm_counter),
        "max_nm_observed": max(nm_counter.keys()) if nm_counter else None,
        "max_mapq_observed": max(mapq_counter.keys()) if mapq_counter else None,
    }

    return pair_df, mapq_df, nm_df, stats


def plot_joint(pair_df, mapq_df, nm_df, out_png, sample, max_nm=20, mapq_max=255):
    # plot table uses capped NM for display
    if pair_df.empty:
        fig, ax = plt.subplots(figsize=(8, 6))
        ax.text(0.5, 0.5, "No data to plot", ha="center", va="center", fontsize=14)
        ax.axis("off")
        plt.tight_layout()
        plt.savefig(out_png, dpi=300)
        plt.close()
        return

    plot_df = pair_df.copy()
    plot_df["NM_plot"] = plot_df["NM"].clip(upper=max_nm)
    plot_df["MAPQ_plot"] = plot_df["MAPQ"].clip(upper=mapq_max)

    # merge same capped bins
    plot_df = (
        plot_df.groupby(["MAPQ_plot", "NM_plot"], as_index=False)["count"]
        .sum()
        .rename(columns={"MAPQ_plot": "MAPQ", "NM_plot": "NM"})
    )

    # complete marginal arrays
    mapq_full = pd.DataFrame({"MAPQ": np.arange(0, mapq_max + 1)})
    mapq_plot = mapq_full.merge(mapq_df, on="MAPQ", how="left").fillna(0)

    nm_full = pd.DataFrame({"NM": np.arange(0, max_nm + 1)})
    nm_tmp = nm_df.copy()
    nm_tmp["NM"] = nm_tmp["NM"].clip(upper=max_nm)
    nm_tmp = nm_tmp.groupby("NM", as_index=False)["count"].sum()
    nm_tmp["fraction"] = nm_tmp["count"] / nm_tmp["count"].sum()
    nm_plot = nm_full.merge(nm_tmp, on="NM", how="left").fillna(0)

    mapq_density = moving_density(mapq_plot["count"].values, sigma=1.5)
    nm_density = moving_density(nm_plot["count"].values, sigma=1.2)

    fig = plt.figure(figsize=(12, 10))
    gs = GridSpec(
        2, 2,
        width_ratios=[4.5, 1.5],
        height_ratios=[1.5, 4.5],
        wspace=0.05,
        hspace=0.05
    )

    ax_top = fig.add_subplot(gs[0, 0])
    ax_main = fig.add_subplot(gs[1, 0])
    ax_right = fig.add_subplot(gs[1, 1], sharey=ax_main)

    # main bubble plot
    sizes = 20 + 300 * (np.log10(plot_df["count"] + 1) / np.log10(plot_df["count"].max() + 1))
    sc = ax_main.scatter(
        plot_df["MAPQ"],
        plot_df["NM"],
        s=sizes,
        c=np.log10(plot_df["count"] + 1),
        cmap="viridis",
        alpha=0.75,
        edgecolors="black",
        linewidths=0.3
    )
    cbar = plt.colorbar(sc, ax=ax_main, fraction=0.046, pad=0.04)
    cbar.set_label("log10(count + 1)")

    ax_main.set_xlabel("MAPQ")
    ax_main.set_ylabel("NM")
    ax_main.set_title(f"{sample}: MAPQ vs NM")
    ax_main.set_xlim(-1, mapq_max + 1)
    ax_main.set_ylim(-0.5, max_nm + 0.5)
    ax_main.grid(alpha=0.25, linestyle="--")

    # top marginal: MAPQ histogram + density
    ax_top.bar(mapq_plot["MAPQ"], mapq_plot["count"], color="#8ecae6", width=1.0, alpha=0.8)
    ax_top2 = ax_top.twinx()
    ax_top2.plot(mapq_plot["MAPQ"], mapq_density, color="#d62828", linewidth=2)

    ax_top.set_ylabel("Count")
    ax_top2.set_ylabel("Density", color="#d62828")
    ax_top.tick_params(axis="x", labelbottom=False)
    ax_top.spines["right"].set_visible(False)
    ax_top.spines["top"].set_visible(False)

    # right marginal: NM histogram + density
    ax_right.barh(nm_plot["NM"], nm_plot["count"], color="#ffb703", alpha=0.8)
    ax_right2 = ax_right.twiny()
    ax_right2.plot(nm_density, nm_plot["NM"], color="#d62828", linewidth=2)

    ax_right.set_xlabel("Count")
    ax_right2.set_xlabel("Density", color="#d62828")
    ax_right.tick_params(axis="y", labelleft=False)
    ax_right.spines["right"].set_visible(False)
    ax_right.spines["top"].set_visible(False)

    plt.tight_layout()
    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close()


def main():
    args = parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    sample = args.sample if args.sample else infer_sample_name(args.bam)
    sample_outdir = os.path.join(args.outdir, sample)
    os.makedirs(sample_outdir, exist_ok=True)

    pair_df, mapq_df, nm_df, stats = analyze_bam(
        bam_path=args.bam,
        max_nm=args.max_nm,
        mapq_max=args.mapq_max,
        skip_secondary=args.skip_secondary,
        skip_supplementary=args.skip_supplementary,
        skip_duplicate=args.skip_duplicate,
        skip_qcfail=args.skip_qcfail,
        skip_unmapped=args.skip_unmapped
    )

    pair_tsv = os.path.join(sample_outdir, f"{sample}.mapq_nm_pair_counts.tsv")
    mapq_tsv = os.path.join(sample_outdir, f"{sample}.mapq_distribution.tsv")
    nm_tsv = os.path.join(sample_outdir, f"{sample}.nm_distribution.tsv")
    stats_txt = os.path.join(sample_outdir, f"{sample}.summary.txt")
    plot_png = os.path.join(sample_outdir, f"{sample}.mapq_nm_jointplot.png")

    pair_df.to_csv(pair_tsv, sep="\t", index=False)
    mapq_df.to_csv(mapq_tsv, sep="\t", index=False)
    nm_df.to_csv(nm_tsv, sep="\t", index=False)

    with open(stats_txt, "w") as f:
        for k, v in stats.items():
            f.write(f"{k}\t{v}\n")

    plot_joint(
        pair_df=pair_df,
        mapq_df=mapq_df,
        nm_df=nm_df,
        out_png=plot_png,
        sample=sample,
        max_nm=args.max_nm,
        mapq_max=args.mapq_max
    )

    print(f"[INFO] Sample: {sample}")
    print(f"[INFO] Pair counts: {pair_tsv}")
    print(f"[INFO] MAPQ distribution: {mapq_tsv}")
    print(f"[INFO] NM distribution: {nm_tsv}")
    print(f"[INFO] Summary: {stats_txt}")
    print(f"[INFO] Plot: {plot_png}")


if __name__ == "__main__":
    main()

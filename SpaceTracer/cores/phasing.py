from collections import defaultdict, Counter
from dataclasses import dataclass
from functools import partial
from math import ceil
import random
from typing import Iterator, List, Optional, Tuple
import gc
import os

import numpy as np
import pandas as pd
import pyranges as pr
import pysam
import scipy.stats

from SpaceTracer.utils.handle_UMI_combine import (
    calculate_UMI_combine_phred,
    calculate_UMI_combine_phred_count_dict,
    get_most_candidate_allele, 
    handle_cigar,
    handle_seq,
    handle_pos,
    handle_quality_matrix,
    handle_seq_type
)

from SpaceTracer.utils.logger import get_logger
from SpaceTracer.utils.parallel import parallel_imap
model_name=__name__
logger = get_logger("<core function>: "+model_name)


@dataclass
class PhaseConfig:
    fasta: str
    bam: str
    germline: str
    indgeno: str
    out_phasing_file: str
    out_cluster_file: str
    gene_bed : str
    phasing_chromosomes: List[str]
    seq_type: str = "visium"
    bin_size: int = 100
    minprior: float = 0.01
    thread: int = 4
    species: str = "human"
    min_total_dp: int = 50
    min_dp: int = 20
    alpha: float = 0.05
    max_dist: int = 10
    phasing_pad: int = 1000
    merge_gap: int = 200
    max_target: int = 100000
    seed: int = 42
    max_region_span: int = 10000
    memory_limit_bytes: Optional[int] = None
    phasing_max_workers: Optional[int] = None

def _dedup_list_keep_order(values):
    return list(dict.fromkeys(values))


def _intervals_close_or_overlap(a, b, gap=0):
    return not (a[1] + gap < b[0] or b[1] + gap < a[0])
    
def read_table_auto(path: str, header=None, comment="#"):
    if path.endswith(".gz"):
        return pd.read_csv(path, sep="\t", header=header, comment=comment, compression="gzip")
    return pd.read_csv(path, sep="\t", header=header, comment=comment)


def parse_germline_to_df(
    germline_file: str, 
    species: str, 
    minprior: float, 
    min_dp: int, 
    min_total_dp: int) -> pd.DataFrame:
    """
    Input expected columns:
    #chrom  site    genotype    allele  allele_count prior
    """
    df = pd.read_csv(germline_file, sep="\t", header=0)

    expected = df.columns.tolist()
    
    if len(expected) >= 6:
        df = df.iloc[:, :6].copy()
        df.columns = ["Chromosome", "Position", "genotype", "allele", "allele_count", "prior"]
    else:
        raise ValueError(f"Unexpected germline format: {germline_file}")

    df = df[df["genotype"] == "het"].copy()
    if df.empty:
        return pd.DataFrame()

    allele_counts = df["allele_count"].str.split(",", expand=True)
    priors = df["prior"].str.split(",", expand=True)
    alleles = df["allele"].str.split(",", expand=True)

    df["allele1"] = alleles[0]
    df["allele2"] = alleles[1]
    df["count1"] = pd.to_numeric(allele_counts[0], errors="coerce")
    df["count2"] = pd.to_numeric(allele_counts[1], errors="coerce")
    df["prior1"] = pd.to_numeric(priors[0], errors="coerce")
    df["prior2"] = pd.to_numeric(priors[1], errors="coerce")
    df["total"] = df["count1"] + df["count2"]

    if species == "human":
        df = df[
            (df["count1"] > min_dp) &
            (df["count2"] > min_dp) &
            (df["prior1"] > minprior) &
            (df["prior2"] > minprior) &
            (df["total"] > min_total_dp)
        ].copy()
    else:
        df = df[
            (df["count1"] > min_dp) &
            (df["count2"] > min_dp) &
            (df["prior1"] > 0) &
            (df["prior2"] > 0) &
            (df["total"] > min_total_dp)
        ].copy()

    out = pd.DataFrame({
        "Chromosome": df["Chromosome"],
        "Start": df["Position"].astype(int),
        "End": df["Position"].astype(int),
        "type": df["genotype"],
        "allele1": df["allele1"],
        "allele2": df["allele2"],
        "allele_count": df["allele_count"],
        "prior": df["prior"],
    })
    return out


def parse_indgeno_mosaic_df(indgeno_file: str) -> pd.DataFrame:
    """
    Expected columns similar to:
    #chrom site ID germline mutant cluster spot_number consensus_read_count genotype p_mosaic Gi vaf_hat
    """
    df= pd.read_csv(indgeno_file, sep="\t", header=0)
    # cols = df.columns.tolist()

    # if len(cols) < 12:
    #     raise ValueError(f"Unexpected indgeno format: {indgeno_file}")
    df = df[df["genotype"] == "mosaic"].copy()

    df["total_umi"] = df["consensus_read_count"].astype(str).apply(
            lambda x: sum(int(i) for i in x.split(",") if i != "")
        )
    df["vaf_hat"] = df["vaf_hat"].astype(float)
    df["mutant_umi"] = df["total_umi"] * df["vaf_hat"]

    out = pd.DataFrame({
        "Chromosome": df["#chrom"],
        "Start": df["pos"].astype(int),
        "End": df["pos"].astype(int),
        "type": df["genotype"],
        "allele1": df["germline"],
        "allele2": df["mutant"],
        "allele_count": df["consensus_read_count"],
        "prior": df["prior_ATCG"],
    })
    return out


def parse_gene_bed_to_df(gene_bed_file: str,phasing_chromosomes:list) -> pd.DataFrame:
    df = read_table_auto(gene_bed_file, header=None, comment=None)
    if df.shape[1] < 4:
        raise ValueError("gene bed file must have at least 4 columns")
    df = df.iloc[:, :4].copy()
    df.columns = ["Chromosome", "Start", "End", "gene_name"]

    df = df[~df["gene_name"].astype(str).str.contains("MT-", na=False)].copy()

    if phasing_chromosomes is not None and len(phasing_chromosomes) > 0:
        keep_set = set(map(str, phasing_chromosomes))
        df = df[df["Chromosome"].astype(str).isin(keep_set)].copy()

    df["Start"] = df["Start"].astype(int)
    df["End"] = df["End"].astype(int)
    return df


def build_gene_and_cluster_df(
    gene_df: pd.DataFrame, 
    het_df: pd.DataFrame, 
    mosaic_df: pd.DataFrame, 
    max_dist: int = 10,
    phasing_pad: int = 1000,
    merge_gap: int = 200
) -> pd.DataFrame:

    het_df = add_site_identifier(het_df)
    mosaic_df = add_site_identifier(mosaic_df)
    cluster_df = prepare_cluster_df(mosaic_df, max_dist=max_dist)
    gene_het_region_df = build_gene_het_region_df(
        gene_df=gene_df,
        het_df=het_df,
        phasing_pad=phasing_pad,
        merge_gap=merge_gap,
    )

    final_df = intersect_cluster_with_gene_het_regions(
        cluster_df=cluster_df,
        gene_het_region_df=gene_het_region_df,
    )
    return final_df


def merge_intervals_with_ids(interval_records, merge_gap=500):
    


    if not interval_records:
        return []

    interval_records = sorted(interval_records, key=lambda x: (x["start"], x["end"]))

    merged = [{
        "start": interval_records[0]["start"],
        "end": interval_records[0]["end"],
        "identifiers": [interval_records[0]["identifier"]],
    }]

    for rec in interval_records[1:]:
        last = merged[-1]

        if _intervals_close_or_overlap(
            (last["start"], last["end"]),
            (rec["start"], rec["end"]),
            gap=merge_gap
        ):
            last["start"] = min(last["start"], rec["start"])
            last["end"] = max(last["end"], rec["end"])
            last["identifiers"].append(rec["identifier"])
        else:
            merged.append({
                "start": rec["start"],
                "end": rec["end"],
                "identifiers": [rec["identifier"]],
            })

    for m in merged:
        m["identifiers"] = _dedup_list_keep_order(m["identifiers"])

    return merged


def add_site_identifier(site_df: pd.DataFrame) -> pd.DataFrame:
    df = site_df.copy()

    if df.empty:
        df["identifier"] = pd.Series(dtype=object)
        return df

    df["Chromosome"] = df["Chromosome"].astype(str)
    df["Start"] = df["Start"].astype(int)
    df["End"] = df["End"].astype(int)

    df["identifier"] = list(zip(
        df["Chromosome"],
        df["Start"],
        df["allele1"].astype(str),
        df["allele2"].astype(str),
    ))
    return df


def prepare_cluster_df(mosaic_df: pd.DataFrame, max_dist: int = 10) -> pd.DataFrame:
    """
    cluster by mosaic_df, and return candidate cluster event 
    return:
    Chromosome, cluster_id, cluster_start, cluster_end, mosaic_sites, n_sites, cluster_span, cluster
    """
    final_cols = [
        "Chromosome", "cluster_id", "cluster_start", "cluster_end",
        "mosaic_sites", "n_sites", "cluster_span", "cluster"
    ]

    if mosaic_df.empty:
        return pd.DataFrame(columns=final_cols)

    df = mosaic_df.copy().sort_values(["Chromosome", "Start"]).reset_index(drop=True)

    df["prev_start"] = df.groupby("Chromosome", observed=False)["Start"].shift(1)
    df["dist_prev"] = df["Start"] - df["prev_start"]

    df["new_cluster"] = (
        df["prev_start"].isna() | (df["dist_prev"] > max_dist)
    ).astype(int)

    df["cluster_id"] = (
        df.groupby("Chromosome", observed=False)["new_cluster"].cumsum()
    )

    cluster_df = (
        df.groupby(["Chromosome", "cluster_id"], observed=False)
        .agg(
            cluster_start=("Start", "min"),
            cluster_end=("Start", "max"),
            mosaic_sites=("identifier", _dedup_list_keep_order),
            n_sites=("Start", "size"),
        )
        .reset_index()
    )

    cluster_df["cluster_span"] = cluster_df["cluster_end"] - cluster_df["cluster_start"]
    cluster_df["cluster"] = cluster_df["n_sites"] > 1

    return cluster_df[final_cols].copy()


def build_gene_het_region_df(
    gene_df: pd.DataFrame,
    het_df: pd.DataFrame,
    phasing_pad: int = 1000,
    merge_gap: int = 500,
) -> pd.DataFrame:
    """
    step 1. intersect the heterzygous and gene
    stpe 2. refine phasing region based on het-in-gene [max(het_pos-phasing_pad, gene_start), min(het_pos+phasing_pad, gene_end)]
    step 3. megre the heterzygous phasing region within one gene, if the distance <= merge_gap

    return:
    Chromosome, gene_name, gene_start, gene_end, region_start, region_end, het_sites
    """
    final_cols = [
        "Chromosome", "gene_name", "gene_start", "gene_end",
        "region_start", "region_end", "het_sites"
    ]

    if gene_df.empty or het_df.empty:
        return pd.DataFrame(columns=final_cols)

    gene_df2 = gene_df.copy()
    gene_df2["Chromosome"] = gene_df2["Chromosome"].astype(str)
    gene_df2["Start"] = gene_df2["Start"].astype(int)
    gene_df2["End"] = gene_df2["End"].astype(int)
    if "_gene_id" not in gene_df2.columns:
        gene_df2["_gene_id"] = range(len(gene_df2))

    het_df2 = het_df.copy()
    het_df2["Chromosome"] = het_df2["Chromosome"].astype(str)
    het_df2["Start"] = het_df2["Start"].astype(int)
    het_df2["End"] = het_df2["End"].astype(int)

    overlap_df = pr.PyRanges(gene_df2).join(pr.PyRanges(het_df2)).df

    if overlap_df.empty:
        return pd.DataFrame(columns=final_cols)

    gene_name_col = "gene_name" if "gene_name" in overlap_df.columns else "gene_name_b"
    gene_start_col = "Start" if "Start" in overlap_df.columns else "Start_b"
    gene_end_col = "End" if "End" in overlap_df.columns else "End_b"

    het_start_col = "Start_b" if "Start_b" in overlap_df.columns else "Start"
    identifier_col = "identifier_b" if "identifier_b" in overlap_df.columns else "identifier"

    if "_gene_id" not in overlap_df.columns and "_gene_id_b" in overlap_df.columns:
        overlap_df["_gene_id"] = overlap_df["_gene_id_b"]

    overlap_df["region_start"] = (
        overlap_df[[het_start_col, gene_start_col]]
        .apply(lambda x: max(int(x.iloc[0]) - phasing_pad, int(x.iloc[1])), axis=1)
    )
    overlap_df["region_end"] = (
        overlap_df[[het_start_col, gene_end_col]]
        .apply(lambda x: min(int(x.iloc[0]) + phasing_pad, int(x.iloc[1])), axis=1)
    )

    results = []

    group_cols = ["Chromosome", "_gene_id", gene_name_col, gene_start_col, gene_end_col]

    for keys, sub_df in overlap_df.groupby(group_cols, observed=False):
        chrom, gene_id, gene_name, gene_start, gene_end = keys

        interval_records = []
        for _, row in sub_df.iterrows():
            interval_records.append({
                "start": int(row["region_start"]),
                "end": int(row["region_end"]),
                "identifier": row[identifier_col],
            })

        merged_blocks = merge_intervals_with_ids(interval_records, merge_gap=merge_gap)

        for block in merged_blocks:
            results.append({
                "Chromosome": str(chrom),
                "gene_name": gene_name,
                "gene_start": int(gene_start),
                "gene_end": int(gene_end),
                "region_start": int(block["start"]),
                "region_end": int(block["end"]),
                "het_sites": block["identifiers"],
            })

    out = pd.DataFrame(results, columns=final_cols)

    if out.empty:
        return pd.DataFrame(columns=final_cols)

    out = out.sort_values(["Chromosome", "region_start", "region_end", "gene_name"]).reset_index(drop=True)
    return out

def intersect_cluster_with_gene_het_regions(
    cluster_df: pd.DataFrame,
    gene_het_region_df: pd.DataFrame,
) -> pd.DataFrame:
    """
    intersect the cluster_df and gene_het_region_df

    return:
    1. candidate cluster events () and gene_het_region -> phasing=True
    2. for the candidate cluster events without germline -> phasing=False, cluster =True

    filter:
    cluster_df["cluster"]==False, and no overlap with gene_het_region_df
    """
    final_cols = [
        "Chromosome", "Start", "End", "gene_name", "cluster_id",
        "cluster_region_start", "cluster_region_end",
        "het_sites", "mosaic_sites", "phasing", "cluster"
    ]

    if cluster_df.empty:
        return pd.DataFrame(columns=final_cols)

    # -------------------------
    # if gene_het_region is empty, remain cluster events
    # -------------------------
    if gene_het_region_df.empty:
        out = cluster_df[cluster_df["cluster"]].copy()

        if out.empty:
            return pd.DataFrame(columns=final_cols)

        out["Start"] = out["cluster_start"]
        out["End"] = out["cluster_end"]
        out["gene_name"] = (
            out["Chromosome"].astype(str) + "_cluster_" + out["cluster_id"].astype(str)
        )
        out["cluster_region_start"] = out["cluster_start"]
        out["cluster_region_end"] = out["cluster_end"]
        out["het_sites"] = [[] for _ in range(len(out))]
        out["phasing"] = False

        out = out[[
            "Chromosome", "Start", "End", "gene_name", "cluster_id",
            "cluster_region_start", "cluster_region_end",
            "het_sites", "mosaic_sites", "phasing", "cluster"
        ]].copy()

        out = out.sort_values(["Chromosome", "Start", "End", "gene_name", "cluster_id"]).reset_index(drop=True)
        return out

    # -------------------------
    # overlap
    # -------------------------
    cluster_pr_df = cluster_df.rename(columns={
        "cluster_start": "Start",
        "cluster_end": "End",
    }).copy()

    gene_region_pr_df = gene_het_region_df.rename(columns={
        "region_start": "Start",
        "region_end": "End",
    }).copy()

    overlap_df = pr.PyRanges(cluster_pr_df).join(pr.PyRanges(gene_region_pr_df)).df

    overlap_results = []

    if not overlap_df.empty:
        for _, row in overlap_df.iterrows():
            cluster_start = int(row["Start"])
            cluster_end = int(row["End"])

            gene_region_start = int(row["Start_b"])
            gene_region_end = int(row["End_b"])

            final_start = min(cluster_start, gene_region_start)
            final_end = max(cluster_end, gene_region_end)

            overlap_results.append({
                "Chromosome": str(row["Chromosome"]),
                "Start": final_start,
                "End": final_end,
                "gene_name": row["gene_name"],
                "cluster_id": row["cluster_id"],
                "cluster_region_start": cluster_start,
                "cluster_region_end": cluster_end,
                "het_sites": row["het_sites"],
                "mosaic_sites": row["mosaic_sites"],
                "phasing": True,
                "cluster": bool(row["cluster"]),
            })

    overlap_out = pd.DataFrame(overlap_results, columns=final_cols)

    # -------------------------
    # non-overlap cluster events only
    # -------------------------
    if overlap_df.empty:
        overlapped_cluster_keys = set()
    else:
        overlapped_cluster_keys = set(
            zip(
                overlap_df["Chromosome"].astype(str),
                overlap_df["cluster_id"]
            )
        )

    cluster_keys = list(zip(cluster_df["Chromosome"].astype(str), cluster_df["cluster_id"]))
    non_overlap_mask = [key not in overlapped_cluster_keys for key in cluster_keys]

    non_overlap_cluster_df = cluster_df.loc[non_overlap_mask].copy()

    # keep those cluster event with ['cluster']==True
    non_overlap_cluster_df = non_overlap_cluster_df[non_overlap_cluster_df["cluster"]].copy()

    if non_overlap_cluster_df.empty:
        non_overlap_out = pd.DataFrame(columns=final_cols)
    else:
        non_overlap_cluster_df["Start"] = non_overlap_cluster_df["cluster_start"]
        non_overlap_cluster_df["End"] = non_overlap_cluster_df["cluster_end"]
        non_overlap_cluster_df["gene_name"] = (
            non_overlap_cluster_df["Chromosome"].astype(str)
            + "_cluster_"
            + non_overlap_cluster_df["cluster_id"].astype(str)
        )
        non_overlap_cluster_df["cluster_region_start"] = non_overlap_cluster_df["cluster_start"]
        non_overlap_cluster_df["cluster_region_end"] = non_overlap_cluster_df["cluster_end"]
        non_overlap_cluster_df["het_sites"] = [[] for _ in range(len(non_overlap_cluster_df))]
        non_overlap_cluster_df["phasing"] = False

        non_overlap_out = non_overlap_cluster_df[[
            "Chromosome", "Start", "End", "gene_name", "cluster_id",
            "cluster_region_start", "cluster_region_end",
            "het_sites", "mosaic_sites", "phasing", "cluster"
        ]].copy()

    # -------------------------
    # merge all
    # -------------------------
    out = pd.concat([overlap_out, non_overlap_out], ignore_index=True)

    if out.empty:
        return pd.DataFrame(columns=final_cols)

    out = out.sort_values(["Chromosome", "Start", "End", "gene_name", "cluster_id"]).reset_index(drop=True)
    return out



############# func for phasing #############
def filter_geno_dict(count_result):
    hSNP_dict = defaultdict(int)
    for k, v in count_result.items():
        new_k = k[-1]
        hSNP_dict[new_k] += v

    hSNP_rank_list = sorted(hSNP_dict.items(), key=lambda item: item[1], reverse=True)
    h1, h2, h3, h4 = 0, 0, 0, 0
    try:
        h1 = int(hSNP_rank_list[0][1]); h1_g = hSNP_rank_list[0][0]
    except:
        h1_g = "NONE"

    try:
        h2 = int(hSNP_rank_list[1][1]); h2_g = hSNP_rank_list[1][0]
    except:
        h2_g = "H2none"

    try:
        h3 = int(hSNP_rank_list[2][1]); h3_g = hSNP_rank_list[2][0]
    except:
        h3_g = ""

    try:
        h4 = int(hSNP_rank_list[3][1]); h4_g = hSNP_rank_list[3][0]
    except:
        h4_g = ""

    dp = h1 + h2 + h3 + h4
    if dp == 0:
        return h1_g, h2_g, {}

    binom_test_p_h1 = scipy.stats.binomtest(ceil(100 * h1 / dp), 100, p=0.5, alternative="two-sided").pvalue
    binom_test_p_h2 = scipy.stats.binomtest(ceil(100 * h2 / dp), 100, p=0.5, alternative="two-sided").pvalue

    if binom_test_p_h1 > 0.05 and binom_test_p_h2 > 0.05:
        short_count_result = {k: v for k, v in count_result.items() if k[-1] in [h1_g, h2_g]}
        geno_rank_list = sorted(short_count_result.items(), key=lambda item: item[1], reverse=True)
        geno_count_dict = dict(geno_rank_list)
    else:
        geno_count_dict = {}

    return h1_g, h2_g, geno_count_dict


def calculate_phased_haplo(geno_count_dict, germline, mutant, h1, h2):
    germline = germline.split(",")
    germline_list = list(set(germline))
    haplo = ""
    annotated_type = ""
    scale_ratio = 5

    if len(germline_list) == 1:
        ref_h1, ref_h2, alt_h1, alt_h2 = 0, 0, 0, 0
        ref_h1 = geno_count_dict.get(germline_list[0] + "," + h1, 0)
        ref_h2 = geno_count_dict.get(germline_list[0] + "," + h2, 0)
        alt_h1 = geno_count_dict.get(mutant + "," + h1, 0)
        alt_h2 = geno_count_dict.get(mutant + "," + h2, 0)

        dp = ref_h1 + ref_h2 + alt_h1 + alt_h2
        mut_allele = alt_h1 + alt_h2
        ref_alt_phased = (ref_h1 > scale_ratio * ref_h2 and alt_h2 > scale_ratio * alt_h1) or \
                         (ref_h2 > scale_ratio * ref_h1 and alt_h1 > scale_ratio * alt_h2)
        ref_differ = (ref_h1 > scale_ratio * ref_h2 or ref_h2 > scale_ratio * ref_h1)
        alt_differ = (alt_h2 > scale_ratio * alt_h1 or alt_h1 > scale_ratio * alt_h2)

        detail_count_list = [
            germline_list[0] + "," + h1 + ":" + str(ref_h1),
            germline_list[0] + "," + h2 + ":" + str(ref_h2),
            mutant + "," + h1 + ":" + str(alt_h1),
            mutant + "," + h2 + ":" + str(alt_h2),
        ]

        if ref_alt_phased and dp > 4:
            haplo = "haplo=2"
            annotated_type = "heterzygous"
        elif (not ref_differ) and alt_differ and alt_h1 + alt_h2 >= 2 and ((ref_h1 > ref_h2 and alt_h1 < alt_h2) or (ref_h2 > ref_h1 and alt_h2 < alt_h1)):
            haplo = "haplo=3"
            annotated_type = "mut_mosaic"
        else:
            haplo = "haplo>3"
            annotated_type = "artifacts"

    elif len(germline_list) == 2:
        germline_list = germline
        ref_h1 = geno_count_dict.get(germline_list[0] + "," + h1, 0)
        ref_h2 = geno_count_dict.get(germline_list[0] + "," + h2, 0)
        alt_h1 = geno_count_dict.get(germline_list[1] + "," + h1, 0)
        alt_h2 = geno_count_dict.get(germline_list[1] + "," + h2, 0)
        mut_h1 = geno_count_dict.get(mutant + "," + h1, 0)
        mut_h2 = geno_count_dict.get(mutant + "," + h2, 0)

        detail_count_list = [
            germline_list[0] + "," + h1 + ":" + str(ref_h1),
            germline_list[0] + "," + h2 + ":" + str(ref_h2),
            germline_list[1] + "," + h1 + ":" + str(alt_h1),
            germline_list[1] + "," + h2 + ":" + str(alt_h2),
            mutant + "," + h1 + ":" + str(mut_h1),
            mutant + "," + h2 + ":" + str(mut_h2),
        ]

        dp = ref_h1 + ref_h2 + alt_h1 + alt_h2 + mut_h1 + mut_h2
        mut_allele = mut_h1 + mut_h2
        ref_alt_phased = (ref_h1 > scale_ratio * ref_h2 and alt_h2 > scale_ratio * alt_h1) or \
                         (ref_h2 > scale_ratio * ref_h1 and alt_h1 > scale_ratio * alt_h2)
        ref_differ = (ref_h1 > scale_ratio * ref_h2 or ref_h2 > scale_ratio * ref_h1)
        alt_differ = (alt_h2 > scale_ratio * alt_h1 or alt_h1 > scale_ratio * alt_h2)
        mut_differ = (mut_h1 > scale_ratio * mut_h2 or mut_h2 > scale_ratio * mut_h1)

        if ref_alt_phased and mut_h1 + mut_h2 < 2:
            haplo = "haplo>3"
            annotated_type = "artifacts"
        elif ref_alt_phased and mut_differ:
            haplo = "haplo=3"
            annotated_type = "mut_mosaic"
        else:
            haplo = "haplo>3"
            annotated_type = "artifacts"

    return dp, mut_allele, haplo, annotated_type, detail_count_list


def judge_origin(germline, annotated_type, detail_count_list):
    def get_count(detail_count_list, index):
        return int(detail_count_list[index].split(":")[-1])

    mut_origin = "NA"
    if len(germline.split(",")) == 2 and annotated_type == "alt_mosaic" and len(detail_count_list) == 6:
        mut_origin = germline
    elif len(germline.split(",")) == 2 and annotated_type == "mut_mosaic" and len(detail_count_list) == 6:
        ref = germline.split(",")[0]
        alt = germline.split(",")[1]
        mut_h1, mut_h2 = get_count(detail_count_list, 4), get_count(detail_count_list, 5)
        if mut_h1 > mut_h2:
            mut_origin = ref if get_count(detail_count_list, 0) > get_count(detail_count_list, 2) else alt
        elif mut_h1 < mut_h2:
            mut_origin = ref if get_count(detail_count_list, 1) > get_count(detail_count_list, 3) else alt
    elif len(germline.split(",")) == 1 and len(detail_count_list) == 4:
        mut_origin = germline

    return mut_origin


def _iter_subspans(start: int, end: int, max_span: int) -> Iterator[Tuple[int, int]]:
    """Split [start, end) into windows of at most max_span bp for bounded BAM fetch."""
    if max_span <= 0 or end - start <= max_span:
        yield start, end
        return
    pos = start
    while pos < end:
        sub_end = min(pos + max_span, end)
        yield pos, sub_end
        pos = sub_end


def _fetch_reads_reservoir(
    bam: pysam.AlignmentFile,
    chrom: str,
    start: int,
    end: int,
    k: int,
    seed: int,
    max_span: int = 10000,
) -> Tuple[List, int, bool]:
    """
    Stream-fetch reads with reservoir sampling so memory stays O(k), not O(region depth).
    Returns (sampled_reads, total_seen, adjusted).
    """
    if k <= 0:
        total_seen = 0
        for sub_start, sub_end in _iter_subspans(start, end, max_span):
            for _ in bam.fetch(chrom, sub_start, sub_end):
                total_seen += 1
        return [], total_seen, total_seen > 0

    rng = random.Random(seed)
    reservoir: List = []
    total_seen = 0

    for sub_start, sub_end in _iter_subspans(start, end, max_span):
        for read in bam.fetch(chrom, sub_start, sub_end):
            total_seen += 1
            if len(reservoir) < k:
                reservoir.append(read)
            else:
                j = rng.randint(0, total_seen - 1)
                if j < k:
                    reservoir[j] = read

    adjusted = total_seen > k
    return reservoir, total_seen, adjusted


def _accumulate_reads_into_dict(
    region_reads,
    run_type,
    bins,
    candidate_allele_info,
    per_read_dict,
    allele_total_count,
    used_barcode,
):
    for reads in region_reads:
        barcode_name, UMI_name = handle_seq_type(reads, run_type, bins)
        if barcode_name is None or UMI_name is None:
            continue

        seq_cut, pos_cut = handle_cigar(reads.cigar)
        cut_seq = handle_seq(reads.seq, seq_cut)
        cut_pos = handle_pos(reads.get_reference_positions(), pos_cut)
        barcode_name = f"{barcode_name}_{UMI_name}"

        for candidate_allele_tuple in candidate_allele_info:
            candidate_allele, ref = candidate_allele_tuple
            pos_index = int(candidate_allele) - 1

            if pos_index not in cut_pos:
                continue

            raw_index = handle_quality_matrix(cut_pos.index(pos_index), reads.seq, cut_seq)
            quality = reads.get_forward_qualities()[raw_index]
            geno = cut_seq[cut_pos.index(pos_index)]

            if geno not in "ATCG":
                continue

            if candidate_allele not in per_read_dict[barcode_name]:
                per_read_dict[barcode_name][candidate_allele] = {
                    "count": defaultdict(int),
                    "quality": {
                        "A": defaultdict(int),
                        "T": defaultdict(int),
                        "C": defaultdict(int),
                        "G": defaultdict(int),
                    },
                }

            per_read_dict[barcode_name][candidate_allele]["count"][geno] += 1
            per_read_dict[barcode_name][candidate_allele]["quality"][geno][quality] += 1

            if barcode_name not in used_barcode[candidate_allele]:
                allele_total_count[candidate_allele] += 1
                used_barcode[candidate_allele].append(barcode_name)


def get_candidate_germline(
    ref_fasta,
    in_bam_name,
    phasing_chromosomes,
    gene_name_col,
    run_type,
    bins,
    min_total_dp,
    alpha,
    max_target,
    seed,
    max_region_span,
    line_dict
):
    effective_max_target = line_dict.pop("_effective_max_target", max_target)

    out_list = []
    cluster_event_list = []

    chr_ = str(line_dict["Chromosome"])

    phasing_flag = bool(line_dict.get("phasing", False))
    cluster_flag = bool(line_dict.get("cluster", False))

    if chr_ not in phasing_chromosomes:
        phasing_flag = False

    if not phasing_flag and not cluster_flag:
        return out_list, cluster_event_list

    pos_s = int(line_dict["Start"])
    pos_e = int(line_dict["End"])
    gene_name = str(line_dict[gene_name_col]).replace('"', "")

    het_sites = line_dict.get("het_sites", [])
    mosaic_sites = line_dict.get("mosaic_sites", [])

    pos_candidate_dict = {"mosaic_pos": [], "informative_SNP": []}
    count_allele = []

    for site in het_sites:
        # site: (Chromosome, Start, allele1, allele2)
        site_chr, site_pos, allele1, allele2 = site
        pos_candidate_dict["informative_SNP"].append(site)
        count_allele.append((str(site_pos), str(allele1)[0]))

    for site in mosaic_sites:
        site_chr, site_pos, allele1, allele2 = site
        pos_candidate_dict["mosaic_pos"].append(site)
        count_allele.append((str(site_pos), str(allele1)[0]))

    if len(pos_candidate_dict["mosaic_pos"]) == 0 and len(pos_candidate_dict["informative_SNP"]) == 0:
        return [],[]

    candidate_allele_info = list(set(count_allele))
    per_read_dict = defaultdict(dict)
    allele_total_count = defaultdict(int)
    used_barcode = defaultdict(list)

    adjusted = False
    try:
        task_seed = seed ^ hash((chr_, pos_s, pos_e, gene_name)) & 0xFFFFFFFF
        with pysam.AlignmentFile(in_bam_name, "rb", reference_filename=ref_fasta) as bam:
            region_reads, original_region_depth, adjusted = _fetch_reads_reservoir(
                bam=bam,
                chrom=chr_,
                start=pos_s,
                end=pos_e,
                k=effective_max_target,
                seed=task_seed,
                max_span=max_region_span,
            )

        _accumulate_reads_into_dict(
            region_reads=region_reads,
            run_type=run_type,
            bins=bins,
            candidate_allele_info=candidate_allele_info,
            per_read_dict=per_read_dict,
            allele_total_count=allele_total_count,
            used_barcode=used_barcode,
        )
        del region_reads

    except Exception as e:
        logger.error(f"Error in fetch chr={chr_} {pos_s}-{pos_e}: {e}")

    per_read_genotypes_count = defaultdict(list)

    for barcode in per_read_dict.keys():
        for candidate_allele, ref in candidate_allele_info:
            if candidate_allele not in per_read_dict[barcode]:
                most_proved_allele = "."
            else:
                count_dict = per_read_dict[barcode][candidate_allele]["count"]
                quality_dict = per_read_dict[barcode][candidate_allele]["quality"]
                if count_dict:
                    phred_dict = calculate_UMI_combine_phred_count_dict(count_dict, quality_dict, weigh=0.5)
                    most_proved_allele, _ = get_most_candidate_allele(phred_dict, ref)
                else:
                    most_proved_allele = "."
            per_read_genotypes_count[barcode].append(most_proved_allele)

    candidate_allele_info_only_pos = [s[0] for s in candidate_allele_info]

    if phasing_flag:
        for mosaic_pos in pos_candidate_dict["mosaic_pos"]:
            chrom, pos, germline, mutant = mosaic_pos

            for info_SNP in pos_candidate_dict["informative_SNP"]:
                chrom_germ, pos_germ, _, _ = info_SNP
                if pos_germ == pos and chrom_germ == chrom:
                    continue

                index_mosaic = candidate_allele_info_only_pos.index(str(pos))
                index_germ = candidate_allele_info_only_pos.index(str(pos_germ))

                short_list = []
                for barcode in per_read_genotypes_count.keys():
                    values = per_read_genotypes_count[barcode]
                    bases = ",".join([values[index_mosaic], values[index_germ]])
                    short_list.append(bases)

                count_result = dict(Counter([bases for bases in short_list if "." not in bases]))
                h1, h2, geno_count_dict = filter_geno_dict(count_result)
                if h2 == "H2none":
                    continue

                if geno_count_dict:
                    new_mut_name = "_".join([str(chrom_germ), str(pos_germ), h1, h2])
                    total_count = allele_total_count[str(pos)]
                    dp, mut_allele, haplo, annotated_type, detail_count_list = calculate_phased_haplo(
                        geno_count_dict, germline, mutant, h1, h2
                    )
                    if mut_allele == 0:
                        continue

                    if annotated_type not in ["artifacts", "heterzygous"]:
                        mut_origin = judge_origin(germline, annotated_type, detail_count_list)
                    else:
                        mut_origin = "NA"

                    detail_count = ";".join(detail_count_list)
                    out_list.append([
                        chrom, str(pos), germline, mutant, mut_origin, new_mut_name, gene_name,
                        total_count, dp, mut_allele, haplo, annotated_type, detail_count,adjusted
                    ])

    if cluster_flag:
        mosaic_list = pos_candidate_dict["mosaic_pos"]

        # sorted by positions
        mosaic_list = sorted(mosaic_list, key=lambda x: (x[0], int(x[1])))

        significant_cluster_sites = set()

        for i in range(len(mosaic_list)):
            chrom1, pos1, germline1, mutant1 = mosaic_list[i]
            pos1_int = int(pos1)

            for j in range(i + 1, len(mosaic_list)):
                chrom2, pos2, germline2, mutant2 = mosaic_list[j]
                pos2_int = int(pos2)

                if chrom1 != chrom2:
                    continue

                if pos2_int - pos1_int > 10:
                    break

                index1 = candidate_allele_info_only_pos.index(str(pos1))
                index2 = candidate_allele_info_only_pos.index(str(pos2))

                short_list = []
                for barcode in per_read_genotypes_count.keys():
                    values = per_read_genotypes_count[barcode]
                    bases = ",".join([values[index1], values[index2]])
                    short_list.append(bases)

                count_result = Counter([bases for bases in short_list if "." not in bases])

                counts = {
                    "ref_ref": count_result.get(f"{germline1},{germline2}", 0),
                    "ref_alt": count_result.get(f"{germline1},{mutant2}", 0),
                    "alt_ref": count_result.get(f"{mutant1},{germline2}", 0),
                    "alt_alt": count_result.get(f"{mutant1},{mutant2}", 0),
                }

                matrix = np.array([
                    [counts["ref_ref"], counts["ref_alt"]],
                    [counts["alt_ref"], counts["alt_alt"]]
                ])

                if matrix.sum() < min_total_dp:
                    continue

                if (matrix < 5).any():
                    stat, p = scipy.stats.fisher_exact(matrix)
                else:
                    stat, p, _, _ = scipy.stats.chi2_contingency(matrix)

                if p < alpha:
                    significant_cluster_sites.add((chrom1, str(pos1), germline1, mutant1, True))
                    significant_cluster_sites.add((chrom2, str(pos2), germline2, mutant2, True))

        cluster_event_list = [
            list(x) for x in sorted(significant_cluster_sites, key=lambda v: (v[0], int(v[1]), v[2], v[3]))
        ]

    del per_read_dict
    del per_read_genotypes_count
    gc.collect()
    return out_list, cluster_event_list


def _pick_effective_max_target(span_bp: int, max_target: int) -> int:
    """Tighten read cap for large phasing regions to limit per-worker memory."""
    if span_bp > 50000:
        return min(max_target, 5000)
    if span_bp > 20000:
        return min(max_target, 10000)
    if span_bp > 10000:
        return min(max_target, 20000)
    return max_target


def _order_phasing_tasks(lines: List[dict]) -> List[dict]:
    """Run heavy (large-span) regions first so workers finish big jobs early."""
    return sorted(
        lines,
        key=lambda x: (-(int(x["End"]) - int(x["Start"])), str(x["Chromosome"]), int(x["Start"])),
    )


def run_phase_mode(config: PhaseConfig) -> str:

    germ_df = parse_germline_to_df(config.germline, config.species,config.minprior,config.min_dp, config.min_total_dp)

    mosaic_df = parse_indgeno_mosaic_df(config.indgeno)
    gene_df = parse_gene_bed_to_df(config.gene_bed,config.phasing_chromosomes)

    final_df = build_gene_and_cluster_df(
                    gene_df=gene_df,
                    het_df=germ_df,
                    mosaic_df=mosaic_df,
                    max_dist=config.max_dist,
                    phasing_pad=config.phasing_pad,
                    merge_gap=config.merge_gap
                )

    del germ_df, mosaic_df, gene_df
    gc.collect()
    
    gene_name_col = "gene_name"
    final_df_tmp_file = config.out_phasing_file.removesuffix(".txt") + ".raw_info.txt"
    final_df.to_csv(final_df_tmp_file,sep="\t",index=False)

    lines = final_df.to_dict("records")
    del final_df
    gc.collect()

    # Per-task effective max_target by region span
    for line in lines:
        span_bp = int(line["End"]) - int(line["Start"])
        line["_effective_max_target"] = _pick_effective_max_target(span_bp, config.max_target)

    lines = _order_phasing_tasks(lines)

    phasing_workers = config.phasing_max_workers or min(config.thread, 8)
    phasing_workers = max(1, min(phasing_workers, config.thread))

    if config.memory_limit_bytes:
        logger.info(
            f"[phasing] tasks={len(lines)}, workers={phasing_workers}, "
            f"max_target={config.max_target}, max_region_span={config.max_region_span}, "
            f"memory_limit={config.memory_limit_bytes / 1024**3:.1f}GB"
        )
    else:
        logger.info(
            f"[phasing] tasks={len(lines)}, workers={phasing_workers}, "
            f"max_target={config.max_target}, max_region_span={config.max_region_span}"
        )

    partial_func = partial(
        get_candidate_germline,
        config.fasta,
        config.bam,
        config.phasing_chromosomes,
        gene_name_col,
        config.seq_type,
        config.bin_size,
        config.min_total_dp,
        config.alpha,
        config.max_target,
        config.seed,
        config.max_region_span,
    )

    out_phasing_file = config.out_phasing_file
    out_cluster_file = config.out_cluster_file

    COLUMNS_phasing = [
        '#chrom', 'pos', 'germline', 'mutant', 'mut_origin', 'new_mut_name',
        'gene_name', 'total_count', 'dp', 'mut_allele', 'haplo',
        'annotated_type', 'detail_count','adjusted'
    ]

    COLUMNS_cluster = ['#chr', 'pos', 'germline', 'mutant', 'cluster_event']

    with open(out_phasing_file, "w") as f_phase, open(out_cluster_file, "w") as f_cluster:
        pd.DataFrame(columns=COLUMNS_phasing).to_csv(
            f_phase, header=True, index=False, sep='\t'
        )
        pd.DataFrame(columns=COLUMNS_cluster).to_csv(
            f_cluster, header=True, index=False, sep='\t'
        )

        n_done = 0
        for _, result in parallel_imap(
            items=lines,
            worker_fn=partial_func,
            max_workers=phasing_workers,
            desc="candidate_germline",
            raise_on_error=True,
            backend="process",
            progress_interval=0.05,
            max_in_flight=phasing_workers,
            memory_limit_bytes=config.memory_limit_bytes,
            logger=logger,
        ):
            if result is None:
                continue

            phasing_result, cluster_result = result

            if phasing_result:
                pd.DataFrame(phasing_result, columns=COLUMNS_phasing).to_csv(
                    f_phase, header=False, index=False, sep='\t'
                )

            if cluster_result:
                pd.DataFrame(cluster_result, columns=COLUMNS_cluster).to_csv(
                    f_cluster, header=False, index=False, sep='\t'
                )

            n_done += 1
            if n_done % 100 == 0:
                gc.collect()
                


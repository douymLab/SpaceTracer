import math
from pathlib import Path
from operator import add

import numpy as np
import pandas as pd
from scipy.special import comb

from SpaceTracer.utils.read_files import load_spot_genotypes_data, load_spot_count_data
from SpaceTracer.utils.utils import str2dict

def rename_df(df):
    if '#chrom' in df.columns:
        df=df.rename(columns={'#chrom':'chrom'})
    return df

class SpotGenoCalculator:
    def __init__(self, bins, epsQ, thr_dp, pop_vaf, cell_num):
        self.bins = bins
        self.epsQ = epsQ
        self.thr_dp = thr_dp
        self.pop_vaf = pop_vaf
        self.cell_num = cell_num

    def run_from_df(self, spot_count_df, ind_geno_df, cluster_df, cluster_vaf_df, output_file):
        raw_spot_count_df = spot_count_df.copy()
        raw_spot_count_df = raw_spot_count_df.fillna('NA')

        cluster_vaf_df = cluster_vaf_df.copy()
        cluster_vaf_df = cluster_vaf_df.drop(["spot_number", "consensus_read_count"], axis=1, errors='ignore')

        if self.bins is None:
            ind_geno_df = ind_geno_df.drop(["cluster", "spot_number", "consensus_read_count", "prior_ATCG"], axis=1, errors='ignore')
            spot_count_df = raw_spot_count_df
        else:
            bins = self.bins
            ind_geno_df = ind_geno_df.drop(["cluster", "spot_number", "consensus_read_count", "prior_ATCG"], axis=1, errors='ignore')
            spot_count_df = bin_df(raw_spot_count_df, bins)

        spot_count_df=rename_df(spot_count_df)
        cluster_df=rename_df(cluster_df)
        ind_geno_df=rename_df(ind_geno_df)
        cluster_vaf_df=rename_df(cluster_vaf_df)

        count_cluster_df = pd.merge(spot_count_df, cluster_df, on='barcode', how='left')
        count_geno_join = pd.merge(count_cluster_df, ind_geno_df, on=['chrom', 'pos', 'strand'], how='inner')
        count_geno_join = count_geno_join.rename(columns={'vaf': 'ind_vaf'})
        count_geno_vaf_join = pd.merge(
            count_geno_join,
            cluster_vaf_df,
            on=['chrom', 'pos', 'strand', 'germline', 'mutant', 'cluster'],
            how='left'
        )
        count_geno_vaf_join = count_geno_vaf_join.rename(columns={'vaf': 'cluster_vaf'})

        if not isinstance(self.cell_num, int) and Path(self.cell_num).exists():
            cell_num_file = pd.read_csv(self.cell_num, sep="\t", header=0)
            cell_num_file = cell_num_file.drop(["cluster", "nUMI", "nREAD"], axis=1, errors='ignore')
            count_geno_vaf_join = pd.merge(count_geno_vaf_join, cell_num_file, on='barcode', how='left')

        out_colname = "\t".join([
            "#chrom", "pos", "strand", "germline", "mutant", "cluster", "spot_barcode", "consensus_read_count",
            "l_germline", "l_mosaic", "max_spot_geno", "G_spot_max", "depth", "vaf", "p_mosaic"
        ])

        count_geno_vaf_join.to_csv(output_file+".tmp",index=False)

        mosaic_mask = count_geno_vaf_join["genotype"].eq("mosaic")
        mosaic_df = count_geno_vaf_join[mosaic_mask]
        data = mosaic_df.values

        if not isinstance(self.cell_num, int) and Path(self.cell_num).exists():
            cell_nums = mosaic_df.iloc[:, -1].to_numpy()
            mosaic_df_no_cell = mosaic_df.iloc[:, :-1]
        else:
            cell_nums = np.full(len(mosaic_df), self.cell_num)
            mosaic_df_no_cell = mosaic_df

        with open(output_file, 'w') as f:
            f.write(out_colname + "\n")
            for i, (_, row) in enumerate(mosaic_df_no_cell.iterrows()):
                spot_geno_info = spot_genotype(
                    join_info=row,   # 直接传 Series
                    cell_num=cell_nums[i],
                    epsQ=self.epsQ,
                    thr_dp=self.thr_dp,
                    pop_vaf=self.pop_vaf
                )
                if spot_geno_info[10] != "NA":
                    f.write("\t".join(map(str, spot_geno_info)) + "\n")
                    
        df = load_spot_genotypes_data(output_file,prefer_parquet=False)
        parquet_file = str(output_file).replace('.out', '.parquet')
        df.to_parquet(parquet_file, index=True, compression='snappy')

    def run(self, spot_count_file, ind_geno_file, cluster_df, cluster_vaf_file, output_file):
        ind_geno_df = pd.read_csv(ind_geno_file, sep="\t")
        raw_spot_count_df = load_spot_count_data(spot_count_file)
        cluster_vaf_df = pd.read_csv(cluster_vaf_file, sep="\t")
        self.run_from_df(raw_spot_count_df, ind_geno_df, cluster_df, cluster_vaf_df, output_file)


def bin_df(df, bins):
    def round_to_nearest_bin(x):
        return int(np.ceil(x / bins) * bins)

    df = df.copy()
    df[['x', 'y']] = df['spot_barcode'].str.split('_', expand=True)
    df['x'] = df['x'].astype(int)
    df['y'] = df['y'].astype(int)
    df['new_x'] = df['x'].apply(round_to_nearest_bin).astype(str)
    df['new_y'] = df['y'].apply(round_to_nearest_bin).astype(str)
    df['spot_barcode'] = df['new_x'] + "_" + df['new_y']
    df = df.drop(columns=['x', 'y', 'new_x', 'new_y'])
    return df


def spot_posterior(germline, mutant, cluster_vaf, qA, qT, qC, qG, cell_num=20, epsQ=20, thr_dp=1000, pop_vaf=1e-5):
    count_filter, q_filter = spot_filter(qA, qT, qC, qG, epsQ=epsQ, thr_dp=thr_dp)
    count_filter_nozero = {k: v for k, v in count_filter.items() if v > 0}
    allele_list_spot = list(count_filter_nozero.keys())
    depth = sum(count_filter.values())

    ind_ref_list = germline.split(",")
    ref_list = [ref for ref in ind_ref_list if ref in allele_list_spot]
    alt = mutant

    if mutant != ".":
        ind_allele_list = ind_ref_list + [mutant]
    else:
        ind_allele_list = ind_ref_list
    allele_list = [allele for allele in allele_list_spot if allele in ind_allele_list]
    n_allele = len(allele_list)

    if n_allele == 0:
        return "NA", ["NA"] * 2, ["NA", "NA", depth, "NA", "NA"]
    elif len(ref_list) == 0:
        ref_list = [ind_ref_list[0]]

    genotype_list = ["germline", "mosaic"]
    alt_count = count_filter[alt]
    qalt_dict = q_filter[alt]
    if alt == ".":
        alt_count = 0

    l = np.array([0] * 2)
    for ref in ref_list:
        ref_count = count_filter[ref]
        qref_dict = q_filter[ref]
        l_update = spot_likelihood(ref_count, alt_count, qref_dict, qalt_dict, cluster_vaf, pop_vaf)
        l = np.array(list(map(add, l, l_update)))

    s_likelihood = sum(l)
    s_likelihood = 1e-15 if s_likelihood == 0 else s_likelihood
    l_norm = [j / s_likelihood for j in l]

    if cluster_vaf < 0.5:
        prior_value = [(1 - 2 * cluster_vaf) ** cell_num, 1 - (1 - 2 * cluster_vaf) ** cell_num]
    else:
        prior_value = [pop_vaf, 1 - pop_vaf]

    posterior = prior_value * l
    s_posterior = sum(posterior)
    s_posterior = 1e-15 if s_posterior == 0 else s_posterior
    p = [j / s_posterior for j in posterior]
    p = np.nan_to_num(p)
    posterior_dict = dict(zip(genotype_list, p))
    max_spot_geno = max(posterior_dict, key=lambda x: posterior_dict[x])

    p_mosaic = posterior_dict['mosaic']
    vaf = alt_count / depth
    G_list = ["0/0", "0/1"]
    G_dict = dict(zip(genotype_list, G_list))
    G_spot_max = G_dict[max_spot_geno]

    if alt_count == 0:
        if p_mosaic > 0.5:
            p_mosaic = 0.5
        if l_norm[1] > l_norm[0]:
            l_norm = [0.5, 0.5]
    elif alt_count > 0 and l_norm[1] < 0.5:
        l_norm = [0.5, 0.5]

    return p, l_norm, [max_spot_geno, G_spot_max, depth, vaf, p_mosaic]


def spot_filter(qA, qT, qC, qG, epsQ=20, thr_dp=1000):
    nucleotide_list = ["A", "T", "C", "G"]
    q_list = [str2dict(q) for q in [qA, qT, qC, qG]]
    q_filter_list = [{k: v for k, v in q.items() if k >= epsQ} for q in q_list]

    depth = sum([sum(q.values()) for q in q_filter_list])
    if depth > thr_dp:
        scaling_factor = thr_dp / depth
        downsampled_q_list = [
            {key: int(round(value * scaling_factor)) for key, value in q.items()}
            for q in q_filter_list
        ]
        q_filter_list = downsampled_q_list

    q_filter = dict(zip(nucleotide_list, q_filter_list))
    count_list = [sum(q.values()) for q in q_filter_list]
    count_dict = dict(zip(nucleotide_list, count_list))
    return count_dict, q_filter


def spot_likelihood(ref_count, alt_count, qref_dict, qalt_dict, cluster_vaf, pop_vaf=1e-5):
    depth = ref_count + alt_count

    q_refhom = math.log10(1)
    if ref_count != 0:
        q_refhom = q_refhom + sum(math.log10(1 - 0.1 ** (i / 10)) * qref_dict[i] for i in qref_dict.keys())
    if alt_count != 0:
        q_refhom = q_refhom + sum(math.log10(0.1 ** (i / 10)) * qalt_dict[i] for i in qalt_dict.keys())
    l_refhom = math.log10(comb(depth, alt_count, exact=True)) + q_refhom
    l_refhom = 10 ** l_refhom

    cluster_vaf = pop_vaf if cluster_vaf == 0 else cluster_vaf
    cluster_vaf = 1 - pop_vaf if cluster_vaf == 1 else cluster_vaf

    r = 0
    if ref_count != 0:
        r = r + sum([0.1 ** (float(i) / 10) * qref_dict[i] for i in qref_dict.keys()])
    if alt_count != 0:
        r = r + sum([(1 - 0.1 ** (float(i) / 10)) * qalt_dict[i] for i in qalt_dict.keys()])
    l_mosaic = math.log10(comb(depth, alt_count, exact=True)) + r * math.log10(cluster_vaf) + (depth - r) * math.log10(1 - cluster_vaf)
    l_mosaic = 10 ** l_mosaic

    l = [l_refhom, l_mosaic]
    return l


# def spot_genotype(join_info, cell_num=20, epsQ=20, thr_dp=1000, pop_vaf=1e-5):
#     qA = join_info[8]
#     qT = join_info[9]
#     qC = join_info[10]
#     qG = join_info[11]
#     cluster = join_info[12]
#     germline = join_info[13]
#     mutant = join_info[14]
#     cluster_vaf = 0 if join_info[19] == "NA" else float(join_info[19])

#     _, l_norm, spot_geno = spot_posterior(
#         germline, mutant, cluster_vaf, qA, qT, qC, qG,
#         cell_num=cell_num, epsQ=epsQ, thr_dp=thr_dp, pop_vaf=pop_vaf
#     )

#     output = join_info[0:3] + [germline, mutant, str(cluster)] + [join_info[5] , join_info[7]] + l_norm + [str(i) for i in spot_geno]
#     return output

def spot_genotype(join_info, cell_num=20, epsQ=20, thr_dp=1000, pop_vaf=1e-5):
    """
    must contain: chrom,pos,strand,barcode,consensus_read_count,qA,qT,qC,qG,cluster,germline,mutant,cluster_vaf
    """
    # 支持传入 pd.Series / dict
    if hasattr(join_info, "to_dict"):
        row = join_info.to_dict()
    else:
        row = dict(join_info)

    qA = row.get("qA", "NA")
    qT = row.get("qT", "NA")
    qC = row.get("qC", "NA")
    qG = row.get("qG", "NA")

    cluster = row.get("cluster", "NA")
    germline = row.get("germline", "NA")
    mutant = row.get("mutant", "NA")

    cluster_vaf_raw = row.get("cluster_vaf", "NA")
    cluster_vaf = 0 if cluster_vaf_raw in ("NA", None, "") else float(cluster_vaf_raw)

    _, l_norm, spot_geno = spot_posterior(
        germline, mutant, cluster_vaf, qA, qT, qC, qG,
        cell_num=cell_num, epsQ=epsQ, thr_dp=thr_dp, pop_vaf=pop_vaf
    )

    # 输出字段顺序保持你原来的定义：
    # join_info[0:3] + [germline, mutant, cluster] + [barcode, consensus_read_count] + l_norm + spot_geno
    output = [
        str(row.get("chrom", "NA")),
        str(row.get("pos", "NA")),
        str(row.get("strand", "NA")),
        str(germline),
        str(mutant),
        str(cluster),
        str(row.get("barcode", "NA")),
        str(row.get("consensus_read_count", "NA")),
    ] + l_norm + [str(i) for i in spot_geno]

    return output
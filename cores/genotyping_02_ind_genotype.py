import gc
from operator import add

import math
import numpy as np
import pandas as pd
from scipy.special import comb, beta

from SpaceTracer.utils.utils import str2dict


def individual_posterior(ref, alt, qA, qT, qC, qG, fA, fT, fC, fG,
                         mu=1e-7, thr_dp=1000, pop_vaf=1e-5, filter_oneallele=True):
    count_dict, q_dicts, count_dict_original = quality_format(qA, qT, qC, qG, thr_dp=thr_dp)
    count_dict_nozero = {k: v for k, v in count_dict.items() if v > 0}
    n_allele = len(count_dict_nozero)
    allele_list = list(count_dict_nozero.keys())

    nucleotide_list = ["A", "T", "C", "G"]
    f_list = [pop_vaf if i == 0 else i for i in [float(fA), float(fT), float(fC), float(fG)]]
    f = dict(zip(nucleotide_list, f_list))

    germline_geno = None
    germline_allele = None
    germ_prior = None

    if n_allele == 0:
        return [["NA", "NA", "NA", ".", ".", "NA"]], ["NA", "NA", "NA", "NA"]

    elif n_allele == 1:
        allele = allele_list[0]
        candidate_list = ["A", "T", "C", "G"]
        candidate_list.remove(allele)
        tmp = candidate_list[0]
        if allele == ref:
            if filter_oneallele:
                return [["NA", "NA", "NA", ref, ".", "NA"]], ["NA", "NA", "NA", "NA"]
            ref_list = [ref]
            alt_list = [tmp]
            ind_ref = ref
            germline_geno = "refhom"
            germline_allele = ref
            germ_prior = str(f[ref])
        else:
            if filter_oneallele or ref == "N":
                return [["NA", "NA", "NA", ref, allele, "NA"]], ["NA", "NA", "NA", "NA"]
            ref_list = [tmp]
            alt_list = [allele]
            ind_ref = ref
            germline_geno = "althom"
            germline_allele = ",".join([ref, allele])
            try:
                germ_prior = ",".join([str(f[i]) for i in [ref, allele]])
            except Exception:
                print("Opps! The germ prior cannot be detected, in 'individual posterior'")

    else:
        if n_allele == 2 and ref in allele_list:
            ref_list = [ref]
            ind_ref = ref
            allele_list.remove(ref)
            alt_list = [allele_list[0]]
            alt = alt_list[0]
            marker_withref = True
        else:
            count_sort = sorted(count_dict_nozero.items(), key=lambda x: x[1], reverse=True)
            allele_list = list(dict(count_sort).keys())
            if ref in allele_list:
                allele_list.remove(ref)
            alt = allele_list[0]
            marker_withref = False

        ref_count = count_dict[ref]
        qref_dict = q_dicts[ref]
        alt_count = count_dict[alt]
        qalt_dict = q_dicts[alt]
        l_germline = genotype_likelihood(ref_count, alt_count, qref_dict, qalt_dict)

        prior_refhom = f[ref] ** 2
        prior_het = 2 * f[ref] * f[alt]
        prior_althom = f[alt] ** 2
        prior_mosaic = mu
        prior = np.array([prior_refhom, prior_althom, prior_het, prior_mosaic])
        p_germline = prior * l_germline
        p_germline = p_germline / sum(p_germline)

        germline_genotype_list = ["refhom", "althom", "het", "mosaic"]
        germline_dict = dict(zip(germline_genotype_list, p_germline))
        germline_geno = max(germline_dict, key=lambda x: germline_dict[x])
        if germline_geno == "mosaic":
            germline_geno = "refhom"

        if n_allele == 2 and marker_withref:
            if germline_geno == "refhom":
                germline_allele = ref
                germ_count = str(count_dict_original[ref])
                germ_prior = str(f[ref])
            elif germline_geno == "het" or germline_geno == "althom":
                germline_geno = "het"
                germline_allele = ",".join([ref, alt])
                germ_count = ",".join([str(count_dict_original[i]) for i in [ref, alt]])
                germ_prior = ",".join([str(f[i]) for i in [ref, alt]])
        else:
            if germline_geno == "refhom":
                ref_list = [ref]
            else:
                ref_list = [ref, alt]
                allele_list.remove(alt)
            alt_list = allele_list
            ind_ref = ",".join(ref_list)
            germline_allele = ind_ref
            germ_count = ",".join([str(count_dict_original[i]) for i in ref_list])
            germ_prior = ",".join([str(f[i]) for i in ref_list])

    res_list = []
    if len(ref_list) > 1:
        potential_geno = "refhom"
    else:
        potential_geno = germline_geno

    for alt in alt_list:
        alt_count = count_dict[alt]
        qalt_dict = q_dicts[alt]
        prior_althom = f[alt] ** 2
        prior_mosaic = mu
        p = np.array([0] * 4)

        for ref in ref_list:
            prior_refhom = f[ref] ** 2
            prior_het = 2 * f[ref] * f[alt]
            ref_count = count_dict[ref]
            qref_dict = q_dicts[ref]
            l_update = genotype_likelihood(ref_count, alt_count, qref_dict, qalt_dict)
            prior = np.array([prior_refhom, prior_althom, prior_het, prior_mosaic])
            p_update = prior * l_update
            p = np.array(list(map(add, p, p_update)))

        genotype_list = ["refhom", "althom", "het", "mosaic"]
        p_dict = dict(zip(genotype_list, p))
        p_ind = [p_dict[potential_geno], p_dict["mosaic"]]
        p_ind = p_ind / sum(p_ind)

        indgeno_list = ["germline", "mosaic"]
        pind_dict = dict(zip(indgeno_list, p_ind))
        geno_max = max(pind_dict, key=lambda x: pind_dict[x])
        p_mosaic = pind_dict["mosaic"]

        G_list = ["0/0", "1/1", "0/1", "0/1"]
        G_dict = dict(zip(genotype_list, G_list))
        geno_ind = germline_geno if geno_max == "germline" else "mosaic"
        G_max = G_dict[geno_ind]

        depth = sum(count_dict.values())
        vaf = alt_count / depth

        res = [geno_max, p_mosaic, G_max, ind_ref, alt, vaf]
        res_list.append(res)

    germline_list = [germline_geno, germline_allele, germ_count, germ_prior]
    return res_list, germline_list


def quality_format(qA, qT, qC, qG, thr_dp=1000):
    nucleotide_list = ["A", "T", "C", "G"]
    q_list = [str2dict(q) for q in [qA, qT, qC, qG]]
    count_list = [sum(q.values()) for q in q_list]
    depth = sum(count_list)
    count_dict_original = dict(zip(nucleotide_list, count_list))

    if depth > thr_dp:
        scaling_factor = thr_dp / depth
        downsampled_q_list = [
            {key: int(round(value * scaling_factor)) for key, value in q.items()}
            for q in q_list
        ]
        q_list = downsampled_q_list

    q_dicts = dict(zip(nucleotide_list, q_list))
    count_list = [sum(q.values()) for q in q_list]
    count_dict = dict(zip(nucleotide_list, count_list))
    return count_dict, q_dicts, count_dict_original


def genotype_likelihood(ref_count, alt_count, qref_dict, qalt_dict):
    depth = ref_count + alt_count

    l_het = math.log10(comb(depth, alt_count, exact=True)) + math.log10(0.5) * depth
    l_het = 10 ** l_het

    q_refhom = math.log10(1)
    if ref_count != 0:
        q_refhom = q_refhom + sum(math.log10(1 - 0.1 ** (i / 10)) * qref_dict[i] for i in qref_dict.keys())
    if alt_count != 0:
        q_refhom = q_refhom + sum(math.log10(0.1 ** (i / 10)) * qalt_dict[i] for i in qalt_dict.keys())
    l_refhom = math.log10(comb(depth, alt_count, exact=True)) + q_refhom
    l_refhom = 10 ** l_refhom

    q_althom = math.log10(1)
    if ref_count != 0:
        q_althom = q_althom + sum(math.log10(0.1 ** (i / 10)) * qref_dict[i] for i in qref_dict.keys())
    if alt_count != 0:
        q_althom = q_althom + sum(math.log10(1 - 0.1 ** (i / 10)) * qalt_dict[i] for i in qalt_dict.keys())
    l_althom = math.log10(comb(depth, alt_count, exact=True)) + q_althom
    l_althom = 10 ** l_althom

    r = 0
    if ref_count != 0:
        r = r + sum([0.1 ** (float(i) / 10) * qref_dict[i] for i in qref_dict.keys()])
    if alt_count != 0:
        r = r + sum([(1 - 0.1 ** (float(i) / 10)) * qalt_dict[i] for i in qalt_dict.keys()])
    l_mosaic = beta(r + 1, depth - r + 1)
    if l_mosaic > 0:
        l_mosaic = math.log10(l_mosaic) + math.log10(comb(depth, alt_count, exact=True))
        l_mosaic = 10 ** l_mosaic

    l = [l_refhom, l_althom, l_het, l_mosaic]
    return l


def individual_genotype(row, mu=1e-7, thr_dp=1000, pop_vaf=1e-5, filter_oneallele=True):
    ref = row['ref']
    alt = row['alt']
    qA, qT, qC, qG = row[['qA', 'qT', 'qC', 'qG']]

    fA = float(row['fA'])
    fT = float(row['fT'])
    fC = float(row['fC'])
    fG = float(row['fG'])
    prior_string = ",".join([str(i) for i in [fA, fT, fC, fG]])

    res_list, germline_list = individual_posterior(
        ref, alt, qA, qT, qC, qG, fA, fT, fC, fG,
        mu=mu, thr_dp=thr_dp, pop_vaf=pop_vaf,
        filter_oneallele=filter_oneallele
    )

    output_list = []
    write_germline = False
    for res in res_list:
        [geno_max, p_mosaic, G_max, ind_ref, alt, vaf] = res
        if geno_max != "NA":
            output = [
                row['chrom'], row['pos'], row['strand'],
                ind_ref, alt,
                row['cluster'], row['spot_number'], row['consensus_read_count'], prior_string,
                geno_max, p_mosaic, G_max, vaf
            ]
            output_list.append(output)
            write_germline = True

    if write_germline:
        germline_output = [row['chrom'], row['pos']] + germline_list
    else:
        germline_output = []

    return output_list, germline_output


class IndGenoCalculator:
    def _load_count_file(self, file):
        df = pd.read_csv(file, sep="\t", header=None, comment="#", keep_default_na=False)
        columns = ['chrom', 'pos', 'strand', 'ref', 'alt', 'cluster', 'spot_number',
                   'consensus_read_count', 'qA', 'qT', 'qC', 'qG']
        df.columns = columns
        if '#chrom' in df.columns and 'chrom' not in df.columns:
            df = df.rename(columns={'#chrom': 'chrom'})
        df['identifier'] = df['chrom'] + "_" + df['pos'].astype(str)
        return df

    def _merge_inputs_for_ind_geno_df(self, ind_count_df, prior_df):
        df = ind_count_df.copy()
        if '#chrom' in df.columns and 'chrom' not in df.columns:
            df = df.rename(columns={'#chrom': 'chrom'})
        if 'identifier' not in df.columns:
            df['identifier'] = df['chrom'] + "_" + df['pos'].astype(str)

        merged = df.merge(
            prior_df[['identifier', 'fA', 'fT', 'fC', 'fG']],
            on='identifier',
            how='left'
        )
        merged.fillna(0, inplace=True)
        return merged

    def _merge_inputs_for_ind_geno(self, ind_count_file, prior_df):
        ind_count_df = self._load_count_file(ind_count_file)
        return self._merge_inputs_for_ind_geno_df(ind_count_df, prior_df)

    def calculate_individual_genotype_df(
        self,
        ind_count_df,
        prior_df,
        geno_file=None,
        geno_filter_file=None,
        germline_file=None,
        mu=1e-7,
        thr_dp=1000,
        pop_vaf=1e-5,
        filter_oneallele=True
    ):
        df = self._merge_inputs_for_ind_geno_df(ind_count_df, prior_df)

        results = df.apply(
            lambda row: individual_genotype(
                row, mu, thr_dp, pop_vaf, filter_oneallele
            ),
            axis=1
        )

        del df
        gc.collect()

        all_genotypes = []
        all_germlines = []

        for genotype_list, germline in results:
            if genotype_list:
                for geno in genotype_list:
                    all_genotypes.append(geno)
            if germline:
                all_germlines.append(germline)

        geno_columns = ['#chrom', 'pos', 'strand', 'germline', 'mutant',
                        'cluster', 'spot_number', 'consensus_read_count', 'prior_ATCG',
                        'genotype', 'p_mosaic', 'Gi', 'vaf']

        germ_columns = ['#chrom', 'pos', 'germline_geno', 'germline_allele',
                        'germ_count', 'germ_prior']

        geno_df = pd.DataFrame(all_genotypes, columns=geno_columns)
        germ_df = pd.DataFrame(all_germlines, columns=germ_columns)

        if geno_file is not None:
            geno_df.to_csv(geno_file, sep="\t", index=None, na_rep='NA')
        if germline_file is not None:
            germ_df.to_csv(germline_file, sep="\t", index=None, na_rep='NA')

        if geno_df.empty:
            geno_filter_df = pd.DataFrame(columns=geno_columns)
        else:
            geno_df = geno_df.copy()
            geno_df['spot_number'] = geno_df['spot_number'].astype(int)
            geno_filter_df = geno_df[
                (geno_df['spot_number'] >= 30) &
                (geno_df['genotype'] == 'mosaic')
            ].copy()

            if not geno_filter_df.empty:
                geno_filter_df['vaf'] = geno_filter_df['vaf'].astype(float)
                geno_filter_df = geno_filter_df[geno_filter_df["vaf"] <= 0.3]

        if geno_filter_file is not None:
            geno_filter_df.to_csv(geno_filter_file, sep="\t", index=None, na_rep='NA')

            new_file = pd.DataFrame({
                'mutation_id': geno_filter_df['#chrom'].astype(str) + '_' +
                               geno_filter_df['pos'].astype(str) + '_' +
                               geno_filter_df['germline'] + '_' +
                               geno_filter_df['mutant']
            })
            new_file.to_csv(geno_filter_file + ".mutation.list", sep='\t', index=False, header=False)

        return geno_df, geno_filter_df, germ_df

    def calculate_individual_genotype(
        self,
        ind_count_file,
        prior_df,
        geno_file,
        geno_filter_file,
        germline_file,
        mu,
        thr_dp,
        pop_vaf,
        filter_oneallele
    ):
        ind_count_df = self._load_count_file(ind_count_file)
        return self.calculate_individual_genotype_df(
            ind_count_df=ind_count_df,
            prior_df=prior_df,
            geno_file=geno_file,
            geno_filter_file=geno_filter_file,
            germline_file=germline_file,
            mu=mu,
            thr_dp=thr_dp,
            pop_vaf=pop_vaf,
            filter_oneallele=filter_oneallele
        )


def ClusterVAFCalculator_from_df(ind_geno_df, cluster_count_df, outfile_path=None):
    if ind_geno_df.empty:
        result = pd.DataFrame(columns=[
            '#chrom', 'pos', 'strand', 'germline', 'mutant',
            'cluster', 'spot_number', 'consensus_read_count', 'vaf'
        ])
        if outfile_path is not None:
            result.to_csv(outfile_path, sep='\t', index=False)
        return result

    cluster_df = cluster_count_df.copy()
    if 'chrom' in cluster_df.columns and '#chrom' not in cluster_df.columns:
        cluster_df = cluster_df.rename(columns={'chrom': '#chrom'})

    ind_df = ind_geno_df.copy()

    ind_df['key'] = ind_df['#chrom'].astype(str) + "_" + ind_df['pos'].astype(str)
    cluster_df['key'] = cluster_df['#chrom'].astype(str) + "_" + cluster_df['pos'].astype(str)

    merged = pd.merge(cluster_df, ind_df[['key', 'germline', 'mutant']], on='key')
    counts = merged['consensus_read_count'].str.split(',', expand=True).astype(int)
    counts.columns = ['A_count', 'T_count', 'C_count', 'G_count']

    merged = pd.concat([merged, counts], axis=1)
    merged['depth'] = merged[['A_count', 'T_count', 'C_count', 'G_count']].sum(axis=1)

    mutant_map = {'A': 0, 'T': 1, 'C': 2, 'G': 3}
    col_indices = merged['mutant'].map(mutant_map).values
    counts_matrix = merged[['A_count', 'T_count', 'C_count', 'G_count']].values
    merged['mutant_count'] = counts_matrix[np.arange(len(merged)), col_indices]

    merged['vaf'] = np.where(
        merged['depth'] > 0,
        merged['mutant_count'] / merged['depth'],
        "NA"
    )

    result = merged[[
        '#chrom', 'pos', 'strand', 'germline', 'mutant',
        'cluster', 'spot_number', 'consensus_read_count', 'vaf'
    ]]

    if outfile_path is not None:
        result.to_csv(outfile_path, sep='\t', index=False)

    return result


def ClusterVAFCalculator(ind_geno_file, cluster_count_file, outfile_path):
    ind_geno_df = pd.read_csv(ind_geno_file, sep="\t")
    cluster_count_df = pd.read_csv(cluster_count_file, sep="\t")
    return ClusterVAFCalculator_from_df(ind_geno_df, cluster_count_df, outfile_path=outfile_path)

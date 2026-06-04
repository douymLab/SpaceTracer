
import os
import pyranges as pr
import pandas as pd
from scipy.stats import binom,binomtest
import statsmodels.stats.multitest as smm
import pandas as pd
import numpy as np

def binomial_test(row):
    h1, h2 = map(int, row['count'].split(','))
    p_value = binomtest(h1, h1 + h2, p=0.5, alternative='two-sided').pvalue
    return p_value


def read_and_filter_germline_file(germline_file, count_threshold, prior_threshold, p_threshold):
    germline_sites = pd.read_csv(germline_file, sep="\t", header=None,
                                  names=["#chrom", "site", "genotype", "allele", "count", "prior"])

    germline_sites = germline_sites[germline_sites["genotype"] == "het"]
    if germline_sites.empty:
        return pd.DataFrame()

    counts = germline_sites['count'].str.split(',', expand=True).astype(int)
    counts.columns = ['ref_count', 'alt_count']
    
    priors = germline_sites['prior'].str.split(',', expand=True).astype(float)
    priors.columns = ['ref_prior', 'alt_prior']
    
    germline_sites = pd.concat([germline_sites, counts, priors], axis=1)
    
    germline_sites['count_sum'] = counts.sum(axis=1)
    germline_sites['prior_valid'] = (priors > prior_threshold).all(axis=1)
    
    filtered_df = germline_sites[(germline_sites['count_sum'] > count_threshold) & (germline_sites['prior_valid'])]
    if filtered_df.empty:
        return filtered_df
    
    else:
        filtered_df['binomial_p_value'] = filtered_df.apply(binomial_test, axis=1)
        filtered_df = filtered_df[filtered_df['binomial_p_value'] <= p_threshold]
        
        if filtered_df.empty:
            return filtered_df
        else:
            filtered_df = filtered_df.rename(columns={"#chrom": "chrom", "site": "pos"})
            split_result = filtered_df['allele'].str.split(',', expand=True)
            filtered_df[['ref', 'alt']] = filtered_df['allele'].str.split(',', expand=True)

    return filtered_df


def annotate_dbsnp_tabix(filtered_df, dbsnp_vcf_gz):
    import pysam
    filtered_df = filtered_df.copy()
    tb = pysam.TabixFile(dbsnp_vcf_gz)

    ids = []
    for chrom, pos, ref, alt in filtered_df[["chrom","pos","ref","alt"]].itertuples(index=False, name=None):
        chrom = str(chrom)
        pos = int(pos)
        hit = "NA"
        try:
            # tabix is 0-based half-open；VCF POS is 1-based
            for rec in tb.fetch(chrom, pos-1, pos):
                if rec.startswith("#"):
                    continue
                fields = rec.split("\t")
                v_chrom, v_pos, v_id, v_ref, v_alt = fields[0], int(fields[1]), fields[2], fields[3], fields[4]
                if v_pos == pos and v_ref == ref and alt in v_alt.split(","):
                    hit = v_id
                    break
        except ValueError:
            hit = "NA"
        ids.append(hit)

    filtered_df["dbSNPID"] = ids
    return filtered_df


def annotate_dbsnp(filtered_df, dbsnp_vcf_file):
    #read dbSNP VCF
    dbsnp_cols = ['chrom', 'pos', 'dbSNPID', 'ref', 'alt', 'qual', 'filter', 'info']
    dbsnp_df = pd.read_csv(dbsnp_vcf_file, 
                           sep='\t', 
                           comment='#',
                           names=dbsnp_cols,
                           usecols=[0, 1, 2, 3, 4])  
    
    # build a key to speed up merging
    filtered_df['pos_key'] = filtered_df['pos'].astype(str)
    filtered_df['ref_key'] = filtered_df['ref'].astype(str)
    filtered_df['alt_key'] = filtered_df['alt'].astype(str)
    filtered_df['match_key'] = (filtered_df['chrom'].astype(str) + ':' + 
                                filtered_df['pos_key'] + ':' + 
                                filtered_df['ref_key'] + ':' + 
                                filtered_df['alt_key'])
    
    dbsnp_df['match_key'] = (dbsnp_df['chrom'].astype(str) + ':' + 
                             dbsnp_df['pos'].astype(str) + ':' + 
                             dbsnp_df['ref'].astype(str) + ':' + 
                             dbsnp_df['alt'].astype(str))
    
    dbsnp_dict = dbsnp_df.set_index('match_key')['dbSNPID'].to_dict()
    
    filtered_df['dbSNPID'] = filtered_df['match_key'].map(dbsnp_dict).fillna('NA')
    
    filtered_df = filtered_df.drop(['pos_key', 'ref_key', 'alt_key', 'match_key'], axis=1)
    
    return filtered_df

def annotate_genes(filtered_df, gene_bed, default_range=150):
    filtered_df = filtered_df.copy()

    filtered_df["chrom"] = filtered_df["chrom"].astype(str)
    filtered_df["pos"] = pd.to_numeric(filtered_df["pos"], errors="coerce")
    filtered_df = filtered_df.dropna(subset=["pos"])
    filtered_df["pos"] = filtered_df["pos"].astype(int)

    if not gene_bed:
        filtered_df["gene"] = "NA"
        filtered_df["region_chrom"] = filtered_df["chrom"]
        filtered_df["region_start"] = filtered_df["pos"] - default_range
        filtered_df["region_end"] = filtered_df["pos"] + default_range
        return filtered_df

    try:
        genes = pr.read_bed(gene_bed)
        if len(genes) == 0:
            filtered_df["gene"] = "NA"
            filtered_df["region_chrom"] = filtered_df["chrom"]
            filtered_df["region_start"] = filtered_df["pos"] - default_range
            filtered_df["region_end"] = filtered_df["pos"] + default_range
            return filtered_df

        sites_pr = pr.PyRanges(
            chromosomes=filtered_df["chrom"].values,
            starts=filtered_df["pos"].values,
            ends=(filtered_df["pos"] + 1).values,
        )
        sites_pr.Index = filtered_df.index.values

        intersect = sites_pr.join(genes, how="left", apply_strand_suffix=False)

        filtered_df["gene"] = "NA"
        filtered_df["region_start"] = filtered_df["pos"] - default_range
        filtered_df["region_end"] = filtered_df["pos"] + default_range
        filtered_df["region_chrom"] = filtered_df["chrom"]

        if len(intersect) > 0:
            dfj = intersect.df

            if "Name" in dfj.columns:
                dfj = dfj.dropna(subset=["Name"])

            if len(dfj) > 0:
                gene_info = (
                    dfj.groupby("Index")
                    .agg({
                        "Name": lambda x: ",".join(map(str, x)),
                        "Start": "min",
                        "End": "max",
                        "Chromosome": "first",
                    })
                    .reset_index()
                )

                for _, row in gene_info.iterrows():
                    idx = row["Index"]
                    filtered_df.loc[idx, "gene"] = row["Name"]
                    filtered_df.loc[idx, "region_start"] = int(row["Start"])
                    filtered_df.loc[idx, "region_end"] = int(row["End"])
                    filtered_df.loc[idx, "region_chrom"] = row["Chromosome"]

    except Exception as e:
        filtered_df["gene"] = "NA"
        filtered_df["region_chrom"] = filtered_df["chrom"]
        filtered_df["region_start"] = filtered_df["pos"] - default_range
        filtered_df["region_end"] = filtered_df["pos"] + default_range

    return filtered_df

# main get ase
def get_ase_germline_sites(germline_file,dbsnp_vcf_file,ase_germline_file,gene_bed="",
                            count_threshold=50,prior_threshold=0.0001,p_threshold=0.05,default_range=150):
    # step1: read germline file
    filtered_df=read_and_filter_germline_file(germline_file,count_threshold,prior_threshold,p_threshold)

    if filtered_df.empty:
        return pd.DataFrame()
    # step2: add dbSNP info
    filtered_df=annotate_dbsnp_tabix(filtered_df, dbsnp_vcf_file)
    
    # step3: add gene region info
    filtered_df=annotate_genes(filtered_df, gene_bed,default_range)

    short_df=filtered_df[['chrom', 'pos', 'ref', 'alt','ref_count', 'alt_count', 'ref_prior', 'alt_prior','binomial_p_value','dbSNPID','gene','region_chrom','region_start','region_end']]
    # short_df.to_csv(ase_germline_file,sep="\t",index=False)
    return short_df


def intersect_somatic_with_ase(somatic_df, ase_germline_df,p_threshold=0.05):
    ase_series = pd.Series('Unknown', index=somatic_df.index)

    if ase_germline_df.empty:
        return ase_series 

    somatic_pr = pr.PyRanges(pd.DataFrame({
        'Chromosome': somatic_df['chrom'],
        'Start': somatic_df['pos'] - 1,
        'End': somatic_df['pos'],
        'somatic_alt_count': somatic_df['alt_count'],
        'somatic_total': somatic_df['count']
    }))

    ase_pr = pr.PyRanges(pd.DataFrame({
        'Chromosome': ase_germline_df['region_chrom'],
        'Start': ase_germline_df['region_start'],
        'End': ase_germline_df['region_end'],
        'germ_ref_count': ase_germline_df['ref_count'],
        'germ_alt_count': ase_germline_df['alt_count'],
    }))

    intersect = somatic_pr.join(ase_pr, how=None)

    # result_df['ase_direction_consistent'] = np.nan
    if len(intersect.df) == 0:
        return ase_series

    intersect_df = intersect.df.copy()
    intersect_df['germ_total'] = intersect_df['germ_ref_count'] + intersect_df['germ_alt_count']
    intersect_df['germ_vaf'] = intersect_df['germ_ref_count'] / intersect_df['germ_total']
    # intersect_df['somatic_key'] = intersect_df['Chromosome'] + '_' + intersect_df['End'].astype(str)
    intersect_df['somatic_key'] = intersect_df['Chromosome'].astype(str) + '_' + intersect_df['End'].astype(str)

    somatic_key  = somatic_df['chrom'] + '_' + somatic_df['pos'].astype(str)


    for key, group in intersect_df.groupby('somatic_key'):
        mask = somatic_key == key
        if mask.sum() == 0:
            continue

        somatic_row = somatic_df[mask].iloc[0]
        alt_count = int(somatic_row['alt_count'])
        count_sum = int(somatic_row['count'])

        ase_list = []
        # directions = []

        for _, row in group.iterrows():
            germ_vaf = row['germ_vaf']
            p = min(germ_vaf, 1 - germ_vaf)
            mosaic_p = binom.cdf(alt_count, count_sum, p=p)
            ase_list.append(mosaic_p)

            # if germ_vaf > 0.65:
            #     directions.append('ref_bias')
            # elif germ_vaf < 0.35:
            #     directions.append('alt_bias')
            # else:
            #     directions.append('balance')

        if len(ase_list) == 1:
            p_adj = ase_list
        else:
            _, p_adj, _, _ = smm.multipletests(ase_list, method='fdr_bh')

        # ase = 'True' if max(p_adj) >= p_threshold else 'False'
        ase_series[mask] = 'True' if max(p_adj) >= p_threshold else 'False'
       
    return ase_series

 
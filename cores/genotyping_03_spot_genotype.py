
######## spot genotype

import math
from pathlib import Path
from typing import List

import numpy as np
import pandas as pd
from operator import add
from scipy.special import comb

from SpaceTracer.utils.read_files import load_spot_genotypes_data,load_spot_count_data
from SpaceTracer.utils.utils import str2dict


class SpotGenoCalculator:
    def __init__(self,bins,epsQ,thr_dp,pop_vaf,cell_num):
        self.bins=bins
        self.epsQ=epsQ
        self.thr_dp=thr_dp
        self.pop_vaf=pop_vaf
        self.cell_num=cell_num
    
    def run(self,spot_count_file,ind_geno_file,cluster_df,cluster_vaf_file,output_file):
        ind_geno_df=pd.read_csv(ind_geno_file,sep="\t")
        raw_spot_count_df = load_spot_count_data(spot_count_file)
        raw_spot_count_df = raw_spot_count_df.fillna('NA')
        cluster_vaf_df=pd.read_csv(cluster_vaf_file,sep="\t")
        cluster_vaf_df = cluster_vaf_df.drop(["spot_number","consensus_read_count"], axis=1)

        if self.bins is None:
        # read in inputs
            ind_geno_df = ind_geno_df.drop(["cluster","spot_number","consensus_read_count","prior_ATCG"], axis=1)
            spot_count_df=raw_spot_count_df
        
        else:
            bins=self.bins
            ind_geno_df = ind_geno_df.drop(["cluster","spot_number","consensus_read_count","prior_ATCG"], axis=1)
            spot_count_df=bin_df(raw_spot_count_df,bins)

        del raw_spot_count_df

        # join the data
        count_cluster_df = pd.merge(spot_count_df, cluster_df, on='barcode', how='left')
        count_geno_join = pd.merge(count_cluster_df, ind_geno_df, on=['#chrom','pos','strand'], how='inner')
        count_geno_join = count_geno_join.rename(columns={'vaf': 'ind_vaf'})
        count_geno_vaf_join = pd.merge(count_geno_join, cluster_vaf_df, on=['#chrom','pos','strand','germline','mutant','cluster'], how='left')
        count_geno_vaf_join = count_geno_vaf_join.rename(columns={'vaf': 'cluster_vaf'})

        # input the cell number for each spot
        if not isinstance(self.cell_num, int) and Path(self.cell_num).exists():
            cell_num_file = pd.read_csv(self.cell_num,sep="\t",header=0)
            cell_num_file = cell_num_file.drop(["cluster", "nUMI", "nREAD"], axis=1)
            count_geno_vaf_join = pd.merge(count_geno_vaf_join, cell_num_file, on='barcode', how='left')

        out_colname = "\t".join(["#chrom","pos","strand","germline","mutant","cluster","spot_barcode","consensus_read_count", \
                                "l_germline", "l_mosaic", "max_spot_geno","G_spot_max","depth","vaf","p_mosaic"])
        mosaic_mask = count_geno_vaf_join.iloc[:, 14] == "mosaic"
        mosaic_df = count_geno_vaf_join[mosaic_mask]

        data = mosaic_df.values
        if not isinstance(self.cell_num, int) and Path(self.cell_num).exists():
            cell_nums = data[:, -1] 
            data = data[:, :-1]     
        else:
            cell_nums = np.full(len(data), self.cell_num)
        # print("****++++***+++***",output_file)
        with open(output_file, 'w') as f:
            f.write(out_colname+"\n")
            for i in range(len(data)):
                spot_geno_info = spot_genotype(
                    join_info=[str(x) for x in data[i]], 
                    cell_num=cell_nums[i],
                    epsQ=self.epsQ,
                    thr_dp=self.thr_dp, 
                    pop_vaf=self.pop_vaf
                )
                if spot_geno_info[10] != "NA":
                    f.write("\t".join(map(str, spot_geno_info)) + "\n")
        # Read the TSV we just wrote and convert to Parquet with index
        df = load_spot_genotypes_data(output_file) 
        # print(df)
        parquet_file = str(output_file).replace('.out', '.parquet')
        df.to_parquet(parquet_file, index=True, compression='snappy')
        



def bin_df(df,bins):
    def round_to_nearest_bin(x):
        return int(np.ceil(x / bins) * bins)

    df[['x', 'y']] = df['spot_barcode'].str.split('_', expand=True)
    df['x'] = df['x'].astype(int)
    df['y'] = df['y'].astype(int)
    df['new_x'] = df['x'].apply(round_to_nearest_bin).astype(str)
    df['new_y'] = df['y'].apply(round_to_nearest_bin).astype(str)
    df['spot_barcode']=df['new_x']+"_"+df['new_y']
    df.drop(columns=['x', 'y', 'new_x', 'new_y'], inplace=True)
    return df


def spot_posterior(germline, mutant, cluster_vaf, qA, qT, qC, qG, cell_num=20, epsQ=20, thr_dp=1000, pop_vaf=1e-5):
    """
    Calculate the posterior probability of the four spot genotypes:
    refhom, althom, het amd mosaic

    Note: this function also treat following two situations:
    - if the spot does not contain any mutant allele but has high mosaic posterior(>0.5), 
      then mosaic_posterior=0.5, germline_likelihood=mosaic_likelihood=0.5
    - if the spot has mutant allele but the mosaic likelihood is low (<0.5),
      then germline_likelihood=mosaic_likelihood=0.5

    Inputs:
        germline, mutant - the germline alleles and mutant allele on the individual level
        cluster_vaf - the allele frequency (AF) of the mutant (alternative) allele on the cluster level
        qA, qT, qC, qG - the consensus read qualities for each nucleotide
        cell_num - the number of cell of the spot (default=20)
        epsQ - the threshold for the consensus read quality (default=20)
        thr_dp - the threshold for the total depth (default=1000)
                 when depth > thr_dp, downsample the allele numbers to let the depth = thr_dp
        pop_vaf - the probability of mutant allele in population (default=1e-5)

    Outputs: p, l_norm, [max_spot_geno, G_spot_max, depth, vaf, p_mosaic]
        p - the normalized array of the posterior probabilities for all spot genotypes
        l_norm - the normalized array of the likelihoods for all spot genotypes
        max_spot_geno - the spot genotype with the highest posterior probability
        G_spot_max - represent the spot genotype by 0 and 1
                     i.e. refhom(0/0), althom(1/1), het(0/1), mosaic(0/1)
        depth - the UMI count depth for the spot
        vaf - the mutant allele frequency
        p_mosaic - the probability of having mosaic mutation of the spot (mosaic posterior probability)
    """
    # delete the low quality consensus reads and the low AF alternative alleles
    # get the consensus read and quality dictionaries
    count_filter, q_filter = spot_filter(qA, qT, qC, qG, epsQ=epsQ, thr_dp=thr_dp)
    count_filter_nozero = {k:v for k,v  in count_filter.items() if v>0}
    allele_list_spot = list(count_filter_nozero.keys())
    depth = sum(count_filter.values())

    # get the spot reference allele list
    ind_ref_list = germline.split(",")
    ref_list = [ref for ref in ind_ref_list if ref in allele_list_spot]
    alt = mutant
    # get the spot allele list except the errors
    if mutant != ".":
        ind_allele_list = ind_ref_list + [mutant]
    else:
        ind_allele_list = ind_ref_list
    allele_list = [allele for allele in allele_list_spot if allele in ind_allele_list]
    n_allele = len(allele_list)

    # return NAs if the candidate alleles are none
    if n_allele == 0:
        return "NA", ["NA"]*2, ["NA", "NA", depth, "NA", "NA"]
    # set fake ref allele if no actual ref allele
    elif len(ref_list) == 0:
        ref_list = [ind_ref_list[0]]

    # get the alternative allele info
    genotype_list = ["germline", "mosaic"]
    alt_count = count_filter[alt]
    qalt_dict = q_filter[alt]
    # if no alternative allele, give alt allele count as 0
    if alt == ".":
        alt_count = 0

    ## find the highest genotype
    # calcualte the likelihoods for each genotype
    # and combine the likelihoods for each reference candidate
    l = np.array([0]*2)  # initial likelihoods
    for ref in ref_list:
        # likelihood
        ref_count = count_filter[ref]
        qref_dict = q_filter[ref]
        l_update = spot_likelihood(ref_count, alt_count, qref_dict, qalt_dict, cluster_vaf, pop_vaf)
        l = np.array(list(map(add, l, l_update)))

    # normalized likelihood
    s_likelihood = sum(l)
    s_likelihood = 1e-15 if s_likelihood==0 else s_likelihood # avoid all 0
    l_norm = [j/s_likelihood for j in l]
    
    ## prior probability
    if cluster_vaf < 0.5:
        prior_value = [(1-2*cluster_vaf)**cell_num, 1-(1-2*cluster_vaf)**cell_num]
    else:
        prior_value = [pop_vaf, 1-pop_vaf]
    # prior_value = [[1-2*mu, 0, 2*mu, 0], [0, 1-2*mu, 2*mu, 0], [0, 0, 1, 0], [(1-2*cluster_vaf)**cell_num, 0, 0, 1-(1-2*cluster_vaf)**cell_num]]
    # prior_dict = dict(zip(genotype_list, prior_value))
    
    ## calculate the posterior values
    # multiply prior
    posterior = prior_value * l
    # normalise
    s_posterior = sum(posterior)
    s_posterior = 1e-15 if s_posterior==0 else s_posterior # avoid all 0
    p = [j/s_posterior for j in posterior]
    # transfer all nan probability to 0 to avoid 0 ind genotype
    p = np.nan_to_num(p)
    posterior_dict = dict(zip(genotype_list, p))
    ## find the highest genotype
    max_spot_geno = max(posterior_dict, key= lambda x: posterior_dict[x])

    ## probability of mosaic mutation
    p_mosaic = posterior_dict['mosaic']
    ## calculate the mutant allele frequency at the spot level
    vaf = alt_count / depth
    ## represent it by 0 and 1
    G_list = ["0/0", "0/1"]
    G_dict = dict(zip(genotype_list, G_list))
    G_spot_max = G_dict[max_spot_geno]

    # balance the conflict situation
    if alt_count == 0:
        if p_mosaic > 0.5:
            p_mosaic = 0.5
        if l_norm[1] > l_norm[0]:
            l_norm = [0.5, 0.5]
    elif alt_count > 0 and l_norm[1] < 0.5:
        l_norm = [0.5, 0.5]
    
    ## output
    return p, l_norm, [max_spot_geno, G_spot_max, depth, vaf, p_mosaic]



def spot_filter(qA, qT, qC, qG, epsQ=20, thr_dp=1000):
    """
    Delete the low quality (quality < epsQ) consensus reads
    Also downsample the allele numbers if the depth is very high

    Outputs:
        count_dict: a dictionary for the number of consensus reads after filtering
        q_filter: a dictionary for the quality of consensus reads after filtering
    """
    nucleotide_list = ["A", "T", "C", "G"]
    # transfrom the string to dict
    q_list = [str2dict(q) for q in [qA, qT, qC, qG]]
    # delete the low quality consensus reads
    q_filter_list = [{k:v for k,v in q.items() if k>=epsQ} for q in q_list]

    # downsample if the depth is too high
    depth = sum([sum(q.values()) for q in q_filter_list])
    if depth > thr_dp:
        # scale for each quality dict
        scaling_factor = thr_dp / depth
        downsampled_q_list = [{key: int(round(value * scaling_factor)) for key, value in q.items()} for q in q_filter_list]
        q_filter_list = downsampled_q_list
    q_filter = dict(zip(nucleotide_list, q_filter_list))
    
    # calculate the numbers of the nucleotides
    count_list = [sum(q.values()) for q in q_filter_list]
    count_dict = dict(zip(nucleotide_list, count_list))

    return count_dict, q_filter



def spot_likelihood(ref_count, alt_count, qref_dict, qalt_dict, cluster_vaf, pop_vaf=1e-5):
    """
    Calculate the four genotype likelihoods with one ref allele for the spots
    genotype: refhom, althom, het amd mosaic
    here only outputs two likelihoods: refhom and mosaic, 
    cause we could see these two values as the likelihoods for germline genotype or mosaic mutation

    Input:
        ref_count, alt_count - the number of the consensus reads for the reference and alternative alleles
        qref_dict, qalt_dict - the dictionaries for the qualities of the consensus reads on the ref and alt alleles
        cluster_vaf - the allele frequency (AF) of the mutant (alternative) allele on the cluster level
        pop_vaf - the probability of mutant allele in population (default=1e-5)

    Output:
        the likelihood array for the four / two genotypes
    """
    ## depth
    depth = ref_count + alt_count
    
    # ## het likelihood
    # l_het = math.log10(comb(depth,alt_count,exact=True)) + math.log10(0.5)*depth
    # l_het = 10 ** l_het

    ## refhom likelihood
    q_refhom = math.log10(1)
    if ref_count != 0:
        q_refhom = q_refhom + sum(math.log10(1-0.1**(i/10)) * qref_dict[i] for i in qref_dict.keys())
    if alt_count != 0:
        q_refhom = q_refhom + sum(math.log10(0.1**(i/10)) * qalt_dict[i] for i in qalt_dict.keys())
    l_refhom = math.log10(comb(depth,alt_count,exact=True)) + q_refhom
    l_refhom = 10 ** l_refhom

    # ## althom likelihood
    # q_althom = math.log10(1)
    # if ref_count != 0:
    #     q_althom = q_althom + sum(math.log10(0.1**(i/10)) * qref_dict[i] for i in qref_dict.keys())
    # if alt_count != 0:
    #     q_althom = q_althom + sum(math.log10(1-0.1**(i/10)) * qalt_dict[i] for i in qalt_dict.keys())
    # l_althom = math.log10(comb(depth,alt_count,exact=True)) + q_althom
    # l_althom = 10 ** l_althom
    
    ## mosaic likelihood
    # avoid cluster_vaf=0/1
    cluster_vaf = pop_vaf if cluster_vaf==0 else cluster_vaf
    cluster_vaf = 1-pop_vaf if cluster_vaf==1 else cluster_vaf
    # calculate
    r = 0
    if ref_count != 0:
        r = r + sum([0.1**(float(i)/10) * qref_dict[i] for i in qref_dict.keys()])
    if alt_count != 0:
        r = r + sum([(1-0.1**(float(i)/10)) * qalt_dict[i] for i in qalt_dict.keys()])
    l_mosaic = math.log10(comb(depth,alt_count,exact=True)) + r*math.log10(cluster_vaf) + (depth-r)*math.log10(1-cluster_vaf)
    l_mosaic = 10 ** l_mosaic   

    # combine the likelihoods
    # l = [l_refhom, l_althom, l_het, l_mosaic]
    l = [l_refhom, l_mosaic]
    # l_round = [100 if i > 100 else i for i in l]
    # l_scale = np.array([10**i for i in l_round])
    # l_scale = np.array([10**i for i in l])
    
    return l


def spot_genotype(join_info, cell_num=20, epsQ=20, thr_dp=1000, pop_vaf=1e-5):
    """
    Run the spot genotype posterior function

    Input:
        join_info - the join list of spot_site_info and initial_geno
        ['chrom', 'pos', 'strand', 'ref', 'alt', 'spot_barcode',
       'consensus_read_count', 'qA', 'qT', 'qC', 'qG', 'cluster', 'germline',
       'mutant', 'genotype', 'p_mosaic', 'Gi', 'ind_vaf', 'cluster_vaf']
        cell_num - an integer denotes the number of cells in the spot
    Output:
        spot_geno - a dataframe with only one row containing the following information
        #chrom	pos	strand	germline	mutant	cluster spot_barcode    consensus_read_count    l_refhom   l_somatic   max_spot_geno   G_spot_max  depth   vaf    p_mosaic 
    """
    # get info from input
    qA = join_info[7]
    qT = join_info[8]
    qC = join_info[9]
    qG = join_info[10]
    cluster = join_info[11]
    germline = join_info[12]
    mutant = join_info[13]
    cluster_vaf = 0 if join_info[18]=="NA" else float(join_info[18])
    # run
    _, l_norm, spot_geno = spot_posterior(germline, mutant, cluster_vaf, qA, qT, qC, qG, cell_num=cell_num, \
                                            epsQ=epsQ, thr_dp=thr_dp, pop_vaf=pop_vaf)
    # format output as an array

    output = join_info[0:3] + [germline, mutant, str(cluster)] + join_info[5:7] + l_norm + [str(i) for i in spot_geno]

    return output

def convert_tsv_to_parquet(tsv_file: str, output_parquet: str = None):
    if output_parquet is None:
        output_parquet = tsv_file.replace('.txt', '.parquet')
    
    print(f"Converting {tsv_file} → {output_parquet}")
    
    df = load_spot_genotypes_data(tsv_file, prefer_parquet=False)
    df.to_parquet(output_parquet, index=True, compression='snappy')
    

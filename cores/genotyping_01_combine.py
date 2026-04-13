import pandas as pd
from collections import Counter
from SpaceTracer.utils.read_files import load_spot_count_data
from SpaceTracer.utils.utils import str2dict
from scipy.stats import binomtest


# ── combine spot count ────────────────────────────────────────────────────────
class UMICombiner_from_spot:
    """UMI combine to cluster/ind level, from spot/bin level"""
    
    def __init__(self, epsQ=20):
        self.epsQ = epsQ
        self.BASES = ['A', 'T', 'C', 'G']
    
    def _save(self,grouped,output_file):
        columns = ['chrom', 'pos', 'strand', 'ref', 'alt', 'cluster', 'spot_number', 'consensus_read_count', 'qA', 'qT', 'qC', 'qG']
        grouped = grouped[grouped['spot_number'] > 0]
        grouped = grouped[columns]
        grouped.rename(columns={'chrom': '#chrom'}, inplace=True)

        grouped.to_csv(output_file, sep="\t", index=False, na_rep='NA')


    def combine_cluster(self, umi_count_file, output_file, cluster_df):
        """combine UMI to cluster/ind level"""
        df = self._load_data(umi_count_file)
        df.rename(columns={'#chrom': 'chrom'}, inplace=True)
        # print("******************",df)
        # print("******************-------",cluster_df)
        if cluster_df.empty:
            df['cluster'] = "bulk"
        else:
            df = self._merge_clusters(df, cluster_df)
        
        result = self._aggregate(df)
        # print("******************result",result)

        self._save(result, output_file)

    
    def combine_bins(self, umi_count_file, output_file, cluster_df, bins=100):
        """combine UMI from bin level to cluster level(only for Stereo-seq)"""
        df = self._load_data(umi_count_file)
        df = self._rebin_spots(df, bins)
        if cluster_df.empty:
            df['cluster'] = "bulk"
        else:
            df = self._merge_clusters(df, cluster_df)

        result = self._aggregate(df)
        self._save(result, output_file)

    
    def _load_data(self, file_path):
        """load data"""
        df= load_spot_count_data(file_path, sep="\t", header=None,
                            names=['chrom', 'pos', 'strand', 'ref', 'alt', 
                                'barcode', 'consensus_read_count', 'qA', 'qT', 'qC', 'qG'],
                            comment="#")

        # add the spot number column
        df.insert(6, "spot_number", 1)
        return df
    
    def _merge_clusters(self,df,cluster_df):
        """ 
        the cluster file must be formatted as :
        barcode    cluster

        # Note: if for stereo-seq data, the barcode is "x_y", which should be corresponding with bin size,
        for example, the raw (x,y) is (111, 222), the bin 50 level (x,y) is (100,200), 
        so the cluster file may be :
        barcode    cluster
        100_200 cluster1
        """
        df = df.merge(cluster_df, on=['barcode'])
        
        return df

    def _aggregate(self, df):
        """aggregate data"""
        return df.groupby(['chrom', 'pos', 'strand', 'ref', 'cluster']).agg({
            'alt': combine_alt,
            'spot_number': 'sum',
            'consensus_read_count': combine_UMI_count,
            **{f'q{base}': lambda x: combine_q_columns(x,self.epsQ) for base in self.BASES}
        }).reset_index()
    
    def _rebin_spots(self,df,bins):
        count=df
        count[['x', 'y']] = count['barcode'].str.split('_', expand=True)
        count['x'] = count['x'].astype('int32')
        count['y'] = count['y'].astype('int32')
        count['new_x'] = ((count['x'] // bins) * bins).astype('int32') + bins
        count['new_y'] = ((count['y'] // bins) * bins).astype('int32') + bins
        count['barcode'] = count['new_x'].astype(str) + "_" + count['new_y'].astype(str)
        count.drop(columns=['x', 'y', 'new_x', 'new_y'], inplace=True)
        return count

# def UMI_count_cluster(umi_count_file, output_file, type="ind", cluster_file=None, epsQ=20):
#     """
#     Generate the new UMI count file at cluster or individual level
#     by combining the UMI count file at the spot level according to the clusters or combine all.
#     Also delete the low quality (quality < epsQ) consensus reads.

#     Inputs:
#         umi_count_file - the file contain the UMI count info at the spot level
#         output_file - the path of the output file
#         type - "cluster" if combine the info at the cluster level
#                 "ind" if combine the info at the individual level
#                 (default="ind")
#         cluster_file - the file of the clusters (default="None")
#         epsQ - the threshold for the consensus read quality (default=20)
#     """

#     # read in data: the spot UMI count file
#     count = pd.read_csv(umi_count_file, sep="\t", header=None, names=['chr', 'pos', 'ID', 'ref', 'alt', 'barcode', 'consensus_read_count', 'qA', 'qT', 'qC', 'qG'], \
#                         keep_default_na=False, comment = "#")
#     # add the spot number column
#     count.insert(6, "spot_number", 1)
#     if type == "cluster":
#         # cluster info
#         cluster = pd.read_csv(cluster_file, sep="\t", header=None, names=['barcode', 'cluster'], na_values=[])
#         cluster['cluster'] = cluster['cluster'].apply(lambda x: str(int(x)) if isinstance(x, float) and x.is_integer() 
#                                                         else str(x) if pd.notnull(x) else "NA")
#         # merge
#         df = count.merge(cluster, on=['barcode'])
#     else:
#         df = count
#         df['cluster'] = "bulk"
    
#     # combine the columns by the clusters
#     grouped = df.groupby(['chr', 'pos', 'ID', 'ref', 'cluster']).agg({
#         'alt': combine_alt,
#         'spot_number': 'sum',
#         'consensus_read_count': combine_UMI_count,
#         'qA': lambda series: combine_q_columns(series, epsQ=epsQ),
#         'qT': lambda series: combine_q_columns(series, epsQ=epsQ),
#         'qC': lambda series: combine_q_columns(series, epsQ=epsQ),
#         'qG': lambda series: combine_q_columns(series, epsQ=epsQ),
#     }).reset_index()
#     columns = ['chr', 'pos', 'ID', 'ref', 'alt', 'cluster', 'spot_number', 'consensus_read_count', 'qA', 'qT', 'qC', 'qG']
#     grouped = grouped[columns]
#     # write the file
#     grouped.rename(columns={'chr': '#chrom'}, inplace=True)
#     grouped.to_csv(output_file, sep="\t", index=None, na_rep='NA')
    

# ##### Note, the following function is build for stereo-seq, which is clustered in bin N level, not bin1 
# def UMI_count_bins(umi_count_file, output_file, bins=100,cluster_file=None, epsQ=20):
#     """
#     Generate the new UMI count file at cluster or individual level
#     by combining the UMI count file at the spot level according to the clusters or combine all.
#     Also delete the low quality (quality < epsQ) consensus reads.

#     Inputs:
#         umi_count_file - the file contain the UMI count info at the spot level
#         output_file - the path of the output file
#         in_type - "bin" only for stereo_seq, to combine bin into bin50/100/200 level
#         bins - the bin levels to combine
#         cluster_file - the file of the clusters (default="None"). 
#                         both cells with cluster and bin with cluster are supported.
#                         such as:
#                         16600	31100	telencephalon
#         epsQ - the threshold for the consensus read quality (default=20)
#     """

#     count_dtypes = {
#         'chr': 'category',
#         'pos': 'int32',
#         'ID': 'category',
#         'ref': 'category',
#         'alt': 'category',
#         'barcode': 'category',
#         'consensus_read_count': 'category',
#         'qA': 'category',
#         'qT': 'category',
#         'qC': 'category',
#         'qG': 'category'
#     }
    
#     count = pd.read_csv(
#         umi_count_file, sep="\t", header=None, 
#         names=['chr', 'pos', 'ID', 'ref', 'alt', 'barcode', 'consensus_read_count', 'qA', 'qT', 'qC', 'qG'],
#         dtype=count_dtypes,       # 显式指定类型
#         usecols=range(11),        # 防止加载冗余列
#         comment="#",
#         keep_default_na=True
#     )
#     count.insert(6, "spot_number", np.int8(1))  # 小类型优化

#     count[['x', 'y']] = count['barcode'].str.split('_', expand=True)
#     count['x'] = count['x'].astype('int32')
#     count['y'] = count['y'].astype('int32')
#     count['new_x'] = ((count['x'] // bins) * bins).astype('int32') + bins
#     count['new_y'] = ((count['y'] // bins) * bins).astype('int32') + bins
#     count['barcode'] = count['new_x'].astype(str) + "_" + count['new_y'].astype(str)
#     count.drop(columns=['x', 'y', 'new_x', 'new_y'], inplace=True)
#     cell_info=pd.read_csv(cluster_file, sep="\t", header=None, names=['x','y','cluster'], na_values=[])
#     cell_info['x']=cell_info['x'].astype(str)
#     cell_info['y']=cell_info['y'].astype(str)
#     cell_info['barcode']=cell_info['x']+"_"+cell_info['y']
#     cell_info.drop(columns=['x', 'y'], inplace=True)
#     df = count.merge(cell_info, on=['barcode'])

#     # combine the columns by the clusters
#     grouped = df.groupby(['chr', 'pos', 'ID', 'ref', 'cluster']).agg({
#         'alt': combine_alt,
#         'spot_number': 'sum',
#         'consensus_read_count': combine_UMI_count,
#         'qA': lambda series: combine_q_columns(series, epsQ=epsQ),
#         'qT': lambda series: combine_q_columns(series, epsQ=epsQ),
#         'qC': lambda series: combine_q_columns(series, epsQ=epsQ),
#         'qG': lambda series: combine_q_columns(series, epsQ=epsQ),
#     }).reset_index()
#     columns = ['chr', 'pos', 'ID', 'ref', 'alt', 'cluster','spot_number','consensus_read_count', 'qA', 'qT', 'qC', 'qG']
#     grouped = grouped[columns]
#     # write the file
#     grouped = grouped[grouped['spot_number'] > 0]
#     grouped.rename(columns={'chr': '#chrom'}, inplace=True)
#     grouped.to_csv(output_file, sep="\t", index=None, na_rep='NA')
    
# #### utils functions
# def combine_alt(series):
#     """Combine the alt column"""
#     combined = []
#     for alt in series:
#         if alt != ".":
#             alt_list = alt.split(',')
#             combined = list(set(combined) | set(alt_list))
#     if not bool(combined):
#         result = "."
#     else:
#         result = ",".join(combined)
#     return result

# def combine_UMI_count(series):
#     """Combine UMI count column"""
#     combined = [0, 0, 0, 0]
#     for consensus_read_count in series:
#         if isinstance(consensus_read_count, str) and consensus_read_count!="":  # Check if consensus_read_count is a string
#             count = [int(i) for i in consensus_read_count.split(',')]
#             combined = [(a + b) for a, b in zip(combined, count)]
#     # convert to string
#     result = ','.join([str(i) for i in combined])
#     return result


# def combine_q_columns(series, epsQ=20):
#     """Combine quality columns"""
#     combined = {}
#     for q in series:
#         if isinstance(q, str) and q != "":  # Check if q_column is a string
#             q_dict = str2dict(q)
#             q_filter = {k:v for k,v in q_dict.items() if k>=epsQ}
#             combined = dict(Counter(combined) + Counter(q_filter))
#     # convert to string
#     if not bool(combined) or combined=={}:
#         result = "NA"
#     else:
#         result = ','.join(f'{int(key)}:{value}' for key, value in combined.items())
#     return result


# ── filter count in cluster level ────────────────────────────────────────────────────────
class ClusterAlleleFilter:
    """filter allele"""
    
    def __init__(self, alpha=0.05, epsAF=0.003):
        self.alpha = alpha
        self.epsAF = epsAF
        self.BASES = ['A', 'T', 'C', 'G']
    
    def filter(self, cluster_count_file, output_file):
        # cluster_allele_filter(cluster_count_file, output_file, alpha=0.05, epsAF=0.003)
        df = pd.read_csv(cluster_count_file, sep="\t", header=None, \
                    names=['chrom', 'pos', 'strand', 'ref', 'alt', 'cluster', 'spot_number', 'consensus_read_count', 'qA', 'qT', 'qC', 'qG'], \
                    keep_default_na=False, comment = "#")
        # df=cluster_count_df
        if df.empty:
            return
        
        df[['keepA', 'keepT', 'keepC', 'keepG']] = df.apply(
            self._allele_filter_row, axis=1, result_type='expand'
        )
        
        allele_cluster = self._aggregate_by_position(df)
        df = df.merge(allele_cluster, on=['chrom', 'pos'])
        df[['qA_final', 'qT_final', 'qC_final', 'qG_final']] = df.apply(self._quality_choose, axis=1, result_type='expand')
        df = df[['chrom', 'pos', 'strand', 'ref', 'alt', 'cluster', 'spot_number', 'consensus_read_count', 'qA_final', 'qT_final', 'qC_final', 'qG_final']]
        # rename the keep columns
        df = df.rename(columns={
            'qA_final': 'qA',
            'qT_final': 'qT',
            'qC_final': 'qC',
            'qG_final': 'qG'
        })
        # calculate the consensus read counts
        df['consensus_read_count'] = df.apply(self._count_umi, axis=1, result_type='expand')
        df['alt'] = df.apply(self._check_alt, axis=1)
        df.rename(columns={'chrom': '#chrom'}, inplace=True)
        
        df.to_csv(output_file, sep="\t", index=None, na_rep='NA')
    
    def _allele_filter_row(self, row):
        """
        Delete the alternative alleles with low allele frequency by performing a binomial test for each cluster
        (AF < epsAF at alpha significance level)

        Outputs: columns show whether the allele is kept or not
        """
        nucleotide_list = ["A", "T", "C", "G"]
        # transform the string to dict
        try:
            q_dicts = [str2dict(row[x]) for x in ["qA", "qT", "qC", "qG"]]
        except:
            print(f"Cannot find the q_dicts in {row}")
        q_filter = dict(zip(nucleotide_list, q_dicts))
        
        # calculate the numbers of the nucleotides
        count_list = [sum(q_filter[x].values()) for x in nucleotide_list]
        depth = sum(count_list)
        count_dict = dict(zip(nucleotide_list, count_list))

        # binomial test for allele frequency
        count_filter = count_dict.copy()
        if depth != 0:
            for allele, count in count_dict.items():
                result = binomtest(count, n=depth, p=self.epsAF, alternative='greater')
                if result.pvalue >= self.alpha:
                    count_filter[allele] = 0
                    q_filter[allele] = {}

        allele_pass = [bool(count_filter[x]) for x in nucleotide_list]
        # true if keep the allele, false if not
        return pd.Series([allele_pass[0], allele_pass[1], allele_pass[2], allele_pass[3]])

    def _aggregate_by_position(self, df):
        return df.groupby(['chrom', 'pos']).agg({
            f'keep{base}': 'any' for base in self.BASES
        }).rename(columns={f'keep{base}': f'tot{base}' for base in self.BASES})
    

    def _quality_choose(self,row):
        """
        keep the allele if at least one cluster keeps the allele
        otherwise delete the allele below the background error
        """
        qA = row['qA'] if row['totA'] else "NA"
        qT = row['qT'] if row['totT'] else "NA"
        qC = row['qC'] if row['totC'] else "NA"
        qG = row['qG'] if row['totG'] else "NA"
        return pd.Series([qA, qT, qC, qG])

    def _count_umi(self,row):
        """
        Count the consensus reads after filtering
        """
        nucleotide_list = ["A", "T", "C", "G"]
        # transform the string to dict
        q_dicts = [str2dict(row[x]) for x in ["qA", "qT", "qC", "qG"]]
        q_all = dict(zip(nucleotide_list, q_dicts))
        # calculate the numbers of the nucleotides
        count_list = [sum(q_all[x].values()) for x in nucleotide_list]

        # transform format
        count = ','.join([str(i) for i in count_list])
        # true if keep the allele, false if not
        return count

    def _check_alt(self,row):
        """
        Check the remaining alternative alleles after filtering
        """
        nucleotide_list = ["A", "T", "C", "G"]
        # get the number of consensus reads for each allele
        counts = list(map(int, row['consensus_read_count'].split(',')))
        # reference allele
        ref = row['ref']
        # create a list of (allele, count) pairs excluding the ref allele and counts of zero
        alt_counts = [(nucleotide_list[i], counts[i]) for i in range(4) if counts[i] > 0 and nucleotide_list[i] != ref]
        # sort the list by count in descending order
        alt_counts_sorted = sorted(alt_counts, key=lambda x: x[1], reverse=True)
        # join the alleles to create the alt column
        alt = ','.join([allele for allele, _ in alt_counts_sorted]) if len(alt_counts_sorted)>0 else '.'
        return alt


# raw version
# def quality_choose(row):
#     """
#     keep the allele if at least one cluster keeps the allele
#     otherwise delete the allele below the background error
#     """
#     qA = row['qA'] if row['totA'] else "NA"
#     qT = row['qT'] if row['totT'] else "NA"
#     qC = row['qC'] if row['totC'] else "NA"
#     qG = row['qG'] if row['totG'] else "NA"
#     return pd.Series([qA, qT, qC, qG])

# def count_umi(row):
#     """
#     Count the consensus reads after filtering
#     """
#     nucleotide_list = ["A", "T", "C", "G"]
#     # transform the string to dict
#     q_dicts = [str2dict(row[x]) for x in ["qA", "qT", "qC", "qG"]]
#     q_all = dict(zip(nucleotide_list, q_dicts))
#     # calculate the numbers of the nucleotides
#     count_list = [sum(q_all[x].values()) for x in nucleotide_list]

#     # transform format
#     count = ','.join([str(i) for i in count_list])
#     # true if keep the allele, false if not
#     return count

# def check_alt(row):
#     """
#     Check the remaining alternative alleles after filtering
#     """
#     nucleotide_list = ["A", "T", "C", "G"]
#     # get the number of consensus reads for each allele
#     counts = list(map(int, row['consensus_read_count'].split(',')))
#     # reference allele
#     ref = row['ref']
#     # create a list of (allele, count) pairs excluding the ref allele and counts of zero
#     alt_counts = [(nucleotide_list[i], counts[i]) for i in range(4) if counts[i] > 0 and nucleotide_list[i] != ref]
#     # sort the list by count in descending order
#     alt_counts_sorted = sorted(alt_counts, key=lambda x: x[1], reverse=True)
#     # join the alleles to create the alt column
#     alt = ','.join([allele for allele, _ in alt_counts_sorted]) if len(alt_counts_sorted)>0 else '.'
#     return alt

# def cluster_allele_filter(cluster_count_file, output_file, alpha=0.05, epsAF=0.003):
#     """
#     Get the filtered consensus read counts and qualities for clusters
#     here we keep the alternative alleles if they appears at least in one cluster

#     Inputs:
#         cluster_count_file - the umi count file at the cluster level after quality filtering
#         output_file - the path for the output file
#         alpha - the significance level for the binomial test when filtering AF (default=0.05)
#         epsAF - the threshold for the alternative allele frequency or called the background error (default=0.003)
#                 ignore the alternative allele with AF < epsAF by binomial test in the case of multiple alternative alleles
#     """
#     df = pd.read_csv(cluster_count_file, sep="\t", header=None, \
#                     names=['chr', 'pos', 'ID', 'ref', 'alt', 'cluster', 'spot_number', 'consensus_read_count', 'qA', 'qT', 'qC', 'qG'], \
#                     keep_default_na=False, comment = "#")
#     if not df.empty:
#         # create the filtered columns
#         df[['keepA', 'keepT', 'keepC', 'keepG']] = df.apply(lambda row: allele_filter(row, alpha=alpha, epsAF=epsAF), axis=1, result_type='expand')

#         # combine the columns by the clusters
#         def boolean_sum(series):
#             return bool(series.sum())
#         allele_cluster = df.groupby(['chr', 'pos']).agg({
#             'keepA': boolean_sum,
#             'keepT': boolean_sum,
#             'keepC': boolean_sum,
#             'keepG': boolean_sum
#         }).reset_index()
#         # rename the keep columns
#         allele_cluster = allele_cluster.rename(columns={
#             'keepA': 'totA',
#             'keepT': 'totT',
#             'keepC': 'totC',
#             'keepG': 'totG'
#         })
#         df = df.merge(allele_cluster, on=['chr', 'pos'])

#         # fill in the quality counts if kept otherwise deleted
#         df[['qA_final', 'qT_final', 'qC_final', 'qG_final']] = df.apply(quality_choose, axis=1, result_type='expand')
#         df = df[['chr', 'pos', 'ID', 'ref', 'alt', 'cluster', 'spot_number', 'consensus_read_count', 'qA_final', 'qT_final', 'qC_final', 'qG_final']]
#         # rename the keep columns
#         df = df.rename(columns={
#             'qA_final': 'qA',
#             'qT_final': 'qT',
#             'qC_final': 'qC',
#             'qG_final': 'qG'
#         })
#         # calculate the consensus read counts
#         df['consensus_read_count'] = df.apply(count_umi, axis=1, result_type='expand')
#         # find the remaining alternative alleles
#         df['alt'] = df.apply(check_alt, axis=1)

#     # write the file
#     df.rename(columns={'chr': '#chrom'}, inplace=True)
#     df.to_csv(output_file, sep="\t", index=None, na_rep='NA')


# def allele_filter(row, alpha=0.05, epsAF=0.01):
#     """
#     Delete the alternative alleles with low allele frequency by performing a binomial test for each cluster
#     (AF < epsAF at alpha significance level)

#     Outputs: columns show whether the allele is kept or not
#     """
#     nucleotide_list = ["A", "T", "C", "G"]
#     # transform the string to dict
#     try:
#         q_dicts = [str2dict(row[x]) for x in ["qA", "qT", "qC", "qG"]]
#     except:
#         print(row)
#     q_filter = dict(zip(nucleotide_list, q_dicts))
    
#     # calculate the numbers of the nucleotides
#     count_list = [sum(q_filter[x].values()) for x in nucleotide_list]
#     depth = sum(count_list)
#     count_dict = dict(zip(nucleotide_list, count_list))

#     # binomial test for allele frequency
#     count_filter = count_dict.copy()
#     if depth != 0:
#         for allele, count in count_dict.items():
#             result = binomtest(count, n=depth, p=epsAF, alternative='greater')
#             if result.pvalue >= alpha:
#                 count_filter[allele] = 0
#                 q_filter[allele] = {}

#     allele_pass = [bool(count_filter[x]) for x in nucleotide_list]
#     # true if keep the allele, false if not
#     return pd.Series([allele_pass[0], allele_pass[1], allele_pass[2], allele_pass[3]])

# def quality_choose(row):
#     """
#     keep the allele if at least one cluster keeps the allele
#     otherwise delete the allele below the background error
#     """
#     qA = row['qA'] if row['totA'] else "NA"
#     qT = row['qT'] if row['totT'] else "NA"
#     qC = row['qC'] if row['totC'] else "NA"
#     qG = row['qG'] if row['totG'] else "NA"
#     return pd.Series([qA, qT, qC, qG])


# ── filter count in cluster level ────────────────────────────────────────────────────────


class UMICombiner_from_cluster:
# UMI_count_ind_from_cluster(cluster_filter_count_file, output_file, epsQ=20):
    """
    Generate the new UMI count file at  individual level from cluster count file
    by combining the UMI count file at the cluster level.
    So that remains the alternative alleles which pass the allele frequency (AF) binomial test 
    at the cluster level.

    Inputs:
        cluster_filter_count_file - the file contain the UMI count info at the spot level
        output_file - the path of the output file
    """
    def __init__(self,epsQ):
        self.epsQ = epsQ
        self.BASES = ['A', 'T', 'C', 'G']
    
    def _load_data(self, cluster_filter_count_file):
        """load data"""
        df= pd.read_csv(cluster_filter_count_file, sep="\t", header=None, 
                    names=['chrom', 'pos', 'strand', 'ref', 'alt', 'cluster', 'spot_number', 'consensus_read_count', 'qA', 'qT', 'qC', 'qG'], 
                    keep_default_na=False, comment = "#")
        df['cluster'] = df['cluster'].apply(lambda x: str(int(x)) if isinstance(x, float) and x.is_integer() 
                                                                else str(x) if pd.notnull(x) else "NA")
        return df

    def _aggregate(self, df):
        """aggregate data"""
        return df.groupby(['chrom', 'pos', 'strand', 'ref']).agg({
            'alt': combine_alt,
            'spot_number': 'sum',
            'consensus_read_count': combine_UMI_count,
            **{f'q{base}': lambda x: combine_q_columns(x, self.epsQ) for base in self.BASES}
        }).reset_index()
    
    def combine_ind(self,cluster_filter_count_file,output_file):
        df=self._load_data(cluster_filter_count_file)
        grouped=self._aggregate(df)
        grouped['cluster'] = "bulk"
        columns = ['chrom', 'pos', 'strand', 'ref', 'alt', 'cluster', 'spot_number', 'consensus_read_count', 'qA', 'qT', 'qC', 'qG']
        grouped = grouped[columns]
        grouped.rename(columns={'chrom': '#chrom'}, inplace=True)
        grouped.to_csv(output_file, sep="\t", index=None, na_rep='NA')

#### utils functions
def combine_alt(series):
    """Combine the alt column"""
    combined = []
    for alt in series:
        if alt != ".":
            alt_list = alt.split(',')
            combined = list(set(combined) | set(alt_list))
    if not bool(combined):
        result = "."
    else:
        result = ",".join(combined)
    return result

def combine_UMI_count(series):
    """Combine UMI count column"""
    combined = [0, 0, 0, 0]
    for consensus_read_count in series:
        if isinstance(consensus_read_count, str):  # Check if consensus_read_count is a string
            count = [int(i) for i in consensus_read_count.split(',')]
            combined = [(a + b) for a, b in zip(combined, count)]
    # convert to string
    result = ','.join([str(i) for i in combined])
    return result

def combine_q_columns(series, epsQ=20):
    """Combine quality columns"""
    combined = {}
    for q in series:
        if isinstance(q, str):  # Check if q_column is a string
            q_dict = str2dict(q)
            q_filter = {k:v for k,v in q_dict.items() if k>=epsQ}
            combined = dict(Counter(combined) + Counter(q_filter))
    # convert to string
    if not bool(combined):
        result = "NA"
    else:
        result = ','.join(f'{int(key)}:{value}' for key, value in combined.items())
    return result

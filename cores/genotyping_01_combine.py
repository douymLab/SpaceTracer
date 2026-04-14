import pandas as pd
from collections import Counter
from scipy.stats import binomtest

from SpaceTracer.utils.read_files import load_spot_count_data
from SpaceTracer.utils.utils import str2dict


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
        if isinstance(consensus_read_count, str):
            count = [int(i) for i in consensus_read_count.split(',')]
            combined = [(a + b) for a, b in zip(combined, count)]
    result = ','.join([str(i) for i in combined])
    return result


def combine_q_columns(series, epsQ=20):
    """Combine quality columns"""
    combined = {}
    for q in series:
        if isinstance(q, str):
            q_dict = str2dict(q)
            q_filter = {k: v for k, v in q_dict.items() if k >= epsQ}
            combined = dict(Counter(combined) + Counter(q_filter))
    if not bool(combined):
        result = "NA"
    else:
        result = ','.join(f'{int(key)}:{value}' for key, value in combined.items())
    return result


class UMICombiner_from_spot:
    """UMI combine to cluster/ind level, from spot/bin level"""

    def __init__(self, epsQ=20):
        self.epsQ = epsQ
        self.BASES = ['A', 'T', 'C', 'G']

    # ----------------------------
    # I/O helpers
    # ----------------------------
    def _load_data(self, file_path):
        df = load_spot_count_data(
            file_path,
            sep="\t",
            header=None,
            names=[
                'chrom', 'pos', 'strand', 'ref', 'alt',
                'barcode', 'consensus_read_count', 'qA', 'qT', 'qC', 'qG'
            ],
            comment="#"
        )
        if '#chrom' in df.columns and 'chrom' not in df.columns:
            df = df.rename(columns={'#chrom': 'chrom'})
        df.insert(6, "spot_number", 1)
        return df

    def _save(self, grouped, output_file):
        columns = [
            'chrom', 'pos', 'strand', 'ref', 'alt', 'cluster',
            'spot_number', 'consensus_read_count', 'qA', 'qT', 'qC', 'qG'
        ]
        grouped = grouped[grouped['spot_number'] > 0]
        grouped = grouped[columns]
        grouped = grouped.rename(columns={'chrom': '#chrom'})
        grouped.to_csv(output_file, sep="\t", index=False, na_rep='NA')

    # ----------------------------
    # DF API
    # ----------------------------
    def load_df(self, file_path):
        return self._load_data(file_path)

    def save_df(self, df, output_file):
        self._save(df, output_file)

    def combine_cluster_df(self, df, cluster_df):
        df = df.copy()
        if '#chrom' in df.columns and 'chrom' not in df.columns:
            df = df.rename(columns={'#chrom': 'chrom'})

        if cluster_df.empty:
            df['cluster'] = "bulk"
        else:
            df = self._merge_clusters(df, cluster_df)

        result = self._aggregate(df)
        return result

    def combine_bins_df(self, df, cluster_df, bins=100):
        df = df.copy()
        if '#chrom' in df.columns and 'chrom' not in df.columns:
            df = df.rename(columns={'#chrom': 'chrom'})

        df = self._rebin_spots(df, bins)
        if cluster_df.empty:
            df['cluster'] = "bulk"
        else:
            df = self._merge_clusters(df, cluster_df)

        result = self._aggregate(df)
        return result

    # ----------------------------
    # File API (backward compatible)
    # ----------------------------
    def combine_cluster(self, umi_count_file, output_file, cluster_df):
        df = self._load_data(umi_count_file)
        result = self.combine_cluster_df(df, cluster_df)
        self._save(result, output_file)

    def combine_bins(self, umi_count_file, output_file, cluster_df, bins=100):
        df = self._load_data(umi_count_file)
        result = self.combine_bins_df(df, cluster_df, bins=bins)
        self._save(result, output_file)

    # ----------------------------
    # internals
    # ----------------------------
    def _merge_clusters(self, df, cluster_df):
        return df.merge(cluster_df, on=['barcode'])

    def _aggregate(self, df):
        return df.groupby(['chrom', 'pos', 'strand', 'ref', 'cluster']).agg({
            'alt': combine_alt,
            'spot_number': 'sum',
            'consensus_read_count': combine_UMI_count,
            **{f'q{base}': (lambda x, _epsQ=self.epsQ: combine_q_columns(x, _epsQ)) for base in self.BASES}
        }).reset_index()

    def _rebin_spots(self, df, bins):
        count = df.copy()
        count[['x', 'y']] = count['barcode'].str.split('_', expand=True)
        count['x'] = count['x'].astype('int32')
        count['y'] = count['y'].astype('int32')
        count['new_x'] = ((count['x'] // bins) * bins).astype('int32') + bins
        count['new_y'] = ((count['y'] // bins) * bins).astype('int32') + bins
        count['barcode'] = count['new_x'].astype(str) + "_" + count['new_y'].astype(str)
        count = count.drop(columns=['x', 'y', 'new_x', 'new_y'])
        return count


class ClusterAlleleFilter:
    """filter allele"""

    def __init__(self, alpha=0.05, epsAF=0.003):
        self.alpha = alpha
        self.epsAF = epsAF
        self.BASES = ['A', 'T', 'C', 'G']

    # ----------------------------
    # I/O helpers
    # ----------------------------
    def _load_df(self, cluster_count_file):
        df = pd.read_csv(
            cluster_count_file,
            sep="\t",
            header=None,
            names=[
                'chrom', 'pos', 'strand', 'ref', 'alt', 'cluster',
                'spot_number', 'consensus_read_count', 'qA', 'qT', 'qC', 'qG'
            ],
            keep_default_na=False,
            comment="#"
        )
        if '#chrom' in df.columns and 'chrom' not in df.columns:
            df = df.rename(columns={'#chrom': 'chrom'})
        return df

    def _save_df(self, df, output_file):
        out = df.copy()
        out = out.rename(columns={'chrom': '#chrom'})
        out.to_csv(output_file, sep="\t", index=None, na_rep='NA')

    # ----------------------------
    # DF API
    # ----------------------------
    def filter_df(self, df):
        df = df.copy()
        if '#chrom' in df.columns and 'chrom' not in df.columns:
            df = df.rename(columns={'#chrom': 'chrom'})

        if df.empty:
            return df

        df[['keepA', 'keepT', 'keepC', 'keepG']] = df.apply(
            self._allele_filter_row, axis=1, result_type='expand'
        )

        allele_cluster = self._aggregate_by_position(df)
        df = df.merge(allele_cluster, on=['chrom', 'pos'])
        df[['qA_final', 'qT_final', 'qC_final', 'qG_final']] = df.apply(
            self._quality_choose, axis=1, result_type='expand'
        )

        df = df[[
            'chrom', 'pos', 'strand', 'ref', 'alt', 'cluster',
            'spot_number', 'consensus_read_count',
            'qA_final', 'qT_final', 'qC_final', 'qG_final'
        ]]

        df = df.rename(columns={
            'qA_final': 'qA',
            'qT_final': 'qT',
            'qC_final': 'qC',
            'qG_final': 'qG'
        })

        df['consensus_read_count'] = df.apply(self._count_umi, axis=1, result_type='expand')
        df['alt'] = df.apply(self._check_alt, axis=1)
        return df

    # ----------------------------
    # File API (backward compatible)
    # ----------------------------
    def filter(self, cluster_count_file, output_file):
        df = self._load_df(cluster_count_file)
        if df.empty:
            return
        out = self.filter_df(df)
        self._save_df(out, output_file)

    # ----------------------------
    # internals
    # ----------------------------
    def _allele_filter_row(self, row):
        nucleotide_list = ["A", "T", "C", "G"]
        try:
            q_dicts = [str2dict(row[x]) for x in ["qA", "qT", "qC", "qG"]]
        except Exception:
            print(f"Cannot find the q_dicts in {row}")
            q_dicts = [{}, {}, {}, {}]

        q_filter = dict(zip(nucleotide_list, q_dicts))
        count_list = [sum(q_filter[x].values()) for x in nucleotide_list]
        depth = sum(count_list)
        count_dict = dict(zip(nucleotide_list, count_list))

        count_filter = count_dict.copy()
        if depth != 0:
            for allele, count in count_dict.items():
                result = binomtest(count, n=depth, p=self.epsAF, alternative='greater')
                if result.pvalue >= self.alpha:
                    count_filter[allele] = 0
                    q_filter[allele] = {}

        allele_pass = [bool(count_filter[x]) for x in nucleotide_list]
        return pd.Series([allele_pass[0], allele_pass[1], allele_pass[2], allele_pass[3]])

    def _aggregate_by_position(self, df):
        return df.groupby(['chrom', 'pos']).agg({
            f'keep{base}': 'any' for base in self.BASES
        }).rename(columns={f'keep{base}': f'tot{base}' for base in self.BASES}).reset_index()

    def _quality_choose(self, row):
        qA = row['qA'] if row['totA'] else "NA"
        qT = row['qT'] if row['totT'] else "NA"
        qC = row['qC'] if row['totC'] else "NA"
        qG = row['qG'] if row['totG'] else "NA"
        return pd.Series([qA, qT, qC, qG])

    def _count_umi(self, row):
        nucleotide_list = ["A", "T", "C", "G"]
        q_dicts = [str2dict(row[x]) for x in ["qA", "qT", "qC", "qG"]]
        q_all = dict(zip(nucleotide_list, q_dicts))
        count_list = [sum(q_all[x].values()) for x in nucleotide_list]
        count = ','.join([str(i) for i in count_list])
        return count

    def _check_alt(self, row):
        nucleotide_list = ["A", "T", "C", "G"]
        counts = list(map(int, row['consensus_read_count'].split(',')))
        ref = row['ref']
        alt_counts = [
            (nucleotide_list[i], counts[i])
            for i in range(4)
            if counts[i] > 0 and nucleotide_list[i] != ref
        ]
        alt_counts_sorted = sorted(alt_counts, key=lambda x: x[1], reverse=True)
        alt = ','.join([allele for allele, _ in alt_counts_sorted]) if len(alt_counts_sorted) > 0 else '.'
        return alt


class UMICombiner_from_cluster:
    def __init__(self, epsQ):
        self.epsQ = epsQ
        self.BASES = ['A', 'T', 'C', 'G']

    def _load_data(self, cluster_filter_count_file):
        df = pd.read_csv(
            cluster_filter_count_file,
            sep="\t",
            header=None,
            names=[
                'chrom', 'pos', 'strand', 'ref', 'alt', 'cluster',
                'spot_number', 'consensus_read_count', 'qA', 'qT', 'qC', 'qG'
            ],
            keep_default_na=False,
            comment="#"
        )
        if '#chrom' in df.columns and 'chrom' not in df.columns:
            df = df.rename(columns={'#chrom': 'chrom'})

        df['cluster'] = df['cluster'].apply(
            lambda x: str(int(x)) if isinstance(x, float) and x.is_integer()
            else str(x) if pd.notnull(x) else "NA"
        )
        return df

    def _aggregate(self, df):
        return df.groupby(['chrom', 'pos', 'strand', 'ref']).agg({
            'alt': combine_alt,
            'spot_number': 'sum',
            'consensus_read_count': combine_UMI_count,
            **{f'q{base}': (lambda x, _epsQ=self.epsQ: combine_q_columns(x, _epsQ)) for base in self.BASES}
        }).reset_index()

    # DF API
    def combine_ind_df(self, df):
        df = df.copy()
        if '#chrom' in df.columns and 'chrom' not in df.columns:
            df = df.rename(columns={'#chrom': 'chrom'})

        if 'cluster' in df.columns:
            df['cluster'] = df['cluster'].apply(
                lambda x: str(int(x)) if isinstance(x, float) and x.is_integer()
                else str(x) if pd.notnull(x) else "NA"
            )

        grouped = self._aggregate(df)
        grouped['cluster'] = "bulk"
        columns = [
            'chrom', 'pos', 'strand', 'ref', 'alt', 'cluster',
            'spot_number', 'consensus_read_count', 'qA', 'qT', 'qC', 'qG'
        ]
        grouped = grouped[columns]
        return grouped

    # File API
    def combine_ind(self, cluster_filter_count_file, output_file):
        df = self._load_data(cluster_filter_count_file)
        grouped = self.combine_ind_df(df)
        grouped = grouped.rename(columns={'chrom': '#chrom'})
        grouped.to_csv(output_file, sep="\t", index=None, na_rep='NA')

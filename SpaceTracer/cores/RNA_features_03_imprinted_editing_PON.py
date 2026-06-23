import os
import pandas as pd
import pyranges as pr
    

def add_col_from_bed(df: pd.DataFrame, file: str) -> pd.Series:
    if file and os.path.exists(file):
        filter_bed = pr.read_bed(file)
        
        temp_df = pd.DataFrame({
            'Chromosome': df['chrom'].values,
            'Start': df['pos'].values - 1,  
            'End': df['pos'].values,
            'Identifier': df.index.values   
        })
        
        sites = pr.PyRanges(temp_df)
        passed = sites.subtract(filter_bed)
        passed_identifiers = set(passed.Identifier)
        
        result = ~df.index.isin(passed_identifiers)
        
        return pd.Series(result, index=df.index)
        
    else:
        return pd.Series("unknown", index=df.index)


def add_col_from_mutant(df: pd.DataFrame, file: str) -> pd.Series:
    if file and os.path.exists(file):
        filter_df = pd.read_csv(file, sep=r'\s+', header=None, usecols=[0, 1], names=['chrom', 'pos'])
        db_index = pd.MultiIndex.from_frame(filter_df[['chrom', 'pos']])
        query_index = pd.MultiIndex.from_frame(df[['chrom', 'pos']])
        
        return pd.Series(query_index.isin(db_index), index=df.index)
    else:
        return pd.Series(False, index=df.index)


def add_col_from_mutant_from_df(df: pd.DataFrame, filter_df: pd.DataFrame) -> pd.Series:
    if not filter_df.empty:
        filter_index = pd.MultiIndex.from_arrays(
            [filter_df['chrom'], filter_df['pos'], filter_df['ref'], filter_df['alt']]
        )

        query_index = pd.MultiIndex.from_arrays(
            [df['chrom'], df['pos'], df['ref'], df['alt']]
        )

        is_in_pon = query_index.isin(filter_index)
        return pd.Series(is_in_pon, index=df.index)
    else:
        return pd.Series(False, index=df.index)
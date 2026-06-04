import os
import pandas as pd
import pyranges as pr
    

def add_col_from_bed(df: pd.DataFrame, file: str) -> pd.Series:
    if file and os.path.exists(file):
        filter_bed = pr.read_bed(file)
        
        temp_df = pd.DataFrame({
            'Chromosome': df['chrom'].values,
            'Start': df['pos'].values - 1,
            'End': df['pos'].values
        })
        
        sites = pr.PyRanges(temp_df)
        
        sites.Identifier = df.index.values
        
        passed = sites.subtract(filter_bed)
        passed_identifiers = set(passed.Identifier)
        
        return pd.Series(
            ~df.index.isin(passed_identifiers),
            index=df.index
        )
    else:
        return pd.Series("unknown", index=df.index)


def add_col_from_mutant(df: pd.DataFrame, file:str) -> pd.Series:
    if file and os.path.exists(file):
        filter_df = pd.read_csv(file, sep="\t", header=None, 
                            names=['chrom', 'pos', 'ref', 'alt'])
        
        filter_df.index = pd.MultiIndex.from_arrays(
            [filter_df['chrom'], filter_df['pos'], filter_df['ref'], filter_df['alt']],
            names=['chrom', 'pos', 'ref', 'alt']
        )
        
        is_in_pon = df.index.isin(filter_df.index)
        
        return pd.Series(is_in_pon, index=df.index) # True means in pon list
    else:
        return pd.Series("unknown", index=df.index)


def add_col_from_mutant_from_df(df: pd.DataFrame, filter_df) -> pd.Series:
    if not filter_df.empty:
        filter_index = pd.MultiIndex.from_arrays(
            [filter_df['chrom'], filter_df['pos'], filter_df['ref'], filter_df['alt']],
            names=['chrom', 'pos', 'ref', 'alt']
        )

        is_in_pon = df.index.isin(filter_index)
        return pd.Series(is_in_pon, index=df.index)
    else:
        return pd.Series("unknown", index=df.index)

import pandas as pd
import pyranges as pr
    

# def filter_imprinted_genes_simple(identifier_file, imprinted_bed_file):
#     with open(identifier_file) as f:
#         lines = [line.strip() for line in f]
    
#     identifiers = []
#     for line in lines:
#         parts = line.split('_')
#         if len(parts) >= 4:
#             identifiers.append({
#                 'identifier': line,
#                 'chrom': parts[0],
#                 'pos': int(parts[1]),
#                 'ref': parts[2],
#                 'alt': parts[3]
#             })
    
#     df = pd.DataFrame(identifiers)
    
#     sites = pr.PyRanges(
#         chromosomes=df['chrom'].values,
#         starts=df['pos'].values,
#         ends=df['pos'].values
#     )
#     sites.Identifier = df['identifier'].values
    
#     imprinted = pr.read_bed(imprinted_bed_file)
    
#     filtered = sites.subtract(imprinted)

#     return filtered


# def make_filter_from_bed(file):
#     filter_bed = pr.read_bed(file)
    
#     def filter_func(df: pd.DataFrame) -> pd.Series:
#         sites = pr.PyRanges(
#             chromosomes=df['chrom'].values,
#             starts=df['pos'].values - 1,  # 0-based
#             ends=df['pos'].values
#         )
#         sites.Identifier = df['identifier'].values
        
#         passed = sites.subtract(filter_bed)
#         passed_identifiers = set(passed.Identifier)
        
#         return df['identifier'].isin(passed_identifiers)
    
#     return filter_func

def add_col_from_bed(df: pd.DataFrame, file: str) -> pd.Series:
    if file:
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
    if file:
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

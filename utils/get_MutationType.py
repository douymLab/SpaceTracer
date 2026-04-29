import itertools
import os
import pandas as pd
from pyfaidx import Fasta


def get_substitution_type(chrom_dict,chrom,pos):
    '''
    line also could be an identifier
    return_type: 
    '''
    # print(chrom,pos)
    up_base=str(chrom_dict[chrom][pos-1-1])
    down_base=str(chrom_dict[chrom][pos+1-1])

    return up_base,down_base


def trans_base(base):
    try:
        reverse_dict={"A":"T","T":"A","C":"G","G":"C"}
    except:
        return ""
    reverse_base=reverse_dict[base]
    return reverse_base


def get_default_labels(choice: str) -> str:
    if choice not in ["DNA", "RNA", "96", "192"]:
        raise ValueError("Choice must be 'DNA'/'96' or 'RNA'/'192")

    if choice=="DNA" or choice=="96":
        mid_list=["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"]
        value=96
    else:
        mid_list=["C>A", "C>G", "C>T", "T>A", "T>C", "T>G","A>C","A>G","A>T","G>A","G>C","G>T"]
        value=192

    first = ["A", "T", "C", "G"]
    inner_bracket = [[x] * 16 for x in mid_list]
    inner_bracket = [item for sublist in inner_bracket for item in sublist]
    outter_bracket = [x for x in list(itertools.product(first, first))]
    result = [
        outter_bracket[f % 16][0]
        + "["
        + inner_bracket[f]
        + "]"
        + outter_bracket[f % 16][1]
        for f in range(0, value)
    ]
    return result


def get_mutation_type(df, fasta_file, mode):
    """
    mode: 'DNA'/'96'/96 or 'RNA'/'192'/192 or 'both'
    """
    if isinstance(fasta_file, str):
        if os.path.exists(fasta_file):
            chromosome_dict = Fasta(fasta_file)
        else:
            raise FileNotFoundError(f"fasta_file not found! please check {fasta_file}!")
    elif isinstance(fasta_file, dict):
        chromosome_dict = fasta_file
    else:
        raise TypeError(f"fasta_file must be str or dict, got {type(fasta_file)}")

    if mode in ['RNA', '192', 192] and 'strand' not in df.columns:
        raise TypeError("RNA mode requires strand column (+/-)")
    
    contexts = df.apply(
        lambda r: get_substitution_type(chromosome_dict, r['chrom'], r['pos']), 
        axis=1
    )
    df = df.copy()
    df['up_base'] = contexts.apply(lambda x: x[0])
    df['down_base'] = contexts.apply(lambda x: x[1])
    
    mask = ~df[['up_base', 'down_base', 'ref', 'alt']].isin(['N']).any(axis=1)
    df = df[mask]
    
    if mode in ['RNA', '192', 192]:
        if not df['strand'].isin(['+', '-']).all():
            raise TypeError("your input is wrong! strand is not +/-.")
    
    def _make_DNAMutationType(row):
        up, down, ref, alt = row['up_base'], row['down_base'], row['ref'], row['alt']
        if ref in 'CT':
            return f"{up}[{ref}>{alt}]{down}"
        else:
            return f"{trans_base(down)}[{trans_base(ref)}>{trans_base(alt)}]{trans_base(up)}"
        
    def _make_RNAMutationType(row):  # RNA mode
        up, down, ref, alt = row['up_base'], row['down_base'], row['ref'], row['alt']
        if row['strand'] == '-':
            return f"{trans_base(down)}[{trans_base(ref)}>{trans_base(alt)}]{trans_base(up)}"
        else:
            return f"{up}[{ref}>{alt}]{down}"

    if mode in ['DNA', '96', 96, 'both']:
        df["DNAMutationType"]=df.apply(_make_DNAMutationType, axis=1)
    if mode in ['RNA', '192', 192,'both']:
        df["RNAMutationType"]=df.apply(_make_RNAMutationType, axis=1)

    return df


def make_mutation_order():
    first = ["A", "C", "G", "T"]
    inner_bracket = [[x] * 16 for x in ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G", "A>C", "A>G", "A>T", "G>A", "G>C", "G>T"]]
    inner_bracket = [item for sublist in inner_bracket for item in sublist]
    outter_bracket = list(itertools.product(first, first))
    
    return [
        outter_bracket[f % 16][0] + "[" + inner_bracket[f] + "]" + outter_bracket[f % 16][1]
        for f in range(192)
    ]


def reorder_mutation_df(mutation_counter):
    order = make_mutation_order()
    return pd.DataFrame({
        'MutationType': order,
        'Count': [mutation_counter.get(m, 0) for m in order]
    })
from pathlib import Path
from typing import List, Optional
import pyranges as pr
import pandas as pd
import pysam
import json
import sqlite3

## barcode file read (for "tissue_positions.csv" file )
def handle_barcode(barcode_file,in_tissue_choose=0):
    '''
    input: 
    barcode_file
    (barcode,in_tissue_or_not,simplified_location_x,simplified_location_y,real_location_x,real_location_y) #in_tissue_or_not:0 means not in tissue, 1 means in tissue
    ACGCCTGACACGCGCT-1,0,0,0,323,308
    TACCGATCCAACACTT-1,0,1,1,334,326
    
    argument:
    in_tissue_choose: only the in tissue barcode will be chosen

    output:
    dict: {barcode: (in_tissue, pos1, pos2),...}
    '''
    barcode_dict={}
    f=open(barcode_file,"r")
    for line in f.readlines():
        s = line.strip().split(",")
        barcode=s[0]
        if barcode!="barcode":
            in_tissue=int(s[1])
            pos1= int(s[2]); pos2=int(s[3])
            if in_tissue_choose==0 and in_tissue==1:
                barcode_dict[barcode]= (pos1, pos2)
            elif in_tissue_choose==1:
                barcode_dict[barcode]= (pos1, pos2)
    return barcode_dict


## bam file read
def read_bam(bam_file):
    f1 = pysam.AlignmentFile(bam_file,"r")
    return f1

## spot count load
def load_spot_count_data(
    file_path,
    sep: Optional[str] = '\t',
    header: Optional[int] = 0,
    names: Optional[List[str]] = None,
    comment: Optional[str] = None,
    prefer_parquet: bool = True,
    create_parquet: bool = False,
    **kwargs,
) -> pd.DataFrame:
    file_path = Path(file_path)
    
    if not file_path.exists():
        raise FileNotFoundError(f"File not found: {file_path}")
    
    # ── Strategy 1: Check for Parquet version ────────────────────────────
    if prefer_parquet and file_path.suffix.lower() != '.parquet':
        parquet_path = file_path.with_suffix('.parquet')
        
        if parquet_path.exists():
            return load_parquet(parquet_path, **kwargs)
    
    # ── Strategy 2: Load original file ───────────────────────────────────
    if file_path.suffix.lower() == '.parquet':
        return load_parquet(file_path, **kwargs)
    else:
        df = load_text_file(
            file_path,
            sep=sep,
            header=header,
            names=names,
            comment=comment,
            **kwargs
        )
        
        # ── Strategy 3: Create Parquet cache (optional) ──────────────────
        if create_parquet:
            parquet_path = file_path.with_suffix('.parquet')
            if not parquet_path.exists():
                df.to_parquet(parquet_path, index=False, compression='snappy')
        
        return df


def load_parquet(file_path: Path, **kwargs) -> pd.DataFrame:
    """Load parquet file with proper kwargs handling."""
    if 'usecols' in kwargs and 'columns' not in kwargs:
        kwargs['columns'] = kwargs.pop('usecols')
    
    for key in ['sep', 'header', 'names', 'comment']:
        kwargs.pop(key, None)
    
    return pd.read_parquet(file_path, **kwargs)


def load_text_file(
    file_path: Path,
    sep: Optional[str] = None,
    header: Optional[int] = 0,
    names: Optional[List[str]] = None,
    comment: Optional[str] = None,
    **kwargs,
) -> pd.DataFrame:
    """Load text file (CSV/TSV/TXT)."""
    if sep is None:
        sep = '\t' if file_path.suffix.lower() in ['.tsv', '.txt'] else ','
    
    return pd.read_csv(
        file_path,
        sep=sep,
        header=header,
        names=names,
        comment=comment,
        **kwargs
    )


## spot genotype load 
def load_spot_genotypes_data(file: str, prefer_parquet: bool = True) -> pd.DataFrame:
    """
    Load spot genotype file, auto-detecting format.
    
    Logic:
    ──────
    1. If .parquet exists and prefer_parquet=True → load Parquet (fast)
    2. Otherwise → load TSV and set_index (slower)
    
    Returns
    ───────
    DataFrame with 'identifier' as index (ready for .loc[] queries)
    """
    file_path = Path(file)
    
    # Check for Parquet version
    if file_path.suffix in ['.txt', '.tsv', '.out']:
        parquet_path = file_path.with_suffix('.parquet')
    else:
        parquet_path = file_path.with_name(file_path.stem + '.parquet')
    
    # parquet
    if prefer_parquet and parquet_path.exists():
        # print(f"Loading Parquet: {parquet_path}")
        df = pd.read_parquet(parquet_path)
        return df
    
    colnames = [
        "chr", "pos", "strand", "germline", "mutant", "cluster",
        "spot_barcode", "consensus_read_count", "l_germline", "l_mosaic",
        "max_spot_geno", "G_spot_max", "depth", "vaf", "p_mosaic"
    ]
    
    # csv
    df = pd.read_csv(file, sep='\t', header=None, names=colnames, comment="#")
    df['pos'] = df['pos'].astype(str)
    df['identifier'] = (
        df['chr'] + '_' + df['pos'] + '_' +
        df['germline'] + '_' + df['mutant'].astype(str)
    )
    df = df.set_index('identifier')
    
    return df


## spatial feature
def load_spatial_features(file_path: str,prefer_parquet: bool = True) -> pd.DataFrame:
    file_path = Path(file_path)
    if not file_path.exists():
        raise FileNotFoundError(f"File not found: {file_path}")
    
    if prefer_parquet and file_path.suffix != '.parquet':
        parquet_path = file_path.with_suffix('.parquet')
        
        if parquet_path.exists():
            return pd.read_parquet(parquet_path)
    
    if file_path.suffix == '.parquet':
        return pd.read_parquet(file_path)
    
    df = pd.read_csv(file_path, sep='\t', comment='#')
    
    # Build identifier index
    df['identifier'] = (
        df['#chrom'] + '_' +
        df['pos'].astype(str) + '_' +
        df['ref'] + '_' +
        df['alt']
    )
    df = df.set_index('identifier')
    
    return df


def load_gtf_file(file: str, save_bed: bool = True):
    file_path = Path(file)
    if file_path.suffix in ['.gtf', '.gtf.gz']:
        bed_path = file_path.with_suffix('.bed')
    else:
        bed_path = file_path.with_name(file_path.stem + '.bed')
    
    # read exist bed
    if save_bed and bed_path.exists():
        df = pr.read_bed(bed_path)
    
    gr = pr.read_gtf(file_path)
    genes = gr[gr.Feature == "gene"]
    try:
        genes.to_bed(bed_path)
    except:
        pass
    
    short_df=df[["Chromosome","Start","End","Strand","gene_name"]]
    short_df = short_df.rename(columns={"gene_name": "GeneName"})
    return short_df



def load_manifest_files_for_chunk_files(manifest_file, selected_groups=None):
    """
    Parameters
    ----------
    manifest_file : str
    selected_groups : list[str] or None
        eg: ["chr10", "chr11"]
    Returns
    -------
    list[str]
    """
    with open(manifest_file, "r") as f:
        manifest = json.load(f)

    chromosome_groups = manifest.get("chromosome_groups", {})
    chunk_files = []

    for group_name, group_info in chromosome_groups.items():
        if group_name=='all':
            files = group_info.get("files", [])
            return files
        else:
            if selected_groups is not None and group_name not in selected_groups:
                continue

        files = group_info.get("files", [])
        chunk_files.extend(files)

    return chunk_files



def load_chunk_files_from_db(db_path):
    chunk_files = []
    conn = sqlite3.connect(db_path)
    try:
        cur = conn.cursor()
        cur.execute("""
            SELECT chunk_file
            FROM chunks
            ORDER BY cost DESC, chrom, chunk_idx
        """)
        rows = cur.fetchall()
        chunk_files = [row[0] for row in rows if row[0]]
    finally:
        conn.close()

    return chunk_files


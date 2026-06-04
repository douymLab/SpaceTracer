#!/usr/bin/env python3

# phylosolid/data_loader.py
"""
Data loader for PhyloSOLID pipeline.

Functions:
 - load_posterior(path)
 - load_reads(path)
 - parse_allele_entry(entry)  # internal, robust parser for allele count formats
 - derive_M_and_C_from_reads(df_reads_parsed)
 - load_features(path)
 - load_likelihoods(path)
 - load_all(inputpath) -> dict with keys: P, M, C, features, ll_mut, ll_unmut

Notes / assumptions:
 - The files in your repo (as in your example) are tab-separated text with the
   first column as row names (index). Typical usage in your notebook/script was:
       pd.read_csv(path, sep='\t', index_col=0).T
   This loader follows the same convention and returns DataFrames with:
       index = cells, columns = mutations (i.e. already transposed).
 - For allele-count file parsing, supported cell entry formats include:
     "ref:alt", "ref,alt", "ref/alt", "ref|alt", "ref;alt", "ref alt", or two numeric cols.
   If a single numeric value is present:
     - if it's an integer -> treated as coverage (alt unknown) -> M set to nan, C set to that int.
     - if it's float between 0 and 1 -> treated as MAF -> M set and C set to 1.
   Missing values (NA, ".", "") -> coverage 0 (uncovered).
 - If your real files use a different format, adjust parse_allele_entry accordingly.
"""
import os
import re
import math
import numpy as np
import pandas as pd
from typing import Tuple, Optional, Dict, Any

# -------------------------
# Low-level parsing helpers
# -------------------------
def _is_integer_string(s: str) -> bool:
    try:
        int(s)
        return True
    except Exception:
        return False

def _is_float_string(s: str) -> bool:
    try:
        float(s)
        return True
    except Exception:
        return False

def parse_allele_entry(entry: Any) -> Tuple[Optional[int], Optional[int], Optional[float]]:
    """
    Parse one cell x mutation entry from allele count matrix.
    Returns (ref_count, alt_count, maf)
      - if ref_count and alt_count both available -> maf = alt / (ref+alt), coverage = ref+alt
      - if single integer -> treat as coverage (ref unknown) -> ref=None, alt=None, maf=None, coverage=that int (handled by caller)
      - if single float 0..1 -> treat as maf -> coverage=1
      - if NA or '.', return (None, None, None) and caller will set coverage=0
    Supported separators: ':', ',', '/', '|', ';', whitespace
    """
    if pd.isna(entry):
        return None, None, None
    s = str(entry).strip()
    if s == "" or s == "." or s.lower() == "nan":
        return None, None, None
    
    # Try common two-number separators (alt/total format)
    for sep in [":", ",", "/", "|", ";", " "]:
        if sep in s:
            parts = [p for p in re.split(re.escape(sep), s) if p != ""]
            if len(parts) >= 2:
                # try parse as alt/total format (VCF style)
                a, b = parts[0], parts[1]
                if _is_integer_string(a) and _is_integer_string(b):
                    alt = int(a)
                    total = int(b)
                    if total >= alt and total >= 0:
                        ref = total - alt
                        maf = alt / total if total > 0 else 0.0
                        return ref, alt, maf
                    else:
                        return None, None, None
                # allow floats too
                if _is_float_string(a) and _is_float_string(b):
                    alt = float(a)
                    total = float(b)
                    if total >= alt and total >= 0:
                        ref = total - alt
                        maf = alt / total if total > 0 else 0.0
                        return int(round(ref)), int(round(alt)), maf
                    else:
                        return None, None, None
    
    # Single integer: treat as coverage
    if _is_integer_string(s):
        return int(s), None, None
    
    # Single float: treat as MAF
    if _is_float_string(s):
        val = float(s)
        if 0.0 <= val <= 1.0:
            return None, None, float(val)
        else:
            # If outside 0-1 range, treat as coverage
            return int(round(val)), None, None
    
    return None, None, None

# -------------------------
# Loaders
# -------------------------
def _read_and_transpose(path: str, sep: str = "\t", dtype=None) -> pd.DataFrame:
    """
    Helper: read TSV-like with index_col=0 and transpose, keeping original row/col names.
    """
    df = pd.read_csv(path, sep=sep, index_col=0, dtype=dtype)
    return df.T

def load_posterior(posterior_path: str) -> pd.DataFrame:
    """
    Load posterior matrix (cells x mutations). Returns DataFrame (index=cells, columns=mutations).
    """
    if not os.path.exists(posterior_path):
        raise FileNotFoundError(posterior_path)
    df = _read_and_transpose(posterior_path)
    # try cast to float
    df = df.astype(float)
    return df

def load_likelihoods(ll_path: str) -> pd.DataFrame:
    """
    Load likelihood matrix (cells x mutations)
    """
    if not os.path.exists(ll_path):
        raise FileNotFoundError(ll_path)
    df = _read_and_transpose(ll_path)
    # numeric
    df = df.apply(pd.to_numeric, errors='coerce')
    return df

def load_features(features_path: str) -> pd.DataFrame:
    """
    Load features.preprocess_items.txt and transpose to (cells x features)
    """
    if not os.path.exists(features_path):
        raise FileNotFoundError(features_path)
    df = _read_and_transpose(features_path)
    return df

def load_reads(reads_path: str) -> pd.DataFrame:
    """
    Load allele count matrix and transpose to (cells x mutations).
    Returns raw DataFrame of strings/numbers (not yet converted to M/C).
    """
    if not os.path.exists(reads_path):
        raise FileNotFoundError(reads_path)
    # read as string to preserve formats like "12:3"
    df = pd.read_csv(reads_path, sep="\t", index_col=0, dtype=str)
    return df.T

def derive_MCA_from_reads(df_reads: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    输入: df_reads (cells x muts)，每个元素可能是:
      - "alt:total", "alt/total", "alt,total" 等 (alt/total格式)
      - 单个 coverage (整数)
      - 单个 maf (0~1 之间浮点数)
      - NA, "." -> 覆盖度=0
    返回:
      - V: MAF(VAF, mutant/variant allele frequency) 矩阵 (float, NaN 表示无法计算)
      - C: coverage 矩阵 (int) = total
      - A: alt count 矩阵 (int) = alt
    """
    cells = df_reads.index
    muts = df_reads.columns
    V = pd.DataFrame(index=cells, columns=muts, dtype=float)
    C = pd.DataFrame(index=cells, columns=muts, dtype=int)
    A = pd.DataFrame(index=cells, columns=muts, dtype=int)
    
    for i in cells:
        for j in muts:
            entry = df_reads.at[i, j]
            ref, alt, maf = parse_allele_entry(entry)
            
            if ref is None and alt is None and maf is None:
                # uncovered
                V.at[i, j] = np.nan
                C.at[i, j] = 0
                A.at[i, j] = 0
            
            elif maf is not None and ref is None and alt is None:
                # direct maf
                V.at[i, j] = float(maf)
                C.at[i, j] = 1
                A.at[i, j] = int(round(maf))
            
            elif ref is not None and alt is not None:
                # alt/total format: ref = total - alt, alt = alt
                total = ref + alt
                if total == 0:
                    V.at[i, j] = np.nan
                    C.at[i, j] = 0
                    A.at[i, j] = 0
                else:
                    V.at[i, j] = alt / total
                    C.at[i, j] = total
                    A.at[i, j] = alt
            
            elif ref is not None and alt is None:
                # only coverage (total)
                C.at[i, j] = int(ref)
                V.at[i, j] = np.nan
                A.at[i, j] = 0
            
            else:
                # fallback
                C.at[i, j] = 0
                V.at[i, j] = np.nan
                A.at[i, j] = 0
    
    return V, C, A

# -------------------------
# High level loader
# -------------------------
def load_all(inputpath: str) -> Dict[str, pd.DataFrame]:
    """
    Read the standard set of input files from a directory (matching your example):
      - data.posterior_matrix.txt
      - data.allele_count.txt
      - features.preprocess_items.txt
      - data.likelihood_mut_matrix.txt
      - data.likelihood_unmut_matrix.txt
    
    Returns a dict with keys: P, V, C, features, ll_mut, ll_unmut
    """
    files = {
        'posterior': os.path.join(inputpath, "data.posterior_matrix.txt"),
        'reads': os.path.join(inputpath, "data.allele_count.txt"),
        'features': os.path.join(inputpath, "features.preprocess_items.txt"),
        'll_mut': os.path.join(inputpath, "data.likelihood_mut_matrix.txt"),
        'll_unmut': os.path.join(inputpath, "data.likelihood_unmut_matrix.txt"),
    }
    
    P = load_posterior(files['posterior'])
    df_reads = load_reads(files['reads'])
    V, C, A = derive_MCA_from_reads(df_reads)
    features = load_features(files['features'])
    ll_mut = load_likelihoods(files['ll_mut'])
    ll_unmut = load_likelihoods(files['ll_unmut'])
    
    # Sanity checks: align indices / columns (cells x muts)
    # We will intersect cell sets and mutation sets to get consistent matrices
    cells = sorted(set(P.index) & set(V.index) & set(C.index) & set(features.index) & set(ll_mut.index) & set(ll_unmut.index))
    muts = sorted(set(P.columns) & set(V.columns) & set(C.columns) & set(ll_mut.columns) & set(ll_unmut.columns))
    
    if len(cells) == 0 or len(muts) == 0:
        # if full intersection yields empty, fallback to using P's axes and try to reindex others with union (keeping NaNs)
        cells = list(P.index)
        muts = list(P.columns)
    
    # reindex consistently (this may introduce NaNs if some matrices missing entries)
    P = P.reindex(index=cells, columns=muts)
    V = V.reindex(index=cells, columns=muts)
    A = A.reindex(index=cells, columns=muts).fillna(0).astype(int)
    C = C.reindex(index=cells, columns=muts).fillna(0).astype(int)
    # features = features.reindex(index=cells)  # features can be fewer cols; keep as-is
    ll_mut = ll_mut.reindex(index=cells, columns=muts)
    ll_unmut = ll_unmut.reindex(index=cells, columns=muts)
    
    return {
        "P": P,
        "V": V,
        "A": A,
        "C": C,
        "df_reads": df_reads,
        "features": features,
        "ll_mut": ll_mut,
        "ll_unmut": ll_unmut
    }

# -------------------------
# CLI compatibility (optional)
# -------------------------
def main():
    """
    Minimal CLI compatible with the arguments in your snippet.
    You can also import load_all from this module.
    """
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--inputpath", default="./data", type=str)
    parser.add_argument("-o", "--outputpath", default="./results", type=str)
    args = parser.parse_args()
    
    print(f"Loading from: {args.inputpath}")
    out = load_all(args.inputpath)
    print("Loaded matrices:")
    for k, v in out.items():
        if isinstance(v, pd.DataFrame):
            print(f" - {k}: shape {v.shape}")
    # save a small sanity-check preview
    if not os.path.exists(args.outputpath):
        os.makedirs(args.outputpath)
    out['P'].iloc[:5, :5].to_csv(os.path.join(args.outputpath, "preview.P_posterior.csv"))
    out['V'].iloc[:5, :5].to_csv(os.path.join(args.outputpath, "preview.M_mutantAlleleFrequency.csv"))
    out['A'].iloc[:5, :5].to_csv(os.path.join(args.outputpath, "preview.A_mutantAlleleCount.csv"))
    out['C'].iloc[:5, :5].to_csv(os.path.join(args.outputpath, "preview.C_coverage.csv"))
    print(f"Previews written to {args.outputpath}")

if __name__ == "__main__":
    main()

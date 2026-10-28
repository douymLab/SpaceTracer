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
import logging
import os
import math
from concurrent.futures import ThreadPoolExecutor, as_completed
from functools import lru_cache
import numpy as np
import pandas as pd
from typing import Tuple, Optional, Dict, Any

logger = logging.getLogger(__name__)

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

def _parse_alt_total_pair(left: str, right: str) -> Tuple[Optional[int], Optional[int], Optional[float]]:
    if _is_integer_string(left) and _is_integer_string(right):
        alt = int(left)
        total = int(right)
        if total >= alt and total >= 0:
            ref = total - alt
            maf = alt / total if total > 0 else 0.0
            return ref, alt, maf
        return None, None, None

    if _is_float_string(left) and _is_float_string(right):
        alt = float(left)
        total = float(right)
        if total >= alt and total >= 0:
            ref = total - alt
            maf = alt / total if total > 0 else 0.0
            return int(round(ref)), int(round(alt)), maf
        return None, None, None

    return None, None, None

@lru_cache(maxsize=65536)
def _parse_allele_string(s: str) -> Tuple[Optional[int], Optional[int], Optional[float]]:
    s = s.strip()
    if not s:
        return None, None, None

    lowered = s.lower()
    if s == "." or lowered == "na" or lowered == "nan":
        return None, None, None

    for sep in ["/", ":", ",", "|", ";"]:
        if sep in s:
            left, _, right = s.partition(sep)
            if left and right:
                parsed = _parse_alt_total_pair(left.strip(), right.strip())
                if parsed != (None, None, None):
                    return parsed

    parts = s.split()
    if len(parts) >= 2:
        parsed = _parse_alt_total_pair(parts[0], parts[1])
        if parsed != (None, None, None):
            return parsed

    if _is_integer_string(s):
        return int(s), None, None

    if _is_float_string(s):
        val = float(s)
        if 0.0 <= val <= 1.0:
            return None, None, float(val)
        return int(round(val)), None, None

    return None, None, None

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
    if entry is None or entry is pd.NA:
        return None, None, None
    if isinstance(entry, float) and math.isnan(entry):
        return None, None, None
    return _parse_allele_string(str(entry))

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

def _resolve_load_workers(num_workers: Optional[int], num_columns: int) -> int:
    if num_columns <= 1:
        return 1

    if num_workers is not None:
        requested = int(num_workers)
        if requested > 0:
            return min(requested, num_columns)

    cpu_count = os.cpu_count() or 1
    if cpu_count <= 1 or num_columns < 8:
        return 1
    return min(num_columns, min(cpu_count, 8))

def _parse_reads_block(values: np.ndarray, start_col: int, end_col: int) -> Tuple[int, np.ndarray, np.ndarray, np.ndarray]:
    row_count = values.shape[0]
    width = end_col - start_col
    v_block = np.full((row_count, width), np.nan, dtype=float)
    c_block = np.zeros((row_count, width), dtype=np.int32)
    a_block = np.zeros((row_count, width), dtype=np.int32)

    for local_col, col_idx in enumerate(range(start_col, end_col)):
        v_col = v_block[:, local_col]
        c_col = c_block[:, local_col]
        a_col = a_block[:, local_col]

        for row_idx, entry in enumerate(values[:, col_idx]):
            ref, alt, maf = parse_allele_entry(entry)

            if ref is None and alt is None:
                if maf is not None:
                    v_col[row_idx] = float(maf)
                    c_col[row_idx] = 1
                    a_col[row_idx] = int(round(maf))
                continue

            if ref is not None and alt is not None:
                total = ref + alt
                if total > 0:
                    v_col[row_idx] = alt / total
                    c_col[row_idx] = total
                    a_col[row_idx] = alt
                continue

            if ref is not None:
                c_col[row_idx] = int(ref)

    return start_col, v_block, c_block, a_block

def derive_MCA_from_reads(
    df_reads: pd.DataFrame,
    num_workers: Optional[int] = None,
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
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
    values = df_reads.to_numpy(dtype=object, copy=False)
    worker_count = _resolve_load_workers(num_workers, len(muts))

    if worker_count == 1:
        blocks = [_parse_reads_block(values, 0, len(muts))]
    else:
        chunk_size = math.ceil(len(muts) / worker_count)
        ranges = [
            (start_col, min(start_col + chunk_size, len(muts)))
            for start_col in range(0, len(muts), chunk_size)
        ]
        logger.info(
            "Parsing allele count matrix with %d worker threads across %d mutation columns",
            worker_count,
            len(muts),
        )
        blocks = []
        with ThreadPoolExecutor(max_workers=worker_count) as executor:
            futures = [
                executor.submit(_parse_reads_block, values, start_col, end_col)
                for start_col, end_col in ranges
            ]
            for future in as_completed(futures):
                blocks.append(future.result())

        blocks.sort(key=lambda item: item[0])

    v_values = np.concatenate([block[1] for block in blocks], axis=1)
    c_values = np.concatenate([block[2] for block in blocks], axis=1)
    a_values = np.concatenate([block[3] for block in blocks], axis=1)

    V = pd.DataFrame(v_values, index=cells, columns=muts)
    C = pd.DataFrame(c_values, index=cells, columns=muts)
    A = pd.DataFrame(a_values, index=cells, columns=muts)
    return V, C, A

# -------------------------
# High level loader
# -------------------------
def load_all(inputpath: str, load_workers: Optional[int] = None) -> Dict[str, pd.DataFrame]:
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
    V, C, A = derive_MCA_from_reads(df_reads, num_workers=load_workers)
    V_raw_reads = V.copy()
    C_raw_reads = C.copy()
    A_raw_reads = A.copy()
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
        "V_raw_reads": V_raw_reads,
        "A_raw_reads": A_raw_reads,
        "C_raw_reads": C_raw_reads,
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
    parser.add_argument("--load_workers", default=None, type=int)
    args = parser.parse_args()
    
    print(f"Loading from: {args.inputpath}")
    out = load_all(args.inputpath, load_workers=args.load_workers)
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

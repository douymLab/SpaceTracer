from collections import defaultdict
from math import log10
import os
import statistics
from typing import Dict, List, Tuple
import numpy as np
import csv
import os
import psutil


def check_dir(dir):
    if os.path.exists(dir):
        pass
    else:
        os.makedirs(dir,exist_ok=True)
    

def check_file(file,delete=False):
    if os.path.exists(file) and delete==True:
        os.remove(file)
        return False
    elif os.path.exists(file) and delete==False:
        return True
    else:
        return False


def parse_identifier(pos_name):
    s_item=pos_name.strip().split("_")
    chrom=str(s_item[0])
    pos=int(s_item[1])
    ref=str(s_item[2])
    alt=str(s_item[3])

    return chrom,pos,ref,alt


def round_to_nearest_bin(x,bins):
    return int(np.ceil(x / bins) * bins)


def str2dict(q):
    """Convert quality string to dictionary"""
    dict = {int(i.split(':')[0]):int(i.split(':')[1]) for i in q.split(',') if q !="NA"}
    return dict
    

def str2bool(v):
    """ensure boolean input"""
    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise ValueError('Boolean value expected.')


def refine_mean(in_list):
    try:
        return statistics.mean(in_list)
    except:
        # print("wrong in refine mean for", in_list)
        return "NA"
    

def refine_median(in_list):
    try:
        return statistics.median(in_list)
    except:
        # print("wrong in refine median for", in_list)
        return "NA"


def refine_diff(a,b):
    try:
        return a - b
    except:
        # print("wrong in diff for: a:", a,"b:",b)
        return "NA"


def get_regions(identifiers: List[str], 
                max_region_size: int = 20000,
                max_variants_per_region: int = 100) -> List[Dict]:
 
    chrom_variants = defaultdict(list)
    for identifier in identifiers:
        chrom, pos, ref, alt = parse_identifier(identifier)
        chrom_variants[chrom].append((pos, (chrom, pos, ref, alt)))
    
    regions = []
    
    for chrom, variants in chrom_variants.items():
        variants.sort(key=lambda x: x[0])
        
        current_region_variants = [variants[0][1]]
        region_start = variants[0][0]
        
        for pos, variant_tuple in variants[1:]:
            region_span = pos - region_start
            variant_count = len(current_region_variants)
            
            if (region_span > max_region_size or 
                variant_count >= max_variants_per_region):
                
                region_end = current_region_variants[-1][1]  
                
                regions.append({
                    'chrom': chrom,
                    'start': max(0, region_start - 1),
                    'end': region_end,
                    'variants': current_region_variants 
                })
                
                current_region_variants = [variant_tuple]
                region_start = pos
            else:
                current_region_variants.append(variant_tuple)
        
        if current_region_variants:
            region_end = current_region_variants[-1][1]  
            regions.append({
                'chrom': chrom,
                'start': max(0, region_start - 1),
                'end': region_end,
                'variants': current_region_variants  
            })
    
    return regions


def build_region_tasks_for_UMI_combine(
    df,
    max_region_size: int = 20000,
    max_variants_per_region: int = 100,
) -> List[Tuple[str, int, int, List[Tuple[int, str, str, bool, bool]]]]:
    """
    return:
        [
            (
                chrom,
                region_start,
                region_end,
                [
                    (pos, ref, alt, check_mosaic, check_error),
                    ...
                ]
            ),
            ...
        ]
    """

    chrom_sites = defaultdict(list)

    needed_cols = ["#chrom", "pos", "ref", "alt1", "type"]
    sub_df = df[needed_cols].copy()

    for _, row in sub_df.iterrows():
        chrom = str(row["#chrom"])
        pos = int(row["pos"])
        ref = str(row["ref"])
        alt = str(row["alt1"])

        types = {t.strip() for t in str(row["type"]).split(",")}
        check_mosaic = "candidate_somatic" in types
        check_error = "candidate_error" in types

        chrom_sites[chrom].append((pos, ref, alt, check_mosaic, check_error))

    region_tasks = []

    for chrom, sites in chrom_sites.items():
        sites.sort(key=lambda x: x[0])
        current_sites = [sites[0]]
        region_start_pos = sites[0][0]

        for site in sites[1:]:
            pos = site[0]
            region_span = pos - region_start_pos
            variant_count = len(current_sites)

            if (
                region_span > max_region_size
                or variant_count >= max_variants_per_region
            ):
                region_end = current_sites[-1][0]

                region_tasks.append((
                    chrom,
                    max(0, region_start_pos - 1),
                    region_end,
                    current_sites
                ))

                current_sites = [site]
                region_start_pos = pos
            else:
                current_sites.append(site)

        if current_sites:
            region_end = current_sites[-1][0]
            region_tasks.append((
                chrom,
                max(0, region_start_pos - 1),
                region_end,
                current_sites
            ))

    return region_tasks


def handle_p_value_log10(p_value):
    try:
        if not p_value:
            return "no"

        in_type=type(p_value)
        if in_type==list:
            p_list=[]
            for p_val in p_value:
                if p_val==0 or p_val=="0":
                    p_val=1e-300
                elif p_val=="NA":
                    p_val=1
                p_list.append(log10(float(p_val)))
            return p_list
            
        elif in_type==str:
            if p_value=="0":
                p_value=1e-300
            elif p_value=="NA":
                p_value=1
            p=log10(float(p_value))
            return p
        
        elif in_type==float or in_type==np.float64:
            if p_value == 0.0:
                p_value=1e-300
            p=log10(float(p_value))
            return p
    except:
        # print("p_val:",p_value)
        return "no"
        

def barcode_cell_mapping(mapping_file):
    import pandas as pd
    if mapping_file in ['','None',None]:
        return {}
    
    elif os.path.exists(mapping_file):
        df = pd.read_csv(mapping_file, sep='\t', header=None, names=["CB", "cell"])  
        return dict(zip(df['CB'], df['cell']))
    else:
        raise FileNotFoundError(f"Mapping file '{mapping_file}' does not exist")


def list2min(value):
    """find the minimum value in the comma-separated string"""
    # check if the input is a float number
    if isinstance(value, float):
        return value
    elif isinstance(value, int):
        return value
    else:
        return min(float(num) if num!="no" else np.nan for num in value.split(','))


def as_int(x):
    return None if x is None or x == "None" else int(x)

def as_float(x):
    return None if x is None or x == "None" else float(x)

def as_str(x):
    return None if x is None or x == "None" else str(x)




def load_manifest_tsv(manifest_file: str) -> List[Dict[str, str]]:
    rows = []
    with open(manifest_file, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            rows.append(row)
    return rows


def save_manifest_tsv(rows: List[Dict[str, str]], output_file: str):
    valid_rows = [row for row in rows if isinstance(row, dict) and row]

    out_dir = os.path.dirname(output_file)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    if not valid_rows:
        print(f"[manifest] no valid rows, wrote empty manifest: {output_file}")
        with open(output_file, "w", newline="", encoding="utf-8") as f:
            f.write("")
        return

    fieldnames = []
    seen = set()
    for row in valid_rows:
        for key in row:
            if key not in seen:
                seen.add(key)
                fieldnames.append(key)

    with open(output_file, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(valid_rows)


try: #************************
    from pympler import asizeof
    HAS_PYMPLER = True
except Exception:
    HAS_PYMPLER = False

_PROC = psutil.Process(os.getpid())



def _proc_mem_mb():
    rss = 0.0
    uss = None
    vms = None
    try:
        proc = psutil.Process(os.getpid())
        mi = proc.memory_info()
        rss = mi.rss / 1024 / 1024
        vms = mi.vms / 1024 / 1024
    except Exception:
        pass

    try:
        proc = psutil.Process(os.getpid())
        full = proc.memory_full_info()
        uss_val = getattr(full, "uss", None)
        if uss_val is not None:
            uss = uss_val / 1024 / 1024
    except Exception:
        pass

    return rss, uss, vms

def log_worker_mem(logger, label, extra=""):
    rss, uss, vms = _proc_mem_mb()
    uss_text = f"{uss:.1f}MB" if uss is not None else "N/A"
    vms_text = f"{vms:.1f}MB" if vms is not None else "N/A"
    logger.info(f"{label}... RSS={rss:.1f}MB USS={uss_text} VMS={vms_text} extra={extra}...")

def obj_size_mb(obj):
    if not HAS_PYMPLER:
        return None
    try:
        return asizeof.asizeof(obj) / 1024 / 1024
    except Exception:
        return None

def obj_size_text(obj, name):
    size = obj_size_mb(obj)
    if size is None:
        return f"{name}=N/A"
    return f"{name}={size:.2f}MB"


def handel_hd_barcode_mapping(raw_barcode_mappings_parquet,bin_size):
    import pandas as pd 
    barcode_mapping=pd.read_parquet(raw_barcode_mappings_parquet)
    barcode_mapping.to_csv()

    
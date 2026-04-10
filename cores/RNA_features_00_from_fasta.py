
from SpaceTracer.utils.get_MutationType import *




from functools import partial
import multiprocessing
import os
import subprocess
import pandas as pd
from pyfaidx import Fasta
from SpaceTracer.utils.get_genome_info import get_chr_size


def find_nearest_polymer_with_boundary_check(context_20bp, mutation_pos_in_context, min_repeat):
    bases = ["A","T","C","G"]
    nearest_polymer = None
    min_distance = float('inf')
    polymer_start = polymer_end = -1  # 记录polymer的位置
    
    for base in bases:
        poly = base * min_repeat
        start_idx = 0
        while True:
            start_idx = context_20bp.find(poly, start_idx)
            if start_idx == -1:
                break
                
            end_idx = start_idx + len(poly) - 1
            
            if start_idx <= mutation_pos_in_context <= end_idx:
                # 突变在polymer内部，记录polymer位置并立即返回
                polymer_start, polymer_end = start_idx, end_idx
                return (base, 0, polymer_start, polymer_end)
            
            # 计算外部polymer的距离
            if end_idx < mutation_pos_in_context:  # polymer在左侧
                distance = end_idx - mutation_pos_in_context  # 负值
            else:  # polymer在右侧
                distance = start_idx - mutation_pos_in_context  # 正值
                
            if abs(distance) < abs(min_distance):
                min_distance = distance
                nearest_polymer = base
                polymer_start, polymer_end = start_idx, end_idx
                
            start_idx += 1
    
    if nearest_polymer:
        return (nearest_polymer, min_distance, polymer_start, polymer_end)
    return (None, None, -1, -1)


def add_features_from_fasta(df: pd.DataFrame, reference_fasta: str, previous_base=5) -> pd.DataFrame:
    reference = Fasta(reference_fasta) 
    chr_size = get_chr_size(reference_fasta + ".fai")
    base = {'A':'T', 'T':'A', 'G':'C', 'C':'G'}
    min_repeat=previous_base+1
    def _process_one(chrom, pos, ref, alt, strand):
        context_20bp=str(reference[chrom][max(1,int(pos)-11):min(int(pos)+10,int(chr_size[chrom]))])
        GCcontent=(context_20bp.count('G')+context_20bp.count('C'))/len(context_20bp)
        # context=reference[chrom][int(pos)-2:int(pos)+1]
        # context2=(base[str(reference[chrom][int(pos)-2:int(pos)-1])]+base[str(reference[chrom][int(pos)-1:int(pos)])]+base[str(reference[chrom][int(pos):int(pos)+1])])[::-1]
        # context_10bp=str(reference[chrom][max(1,int(pos)-6):min(int(pos)+5,int(chr_size[chrom]))])
        
        ref=ref.split(",")[0]
        up_base=str(reference[chrom][pos-1-1])
        down_base=str(reference[chrom][pos+1-1])
        if ref in "CT":
            sub_name=ref+">"+alt
            dna_stat_name=up_base+"["+sub_name+"]"+down_base
        else:
            sub_name=base[ref]+">"+base[alt]
            dna_stat_name=base[down_base]+"["+ sub_name + "]"+ base[up_base]
        
        if strand=="+":
            sub_name=ref+">"+alt
            rna_stat_name=up_base+"["+sub_name+"]"+down_base
        else:
            sub_name=base[ref]+">"+base[alt]
            rna_stat_name=base[down_base]+"["+ sub_name + "]"+ base[up_base]
        
        DNAMutationType=dna_stat_name
        RNAMutationType=rna_stat_name

        # is followed by ploy alt reads?
        # if strand =="-":
        #     previous_bases=str(reference[chrom][max(1,int(pos)-1-previous_base):int(pos)-1])
        # else:
        #     previous_bases=str(reference[chrom][int(pos):min(int(pos)+previous_base,int(chr_size[chrom]))])
        # if str(alt)*previous_base==previous_bases:
        #     # print(str(alt)*3)
        #     equal_to_previous_bases=True
        # else:
        #     equal_to_previous_bases=False

        changed_context_bases=str(reference[chrom][max(1,int(pos)-1-previous_base):int(pos)-1])+alt+str(reference[chrom][int(pos):min(int(pos)+previous_base,int(chr_size[chrom]))])
        cause_poly_alt=str(alt)*(previous_base+1) in changed_context_bases

        pos_in_context = 10
        
        polymer_type, distance, poly_start, poly_end = find_nearest_polymer_with_boundary_check(context_20bp, pos_in_context, min_repeat)
        equal_to_flanking=False
        if polymer_type:
            if distance == 0:
                left_dist = pos_in_context - poly_start
                right_dist = poly_end - pos_in_context
                # distance_to_boundary = min(left_dist, right_dist)
        
                if left_dist <= right_dist:
                    boundary_pos = poly_start - 1  # 左边界外
                    # boundary_direction = "left"
                else:
                    boundary_pos = poly_start + 1    # 右边界外
                    # boundary_direction = "right"
        
                # 检查边界外1-2bp的碱基是否与alt相同
                alt_matches = False
                if 0 <= boundary_pos < len(context_20bp):  # 确保边界位置有效
                    boundary_base = context_20bp[boundary_pos]
                    alt_matches = (alt == boundary_base)
                if alt_matches==True:
                    equal_to_flanking=True
                # return (distance_to_boundary, alt_matches, boundary_direction)
            
        # if polymer_type:
        #     if distance == 0:
        #         print(f"Mutation is INSIDE a poly{polymer_type} region, and {equal_to_flanking} equal_to_flanking")
        #     else:
        #         print(f"Nearest poly{polymer_type} region is at distance {distance} ({'left' if distance < 0 else 'right'})")
        # else:
        #     print(f"No polymer region (≥{min_repeat}bp) found nearby")

        return DNAMutationType, RNAMutationType,GCcontent,cause_poly_alt,polymer_type

    results = []
    for row in df.itertuples():
        result = _process_one(row.chrom, row.pos, row.ref, row.alt, row.strand)
        results.append(result)
    
    result_df = pd.DataFrame(results, index=df.index)
    return result_df 


from collections import defaultdict
from functools import partial
import re
from typing import List, Dict
import pysam
import random
import numpy as np
import os
import statsmodels.stats.multitest as smm

# from SpaceTracer.cores.read_feature_02_extract import process_reads_for_variant
from SpaceTracer.utils.get_read_level_feature import check_UMIconsistence_for_each_geno, combine_info_from_cigar, do_wilicox_sum_test, get_indel_info, judge_pos_in_indel
from SpaceTracer.utils.handle_UMI_combine import handle_pos, handle_quality_matrix, handle_seq,calculate_UMI_combine_phred, get_most_candidate_allele, handle_seq_type


def handel_bam_file_for_region(bam_file,region_dict,run_type,bins,cell_dict={},readLen=120,downsample=False,target_depth=2000,seed=42):
    region_reads = []

    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for read in bam.fetch(region_dict['chrom'], region_dict['start'], region_dict['end']):
            region_reads.append(read)

    var_result_dict={}
    for var in region_dict['variants']:
        var_reads = [r for r in region_reads 
                    if r.reference_start <= var[1] - 1 < r.reference_end]
        
        # downsample
        total_count = len(var_reads)
        if downsample and total_count > target_depth:
            random.seed(seed)
            sampled_reads = random.sample(var_reads, target_depth)
        else:
            sampled_reads = var_reads

        var_result_dict[var]=process_reads_for_variant(sampled_reads,var,run_type,bins,cell_dict,readLen)
    
    return var_result_dict



def process_reads_for_variant(sampled_reads,var,run_type,bins,cell_dict={},readLen=120):
    chrom,pos,ref,alt=var

    if "," in ref:
        one_ref=ref[0]
    else:
        one_ref=ref

    result_dict={"A":defaultdict(list), "T":defaultdict(list), "C":defaultdict(list), "G":defaultdict(list), "del":defaultdict(list)}
    for geno in "ATCG":
        result_dict[geno]["dp"]=0
        result_dict[geno]["dp_consensus"]=0
        result_dict[geno]["reverse_dp"]=0
        result_dict[geno]["forward_dp"]=0
        result_dict[geno]["edist"]=[0]*readLen
        result_dict[geno]["GenoSpotNum"]=0

    dp=0
    # in_bam_read=pysam.AlignmentFile(bam_file, "rb") # , reference_filename=ref_fasta)
    pos_index = pos-1
    barcode_name=[]
    site_barcode_UMI_dict={}
    
    # if downsample:
    #     reads, original_depth = fetch_and_downsample_optimized(bam_file, chrom, pos, target_depth, seed)
    #     # print(f"Original depth: {original_depth}, Sampled reads: {len(sampled_reads)}")
    # else:
    #     reads=in_bam_read.fetch(chrom, pos-1, pos)

    for read in sampled_reads:
        barcode_name,UMI_name=handle_seq_type(read,run_type,bins,cell_dict)
        if not barcode_name:
            continue

        pos_index=pos-1
        seq_soft_cut, ins_info, del_info, seq_hard_clip = combine_info_from_cigar(read.cigar)
        cut_seq=handle_seq(read.seq, seq_soft_cut)
        cut_pos=handle_pos(read.get_reference_positions(), ins_info)
        indel_pos_list=judge_pos_in_indel(ins_info,del_info,read.get_reference_positions())

        if pos_index in cut_pos or pos in indel_pos_list:
            dp+=1
            if pos_index in indel_pos_list:
                geno="del"
                result_dict[geno]["is_indel"].append(1)
                result_dict[geno]["baseq"]=[]

            else:
                edist=cut_pos.index(pos_index)  
                geno = cut_seq[edist]  
                if geno not in "ATCG":
                    continue
                # print(geno,CB,UB)
                epos=edist/len(cut_pos)
                result_dict[geno]["epos"].append(epos)
                try:
                    result_dict[geno]["edist"][edist]+=1
                except:
                    print(edist)
                result_dict["is_indel"]=[]
                raw_index = handle_quality_matrix(cut_pos.index(pos_index),read.seq,cut_seq)
                quality=read.get_forward_qualities()[raw_index]
                result_dict[geno]["baseq"].append(quality)
                # if result_dict[geno]["dp"]!=[]:
                result_dict[geno]["dp"]+=1
                # else:
                #     result_dict[geno]["dp"]=0
            # effective_DP += 1
        
            #number_mismatch; is_reverse; mapping_quality
                number_mismatch=read.get_tag("nM"); result_dict[geno]["number_mismatch"].append(number_mismatch)
                is_reverse=read.is_reverse; result_dict[geno]["is_reverse"].append(is_reverse)
                map_q=read.mapq; result_dict[geno]["map_q"].append(map_q)
                
                number_mapper=read.get_tag("NH"); result_dict[geno]["number_mapper"].append(number_mapper)

                #soft_clip_length and hard_clip_length
                left_softclip=0 if seq_soft_cut[0]==None else seq_soft_cut[0]
                right_softclip=0 if seq_soft_cut[1]==None else len(read.seq)-seq_soft_cut[1]
                softclip_length=left_softclip+right_softclip
                result_dict[geno]["left_softclip"].append(left_softclip)
                result_dict[geno]["right_softclip"].append(right_softclip)
                result_dict[geno]["softclip_length"].append(softclip_length)

                left_hardclip,right_hardclip=seq_hard_clip[0],seq_hard_clip[1]
                hardclip_length=left_hardclip+right_hardclip
                result_dict[geno]["left_hardclip"].append(left_hardclip)
                result_dict[geno]["right_hardclip"].append(right_hardclip)
                result_dict[geno]["hardclip_length"].append(hardclip_length)

                #indel information, indel number, indel length, indel distance
                ins_num,ins_length,ins_distance=get_indel_info(ins_info,read.get_reference_positions().index(pos_index))
                del_num,del_length,del_distance=get_indel_info(del_info,read.get_reference_positions().index(pos_index))
                result_dict[geno]["ind_num"].append(ins_num+del_num)
                result_dict[geno]["ins_num"].append(ins_num)
                if ins_num==0:
                    result_dict[geno]["ins_length"].append(["no"]); result_dict[geno]["ins_distance"].append(["no"])
                else: #the ins_length and ins_distance are list format
                    result_dict[geno]["ins_length"].append(ins_length) ## append a list
                    result_dict[geno]["ins_distance"].append(ins_distance) ## append a list
                
                result_dict[geno]["del_num"].append(del_num)
                if del_num==0:
                    result_dict[geno]["del_length"].append(["no"]); result_dict[geno]["del_distance"].append(["no"])
                else: # the del_length and del_distance are list format
                    result_dict[geno]["del_length"].append(del_length)
                    result_dict[geno]["del_distance"].append(del_distance)

                # querypos(querypos_p): the distance between pos and read start (doubt: the more far away from 1st seq pos, the lower quality may have), 
                # seqpos_p cycling length, related with strand (note: next_reference_start is only work for PE); 
                # for visium, all reads are read2, so seqpos may same as the len(querypos)
                # left pos: mapping position for the reference start; 
                left_boundary=edist+left_softclip+left_hardclip
                right_boundary=len(cut_pos)-edist + right_softclip + right_hardclip
                result_dict[geno]["left_read_edist"].append(edist)
                result_dict[geno]["right_read_edist"].append(len(cut_pos)-edist)

                left_boundary_remove_clip=edist
                right_boundary_remove_clip=len(cut_pos)-edist
                result_dict[geno]["querypos"].append(left_boundary)
                result_dict[geno]["seqpos"].append(right_boundary)
                if is_reverse in [True,"TRUE","true","True"]:
                    distance_to_end=right_boundary/readLen
                    distance_to_end_value=min(right_boundary,readLen-right_boundary)
                    result_dict[geno]["reverse_dp"]+=1
                    distance_to_end_remove_clip=right_boundary_remove_clip/len(cut_pos)
                    distance_to_end_remove_clip_value=min(left_boundary_remove_clip,right_boundary_remove_clip)
                    # distance_to_end_remove_clip_save=right_boundary_remove_clip
                else:
                    distance_to_end=left_boundary/readLen
                    distance_to_end_value=min(left_boundary,readLen-left_boundary)
                    result_dict[geno]["forward_dp"]+=1
                    distance_to_end_remove_clip=left_boundary_remove_clip/len(cut_pos)
                    distance_to_end_remove_clip_value=min(left_boundary_remove_clip,right_boundary_remove_clip)
                    # distance_to_end_remove_clip_save=left_boundary_remove_clip
                # print(geno,distance_to_end_remove_clip)
                result_dict[geno]["distance_to_end"].append(distance_to_end)
                result_dict[geno]["distance_to_end_remove_clip"].append(distance_to_end_remove_clip)

                leftpos_p=read.reference_start
                rightpos_p=read.reference_end # same as leftpo, can be deleted 
                result_dict[geno]["leftpos_p"].append(leftpos_p)
                result_dict[geno]["rightpos_p"].append(rightpos_p)

                #baseq1b
                if pos_index+1 in cut_pos:
                    baseq1b=read.get_forward_qualities()[raw_index+1]
                else:
                    baseq1b=""
                result_dict[geno]["baseq1b"].append(baseq1b)
                # print(read)
                #gene information
                try:
                    result_dict[geno]["GeneID_list"].append(read.get_tag("GX"))
                except:
                    result_dict[geno]["GeneID_list"].append("no")
                try:
                    result_dict[geno]["GeneName_list"].append(read.get_tag("GN"))
                except:
                    result_dict[geno]["GeneName_list"].append("no")
                try:
                    for item in read.get_tag("TX").split(";"):
                        transcript_id,_,_=item.split(",")
                        result_dict[geno]["TransID_list"].append(transcript_id)
                except:
                    result_dict[geno]["TransID_list"].append("no")

                if barcode_name not in site_barcode_UMI_dict.keys():
                    site_barcode_UMI_dict[barcode_name]=defaultdict(dict)

                if UMI_name not in site_barcode_UMI_dict[barcode_name].keys():
                    site_barcode_UMI_dict[barcode_name][UMI_name]["count"]=defaultdict(int)
                    site_barcode_UMI_dict[barcode_name][UMI_name]["quality"]={"A":defaultdict(int),"T":defaultdict(int),"C":defaultdict(int),"G":defaultdict(int)}
                    # site_barcode_UMI_dict[barcode_name][UMI_name]["context"]=[]
                    site_barcode_UMI_dict[barcode_name][UMI_name]["end"]=[]
                    site_barcode_UMI_dict[barcode_name][UMI_name]["end_remove_clip"]=[]
                    site_barcode_UMI_dict[barcode_name][UMI_name]["end_value"]=[]
                    site_barcode_UMI_dict[barcode_name][UMI_name]["end_remove_clip_value"]=[]

                site_barcode_UMI_dict[barcode_name][UMI_name]["count"][geno]+=1
                site_barcode_UMI_dict[barcode_name][UMI_name]["quality"][geno][quality]+=1
                # site_barcode_UMI_dict[barcode_name][UMI_name]["context"].append(cut_seq[max(0,read_index-4):min(read_index+5,len(cut_seq))])
                site_barcode_UMI_dict[barcode_name][UMI_name]["end"].append(distance_to_end)
                site_barcode_UMI_dict[barcode_name][UMI_name]["end_remove_clip"].append(distance_to_end_remove_clip)
                site_barcode_UMI_dict[barcode_name][UMI_name]["end_value"].append(distance_to_end_value)
                site_barcode_UMI_dict[barcode_name][UMI_name]["end_remove_clip_value"].append(distance_to_end_remove_clip_value)
                
    for barcode in site_barcode_UMI_dict.keys():
        read_have_alt=False
        read_number_per_spot=0
        UMI_number_per_spot=0
        alt_UMI_number_per_spot=0
        UMI_count_by_allele=[0,0,0,0]
        # UMI_dp+=len(site_barcode_UMI_dict[barcode].keys())
        for UMI in site_barcode_UMI_dict[barcode]:             
            count_dict=site_barcode_UMI_dict[barcode][UMI]["count"]
            quality_dict=site_barcode_UMI_dict[barcode][UMI]["quality"]
            phred_dict=calculate_UMI_combine_phred(count_dict,quality_dict,weigh=0.5)
            candidate_allele,phred=get_most_candidate_allele(phred_dict,one_ref)
            result_dict[candidate_allele]["dp_consensus"]+=1
            UMI_count_by_allele["ATCG".index(candidate_allele)]+=1
            threshold=1
            norm_count=check_UMIconsistence_for_each_geno(count_dict,threshold)
            norm_count_remove_single_read=check_UMIconsistence_for_each_geno(count_dict,2)

            if norm_count!=[]:
                for geno,prop in zip("ATCG",norm_count):
                    result_dict[geno]["UMI_consistence_prop"].append(prop)

            if norm_count_remove_single_read!=[]:
                for geno,prop in zip("ATCG",norm_count_remove_single_read):
                    result_dict[geno]["UMI_consistence_prop_remove_single_read"].append(prop)

            for geno in "ATCG":
                if site_barcode_UMI_dict[barcode][UMI]["count"][geno]!=0:
                    # print([count_dict["A"],count_dict["T"],count_dict["C"],count_dict["G"]])
                    result_dict[geno]["read_number_per_UMI"].append(count_dict[geno])
                    read_number_per_spot+=site_barcode_UMI_dict[barcode][UMI]["count"][geno]
                    result_dict[geno]["base_proportion_per_UMI"].append(count_dict[geno]/sum(count_dict.values()))
                    
            UMI_number_per_spot+=1
            end=np.median(site_barcode_UMI_dict[barcode][UMI]["end"])
            end_remove_clip=np.median(site_barcode_UMI_dict[barcode][UMI]["end_remove_clip"])
            end_value=np.median(site_barcode_UMI_dict[barcode][UMI]["end_value"])
            end_remove_clip_value=np.median(site_barcode_UMI_dict[barcode][UMI]["end_remove_clip_value"])

            if candidate_allele==alt:
                read_have_alt=True
                alt_UMI_number_per_spot+=1

            result_dict[candidate_allele]["per_UMI_end"].append(end)
            result_dict[candidate_allele]["per_UMI_end_remove_clip"].append(end_remove_clip)
            result_dict[candidate_allele]["per_UMI_end_value"].append(end_value)
            result_dict[candidate_allele]["per_UMI_end_remove_clip_value"].append(end_remove_clip_value)
        
        if read_have_alt==True:
            # print(barcode)
            result_dict[alt]["GenoSpotNum"]+=1 
            result_dict[alt]["total_read_number_per_spot"].append(read_number_per_spot)
            result_dict[alt]["total_UMI_number_per_spot"].append(UMI_number_per_spot)
            result_dict[alt]["UMI_end"].append(end)
            result_dict[alt]["UMI_end_remove_clip"].append(end_remove_clip)
            result_dict[alt]["UMI_end_value"].append(end_value)
            result_dict[alt]["UMI_end_remove_clip_value"].append(end_remove_clip_value)
            result_dict[alt]["vaf_spot"].append(alt_UMI_number_per_spot/UMI_number_per_spot)

        else:
            result_dict[one_ref]["GenoSpotNum"]+=1
            result_dict[one_ref]["total_read_number_per_spot"].append(read_number_per_spot)
            result_dict[one_ref]["total_UMI_number_per_spot"].append(UMI_number_per_spot)
            result_dict[one_ref]["UMI_end"].append(end)
            result_dict[one_ref]["UMI_end_remove_clip"].append(end_remove_clip)
            result_dict[one_ref]["UMI_end_value"].append(end_value)
            result_dict[one_ref]["UMI_end_remove_clip_value"].append(end_remove_clip_value)
        

        # print(end_value,end_remove_clip_value)
        for geno,count in zip("ATCG",UMI_count_by_allele):
            result_dict[geno]["UMI_number_per_spot"].append(count)

    result_dict['dp']=dp
    return result_dict




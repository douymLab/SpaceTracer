from collections import Counter, defaultdict
from functools import partial
import re
from typing import List, Dict, Optional
import pysam
import random
import numpy as np
from dataclasses import dataclass, field, fields
import os
from typing import Any, Dict
import statsmodels.stats.multitest as smm
import scipy

from SpaceTracer.utils.get_read_level_feature import calculate_rbc_for_paired_wilcoxon, check_UMIconsistence_for_each_geno, combine_info_from_cigar, do_wilicox_sum_test, get_indel_info, judge_pos_in_indel, wilcoxon_with_rbc
from SpaceTracer.utils.handle_UMI_combine import handle_pos, handle_quality_matrix, handle_seq,calculate_UMI_combine_phred, get_most_candidate_allele, handle_seq_type
from SpaceTracer.utils.utils import handle_p_value_log10, parse_identifier, refine_diff, refine_mean, refine_median, round_to_nearest_bin


@dataclass
class readLevelFeatures:
    identifier: str
    
    chrom: str = field(init=False)
    pos: int = field(init=False)
    ref: str = field(init=False)
    alt: str = field(init=False)
    # mut_origin: str = "NA"

    baseq_p: float = None
    baseq_p_adj: float = None
    baseq_rbc: float = None

    alt_baseq1b_p: float = None
    alt_baseq1b_p_adj: float = None
    alt_baseq1b_rbc: float = None

    querypos_p: float = None
    querypos_p_adj: float = None
    querypos_rbc: float = None

    leftpos_p: float = None
    leftpos_p_adj: float = None
    leftpos_rbc: float = None

    seqpos_p: float = None
    seqpos_p_adj: float = None
    seqpos_rbc: float = None

    alt_querypos_num: int = None
    per_alt_UMI_end_remove_clip_mean: int = None
    per_alt_UMI_end_remove_clip_median: int = None
    per_UMI_end_remove_clip_p: float = None
    per_UMI_end_remove_clip_p_adj: float = None
    per_UMI_end_remove_clip_rbc: float = None

    ref_mismatches_mean: float = None
    alt_mismatches_mean: float = None
    mismatches_p: float = None
    mismatches_p_adj: float = None
    mismatches_rbc: float = None

    mapq_p: float = None
    mapq_p_adj: float = None
    mapq_rbc: float = None
    ref_mapq_mean: float = None
    alt_mapq_mean: float = None
    mapq_mean: float = None

    alt_UMI_consistence_prop: float = None
    alt_UMI_avg_consistence_remove_single_read: float = None
    alt_consistence_hard_prop: float = None
    alt_consistent_UMI_prop_strict_remove_single_read: float = None
    alt_consistence_soft_prop: float = None
    alt_consistent_UMI_prop_relaxed_remove_single_read: float = None

    read_number_p: float = None
    read_number_p_adj: float = None
    read_number_rbc: float = None
    alt_read_number_perUMI_max: int = None
    indel_proportion_for_site: float = None

    reads_with_indel_p: float = None
    reads_with_indel_p_adj: float = None
    alt2_proportion_per_UMI: float = None

    multi_mapper_p: float = None
    multi_mapper_p_adj: float = None
    ref_multi_map_prop: float = None
    alt_multi_map_prop: float = None
    multi_map_prop: float = None
    
    strand_bias_p: float = None
    ref_softclip_prop: float = None
    alt_softclip_prop: float = None
    
    softclip_prop: float = None
    softclip_prop_p: float = None
    softclip_prop_p_adj: float = None

    softclip_length_p: float = None
    softclip_length_p_adj: float = None
    softclip_length_rbc: float = None


    # dp
    # dp_consensus
    # ref_allele_count
    # consensus_ref_allele_count
    # alt_allele_count
    # consensus_alt_allele_count
    # alt2_allele_count
    # consensus_alt2_allele_count
    # alt2_proportion_consensus: float = None
    # af

    # vaf_spot_mean
    # vaf_spot_median
    # AFind
    # mean_AFspot
    # spotNum
    # alt_SpotNum
    # DNAMutationType
    # RNAMutationType
    # cause_ploy_alt
    # mosaic_likelihood
    # combine_nearest_phase_haplotype
    # combine_nearest_info_mutant_prop
    # combine_nearest_discordant_prop
    # combine_nearest_phase_distance
    # combine_most_phase_haplotype
    # combine_most_info_mutant_prop
    # combine_most_discordant_prop
    # combine_most_phase_distance
    # GCcontent
    # fref
    # falt
    # anno
    # anno_gene
    # major_read_strand
    # r2
    # wilcoxon_p
    alt_vs_total_dp_paired_wilcoxon_rbc: float = None
    # KS_p
    # KS_s
    # MI_p
    # MI_s
    # mut_rate_prob
    # mut_rate_vaf
    # sf_test_sig
    
    # 定义不导出的字段
    EXCLUDE_FIELDS = {'chrom', 'pos', 'ref', 'alt'}
    
    def __post_init__(self):
        mut_chrom, mut_pos, mut_ref, mut_alt = parse_identifier(self.identifier)
        self.chrom = mut_chrom
        self.pos = mut_pos
        self.ref = mut_ref
        self.alt = mut_alt
  
    def to_dict(self, mode: str = "normal", exclude: set = None) -> Dict[str, Any]:
        if mode == "test":
            return self.test_values()
        
        if exclude is None:
            exclude = self.EXCLUDE_FIELDS
        
        result = {}
        for f in fields(self):
            if f.name in exclude:
                continue
            
            value = getattr(self, f.name)
            
            if f.name.endswith('_p') or f.name.endswith('_p_adj'):
                value = handle_p_value_log10(value)
            
            result[f.name] = value
        
        return result
    
    def test_values(self) -> Dict[str, Any]:
        return {
            'identifier': self.identifier
        }

    @classmethod
    def from_read_info(cls, identifier: str, read_info_dict: Dict[str, Any]):
        if not read_info_dict:
            return None
        
        if read_info_dict.get('dp', 0) == 0:
            return None
        
        feature = cls(identifier=identifier)
        feature._add_read_info(read_info_dict)
        return feature

    @classmethod
    def from_read_info_to_dict(cls, identifier: str, read_info_dict: Dict[str, Any]) -> Optional[Dict[str, Any]]:
        feature = cls.from_read_info(identifier, read_info_dict)
        return feature.to_dict() if feature is not None else None

    def _run(self, read_info_dict: Dict[str, Any]):
        if not read_info_dict:
            return None
        
        if read_info_dict.get('dp', 0) == 0:
            return None
        
        self._add_read_info(read_info_dict)
        return self
        
        # if os.path.exists(bam_file):
        #     read_info_dicts=handel_bam_file_for_region(region_dict,run_type,bins,cell_dict,readLen,downsample,target_depth,seed)
        #     for read_info_dict in read_info_dicts:
        #         if read_info_dict['dp']!=0:
        #             mutation_features=ReadFeature()
        #             mutation_features.add_read_info(read_info_dict,dp)

            # read_info_dict,dp= handel_bam_file(bam_file,mut_chrom,mut_pos,mut_ref[0],mut_alt,run_type,bins,cell_dict,readLen,downsample,targe_dp,seed)
        # print(read_info_dict)
        # if dp!=0:
        #     mutation_features.add_read_info(read_info_dict,dp)

    def _add_read_info(self,read_info_dict):

        dp=read_info_dict['dp']
        def simple_to_get_list(read_info_dict,ref_allele,alt_allele,var):
            # print(ref_allele,alt_allele, ref_allele.split(","))
            return_ref_list=[float(k) for allele in ref_allele.split(",") for k in read_info_dict[allele][var] if k !=""]
            return_alt_list=[float(k) for allele in alt_allele.split(",") for k in read_info_dict[allele][var] if k !=""]
            return return_ref_list,return_alt_list
            
        get_list=partial(simple_to_get_list,read_info_dict,self.ref,self.alt)

        odds_eps=0

        # # gene id and trancript id
        # GeneID_list=[]
        # GeneName_list=[]
        # TransID_list=[]
        # for allele in "ATCG":
        #     GeneID_list += read_info_dict[allele]["GeneID_list"]
        #     GeneName_list += read_info_dict[allele]["GeneName_list"]
        #     TransID_list += read_info_dict[allele]["TransID_list"]

        # def tidy_result(in_list):
        #     a_list=list(set(in_list))
        #     if "no" in a_list:
        #         a_list.remove("no")
        #     if a_list is None:
        #         a_list=[]
        #     b_list=[]    
        #     for m in a_list:
        #         if ";" in m or "," in m:
        #             m_a=re.split(";|,",m)
        #             b_list.append(m_a)
        #         else:
        #             b_list.append(m)
        #     a_list=b_list
        #     return a_list
        # self.gene_id=tidy_result(GeneID_list)
        # self.gene_name=tidy_result(GeneName_list)
        # self.transcript_name=tidy_result(TransID_list)

        ref_mismatches, alt_mismatches=get_list("number_mismatch")
        refine_alt_mismatches=[int(i)-1 if i!=0 else 0 for i in alt_mismatches]
        
        #number_multi_mapper
        ref_mappers, alt_mappers=get_list("number_mapper")
        refine_ref_mappers_uniq = len([i for i in ref_mappers if i == 1]); refine_ref_mappers_multi=len(ref_mappers)-refine_ref_mappers_uniq
        refine_alt_mappers_uniq = len([i for i in alt_mappers if i == 1]); refine_alt_mappers_multi=len(alt_mappers)-refine_alt_mappers_uniq
        multi_mapper_odds, self.multi_mapper_p = scipy.stats.fisher_exact([[refine_alt_mappers_multi, refine_alt_mappers_uniq ],[refine_ref_mappers_multi, refine_ref_mappers_uniq]])
        self.ref_multi_map_prop=refine_ref_mappers_multi/len(ref_mappers) if len(ref_mappers)!=0 else "NA"
        self.alt_multi_map_prop=refine_alt_mappers_multi/len(alt_mappers) if len(alt_mappers)!=0 else "NA"
        all_mappers=len(ref_mappers)+len(alt_mappers)
        self.multi_map_prop=(refine_ref_mappers_multi+refine_alt_mappers_multi)/all_mappers if all_mappers!=0 else "NA"

        self.ref_mismatches_mean=refine_mean(ref_mismatches)
        self.alt_mismatches_mean=refine_mean(refine_alt_mismatches)
        # self.mismatches_s,self.mismatches_p=do_wilicox_sum_test(ref_mismatches,refine_alt_mismatches,type="list")
        mismatches_s,self.mismatches_p,self.mismatches_rbc=wilcoxon_with_rbc(ref_mismatches,refine_alt_mismatches)
        
        ref_is_reverse, alt_is_reverse=get_list("is_reverse")
        refine_ref_is_reverse = len([i for i in ref_is_reverse if i == 1]); refine_ref_is_forward = len(ref_is_reverse)-refine_ref_is_reverse
        refine_alt_is_reverse = len([i for i in alt_is_reverse if i == 1]); refine_alt_is_forward = len(alt_is_reverse)-refine_alt_is_reverse
        strand_bias_odds, self.strand_bias_p = scipy.stats.fisher_exact([[refine_alt_is_reverse, refine_alt_is_forward ],[refine_ref_is_reverse, refine_ref_is_forward]])

        ref_mapq,alt_mapq=get_list("map_q")
        self.ref_mapq_mean=refine_mean(ref_mapq)
        self.alt_mapq_mean=refine_mean(ref_mapq)
        self.mapq_mean=refine_mean(ref_mapq+alt_mapq)
        mapq_s,self.mapq_p,self.mapq_rbc=wilcoxon_with_rbc(ref_mapq,alt_mapq,alternative="greater")

        # ref_map_q, alt_map_q=get_list("map_q")
        # self.mapq_s,self.mapq_p=do_wilicox_sum_test(ref_map_q,alt_map_q,method="greater",type="list")
        # self.mapq_difference=refine_diff(refine_mean(ref_mapq),refine_mean(alt_mapq))
        ref_ind_num,alt_ind_num=get_list("ind_num")
        ref_read_with_indel = len([i for i in ref_ind_num if i > 0]); ref_read_without_indel = len(ref_ind_num)-ref_read_with_indel
        alt_read_with_indel = len([i for i in alt_ind_num if i > 0]); alt_read_without_indel = len(alt_ind_num)-alt_read_with_indel
        reads_with_indel_odds, self.reads_with_indel_p = scipy.stats.fisher_exact([[alt_read_with_indel+odds_eps, alt_read_without_indel+odds_eps ],[ref_read_with_indel+odds_eps, ref_read_without_indel+odds_eps]])
    
        # self.read_depth=dp
        # alt_dp= int(read_info_dict[self.alt]["dp"])
        # self.alt_allele_count=alt_dp
        # if dp!=0:
        #     self.vaf=int(alt_dp)/int(dp)
        # else:
            # self.vaf=0
        if dp!=0:
            self.indel_proportion_for_site=sum(read_info_dict["del"]["is_indel"])/dp
        else:
            self.indel_proportion_for_site=0

        # self.ref_allele_count=sum([float(read_info_dict[allele]["dp"]) for allele in self.ref.split(",") ])
        # self.alt_allele_count=sum([float(read_info_dict[allele]["dp"]) for allele in self.alt.split(",") ])
        # other_count_dict={base:int(read_info_dict[base]["dp"]) for base in "ATCG" if base not in self.ref+self.alt}
        # self.alt2_allele_count=max(list(other_count_dict.values()),default=0)

        # self.read_depth_reverse=sum([int(read_info_dict[allele]["reverse_dp"]) for allele in "ATCG"])
        # self.read_depth_forward=sum([int(read_info_dict[allele]["forward_dp"]) for allele in "ATCG"])
        # self.alt_allele_count_reverse=int(read_info_dict[self.alt]["reverse_dp"])
        # self.alt_allele_count_forward=int(read_info_dict[self.alt]["forward_dp"])

        # if self.read_depth_reverse>=self.read_depth_forward:
        #     self.major_read_strand="-"
        # elif self.read_depth_reverse<self.read_depth_forward:
        #     self.major_read_strand="+"

        ref_baseq, alt_baseq=get_list("baseq")
        baseq_s,self.baseq_p,self.baseq_rbc=wilcoxon_with_rbc(ref_baseq,alt_baseq,alternative="greater")

        ref_baseq1b,alt_baseq1b=get_list("baseq1b")
        # self.ref_baseq1b_s,self.ref_baseq1b_p,self.ref_baseq1b_rbc=wilcoxon_with_rbc(ref_baseq1b,ref_baseq,alternative="greater")
        alt_baseq1b_s,self.alt_baseq1b_p,self.alt_baseq1b_rbc=wilcoxon_with_rbc(alt_baseq1b,alt_baseq,alternative="greater")

        # self.ref_baseq1b_s, self.ref_baseq1b_p=do_wilicox_sum_test(ref_baseq,ref_baseq1b,method="greater",type="list")
        # self.alt_baseq1b_s, self.alt_baseq1b_p=do_wilicox_sum_test(alt_baseq,alt_baseq1b,method="greater",type="list")
        # print(ref_baseq,ref_baseq1b)
        # if self.read_depth!=0:
        #     self.alt2_proportion=int(max([read_info_dict[allele]["dp"] for allele in "ATCG" if allele not in self.ref+self.alt]))/self.read_depth
        candidates = [allele for allele in "ATCG" if allele not in [self.ref, self.alt]]
        if candidates:
            alt2 = max(candidates, key=lambda x: read_info_dict[x]["dp"])
            self.alt2_proportion_per_UMI = refine_mean(read_info_dict[alt2]["base_proportion_per_UMI"])
        else:
            self.alt2 = None
            self.alt2_proportion_per_UMI = 0.0
        # ref_rightpos,alt_rightpos=get_list("rightpos_p")
        # _,self.rightpos_p=do_wilicox_sum_test(ref_rightpos,alt_rightpos,type="list")

        ref_leftpos,alt_leftpos=get_list("leftpos_p")
        # self.leftpos_s,self.leftpos_p=do_wilicox_sum_test(ref_leftpos,alt_leftpos,type="list")
        leftpos_s,self.leftpos_p,self.leftpos_rbc=wilcoxon_with_rbc(ref_leftpos,alt_leftpos)

        ref_seqpos,alt_seqpos=get_list("seqpos")
        # self.seqpos_s,self.seqpos_p=do_wilicox_sum_test(ref_seqpos,alt_seqpos,type="list")
        seqpos_s,self.seqpos_p,self.seqpos_rbc=wilcoxon_with_rbc(ref_seqpos,alt_seqpos)

        ref_querypos,alt_querypos=get_list("querypos")
        # self.querypos_s,self.querypos_p=do_wilicox_sum_test(ref_querypos,alt_querypos,type="list")
        querypos_s,self.querypos_p,self.querypos_rbc=wilcoxon_with_rbc(ref_querypos,alt_querypos)

        # self.ref_querypos_list = ",".join([str(i) for i in ref_querypos]); self.ref_querypos_num=len(set(ref_querypos))
        # self.alt_querypos_list = ",".join([str(i) for i in alt_querypos]); 
        self.alt_querypos_num=len(set(alt_querypos))
     
        # left_ref_edist,left_alt_edist=get_list("left_read_edist")
        # print("left_ref_edist\n",set(left_ref_edist))
        # print("left_alt_edist\n",set(left_alt_edist))
        # self.left_ref_querypos_list_remove_clip = ",".join([str(i) for i in left_ref_edist])
        # self.left_ref_querypos_num_remove_clip=len(set(left_ref_edist))
        # self.left_alt_querypos_list_remove_clip = ",".join([str(i) for i in left_alt_edist])
        # self.left_alt_querypos_num_remove_clip=len(set(left_alt_edist))

        # right_ref_edist,right_alt_edist=get_list("right_read_edist")
        # print("right_ref_edist\n",set(right_ref_edist))
        # print("right_alt_edist\n",set(right_alt_edist))
        # self.right_ref_querypos_list_remove_clip = ",".join([str(i) for i in right_ref_edist])
        # self.right_ref_querypos_num_remove_clip=len(set(right_ref_edist))
        # self.right_alt_querypos_list_remove_clip = ",".join([str(i) for i in right_alt_edist])
        # self.right_alt_querypos_num_remove_clip=len(set(right_alt_edist))

        # self.queryposNum_odds, self.queryposNum_p = scipy.stats.fisher_exact([[refine_alt_mappers_multi, refine_alt_mappers_uniq ],[refine_ref_mappers_multi, refine_ref_mappers_uniq]])

        # distance_list=[]
        # distance_list_remove_clip=[]
        # for allele in "ATCG":
        #     distance_list+=read_info_dict[allele]["distance_to_end"]
        #     distance_list_remove_clip+=read_info_dict[allele]["distance_to_end_remove_clip"]
        # self.mean_distance_to_end=refine_mean(distance_list)
        # self.median_distance_to_end=refine_median(distance_list)
        # self.mean_distance_to_end_remove_clip=refine_mean(distance_list_remove_clip)
        # self.median_distance_to_end_remove_clip=refine_median(distance_list_remove_clip)

        # ref_distance_to_end,alt_distance_to_end=get_list("distance_to_end")
        # self.distance_to_end_s,self.distance_to_end_p=do_wilicox_sum_test(ref_distance_to_end,alt_distance_to_end,type="list")
        # self.distance_to_end_s,self.distance_to_end_p,self.distance_to_end_rbc=wilcoxon_with_rbc(ref_distance_to_end,alt_distance_to_end)

        # ref_distance_to_end_by_UMI,alt_distance_to_end_by_UMI=get_list("UMI_end")
        # self.UMI_end_s,self.UMI_end_p,self.UMI_end_rbc=wilcoxon_with_rbc(ref_distance_to_end_by_UMI,alt_distance_to_end_by_UMI)
        # self.UMI_end_s,self.UMI_end_p=do_wilicox_sum_test(ref_distance_to_end_by_UMI,alt_distance_to_end_by_UMI,type="list")

        # self.ref_UMI_end_mean=refine_mean(ref_distance_to_end_by_UMI)
        # self.ref_UMI_end_median=refine_median(ref_distance_to_end_by_UMI)
        # self.alt_UMI_end_mean=refine_mean(alt_distance_to_end_by_UMI)
        # self.alt_UMI_end_median=refine_median(alt_distance_to_end_by_UMI)

        # ref_distance_to_end_remove_clip_by_UMI,alt_distance_to_end_remove_clip_by_UMI=get_list("UMI_end_remove_clip")
        # self.ref_UMI_end_mean_remove_clip = refine_mean(ref_distance_to_end_remove_clip_by_UMI)  # 修正变量名
        # self.ref_UMI_end_median_remove_clip = refine_median(ref_distance_to_end_remove_clip_by_UMI)  # 修正变量名
        # self.alt_UMI_end_mean_remove_clip = refine_mean(alt_distance_to_end_remove_clip_by_UMI)  # 修正变量名
        # self.alt_UMI_end_median_remove_clip = refine_median(alt_distance_to_end_remove_clip_by_UMI)  # 修正变量名
        # self.UMI_end_remove_clip_s,self.UMI_end_remove_clip_p,self.UMI_end_remove_clip_rbc=wilcoxon_with_rbc(ref_distance_to_end_remove_clip_by_UMI,alt_distance_to_end_remove_clip_by_UMI)

        # ref_distance_to_end_remove_clip,alt_distance_to_end_remove_clip=get_list("distance_to_end_remove_clip")
        # print(ref_distance_to_end_remove_clip,alt_distance_to_end_remove_clip)
        # self.distance_to_end_remove_clip_s,self.distance_to_end_remove_clip_p=do_wilicox_sum_test(ref_distance_to_end_remove_clip,alt_distance_to_end_remove_clip,type="list")
        # print(self.distance_to_end_remove_clip_s,self.distance_to_end_remove_clip_p)
        ## NOTE: save softclip_length and plot them
        
        # ref_distance_to_end_per_UMI, alt_distance_to_end_per_UMI = get_list("per_UMI_end")
        # self.per_ref_UMI_end_mean = refine_mean(ref_distance_to_end_per_UMI)
        # self.per_ref_UMI_end_median = refine_median(ref_distance_to_end_per_UMI)
        # self.per_alt_UMI_end_mean = refine_mean(alt_distance_to_end_per_UMI)
        # self.per_alt_UMI_end_median = refine_median(alt_distance_to_end_per_UMI)
        # self.per_UMI_end_s, self.per_UMI_end_p, self.per_UMI_end_rbc = wilcoxon_with_rbc(ref_distance_to_end_per_UMI, alt_distance_to_end_per_UMI)

        # per_UMI_end_remove_clip 相关计算
        ref_distance_to_end_remove_clip_per_UMI, alt_distance_to_end_remove_clip_per_UMI = get_list("per_UMI_end_remove_clip")
        # self.per_ref_UMI_end_remove_clip_mean = refine_mean(ref_distance_to_end_remove_clip_per_UMI)  # 修正输入变量
        # self.per_ref_UMI_end_remove_clip_median = refine_median(ref_distance_to_end_remove_clip_per_UMI)  # 修正输入变量
        self.per_alt_UMI_end_remove_clip_mean = refine_mean(alt_distance_to_end_remove_clip_per_UMI)  # 修正输入变量
        self.per_alt_UMI_end_remove_clip_median = refine_median(alt_distance_to_end_remove_clip_per_UMI)  # 修正输入变量
        per_UMI_end_remove_clip_s, self.per_UMI_end_remove_clip_p, self.per_UMI_end_remove_clip_rbc = wilcoxon_with_rbc(ref_distance_to_end_remove_clip_per_UMI, alt_distance_to_end_remove_clip_per_UMI)

        # 获取 per_UMI_end_value 数据并计算统计量
        # ref_distance_to_end_value_per_UMI, alt_distance_to_end_value_per_UMI = get_list("per_UMI_end_value")
        # self.per_ref_UMI_end_value_mean = refine_mean(ref_distance_to_end_value_per_UMI)  # 修正变量名
        # self.per_ref_UMI_end_value_median = refine_median(ref_distance_to_end_value_per_UMI)  # 修正变量名
        # self.per_alt_UMI_end_value_mean = refine_mean(alt_distance_to_end_value_per_UMI)  # 修正变量名
        # self.per_alt_UMI_end_value_median = refine_median(alt_distance_to_end_value_per_UMI)  # 修正变量名

        # 获取 per_UMI_end_remove_clip_value 数据并计算统计量
        # ref_distance_to_end_remove_clip_value_per_UMI, alt_distance_to_end_remove_clip_value_per_UMI = get_list("per_UMI_end_remove_clip_value")  # 统一变量名
        # self.per_ref_UMI_end_remove_clip_value_mean = refine_mean(ref_distance_to_end_remove_clip_value_per_UMI)  # 使用正确变量
        # self.per_ref_UMI_end_remove_clip_value_median = refine_median(ref_distance_to_end_remove_clip_value_per_UMI)  # 使用正确变量
        # self.per_alt_UMI_end_remove_clip_value_mean = refine_mean(alt_distance_to_end_remove_clip_value_per_UMI)  # 使用正确变量
        # self.per_alt_UMI_end_remove_clip_value_median = refine_median(alt_distance_to_end_remove_clip_value_per_UMI)  # 使用正确变量

        
        # ref_hardclip_length,alt_hardclip_length=get_list("hardclip_length")
        # self.hardclip_length_s,self.hardclip_length_p=do_wilicox_sum_test(ref_hardclip_length,alt_hardclip_length,type="list")
        # try:
        #     self.ref_hardclip_prop=len([x for x in ref_hardclip_length if x > 10])/len(ref_hardclip_length)
        #     self.alt_hardclip_prop=len([x for x in alt_hardclip_length if x > 10])/len(alt_hardclip_length)
        #     # self.ref_hardclip_length,self.alt_hardclip_length=ref_hardclip_length,alt_hardclip_length
        #     ref_hard_count=len([x for x in ref_hardclip_length if x > 10]); ref_no_hard_count=len(ref_hardclip_length)-ref_hard_count
        #     alt_hard_count=len([x for x in alt_hardclip_length if x > 10]); alt_no_hard_count=len(alt_hardclip_length)-alt_hard_count
        #     self.hardclip_prop_odds, self.hardclip_prop_p=scipy.stats.fisher_exact([[alt_hard_count, alt_no_hard_count ],[ref_hard_count, ref_no_hard_count]])
    
        # except:
        #     print(self.identifier, "dose not have hard_clip_info")

        ref_softclip_length,alt_softclip_length=get_list("softclip_length")
        # self.softclip_length_s,self.softclip_length_p=do_wilicox_sum_test(ref_softclip_length,alt_softclip_length,type="list")
        softclip_length_s,self.softclip_length_p,self.softclip_length_rbc=wilcoxon_with_rbc(ref_softclip_length,alt_softclip_length)
        

        try:    
            self.ref_softclip_prop=len([x for x in ref_softclip_length if x > 10])/len(ref_softclip_length)
            self.alt_softclip_prop=len([x for x in alt_softclip_length if x > 10])/len(alt_softclip_length)
            # self.ref_softclip_length,self.alt_softclip_length=ref_softclip_length,alt_softclip_length
            ref_soft_count=len([x for x in ref_softclip_length if x > 10]); ref_no_soft_count=len(ref_softclip_length)-ref_soft_count
            alt_soft_count=len([x for x in alt_softclip_length if x > 10]); alt_no_soft_count=len(alt_softclip_length)-alt_soft_count
            softclip_prop_odds, self.softclip_prop_p=scipy.stats.fisher_exact([[alt_soft_count, alt_no_soft_count ],[ref_soft_count, ref_no_soft_count]])

            all_softclip=ref_softclip_length+alt_softclip_length
            self.softclip_prop=len([x for x in all_softclip if x > 10])/len(all_softclip) if len(all_softclip)!=0 else "NA"

        except:
            print(self.identifier, "dose not have soft_clip_info")
        # def simple_to_get_list_float(read_info_dict,ref_allele,alt_allele,var):
        #     # print(ref_allele,alt_allele, ref_allele.split(","))
        #     return_ref_list=[k for allele in ref_allele.split(",") for k in read_info_dict[allele][var] if k !=0.0]
        #     return_alt_list=[k for allele in alt_allele.split(",") for k in read_info_dict[allele][var] if k !=0.0]
        #     return return_ref_list,return_alt_list

        # get_list_float=partial(simple_to_get_list_float,read_info_dict,self.ref,self.alt)
        #UMI_consistence_prop
        # _,alt_UMI_consistence=get_list_float("UMI_consistence_prop")
        ##### UMI consistency
        alt_UMI_consistence=0; total_UMI_contain_alt=0; alt_consistence_hard=0; alt_consistence_soft=0
        for allele in self.alt.split(","):
            for k in read_info_dict[allele]["UMI_consistence_prop"]:
                if k !=0.0:
                    alt_UMI_consistence+=k
                    total_UMI_contain_alt+=1
                    if k==1.0:
                        alt_consistence_hard+=1
                        alt_consistence_soft+=1
                    elif k>=0.75:
                        alt_consistence_soft+=1
        try:
            alt_UMI_consistence_prop=alt_UMI_consistence/total_UMI_contain_alt
            alt_consistent_UMI_prop_strict=alt_consistence_hard/total_UMI_contain_alt
            alt_consistent_UMI_prop_relaxed=alt_consistence_soft/total_UMI_contain_alt
            self.alt_UMI_consistence_prop=alt_UMI_consistence_prop
            self.alt_consistent_UMI_prop_strict=alt_consistent_UMI_prop_strict
            self.alt_consistent_UMI_prop_relaxed=alt_consistent_UMI_prop_relaxed
        except:
            # print(self.identifier,"is wrong in consistence.")
            pass

        ##### UMI consistency remove single read
        alt_UMI_consistence_remove_single_read=0; total_UMI_contain_alt_remove_single_read=0
        alt_consistence_hard_remove_single_read=0; alt_consistence_soft_remove_single_read=0
        for allele in self.alt.split(","):
            for k in read_info_dict[allele]["UMI_consistence_prop_remove_single_read"]:
                if k !=0.0:
                    alt_UMI_consistence_remove_single_read+=k
                    total_UMI_contain_alt_remove_single_read+=1
                    if k==1.0:
                        alt_consistence_hard_remove_single_read+=1
                        alt_consistence_soft_remove_single_read+=1
                    elif k>=0.75:
                        alt_consistence_soft_remove_single_read+=1
        try:
            alt_UMI_avg_consistence_remove_single_read=alt_UMI_consistence_remove_single_read/total_UMI_contain_alt_remove_single_read
            alt_consistent_UMI_prop_relaxed_remove_single_read=alt_consistence_soft_remove_single_read/total_UMI_contain_alt_remove_single_read
            alt_consistent_UMI_prop_strict_remove_single_read=alt_consistence_hard_remove_single_read/total_UMI_contain_alt_remove_single_read
            self.alt_UMI_avg_consistence_remove_single_read=alt_UMI_avg_consistence_remove_single_read
            self.alt_consistent_UMI_prop_strict_remove_single_read=alt_consistent_UMI_prop_strict_remove_single_read
            self.alt_consistent_UMI_prop_relaxed_remove_single_read=alt_consistent_UMI_prop_relaxed_remove_single_read
        except:
            # print(self.identifier,"is wrong in consistence remove single read.")
            pass

        try:
            ref_read_number_perUMI,alt_read_number_perUMI=get_list("read_number_per_UMI")
            read_number_s,self.read_number_p,self.read_number_rbc=wilcoxon_with_rbc(ref_read_number_perUMI,alt_read_number_perUMI)
            # self.read_number_s,self.read_number_p=do_wilicox_sum_test(ref_read_number_perUMI,alt_read_number_perUMI,type="list")
            # self.ref_read_number_perUMI_median=refine_median(ref_read_number_perUMI)
            # self.alt_read_number_perUMI_median=refine_median(alt_read_number_perUMI)
            # self.ref_read_number_perUMI_max=max(ref_read_number_perUMI)
            self.alt_read_number_perUMI_max=max(alt_read_number_perUMI)
        except:
            print(self.identifier,"does not have perUMI info")
        
        # try:
        #     ref_read_number_per_spot,alt_read_number_per_spot=get_list("total_read_number_per_spot")
        #     self.read_number_per_spot_s,self.read_number_per_spot_p = do_wilicox_sum_test(ref_read_number_per_spot, alt_read_number_per_spot,method="two-sided",type="list")
        # except:
        #     print(self.identifier,"does not have total_read_number_per_spot info")

        # try:
        #     ref_spot_dp_list, alt_spot_dp_list=get_list("total_UMI_number_per_spot")
        #     # self.UMI_number_per_spot_s,self.UMI_number_per_spot_p = do_wilicox_sum_test(ref_spot_dp_list, alt_spot_dp_list,method="two-sided",type="list")
        #     # self.ref_UMI_number_per_spot_median,self.alt_UMI_number_per_spot_median = np.median(ref_spot_dp_list),np.median(alt_spot_dp_list)
        #     self.ref_UMI_number_per_spot_max,self.alt_UMI_number_per_spot_max = max([0]+ref_spot_dp_list),max([0]+alt_spot_dp_list)
        # except:
        #     print(self.identifier,"does not have total_UMI_number_per_spot info")
        
        alt_umi_number=[float(k) for allele in self.alt.split(",") for k in read_info_dict[allele]["UMI_number_per_spot"]]
        # ref_umi_number=[float(k) for allele in self.ref.split(",") for k in read_info_dict[allele]["UMI_number_per_spot"]]

        umi_depth=[a+t+c+g for a,t,c,g in zip(read_info_dict["A"]["UMI_number_per_spot"], \
                                                read_info_dict["T"]["UMI_number_per_spot"], \
                                                read_info_dict["C"]["UMI_number_per_spot"], \
                                                read_info_dict["G"]["UMI_number_per_spot"])]
        # self.UMI_number_per_spot_median=refine_median(umi_depth)
        # self.UMI_number_per_spot_mean=refine_mean(umi_depth)
        self.alt_vs_total_dp_paired_wilcoxon_rbc=calculate_rbc_for_paired_wilcoxon(umi_depth,alt_umi_number)
        # self.downsample_consensus_UMI_count=sum([float(k) for allele in "ATCG" for k in read_info_dict[allele]["UMI_number_per_spot"]])
        # self.downsample_consensus_alt_allele_count=sum(alt_umi_number)
        
        # alt_consistence_prop=sum(alt_UMI_consistence)/sum([1 for i in read_info_dict[self.alt]["UMI_consistence_prop"] if i !=0.0])
        # total_UMI_contain_alt=[i for i in read_info_dict[self.alt]["UMI_consistence_prop"] if i !=0.0]
        # alt_consistence_hard=sum([1 for i in total_UMI_contain_alt if i ==1.0])/len(total_UMI_contain_alt)
        # alt_consistence_soft=sum([1 for i in total_UMI_contain_alt if i >=0.75])/len(total_UMI_contain_alt)       
        
        # "del_length", "del_distance", "del_num", "ins_distance", "ins_length","ins_num",

        # def simple_to_get_list_facing_list(read_info_dict,ref_allele,alt_allele,var):
        #     # print(read_info_dict[ref_allele][var])
        #     return_ref_list=[i for allele in ref_allele.split(",") for k in read_info_dict[allele][var] for i in k]
        #     return_alt_list=[i for allele in alt_allele.split(",") for k in read_info_dict[allele][var] for i in k]
        #     # return_ref_list=[str(k) for allele in ref_allele.split(",") for k in read_info_dict[allele][var] if k !=""]
        #     # return_alt_list=[str(k) for allele in alt_allele.split(",") for k in read_info_dict[allele][var] if k !=""]
        #     return return_ref_list,return_alt_list
            
        # get_list_face_list=partial(simple_to_get_list_facing_list,read_info_dict,self.ref,self.alt)
        #ref_ins_distance: [[],["no",1],[2]]
        # ref_ins_num,alt_ins_num=get_list("ins_num")
        # ref_ins_length,alt_ins_length=get_list_face_list("ins_length")
        # ref_ins_distance,alt_ins_distance=get_list_face_list("ins_distance")
        # self.ref_ins_num,self.alt_ins_num=ref_ins_num,alt_ins_num
        # self.ref_ins_length,self.alt_ins_length=ref_ins_length,alt_ins_length
        # self.ref_ins_distance,self.alt_ins_distance=ref_ins_distance,alt_ins_distance

        # ref_del_num,alt_del_num=get_list("del_num")
        # ref_del_length,alt_del_length=get_list_face_list("del_length")
        # ref_del_distance,alt_del_distance=get_list_face_list("del_distance")
        # self.ref_del_num,self.alt_del_num=ref_del_num,alt_del_num
        # self.ref_del_length,self.alt_del_length=ref_del_length,alt_del_length
        # self.ref_del_distance,self.alt_del_distance=ref_del_distance,alt_del_distance

        # def get_sum_except_0(indel_number_list):
        #     count_dict=Counter(indel_number_list)
        #     sum_value=sum([int(count_dict[key]) for key in count_dict.keys() if key not in ["0",0]])
        #     total=sum(count_dict.values())
        #     try:
        #         return sum_value/total
        #     except:
        #         return 0
        
        # def get_major_prop_except_no(ref_del_length):
        #     count_dict=Counter(ref_del_length);del count_dict["no"]
        #     if count_dict!={}:
        #         most_common_key = count_dict.most_common(1)[0][0]
        #         most_common_count=count_dict[most_common_key]
        #         indel_num_except_no=sum(count_dict.values())
        #         most_common_prop=int(most_common_count)/int(indel_num_except_no)
        #         return most_common_key,most_common_prop
        #     else:
        #         return "",0
        
        # sum_ref_ins_num=get_sum_except_0(self.ref_ins_num)
        # sum_ref_del_num=get_sum_except_0(self.ref_del_num)
        # sum_alt_ins_num=get_sum_except_0(self.alt_ins_num)
        # sum_alt_del_num=get_sum_except_0(self.alt_del_num)

        # self.ref_ins_prop=sum_ref_ins_num
        # self.ref_del_prop=sum_ref_del_num
        # _,self.ref_ins_major_prop=get_major_prop_except_no(self.ref_ins_distance)
        # _,self.ref_del_major_prop=get_major_prop_except_no(self.ref_del_distance)
        # self.alt_ins_prop=sum_alt_ins_num
        # self.alt_del_prop=sum_alt_del_num
        # _,self.alt_ins_major_prop=get_major_prop_except_no(self.alt_ins_distance)
        # _,self.alt_del_major_prop=get_major_prop_except_no(self.alt_del_distance)

        # features from bcftools: VDB   
        # self.VDB=calc_vdb(read_info_dict[self.alt]["edist"],self.readLen,self.readLen)
        # ref_pos,alt_pos=get_list("epos")
        # self.RPBZ,_=do_wilicox_sum_test(ref_pos,alt_pos,method="two-sided",type="list")
        # self.MQBZ,_=calc_mwu_biasZ(bca->ref_mq,  bca->alt_mq, bca->nqual,1,1) #mapQ have been done before
        # BQBZ= calc_mwu_biasZ(bca->ref_bq,  bca->alt_bq, bca->nqual,0,1) #baseQ have been done before
        # MQSBZ= calc_mwu_biasZ(bca->fwd_mqs, bca->rev_mqs, bca->nqual,0,1) # this have no sense, cause all RNA-seq reads are same direction
        # ref_dp_list=[read_info_dict[allele]["dp"] for allele in self.ref.split(",")]
        # alt_dp_list=[read_info_dict[allele]["dp"] for allele in self.alt.split(",")]
        # self.SGB=calc_SegBias_for_one_sample(sum(ref_dp_list),sum(alt_dp_list))

        # self.downsample_alt_spot_num=read_info_dict[self.alt]["GenoSpotNum"]
        # if self.downsample_spot_num=="no":
        #     self.downsample_spot_num=int(read_info_dict[self.alt]["GenoSpotNum"])+int(read_info_dict[self.ref[0]]["GenoSpotNum"])
        # # print(read_info_dict[self.alt]["vaf_spot"])
        # self.downsample_mut_spots_vaf_mean=refine_mean(read_info_dict[self.alt]["vaf_spot"])
        # self.downsample_mut_spots_vaf_median=refine_median(read_info_dict[self.alt]["vaf_spot"])
        # vaf_values = np.array(read_info_dict[self.alt]["vaf_spot"])  # 确保是数组
        # self.downsample_only_alt_mutant_spot_prop = (vaf_values == 1).mean()

        # p_list=[self.baseq_p,self.ref_baseq1b_p,self.alt_baseq1b_p, \
        #         self.querypos_p,self.leftpos_p,self.seqpos_p, \
        #         self.distance_to_end_p,self.UMI_end_p,self.UMI_end_remove_clip_p, \
        #         self.mismatches_p,self.mapq_p,self.read_number_p,self.softclip_length_p]
        
        p_list=[self.baseq_p,self.ref_baseq1b_p,self.alt_baseq1b_p, \
                self.querypos_p,self.leftpos_p,self.seqpos_p, \
                self.distance_to_end_p, self.UMI_end_p,self.UMI_end_remove_clip_p, \
                self.per_UMI_end_p,self.per_UMI_end_remove_clip_p, \
                self.mismatches_p,self.mapq_p,self.read_number_p,self.softclip_length_p]
        _, p_adj_list, _, _ = smm.multipletests(p_list, method='fdr_bh')
        # self.baseq_p_adj,self.ref_baseq1b_p_adj,self.alt_baseq1b_p_adj, \
        # self.querypos_p_adj,self.leftpos_p_adj,self.seqpos_p_adj, \
        # self.distance_to_end_p_adj,self.UMI_end_p_adj,self.UMI_end_remove_clip_p_adj, \
        # self.mismatches_p_adj,self.mapq_p_adj,self.read_number_p_adj,self.softclip_length_p_adj=p_adj_list
        
        self.baseq_p_adj,self.ref_baseq1b_p_adj,self.alt_baseq1b_p_adj, \
        self.querypos_p_adj,self.leftpos_p_adj,self.seqpos_p_adj, \
        self.distance_to_end_p_adj,self.UMI_end_p_adj,self.UMI_end_remove_clip_p_adj, \
        self.per_UMI_end_p_adj,self.per_UMI_end_remove_clip_p_adj, \
        self.mismatches_p_adj,self.mapq_p_adj,self.read_number_p_adj,self.softclip_length_p_adj=p_adj_list


    def _expand_features(self):
        """
        after combine the features from different files, some features can be expanded
        such as: 
        exonDesOrder: sort descending of exon location
        
        """

        if self.consensus_read_count!="no" and self.ref!="no" and self.alt!="no":
            other_count_dict={base:int(count) for base,count in zip("ATCG",self.consensus_read_count.split(",")) if base not in self.ref+self.alt}
            # print(self.consensus_read_count,"===", other_count_dict)
            alt2_count_consensus=max(list(other_count_dict.values()),default=0)
            consensus_dp=sum([int(i) for i in self.consensus_read_count.split(",")])
            self.consensus_alt2_proportion=alt2_count_consensus/consensus_dp
            self.consensus_UMI_count=consensus_dp
            self.consensus_alt_allele_count=self.consensus_read_count.split(",")["ATCG".index(self.alt)]
            self.consensus_vaf=int(self.consensus_alt_allele_count)/self.consensus_UMI_count




####################### the following is outside functions #############################

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

        var_result_dict[var],_=process_reads_for_variant(sampled_reads,var,run_type,bins,cell_dict,readLen)
    
    return var_result_dict


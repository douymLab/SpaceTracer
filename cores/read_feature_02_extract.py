from functools import partial
from typing import  Dict, Optional
from dataclasses import dataclass, field, fields
from typing import Any, Dict
import statsmodels.stats.multitest as smm
import scipy
import numpy as np

from SpaceTracer.utils.get_read_level_feature import wilcoxon_with_rbc
from SpaceTracer.utils.utils import refine_mean, refine_median


@dataclass
class readLevelFeatures:
    identifier: tuple
    
    chrom: str = field(init=False)
    pos: int = field(init=False)
    ref: str = field(init=False)
    alt: str = field(init=False)
    # mut_origin: str = "NA"

    baseq_p: float = None
    baseq_p_adj: float = None
    baseq_rbc: float = None
    ref_baseq_mean: float = None
    alt_baseq_mean: float = None
    baseq_mean: float = None

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

    ref_mapq255_prop: float = None
    alt_mapq255_prop: float = None
    mapq255_prop: float = None

    alt_UMI_avg_consistence: float = None
    alt_UMI_avg_consistence_remove_single_read: float = None
    alt_consistent_UMI_prop_strict: float = None
    alt_consistent_UMI_prop_strict_remove_single_read: float = None
    alt_consistent_UMI_prop_relaxed: float = None
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
    strand_bias_p_adj: float = None
    ref_softclip_prop: float = None
    alt_softclip_prop: float = None
    
    softclip_prop: float = None
    softclip_prop_p: float = None
    softclip_prop_p_adj: float = None

    softclip_length_p: float = None
    softclip_length_p_adj: float = None
    softclip_length_rbc: float = None

    # alt_vs_total_dp_paired_wilcoxon_rbc: float = None
   
    # not used feature
    EXCLUDE_FIELDS = {'identifier'} #'chrom', 'pos', 'ref', 'alt', 
    
    def __post_init__(self):
        mut_chrom, mut_pos, mut_ref, mut_alt = self.identifier

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
            if not value:
                value=None
            # if f.name.endswith('_p') or f.name.endswith('_p_adj'):
            #     value = handle_p_value_log10(value)
            
            result[f.name] = value
        
        return result
    
    def test_values(self) -> Dict[str, Any]:
        return {
            'identifier': self.identifier
        }

    @classmethod
    def from_read_info(cls, identifier: tuple, read_info_dict: Dict[str, Any]):
        if not read_info_dict:
            return None
        
        if read_info_dict.get('dp', 0) == 0:
            return None
        
        feature = cls(identifier=identifier)
        feature._add_read_info(read_info_dict)
        return feature

    @classmethod
    def from_read_info_to_dict(cls, identifier: tuple, read_info_dict: Dict[str, Any]) -> Optional[Dict[str, Any]]:
        feature = cls.from_read_info(identifier, read_info_dict)
        return feature.to_dict() if feature is not None else None

    def _run(self, read_info_dict: Dict[str, Any]):
        if not read_info_dict:
            return None
        
        if read_info_dict.get('dp', 0) == 0:
            return None
        
        self._add_read_info(read_info_dict)
        return self

    def _add_read_info(self,read_info_dict):
        dp=read_info_dict['dp']
            
        def simple_to_get_list(read_info_dict,ref_allele,alt_allele,var):
            # print(ref_allele,alt_allele, ref_allele.split(","))
            return_ref_list=[float(k) for allele in ref_allele.split(",") for k in read_info_dict[allele][var] if k !=""]
            return_alt_list=[float(k) for allele in alt_allele.split(",") for k in read_info_dict[allele][var] if k !=""]
            return return_ref_list,return_alt_list
            
        get_list=partial(simple_to_get_list,read_info_dict,self.ref,self.alt)

        odds_eps=0

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
        mismatches_s,self.mismatches_p,self.mismatches_rbc=wilcoxon_with_rbc(ref_mismatches,refine_alt_mismatches)
        
        ref_is_reverse, alt_is_reverse=get_list("is_reverse")
        refine_ref_is_reverse = len([i for i in ref_is_reverse if i == 1]); refine_ref_is_forward = len(ref_is_reverse)-refine_ref_is_reverse
        refine_alt_is_reverse = len([i for i in alt_is_reverse if i == 1]); refine_alt_is_forward = len(alt_is_reverse)-refine_alt_is_reverse
        strand_bias_odds, self.strand_bias_p = scipy.stats.fisher_exact([[refine_alt_is_reverse, refine_alt_is_forward ],[refine_ref_is_reverse, refine_ref_is_forward]])

        ref_mapq,alt_mapq=get_list("map_q")

        def calc_255_prop(seq):
            if not seq:
                return 0.0
            try:
                return sum(1 for i in seq if int(i) == 255) / len(seq)
            except (ValueError, TypeError):
                return 0.0

        self.ref_mapq255_prop = calc_255_prop(ref_mapq)
        self.alt_mapq255_prop = calc_255_prop(alt_mapq)

        total_len = len(ref_mapq) + len(alt_mapq)
        if total_len > 0:
            count = sum(1 for i in ref_mapq if int(i) == 255) + \
                    sum(1 for i in alt_mapq if int(i) == 255)
            self.mapq255_prop = count / total_len
        else:
            self.mapq255_prop = 0.0

        mapq_s,self.mapq_p,self.mapq_rbc=wilcoxon_with_rbc(ref_mapq,alt_mapq,alternative="greater")

        ref_ind_num,alt_ind_num=get_list("ind_num")
        ref_read_with_indel = len([i for i in ref_ind_num if i > 0]); ref_read_without_indel = len(ref_ind_num)-ref_read_with_indel
        alt_read_with_indel = len([i for i in alt_ind_num if i > 0]); alt_read_without_indel = len(alt_ind_num)-alt_read_with_indel
        reads_with_indel_odds, self.reads_with_indel_p = scipy.stats.fisher_exact([[alt_read_with_indel+odds_eps, alt_read_without_indel+odds_eps ],[ref_read_with_indel+odds_eps, ref_read_without_indel+odds_eps]])
    
        if dp!=0:
            self.indel_proportion_for_site=sum(read_info_dict["del"]["is_indel"])/dp
        else:
            self.indel_proportion_for_site=0

        ref_baseq, alt_baseq=get_list("baseq")
        baseq_s,self.baseq_p,self.baseq_rbc=wilcoxon_with_rbc(ref_baseq,alt_baseq,alternative="greater")
        self.ref_baseq_mean=refine_mean(ref_baseq)
        self.alt_baseq_mean=refine_mean(alt_baseq)
        self.baseq_mean=refine_mean(ref_baseq+alt_baseq)

        ref_baseq1b,alt_baseq1b=get_list("baseq1b")
        # self.ref_baseq1b_s,self.ref_baseq1b_p,self.ref_baseq1b_rbc=wilcoxon_with_rbc(ref_baseq1b,ref_baseq,alternative="greater")
        alt_baseq1b_s,self.alt_baseq1b_p,self.alt_baseq1b_rbc=wilcoxon_with_rbc(alt_baseq1b,alt_baseq,alternative="greater")

        candidates = [allele for allele in "ATCG" if allele not in [self.ref, self.alt]]
        if candidates:
            alt2 = max(candidates, key=lambda x: read_info_dict[x]["dp"])
            self.alt2_proportion_per_UMI = refine_mean(read_info_dict[alt2]["base_proportion_per_UMI"])
        else:
            self.alt2 = None
            self.alt2_proportion_per_UMI = 0.0

        ref_leftpos,alt_leftpos=get_list("leftpos_p")
        leftpos_s,self.leftpos_p,self.leftpos_rbc=wilcoxon_with_rbc(ref_leftpos,alt_leftpos)

        ref_seqpos,alt_seqpos=get_list("seqpos")
        seqpos_s,self.seqpos_p,self.seqpos_rbc=wilcoxon_with_rbc(ref_seqpos,alt_seqpos)

        ref_querypos,alt_querypos=get_list("querypos")
        querypos_s,self.querypos_p,self.querypos_rbc=wilcoxon_with_rbc(ref_querypos,alt_querypos)
        self.alt_querypos_num=len(set(alt_querypos))
       
        ref_distance_to_end_remove_clip_per_UMI, alt_distance_to_end_remove_clip_per_UMI = get_list("per_UMI_end_remove_clip")
        self.per_alt_UMI_end_remove_clip_mean = refine_mean(alt_distance_to_end_remove_clip_per_UMI)  # 修正输入变量
        self.per_alt_UMI_end_remove_clip_median = refine_median(alt_distance_to_end_remove_clip_per_UMI)  # 修正输入变量
        # print(ref_distance_to_end_remove_clip_per_UMI)
        # print(alt_distance_to_end_remove_clip_per_UMI)
        per_UMI_end_remove_clip_s, self.per_UMI_end_remove_clip_p, self.per_UMI_end_remove_clip_rbc = wilcoxon_with_rbc(ref_distance_to_end_remove_clip_per_UMI, alt_distance_to_end_remove_clip_per_UMI)

        ref_softclip_length,alt_softclip_length=get_list("softclip_length")
        softclip_length_s,self.softclip_length_p,self.softclip_length_rbc=wilcoxon_with_rbc(ref_softclip_length,alt_softclip_length)
        
        try:    
            self.ref_softclip_prop=len([x for x in ref_softclip_length if x > 10])/len(ref_softclip_length)
            self.alt_softclip_prop=len([x for x in alt_softclip_length if x > 10])/len(alt_softclip_length)
            ref_soft_count=len([x for x in ref_softclip_length if x > 10]); ref_no_soft_count=len(ref_softclip_length)-ref_soft_count
            alt_soft_count=len([x for x in alt_softclip_length if x > 10]); alt_no_soft_count=len(alt_softclip_length)-alt_soft_count
            softclip_prop_odds, self.softclip_prop_p=scipy.stats.fisher_exact([[alt_soft_count, alt_no_soft_count ],[ref_soft_count, ref_no_soft_count]])

            all_softclip=ref_softclip_length+alt_softclip_length
            self.softclip_prop=len([x for x in all_softclip if x > 10])/len(all_softclip) if len(all_softclip)!=0 else "NA"

        except:
            print(self.identifier, "dose not have soft_clip_info")
    
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
            self.alt_UMI_avg_consistence=alt_UMI_consistence_prop
            self.alt_consistent_UMI_prop_strict=alt_consistent_UMI_prop_strict
            self.alt_consistent_UMI_prop_relaxed=alt_consistent_UMI_prop_relaxed
        except:
            pass

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
            pass

        try:
            ref_read_number_perUMI,alt_read_number_perUMI=get_list("read_number_per_UMI")
            read_number_s,self.read_number_p,self.read_number_rbc=wilcoxon_with_rbc(ref_read_number_perUMI,alt_read_number_perUMI)
            self.alt_read_number_perUMI_max=max(alt_read_number_perUMI)
        except:
            print(self.identifier,"does not have perUMI info")
        
        # alt_umi_number=[float(k) for allele in self.alt.split(",") for k in read_info_dict[allele]["UMI_number_per_spot"]]
        # ref_umi_number=[float(k) for allele in self.ref.split(",") for k in read_info_dict[allele]["UMI_number_per_spot"]]

        # umi_depth=[a+t+c+g for a,t,c,g in zip(read_info_dict["A"]["UMI_number_per_spot"], \
        #                                         read_info_dict["T"]["UMI_number_per_spot"], \
        #                                         read_info_dict["C"]["UMI_number_per_spot"], \
        #                                         read_info_dict["G"]["UMI_number_per_spot"])]
        # self.alt_vs_total_dp_paired_wilcoxon_rbc=calculate_rbc_for_paired_wilcoxon(umi_depth,alt_umi_number)
        
        # p_list=[self.baseq_p,self.alt_baseq1b_p,self.querypos_p, \
        #         self.leftpos_p, self.seqpos_p, self.per_UMI_end_remove_clip_p, self.mismatches_p, \
        #         self.mapq_p, self.read_number_p, self.softclip_length_p, self.softclip_prop_p, \
        #         self.reads_with_indel_p, self.multi_mapper_p, self.strand_bias_p]

        # name_list=['baseq_p','alt_baseq1b_p','querypos_p', \
        #         'leftpos_p', 'seqpos_p', 'per_UMI_end_remove_clip_p', 'mismatches_p', \
        #         'mapq_p', 'read_number_p', 'softclip_length_p', 'softclip_prop_p', \
        #         'reads_with_indel_p', 'multi_mapper_p', 'strand_bias_p']

        # # for name,p in zip(name_list,p_list):
        # #     print("#",name,":", p)

        # _, p_adj_list, _, _ = smm.multipletests(p_list, method='fdr_bh')
        
        # self.baseq_p_adj,self.alt_baseq1b_p_adj,self.querypos_p_adj, \
        # self.leftpos_p_adj, self.seqpos_p_adj, self.per_UMI_end_remove_clip_p_adj, self.mismatches_p_adj, \
        # self.mapq_p_adj, self.read_number_p_adj, self.softclip_length_p_adj, self.softclip_prop_p_adj, \
        # self.reads_with_indel_p_adj, self.multi_mapper_p_adj, self.strand_bias_p_adj \
        #     = p_adj_list

        p_list = [self.baseq_p, self.alt_baseq1b_p, self.querypos_p,
                self.leftpos_p, self.seqpos_p, self.per_UMI_end_remove_clip_p, self.mismatches_p,
                self.mapq_p, self.read_number_p, self.softclip_length_p, self.softclip_prop_p,
                self.reads_with_indel_p, self.multi_mapper_p, self.strand_bias_p]

        valid_mask = [p is not None and not np.isnan(p) for p in p_list]
        p_valid = [p for p, valid in zip(p_list, valid_mask) if valid]

        p_adj_list = [np.nan] * len(p_list)

        if p_valid:
            _, p_adj_valid, _, _ = smm.multipletests(p_valid, method='fdr_bh')
            j = 0
            for i, valid in enumerate(valid_mask):
                if valid:
                    p_adj_list[i] = p_adj_valid[j]
                    j += 1

        self.baseq_p_adj, self.alt_baseq1b_p_adj, self.querypos_p_adj, \
        self.leftpos_p_adj, self.seqpos_p_adj, self.per_UMI_end_remove_clip_p_adj, self.mismatches_p_adj, \
        self.mapq_p_adj, self.read_number_p_adj, self.softclip_length_p_adj, self.softclip_prop_p_adj, \
        self.reads_with_indel_p_adj, self.multi_mapper_p_adj, self.strand_bias_p_adj \
            = p_adj_list
        
        # self.baseq_p_adj,self.ref_baseq1b_p_adj,self.alt_baseq1b_p_adj, \
        # self.querypos_p_adj,self.leftpos_p_adj,self.seqpos_p_adj, \
        # self.distance_to_end_p_adj,self.UMI_end_p_adj,self.UMI_end_remove_clip_p_adj, \
        # self.per_UMI_end_p_adj,self.per_UMI_end_remove_clip_p_adj, \
        # self.mismatches_p_adj,self.mapq_p_adj,self.read_number_p_adj,self.softclip_length_p_adj=p_adj_list

        # baseq_p_adj,alt_baseq1b_p_adj,querypos_p_adj,leftpos_p_adj, \
        # seqpos_p_adj,per_UMI_end_remove_clip_p_adj, mismatches_p_adj, \
        # mapq_p_adj,read_number_p_adj,reads_with_indel_p_adj, \
        # multi_mapper_p_adj,strand_bias_p_adj, softclip_prop_p_adj,\
        # softclip_length_p_adj


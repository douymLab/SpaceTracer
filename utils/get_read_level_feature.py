import os
from scipy import stats
import numpy as np
import pysam
from collections import Counter, defaultdict
from typing import Tuple

from SpaceTracer.utils.handle_UMI_combine import calculate_UMI_combine_phred, get_most_candidate_allele, handle_quality_matrix, handle_seq, handle_pos, handle_seq_type


def check_UMIconsistence_for_each_geno(count_dict,threshold=1):
    '''
    This function is used to count the consistence or not for each geno and each dict
    Version1: we want to contaion those info: A:9,T:1. Both geno A and T will be counted as 1 UMI inconsistence
    '''
    UMI_DP=sum(count_dict.values())
    if UMI_DP>=threshold:
        norm_count=[count_dict[geno]/UMI_DP for geno in "ATCG"]
        return norm_count
    else:
        return []
    
def get_indel_info(indel_info_list,pos_index_in_raw_matrix):
    """
    similar as handle_pos, input is [(10,1)]. 10 is the indel position to sequence start(exclude soft/hard slip seq)
    output:
    the plus or minus of indel will display the direction between indel and mutation. The plus means indel in the left of mutation
    """
    indel_num=0;indel_length=[];indel_distance=[]
    if len(indel_info_list) == 0:
        indel_num=0
    elif len(indel_info_list) >1 and pos_index_in_raw_matrix=="in":
        indel_num=len(indel_info_list)
        indel_distance=[]
    else:
        indel_num=len(indel_info_list)
        for item in indel_info_list:
            indel_length.append(item[1])
            if pos_index_in_raw_matrix!="in":
                indel_distance.append(pos_index_in_raw_matrix-item[0])
            elif pos_index_in_raw_matrix=="in":
                indel_distance.append(0)
    # indel_length="/".join([str(i) for i in indel_length])
    # indel_distance="/".join([str(i) for i in indel_distance])
    # if indel_num!=0:
    #     print(indel_info_list,pos_index_in_raw_matrix,indel_num,indel_length,indel_distance)
    return indel_num,indel_length,indel_distance

def judge_pos_in_indel(ins_info_list,del_info_list,get_reference_positions):
    pos_list=[]

    def get_pos_list(info_list,pos_list):
        for item in info_list:
            index_num=item[0];count_num=item[1]
            pos=get_reference_positions[index_num-1]
            for i in range(count_num):
                pos_list.append(pos+i+1)
        return pos_list

    pos_list=get_pos_list(ins_info_list,pos_list)
    pos_list=get_pos_list(del_info_list,pos_list)

    return pos_list

def combine_info_from_cigar(cigar_symbol):
    '''
    ## combine: handel cigar;  get_hard_clip_count
    # [(0, 76), (2, 1), (0, 33), (3, 139241), (0, 11)]
    # '76M1D33M139241N11M'
    # the 1st is symbol; and the 2nd is count
    # 0: Match; 1: Insertion; 2: deletion; 3: N; 4: S; 5: H; 6: P; 7: =; 8: X
    '''
    seq_length_before = 0
    pos_length_before = 0

    left_hardclip,right_hardclip=0,0

    seq_cut_start = None; seq_cut_end = None
    pos_cut=[]
    del_info=[]
    for cigars, i in zip(cigar_symbol,range(1,len(cigar_symbol)+1)):
        symbol = cigars[0]
        count = cigars[1]
        if symbol in [6,7,8]:
            # an api for handeling mapping issues "HP=X"
            print(cigar_symbol)  ## LOG
        elif symbol in [0, 1, 2, 4,5]:
            # measure the seq length 
            seq_length_before += count
            if symbol == 0:
                pos_length_before += count
            elif symbol == 4:
                # whether "S" is in this read
                if i == 1:
                    # whether the "S" is in the head or tail
                    seq_cut_start = seq_length_before
                elif i ==len(cigar_symbol):
                    seq_cut_end = seq_length_before
                else:
                    print(cigar_symbol) ## LOG
            elif symbol == 1:
                # whether the "I" is in the cigar
                pos_cut.append((pos_length_before,count))
            elif symbol == 2:
                del_info.append((pos_length_before,count))
            elif symbol==5:
                print(cigar_symbol)  ## LOG
                if i==1:
                    left_hardclip=count
                elif i==len(cigar_symbol):
                    right_hardclip=count
        else:
            pass
    seq_soft_cut = (seq_cut_start, seq_cut_end)
    seq_hard_clip = (left_hardclip,right_hardclip)

    return seq_soft_cut, pos_cut, del_info, seq_hard_clip



def barcode_cell_mapping(mapping_file):
    import pandas as pd
    if mapping_file == "":
        return {}
    
    elif os.path.exists(mapping_file):
        df = pd.read_csv(mapping_file, sep='\t', header=None, names=["CB", "cell"])  
        return dict(zip(df['CB'], df['cell']))
    else:
        raise FileNotFoundError(f"Mapping file '{mapping_file}' does not exist")



def wilcoxon_with_rbc(x, y, alternative='two-sided'):
    """ 
    Perform a Wilcoxon rank-sum test with calculating the Rank-Biserial Correlation (RBC) value.
    Parameters
    ----------
    x,y : array_like
        The data from the two samples.
    alternative : {'two-sided', 'less', 'greater'}, optional
        Defines the alternative hypothesis. Default is 'two-sided'.
        The following options are available:
        * 'two-sided': one of the distributions (underlying `x` or `y`) is
          stochastically greater than the other.
        * 'less': the distribution underlying `x` is stochastically less
          than the distribution underlying `y`.
        * 'greater': the distribution underlying `x` is stochastically greater
          than the distribution underlying `y`.
    Returns
    -------
    statistic : float
        The test statistic under the large-sample approximation that the
        rank sum statistic is normally distributed.
    pvalue : float
        The p-value of the test.
    rbc: float
        The Rank-Biserial Correlation (RBC) value, which is a measure of
        the strength and direction of the association between the two samples.
        It ranges from -1 to 1, where:
        - 1 indicates a perfect positive association,
        - -1 indicates a perfect negative association,
        - 0 indicates no association. 
    """
    x, y = map(np.asarray, (x, y))
    n1 = len(x)
    n2 = len(y)
    alldata = np.concatenate((x, y))
    ranked = stats.rankdata(alldata)
    x_ranks = ranked[:n1]
    R1 = np.sum(x_ranks)
    # Mann-Whitney U statistic for group x
    U = R1 - (n1 * (n1 + 1)) / 2
    # Rank Biserial Correlation (rbc)
    rbc = (2 * U) / (n1 * n2) - 1
    # statistic and p-value with large number normal approximation
    statistic, p_value = stats.ranksums(x, y, alternative=alternative)
    return statistic, p_value, rbc



def calculate_rbc_for_paired_wilcoxon(x, y):
    """ 
    Calculate the Rank-Biserial Correlation (RBC) value for the paired Wilcoxon signed-rank test.
    Parameters
    ----------
    x,y : array_like
        The data from the two samples.
    Returns
    -------
    rbc: float
        The Rank-Biserial Correlation (RBC) value, which is a measure of
        the strength and direction of the association between the two samples.
        It ranges from -1 to 1, where:
        - 1 indicates a perfect positive association,
        - -1 indicates a perfect negative association,
        - 0 indicates no association. 
    """
    from sklearn.preprocessing import StandardScaler
    # format as numpy arrays
    x, y = map(np.asarray, (x, y))
    # test whether the two arrays are of the same length
    if len(x) != len(y):
        raise ValueError("Input arrays must have the same length.")
    try:
        # standardize
        X = np.array([x, y]).T 
        scaler = StandardScaler()
        X_scaled = scaler.fit_transform(X)  # normalizes each column independently
        x_scaled = X_scaled[:, 0]
        y_scaled = X_scaled[:, 1] 
        # Compute RBC
        diff = x_scaled - y_scaled
        non_zero = diff != 0
        ranks = stats.rankdata(np.abs(diff[non_zero]))
        signed_ranks = ranks * np.sign(diff[non_zero])
        W_plus = signed_ranks[signed_ranks > 0].sum()
        W_minus = -signed_ranks[signed_ranks < 0].sum()
        rbc = (W_plus - W_minus) / (W_plus + W_minus)
        return rbc
    except:
        return 0
    


def do_wilicox_sum_test(input_1,input_2,method="two-sided",type="dict"):
    from scipy import stats
    if type=="dict":
        input_list1=[]
        input_list2=[]
        for i in input_1.keys():
            input_list1.extend(i[0]*int(i[1]))

        for m in input_2.keys():
            input_list2.extend(m[0]*int(m[1]))    
    elif type=="list":
        input_list1=input_1
        input_list2=input_2
    else:
        print(f"wrong input in do_wilicox_sum_test: {type}")
        return "no","no"
    
    z_statistic, p_value=stats.ranksums(input_list1,input_list2,alternative=method)
    z_statistic=float(z_statistic);p_value=float(p_value)
    return z_statistic,p_value


### the main function to get read_info_dict 
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
    del in_bam_read
    return result_dict


def detect_read_length(bam_file: str, 
                       sample_size: int = 100,
                       min_consensus: float = 0.9) -> Tuple[int, dict]:

    read_lengths = []
    
    try:
        with pysam.AlignmentFile(bam_file, "rb") as bam:
            for i, read in enumerate(bam):
                if i >= sample_size:
                    break
                
                if read.is_unmapped:
                    continue
                
                read_length = read.query_length
                if read_length is not None and read_length > 0:
                    read_lengths.append(read_length)
    
    except Exception as e:
        raise
    
    if not read_lengths:
        raise ValueError(f"No valid reads found in {bam_file}")
    
    length_counter = Counter(read_lengths)
    total_reads = len(read_lengths)
    
    most_common_length, most_common_count = length_counter.most_common(1)[0]
    consensus_ratio = most_common_count / total_reads
    
    if consensus_ratio < min_consensus:
        raise ValueError(f"The input bam file is under mixed library!")
        
    
    return most_common_length



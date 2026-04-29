

from collections import defaultdict
from functools import reduce
from math import ceil, log10

import numpy as np
from SpaceTracer.utils.logger import get_logger
from SpaceTracer.utils.utils import round_to_nearest_bin

model_name=__name__
logger = get_logger(model_name)


def handle_cigar(cigar_symbol):
    '''
    ## handle cigar
    # [(0, 76), (2, 1), (0, 33), (3, 139241), (0, 11)]
    # '76M1D33M139241N11M'
    # the 1st is symbol; and the 2nd is count
    # 0: Match; 1: Insertion; 2: deletion; 3: N; 4: S; 5: H; 6: P; 7: =; 8: X
    output:
    set_cut may be a turple, contain the softclip information; pos_cut are the indel information: [(before_length, insert number)]
    '''
    seq_length_before = 0
    pos_length_before = 0

    seq_cut_start = None; seq_cut_end = None
    pos_cut=[]
    for cigars, i in zip(cigar_symbol,range(1,len(cigar_symbol)+1)):
        symbol = cigars[0]
        count = cigars[1]
        if symbol in [5,6,7,8]:
            # an api for handling mapping issues "HP=X"
            logger.debug(cigar_symbol)  ## LOG
        elif symbol in [0, 1, 4]:
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
                    logger.debug(cigar_symbol) ## LOG
            elif symbol == 1:
                # whether the "I" is in the cigar
                pos_cut.append((pos_length_before,count))
        else:
            pass
    seq_cut = (seq_cut_start, seq_cut_end)
    return seq_cut, pos_cut


def handle_seq(seq, seq_cut):
    # only support one time for cut
    cut_seq=seq[seq_cut[0]:seq_cut[1]]
    return cut_seq


def handle_pos(pos_matrix,pos_cut):
    if len(pos_cut) == 0:
        cut_pos_matrix = pos_matrix
    elif len(pos_cut) == 1:
        times = pos_cut[0][1]
        pos = pos_cut[0][0]
        cut_pos_matrix = pos_matrix[0:pos] + [""] * times + pos_matrix[pos:]
    else:
        start = 0
        cut_pos_matrix = []
        for item in range(len(pos_cut)):
            pos = pos_cut[item][0]; times = pos_cut[item][1]
            cut_pos_matrix = cut_pos_matrix + pos_matrix[start:pos] + [""] * times
            start=pos
        last_pos = pos_cut[-1][0]
        cut_pos_matrix = cut_pos_matrix + pos_matrix[last_pos:]
    return cut_pos_matrix


def handle_quality_matrix(mutation_in_cutseq_index,seq,cut_seq):
    if len(cut_seq[mutation_in_cutseq_index:]) >= len(cut_seq[:mutation_in_cutseq_index]):
        query_str = cut_seq[mutation_in_cutseq_index:]
        raw_index = seq.index(query_str)
    else:
        query_str = cut_seq[:mutation_in_cutseq_index]
        raw_index = seq.index(query_str) + len(query_str)
    return raw_index


def phred_2_q(phred):
    try:
        phred=int(phred)
        q=1 - (10 ** -(phred/10))
    except ValueError as e:
        q=0
        print(f"wrong phred format: {phred}")
        print("Error:", e)
    return q


def q_2_phred(q):
    try:
        q=float(q)
        if q==1.0:
            phred=40
        else:
            phred=ceil(-log10(1-q) * 10)
    except ValueError as e:
        print(f"wrong q format: {q}")
        print("Error:", e)
        # sys.exit()
        phred=0
    return phred


def trans(qual_list):
    if qual_list ==[]:
        qual_str="NA"
    else:
        qual_str=",".join(qual_list)
    return qual_str


def handel_barcode_name(cell_dict,barcode_name):
    if cell_dict!={}:
        return cell_dict.get(barcode_name, barcode_name)
    else:
        return str(barcode_name)


def handle_seq_type(read,run_type,bins,cell_dict={}):
    if run_type=="visium":
        try:
            CB=read.get_tag("CB").strip()
            UB=read.get_tag("UB").strip()

            barcode_name=str(CB)
            barcode_name=handel_barcode_name(cell_dict,barcode_name)
            UMI_name=barcode_name+"_"+str(UB)
        except:
            return None,None

    elif run_type=="stereo":
        try:
            Cx_raw=int(read.get_tag("Cx"))
            Cy_raw=int(read.get_tag("Cy"))
            if bins !=1:
                Cx=round_to_nearest_bin(Cx_raw,bins)
                Cy=round_to_nearest_bin(Cy_raw,bins)
            else:
                Cx=Cx_raw
                Cy=Cy_raw
            UR=read.get_tag("UR").strip()

            barcode_name=str(Cx)+"_"+str(Cy)
            barcode_name=handel_barcode_name(cell_dict,barcode_name)
            UMI_name=barcode_name+"_"+str(UR)
        except:
            return None,None


    elif run_type=="ST":
        try:
            CB=str(read.get_tag("B0"))
            UB=str(read.get_tag("B3"))
            barcode_name=str(CB)
            barcode_name=handel_barcode_name(cell_dict,barcode_name)
            UMI_name=barcode_name+"_"+str(UB)
        except:
            return None,None

    else:
        raise ValueError(f"The input type is not recognized {run_type}")

    return barcode_name, UMI_name
    

def calculate_UMI_combine_phred_count_dict(count_dict, quality_dict,weigh=0.5):
    """
    The function is used to get all candidate allele and their phred score,
    based on count and quality dict per UMI.
    """
    all_genos=["A","T","C","G"]
    pcr_error = 1e-6
    #no_pcr_error = 1.0 - 3e-5 the reference from smcount
    no_pcr_error = (1.0 - pcr_error) ** 100 # median cycle in RNA-seq is 100 (50-150)
    rightP = 1.0
    sumP = 0.0
    dp=sum(count_dict.values())
    proP_dict=defaultdict(lambda : 1.0)
    pcrP_dict=defaultdict(float)
    likelihood_dict=defaultdict(float)
    phred_dict=defaultdict(float)
    for geno in count_dict.keys():
        ## proP_value means no sequencing error for each geno
        # the likelihood whose allele equal to geno, here the quality is the right prob for one base
        qual_geno_list=[phred_2_q(key)**int(quality_dict[geno][key]) for key in quality_dict[geno].keys()]
        qual_geno=reduce(lambda x, y: x*y, qual_geno_list)
        proP_dict[geno]*=qual_geno
        # the likelihood whose allele not equal to geno
        for other_geno in quality_dict.keys()-set([geno]):
            other_qual_geno_list = [(1-phred_2_q(key))**int(quality_dict[other_geno][key]) for key in quality_dict[other_geno].keys()]
            if other_qual_geno_list == []:
                continue
            other_qual_geno=reduce(lambda x, y: x*y, other_qual_geno_list)
            proP_dict[geno]*=other_qual_geno
        
        ## rightP means no sequencing error, or no base calling error for all base
        rightP = rightP * qual_geno
    
    for geno in all_genos:
        ## pcrP means PCR error
        count_geno = 0 if geno not in count_dict.keys() else count_dict[geno]
        ratio = ( count_geno + 0.5) / (dp + 0.5 * 4)
        pcrP = 10.0 ** (-6.0 * ratio)
        pcrP_dict[geno]=pcrP
    
    # after obtaining [sequencing_error, no_pcr_error, no_sequencing_error, pcr_error], the likelihood of each geno will be calculate
    for geno in all_genos:
        if geno in count_dict.keys():
            base_calling_error = proP_dict[geno]
            no_base_calling_error=rightP
            pcr_error=min([pcrP_dict[char] for char in pcrP_dict.keys() if char != geno])
            likelihood_value = weigh * no_pcr_error * base_calling_error + (1-weigh) * no_base_calling_error * pcr_error 
        else:
            likelihood_value = rightP
            for char in set(all_genos) - set([geno]):
                likelihood_value *= pcrP_dict[char]
                    
        likelihood_dict[geno]=likelihood_value
        sumP += likelihood_value
    
    for geno in likelihood_dict.keys():
        phred_dict[geno] = 0 if sumP <= 0 else q_2_phred(likelihood_dict[geno] / sumP)

    return phred_dict


def calculate_UMI_combine_phred(count_list, quality_dict, weigh=0.5):
    """
    The function is used to get all candidate allele and their phred score,
    based on count and quality dict per UMI.
    """
    all_genos = ["A", "T", "C", "G"]
    pcr_error = 1e-6
    no_pcr_error = (1.0 - pcr_error) ** 100
    rightP = 1.0
    sumP = 0.0
    
    dp = sum(count_list)  
    
    proP_dict = defaultdict(lambda: 1.0)
    pcrP_dict = defaultdict(float)
    likelihood_dict = defaultdict(float)
    phred_dict = defaultdict(float)
    
    for geno, dp_count in zip("ATCG", count_list):
        if dp_count == 0:
            continue
        
        qual_geno_list = [phred_2_q(key) ** int(quality_dict[geno][key]) 
                          for key in quality_dict[geno].keys()]
        qual_geno = reduce(lambda x, y: x * y, qual_geno_list)
        proP_dict[geno] *= qual_geno
        
        for other_geno in quality_dict.keys() - set([geno]):
            other_qual_geno_list = [(1 - phred_2_q(key)) ** int(quality_dict[other_geno][key]) 
                                    for key in quality_dict[other_geno].keys()]
            if other_qual_geno_list:
                other_qual_geno = reduce(lambda x, y: x * y, other_qual_geno_list)
                proP_dict[geno] *= other_qual_geno
        
        rightP = rightP * qual_geno
    
    for index, geno in enumerate("ATCG"):
        count_geno = count_list[index]
        ratio = (count_geno + 0.5) / (dp + 0.5 * 4)
        pcrP = 10.0 ** (-6.0 * ratio)
        pcrP_dict[geno] = pcrP
    

    for index, geno in enumerate("ATCG"):
        if count_list[index] != 0:
            base_calling_error = proP_dict[geno]
            no_base_calling_error = rightP

            other_pcrP = [pcrP_dict[char] for char in pcrP_dict.keys() if char != geno]
            pcr_error_value = min(other_pcrP) if other_pcrP else 1.0
            
            likelihood_value = weigh * no_pcr_error * base_calling_error + \
                              (1 - weigh) * no_base_calling_error * pcr_error_value
        else:
            likelihood_value = rightP
            for char in set(all_genos) - set([geno]):
                likelihood_value *= pcrP_dict[char]
        
        likelihood_dict[geno] = likelihood_value
        sumP += likelihood_value
    
    for geno in likelihood_dict.keys():
        phred_dict[geno] = 0 if sumP <= 0 else q_2_phred(likelihood_dict[geno] / sumP)
    
    return phred_dict


def get_most_candidate_allele(phred_dict,ref_allele):
    """
    To get the most candidate allele and it's phred
    """
    rank_list=sorted(phred_dict.items(), key = lambda item:item[1], reverse=True)
    major_allele=rank_list[0][0]; major_allele_phred=rank_list[0][1]
    ref_allele_phred=phred_dict[ref_allele]

    if  major_allele != ref_allele and ref_allele_phred>=major_allele_phred:
        candidate_allele=ref_allele;phred=ref_allele_phred
    else:
        candidate_allele=major_allele;phred=major_allele_phred

    return candidate_allele,phred

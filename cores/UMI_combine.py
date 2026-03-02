from collections import defaultdict
from functools import reduce
from SpaceTracer.utils.handle_UMI_combine import handle_cigar, handle_pos, handle_quality_matrix, handle_seq, handle_seq_type, phred_2_q, q_2_phred
from SpaceTracer.utils.read_files import read_bam

from SpaceTracer.utils.logger import get_logger
model_name=__name__
logger = get_logger("<core function>: "+model_name)

def combine_UMI_spot_for_both_mosaic_and_error(bam_file,check_mosaic, check_error,run_type,identifier):
    """
    identifier:"chr1_1000_ref", do not provide alt here.
    """
    
    identifier_list=identifier.split("_")
    chrom=identifier_list[0]
    pos=int(identifier_list[1])
    ref=identifier_list[2]
    new_list=[]
    all_genos=["A","T","C","G"]
    error_allele=None
    # print("###",bam_file,check_mosaic, check_error,run_type,identifier)

    handle_bam=read_bam(bam_file)
    site_barcode_UMI_dict,strand=handle_reads_per_pos_read_count_and_strand(handle_bam,chrom,pos,run_type)
    threshold=3
    consensus_read_count, consensus_read_quality, error_allele = UMI_combination_spot_ind_and_judge_error(check_mosaic,check_error,site_barcode_UMI_dict,chrom,pos,ref,threshold)

    if check_mosaic:
        for barcode in consensus_read_count.keys():
            alt=",".join([alt_geno for alt_geno in consensus_read_count[barcode].keys() if alt_geno != ref])
            if alt=="":
                alt="."
            consensus_read_count_list=[]
            consensus_read_qual_dict=defaultdict(str)
            for geno in all_genos:
                geno_count=0 if geno not in consensus_read_count[barcode].keys() else consensus_read_count[barcode][geno]
                consensus_read_count_list.append(str(geno_count))
                
                if geno not in consensus_read_quality[barcode].keys():
                    geno_qual_str="NA" 
                else:
                    geno_qual_list=[str(phred)+":"+str(consensus_read_quality[barcode][geno][phred]) for phred in consensus_read_quality[barcode][geno].keys()]
                    geno_qual_str=",".join(geno_qual_list)
                consensus_read_qual_dict[geno]=geno_qual_str

            consensus_read_count_str=",".join(consensus_read_count_list)

            barcode_list=[chrom,pos,".",ref,alt,barcode,consensus_read_count_str,consensus_read_qual_dict["A"],consensus_read_qual_dict["T"],consensus_read_qual_dict["C"],consensus_read_qual_dict["G"]]
            new_list.append(barcode_list)
    if error_allele:
        error_list=[chrom,pos,ref,error_allele,strand]
    else:
        error_list=[]
    return new_list,error_list


def UMI_combination_spot_ind_and_judge_error(check_mosaic,check_error,site_barcode_UMI_dict,chrom,pos,ref,threshold):
    # mosaics
    spot_number=len(site_barcode_UMI_dict.keys())
    all_genos=["A","T","C","G"]
    consensus_read_count={}
    consensus_read_quality={}
    raw_count_dict=defaultdict(int)
    
    # errors
    pcr_errors,lysis_errors=0,0
    pcr_alts,lysis_alts=[],[]
    UMI_dp=0
    dp=0
    error_allele=None

    for barcode in site_barcode_UMI_dict.keys():
        consensus_read_count[barcode]=defaultdict(int)
        consensus_read_quality[barcode]=defaultdict(dict)
        for UMI in site_barcode_UMI_dict[barcode]:
            count_dict=site_barcode_UMI_dict[barcode][UMI]["count"]

            if check_mosaic: 
                quality_dict=site_barcode_UMI_dict[barcode][UMI]["quality"]
                phred_dict=calculate_UMI_combine_phred(count_dict,quality_dict,weigh=0.5)
                candidate_allele,phred=get_most_candidate_allele(phred_dict,ref)

                consensus_read_count[barcode][candidate_allele]+=1
                if str(phred) not in consensus_read_quality[barcode][candidate_allele].keys():
                    consensus_read_quality[barcode][candidate_allele][str(phred)]=0

                consensus_read_quality[barcode][candidate_allele][str(phred)]+=1
                for geno in count_dict.keys():
                    raw_count_dict[geno]+=count_dict[geno]

            if check_error:
                _, lysis_error,_,lysis_alt=check_errors(count_dict,ref,threshold)
                # pcr_errors += pcr_error
                lysis_errors += lysis_error
                # pcr_alts += pcr_alt
                lysis_alts += lysis_alt

    if check_mosaic:
        alt=",".join([alt_geno for alt_geno in consensus_read_count.keys() if alt_geno != ref])
        if alt =="":
            alt="."

        consensus_read_count_list=[]
        raw_read_count_list=[]

        consensus_read_qual_dict=defaultdict(str)
        for geno in all_genos:
            geno_count=0 if geno not in consensus_read_count.keys() else consensus_read_count[geno]
            consensus_read_count_list.append(str(geno_count))

            raw_count=0 if geno not in raw_count_dict.keys() else raw_count_dict[geno]
            raw_read_count_list.append(str(raw_count))

            if geno not in consensus_read_quality.keys():
                geno_qual_str="NA" 
            else:
                geno_qual_list=[str(phred)+":"+str(consensus_read_quality[geno][phred]) for phred in consensus_read_quality[geno].keys()]
                geno_qual_str=",".join(geno_qual_list)
            consensus_read_qual_dict[geno]=geno_qual_str

        # consensus_read_count_str=",".join(consensus_read_count_list)
        # raw_read_count_str=",".join(raw_read_count_list)
        # site_info=[chrom, pos,".", ref, alt, spot_number, consensus_read_count_str,consensus_read_qual_dict["A"],consensus_read_qual_dict["T"],consensus_read_qual_dict["C"],consensus_read_qual_dict["G"],raw_read_count_str]
    # else:
        # site_info=[]

    if check_error:
        if lysis_errors >1 or lysis_errors==0:
            error_allele=None
        else:
            error_allele=lysis_alts[0]

    return consensus_read_count, consensus_read_quality, error_allele


def check_errors(count_dict,ref,threshold=3):
    count_above_threshold = sum(1 for allele in count_dict if allele != ref and count_dict[allele] >= threshold)
    count_nonzero = sum(1 for allele in count_dict if count_dict[allele] > 0)
    pcr_error = 1 if count_above_threshold >= 1 and count_nonzero >= 2 else 0
    pcr_alt=[]
    if pcr_error:
        pcr_alt = [allele for allele in count_dict if allele != ref and count_dict[allele] >= threshold]
    
    # lysis_error: require all reads in one UMI must same
    lysis_error = 0
    lysis_alt=[]
    for allele in count_dict:
        if allele != ref and count_dict[allele] > threshold:
            if all(count_dict[other_allele] == 0 for other_allele in count_dict if other_allele != allele):
                lysis_error = 1
                lysis_alt.append(allele)
                break

    return pcr_error, lysis_error,pcr_alt,lysis_alt


def handle_reads_per_pos_read_count_and_strand(bam_handle,chrom,pos,run_type):
    '''
    input:
    bam_handle: the bam file handled by pysam
    
    output:
    a dict containing UMI and read information: {barcode_name: {UMI_name: {"count": 1, "quality": {"A":{30:9, 10:1}, "T":{10:1}}}}}
    '''
    reads=bam_handle.fetch(chrom,pos-1,pos,multiple_iterators=True)
    pos_index = pos-1
    site_barcode_UMI_dict={}
    reverse_dp=0
    forward_dp=0
    for item in reads:
        # try:
        #     # a part of reads didn't have the information of "CB", because the "CR" didn't pass QC
        #     if run_type=="visium":
        #         CB=item.get_tag("CB").strip()
        #         UB=item.get_tag("UB").strip()

        #         barcode_name=str(CB)
        #         UMI_name=str(UB)

        #     elif run_type=="stereo":
        #         Cx=str(item.get_tag("Cx"))
        #         Cy=str(item.get_tag("Cy"))
        #         UR=item.get_tag("UR").strip()

        #         barcode_name=Cx+"_"+Cy
        #         UMI_name=str(UR)
        #     elif run_type=="ST":
        #         CB=str(item.get_tag("B0"))
        #         UB=str(item.get_tag("B3"))
        #         barcode_name=str(CB)
        #         UMI_name=str(UB)
        #     else:
        #         continue
                
        # except:
        #     continue
        barcode_name,UMI_name=handle_seq_type(item,run_type,1)
        if barcode_name==None or UMI_name==None:
            continue
        try:
            item.get_reference_positions().index(pos_index)
        except:
            continue

        if item.is_reverse in [True,"TRUE","true","True"]:
            reverse_dp+=1
        else:
            forward_dp+=1        

        seq_cut, pos_cut = handle_cigar(item.cigar)
        cut_seq=handle_seq(item.seq, seq_cut)
        cut_pos=handle_pos(item.get_reference_positions(), pos_cut)

        if pos_index in cut_pos:
            geno = cut_seq[cut_pos.index(pos_index)]
            if geno not in "ATCG":
                logger.debug(f"The site in this read {item.query_name} is not ATCG, pass.")
                continue

            raw_index = handle_quality_matrix(cut_pos.index(pos_index),item.seq,cut_seq)
            try:
                qualities=item.get_forward_qualities()
                quality=qualities[raw_index]
            except:
                logger.debug(f"The site in this read {item.query_name} do not has quality, pass.")
                continue

            if barcode_name not in site_barcode_UMI_dict.keys():
                site_barcode_UMI_dict[barcode_name]=defaultdict(dict)

            if UMI_name not in site_barcode_UMI_dict[barcode_name].keys():
                site_barcode_UMI_dict[barcode_name][UMI_name]["count"]=defaultdict(int)
                site_barcode_UMI_dict[barcode_name][UMI_name]["quality"]={"A":defaultdict(int),"T":defaultdict(int),"C":defaultdict(int),"G":defaultdict(int)}

            site_barcode_UMI_dict[barcode_name][UMI_name]["count"][geno]+=1
            site_barcode_UMI_dict[barcode_name][UMI_name]["quality"][geno][quality]+=1
    
    if reverse_dp>=forward_dp:
        major_read_strand="-"
    elif reverse_dp<forward_dp:
        major_read_strand="+"
    else:
        major_read_strand="unknown"
        
    return site_barcode_UMI_dict, major_read_strand


def calculate_UMI_combine_phred(count_dict, quality_dict,weigh=0.5):
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


# following the last function, 
def get_most_candidate_allele(phred_dict,ref_allele):
    """
    To get the most candidate allele and it's phred
    """
    rank_list=sorted(phred_dict.items(), key = lambda item:item[1], reverse=True)
    major_allele=rank_list[0][0]; major_allele_phred=rank_list[0][1]
    # major_allele_count=count_dict[major_allele]

    # ref_allele_count=count_dict[ref_allele]
    ref_allele_phred=phred_dict[ref_allele]

    if  major_allele != ref_allele and ref_allele_phred>=major_allele_phred:
        candidate_allele=ref_allele;phred=ref_allele_phred
    else:
        candidate_allele=major_allele;phred=major_allele_phred

    return candidate_allele,phred

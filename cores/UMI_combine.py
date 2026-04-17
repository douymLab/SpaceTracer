from collections import defaultdict
from SpaceTracer.utils.handle_UMI_combine import calculate_UMI_combine_phred, get_most_candidate_allele, handle_cigar, handle_pos, handle_quality_matrix, handle_seq, handle_seq_type, phred_2_q, q_2_phred

from SpaceTracer.utils.logger import get_logger
model_name=__name__
logger = get_logger("<core function>: "+model_name)


# def combine_UMI_spot_for_both_mosaic_and_error(reads,check_mosaic, check_error,run_type,identifier):
#     """
#     identifier:"(1000,ref)", do not provide alt here.
#     """
    
#     # identifier_list=identifier.split("_")
#     chrom=identifier[0]
#     pos=int(identifier[1])
#     ref=identifier[2]
#     new_list=[]
#     all_genos=["A","T","C","G"]
#     error_allele=None
#     # print("###",bam_file,check_mosaic, check_error,run_type,identifier)

#     # handle_bam=handel
#     site_barcode_UMI_dict,strand=handle_reads_per_pos_read_count_and_strand(reads,pos,run_type)
#     threshold=3
#     consensus_read_count, consensus_read_quality, error_allele = UMI_combination_spot_ind_and_judge_error(check_mosaic,check_error,site_barcode_UMI_dict,chrom,pos,ref,threshold)

#     if check_mosaic:
#         for barcode in consensus_read_count.keys():
#             alt=",".join([alt_geno for alt_geno in consensus_read_count[barcode].keys() if alt_geno != ref])
#             if alt=="":
#                 alt="."
#             consensus_read_count_list=[]
#             consensus_read_qual_dict=defaultdict(str)
#             for geno in all_genos:
#                 geno_count=0 if geno not in consensus_read_count[barcode].keys() else consensus_read_count[barcode][geno]
#                 consensus_read_count_list.append(str(geno_count))
                
#                 if geno not in consensus_read_quality[barcode].keys():
#                     geno_qual_str="NA" 
#                 else:
#                     geno_qual_list=[str(phred)+":"+str(consensus_read_quality[barcode][geno][phred]) for phred in consensus_read_quality[barcode][geno].keys()]
#                     geno_qual_str=",".join(geno_qual_list)
#                 consensus_read_qual_dict[geno]=geno_qual_str

#             consensus_read_count_str=",".join(consensus_read_count_list)

#             barcode_list=[chrom,pos,strand,ref,alt,barcode,consensus_read_count_str,consensus_read_qual_dict["A"],consensus_read_qual_dict["T"],consensus_read_qual_dict["C"],consensus_read_qual_dict["G"]]
#             new_list.append(barcode_list)
#     if error_allele:
#         error_list=[chrom,pos,ref,error_allele,strand]
#     else:
#         error_list=[]
#     return new_list,error_list


def UMI_combination_spot_ind_and_judge_error(check_mosaic,check_error,site_barcode_UMI_dict,ref,threshold):
    # mosaics
    spot_number=len(site_barcode_UMI_dict.keys())
    all_genos=["A","T","C","G"]
    consensus_read_count={}
    consensus_read_quality={}
    raw_count_list=[0,0,0,0]
    
    # errors
    pcr_errors,lysis_errors=0,0
    pcr_alts,lysis_alts=[],[]
    UMI_dp=0
    dp=0
    error_allele=None

    for barcode in site_barcode_UMI_dict.keys():
        consensus_read_count[barcode]=[0,0,0,0]
        consensus_read_quality[barcode]={}

        for UMI in site_barcode_UMI_dict[barcode]:
            count_list=site_barcode_UMI_dict[barcode][UMI]["count"]

            if check_mosaic: 
                quality_dict=site_barcode_UMI_dict[barcode][UMI]["quality"]
                phred_dict=calculate_UMI_combine_phred(count_list,quality_dict,weigh=0.5)
                candidate_allele,phred=get_most_candidate_allele(phred_dict,ref)
                index="ATCG".index(candidate_allele) # other alleles has been removed before

                consensus_read_count[barcode][index]+=1

                if candidate_allele not in consensus_read_quality[barcode].keys():
                    consensus_read_quality[barcode][candidate_allele] = {}
                    
                if str(phred) not in consensus_read_quality[barcode][candidate_allele].keys():
                    consensus_read_quality[barcode][candidate_allele][str(phred)]=0

                consensus_read_quality[barcode][candidate_allele][str(phred)]+=1

                for index,dp in enumerate(count_list):
                    raw_count_list[index]+=dp

            if check_error:
                _, lysis_error,_,lysis_alt=check_errors(count_list,ref,threshold)
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

        consensus_read_qual_dict={}
        for geno in all_genos:
            geno_count=0 if geno not in consensus_read_count.keys() else consensus_read_count[geno]
            consensus_read_count_list.append(str(geno_count))
            index="ATCG".index(geno)

            raw_count=raw_count_list[index]
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


def check_errors(count_list,ref,threshold=3):
    ref_index="ATCG".index(ref)
    try:
        count_above_threshold = sum(1 for index,dp in enumerate(count_list) if index != ref_index and dp >= threshold)
    except:
        print("*****************",[ index for index,dp in enumerate(count_list)])
        print("*****************",[ dp for index,dp in enumerate(count_list)])
        raise
    
    # count_above_threshold = sum(1 for allele in count_dict if allele != ref and count_dict[allele] >= threshold)

    # count_nonzero = sum(1 for allele in count_dict if count_dict[allele] > 0)
    count_nonzero = sum(1 for dp in count_list if dp > 0)
    pcr_error = 1 if count_above_threshold >= 1 and count_nonzero >= 2 else 0
    pcr_alt=[]
    if pcr_error:
        # pcr_alt = [allele for allele in count_dict if allele != ref and count_dict[allele] >= threshold]
        pcr_alt = ["ATCG"[index] for index,dp in enumerate(count_list) if index != ref_index and dp >= threshold]
    
    # lysis_error: require all reads in one UMI must same
    lysis_error = 0
    lysis_alt=[]
    # for allele in count_dict:
    for index,dp in enumerate(count_list):
        if index != ref_index and dp > threshold:
            if all(count_list[other_allele_index] == 0 for other_allele_index in [0,1,2,3] if other_allele_index != index):
                lysis_error = 1
                lysis_alt.append("ATCG"[index])
                break

    return pcr_error, lysis_error,pcr_alt,lysis_alt


def handle_reads_per_pos_read_count_and_strand(reads,pos,run_type):
    '''
    input:
    bam_handle: the bam file handled by pysam
    
    output:
    a dict containing UMI and read information: {barcode_name: {UMI_name: {"count": 1, "quality": {"A":{30:9, 10:1}, "T":{10:1}}}}}
    '''
    # reads=bam_handle.fetch(chrom,pos-1,pos,multiple_iterators=True)
    pos_index = pos-1
    site_barcode_UMI_dict={}
    reverse_dp=0
    forward_dp=0
    for item in reads:
      
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
                # logger.debug(f"The site in this read {item.query_name} is not ATCG, pass.")
                continue

            raw_index = handle_quality_matrix(cut_pos.index(pos_index),item.seq,cut_seq)
            try:
                qualities=item.get_forward_qualities()
                quality=qualities[raw_index]
            except:
                # logger.debug(f"The site in this read {item.query_name} do not has quality, pass.")
                continue

            if UMI_name not in site_barcode_UMI_dict[barcode_name].keys():
                site_barcode_UMI_dict[barcode_name][UMI_name] = {
                    "count": [0,0,0,0],
                    "quality": {}
                }

            umi_entry = site_barcode_UMI_dict[barcode_name][UMI_name]

            umi_entry["count"]["ATCG".index(geno)] += 1

            if geno not in umi_entry["quality"]:
                umi_entry["quality"][geno] = {}

            umi_entry["quality"][geno][quality] = umi_entry["quality"][geno].get(quality, 0) + 1

    
    if reverse_dp>=forward_dp:
        major_read_strand="-"
    elif reverse_dp<forward_dp:
        major_read_strand="+"
    else:
        major_read_strand="unknown"
        
    return site_barcode_UMI_dict, major_read_strand



def scan_region_reads_once_for_targets(reads, sites, run_type):
    """
    sites:
        [(pos, ref, alt, check_mosaic, check_error), ...]

    返回:
        {
            pos: (site_barcode_UMI_dict, major_read_strand)
        }
    """
    target_pos0_set = {pos - 1 for pos, *_ in sites}

    per_pos_data = {
        pos: {
            "site_barcode_UMI_dict": {},
            "forward_dp": 0,
            "reverse_dp": 0,
        }
        for pos, *_ in sites
    }

    for item in reads:
        barcode_name, UMI_name = handle_seq_type(item, run_type, 1)
        if barcode_name is None or UMI_name is None:
            continue

        try:
            seq_cut, pos_cut = handle_cigar(item.cigar)
            cut_seq = handle_seq(item.seq, seq_cut)
            cut_pos = handle_pos(item.get_reference_positions(), pos_cut)
            qualities = item.get_forward_qualities()
        except Exception:
            continue

        if not cut_pos or not cut_seq or qualities is None:
            continue

        is_reverse = item.is_reverse in [True, "TRUE", "true", "True"]

        for idx, ref_pos0 in enumerate(cut_pos):
            if ref_pos0 not in target_pos0_set:
                continue

            pos1 = ref_pos0 + 1
            geno = cut_seq[idx]

            if geno not in "ATCG":
                continue

            try:
                raw_index = handle_quality_matrix(idx, item.seq, cut_seq)
                quality = qualities[raw_index]
            except Exception:
                continue

            pos_bucket = per_pos_data[pos1]

            if is_reverse:
                pos_bucket["reverse_dp"] += 1
            else:
                pos_bucket["forward_dp"] += 1

            site_barcode_UMI_dict = pos_bucket["site_barcode_UMI_dict"]

            if barcode_name not in site_barcode_UMI_dict:
                site_barcode_UMI_dict[barcode_name] = {}

            if UMI_name not in site_barcode_UMI_dict[barcode_name].keys():
                site_barcode_UMI_dict[barcode_name][UMI_name] = {
                    "count": [0,0,0,0],
                    "quality": {}
                }

            umi_entry = site_barcode_UMI_dict[barcode_name][UMI_name]

            umi_entry["count"]["ATCG".index(geno)] += 1

            if geno not in umi_entry["quality"]:
                umi_entry["quality"][geno] = {}

            umi_entry["quality"][geno][quality] = umi_entry["quality"][geno].get(quality, 0) + 1

    result = {}
    for pos, bucket in per_pos_data.items():
        reverse_dp = bucket["reverse_dp"]
        forward_dp = bucket["forward_dp"]

        if reverse_dp >= forward_dp:
            major_read_strand = "-"
        else:
            major_read_strand = "+"

        result[pos] = (bucket["site_barcode_UMI_dict"], major_read_strand)

    return result


def process_single_region_for_umi_combine(bam_handle, region_info, seq_type):
    chrom, start, end, sites = region_info

    try:
        reads = bam_handle.fetch(chrom, start, end)
        pos_data_map = scan_region_reads_once_for_targets(reads, sites, seq_type)
    except Exception as exc:
        logger.warning(f"[UmiCombine] region failed {chrom}:{start}-{end}: {exc}")
        return []

    results = []
    for pos, ref, alt, check_mosaic, check_error in sites:
        try:
            site_barcode_UMI_dict, strand = pos_data_map.get(pos, ({}, "unknown"))

            mosaic_spot_list, error_list = summarize_UMI_spot_for_both_mosaic_and_error(
                site_barcode_UMI_dict=site_barcode_UMI_dict,
                strand=strand,
                check_mosaic=check_mosaic,
                check_error=check_error,
                identifier=(chrom, pos, ref),
            )
            results.append((mosaic_spot_list, error_list))

        except Exception as exc:
            raise 
            logger.warning(f"[UmiCombine] worker failed for {(chrom, pos, ref)}: {exc}")
            results.append((None, None))

    return results


def summarize_UMI_spot_for_both_mosaic_and_error(site_barcode_UMI_dict,
                strand,
                check_mosaic,
                check_error,
                identifier,
                threshold=3):

    chrom = identifier[0]
    pos = int(identifier[1])
    ref = identifier[2]

    new_list = []
    all_genos = ["A", "T", "C", "G"]
    error_allele = None

    consensus_read_count, consensus_read_quality, error_allele = UMI_combination_spot_ind_and_judge_error(check_mosaic,check_error,site_barcode_UMI_dict,ref,threshold)

    if check_mosaic:
        for barcode in consensus_read_count.keys():
            # alt=",".join([alt_geno for alt_geno in consensus_read_count[barcode].keys() if alt_geno != ref])
            alt=",".join([alt_geno for alt_geno,dp in zip("ATCG",consensus_read_count[barcode]) if alt_geno != ref and dp !=0])
            if alt=="":
                alt="."
            consensus_read_count_list=[str(i) for i in consensus_read_count[barcode]]
            consensus_read_qual_dict=defaultdict(str)
            for geno in all_genos:
                # geno_count=0 if geno not in consensus_read_count[barcode].keys() else consensus_read_count[barcode][geno]
                # consensus_read_count_list.append(str(geno_count))
                
                if geno not in consensus_read_quality[barcode].keys():
                    geno_qual_str="NA" 
                else:
                    geno_qual_list=[str(phred)+":"+str(consensus_read_quality[barcode][geno][phred]) for phred in consensus_read_quality[barcode][geno].keys()]
                    geno_qual_str=",".join(geno_qual_list)
                consensus_read_qual_dict[geno]=geno_qual_str

            consensus_read_count_str=",".join(consensus_read_count_list)

            barcode_list=[chrom,pos,strand,ref,alt,barcode,consensus_read_count_str,consensus_read_qual_dict["A"],consensus_read_qual_dict["T"],consensus_read_qual_dict["C"],consensus_read_qual_dict["G"]]
            new_list.append(barcode_list)

    if error_allele:
        error_list=[chrom,pos,ref,error_allele,strand]
    else:
        error_list=[]

    return new_list,error_list

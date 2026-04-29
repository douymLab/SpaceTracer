from collections import defaultdict, Counter
from functools import partial
import gc
import multiprocessing
import os
# import re
import subprocess
import sys
import uuid
import warnings
# from tabnanny import check
# import time
# from pyfaidx import Fasta
import pandas as pd
import pysam
import argparse
import scipy.stats
import gzip
from pybedtools import BedTool
# import pysamstats
# import sys
from math import log10,ceil
from functools import reduce
from module.UMI_combine import calculate_UMI_combine_phred,get_most_candidate_allele,handle_cigar,handle_seq,handle_pos,handle_quality_matrix,UMI_combination_spot
from utils import get_chr_size, round_to_nearest_bin
# from utils import check_dir
import psutil

# from memory_profiler import memory_usage
# from memory_profiler import profile


    
def print_memory_usage(note=""):
    process = psutil.Process(os.getpid())
    mem = process.memory_info().rss / 1024 / 1024  # 单位MB
    print(f"[内存监控] {note} 当前内存: {mem:.2f} MB")


def handel_candidate_informative_SNP_site(species,germline_file_name, tmp_dir):
    '''
    Input:
    #chrom  site    genotype        allele  prior
    chr12   52487210        refhom  A       1.0
    chr12   52487211        refhom  A       1.0
    chr12   52487212        het     C,G     0.4225,0.5775
    chr12   52487215        het     C,T     0.7311000000000001,0.2689
    

    Output:
    chr12	52487212	52487212    het	C	G	0.4225,0.5775
    chr12	52487215	52487215    het	C	T	0.7311000000000001,0.2689
    '''
    tmp=str(uuid.uuid4())
    
    site_bed_file=os.path.join(tmp_dir,os.path.basename(germline_file_name)+tmp+".tmp.bed")
    # print(site_bed_file)
    # if not os.path.exists(site_bed_file):
    if species=="human":
        run_shell="awk -F'\t' '$3 == \"het\" {split($4, allele, \",\");split($5, allele_count, \",\");if (allele_count[1] > 20 && allele_count[2] > 20) " \
                    "{split($6, prior, \",\");if (prior[1] > 0.01 && prior[2] > 0.01) " \
                    "{total = allele_count[1] + allele_count[2];if (total > 50) {print $1, $2, $2,$3,allele[1], allele[2], $5,$6;}}}}' OFS=\"\\t\" %s > %s" \
                    "" % (germline_file_name, site_bed_file)
    else:
        run_shell="awk -F'\t' '$3 == \"het\" {split($4, allele, \",\");split($5, allele_count, \",\");if (allele_count[1] > 20 && allele_count[2] > 20) " \
                    "{split($6, prior, \",\");if (prior[1] > 0 && prior[2] > 0) " \
                    "{total = allele_count[1] + allele_count[2];if (total > 50) {print $1, $2, $2,$3,allele[1], allele[2], $5,$6;}}}}' OFS=\"\\t\" %s > %s" \
                    "" % (germline_file_name, site_bed_file)

    result=subprocess.run(run_shell,shell=True)
    print_memory_usage("line70")
    if result.returncode!=0:
        print(f'Something wrong when generate the filtered bed file for {germline_file_name}.')
    
    return site_bed_file


def short_bed_file(gene_review_file, site_bed_file, tmp_dir):
    '''
    please make sure the bedtools is availiable. This step is used to get the gene list
    '''
    tmp=str(uuid.uuid4())

    short_gene_review_file=os.path.join(tmp_dir,os.path.basename(gene_review_file)+tmp+".tmp.short.genereview.bed")
    result=subprocess.run("bedtools intersect -a %s -b %s -u | grep -v \"MT-\" > %s" % \
                          (gene_review_file, site_bed_file, short_gene_review_file),shell=True)
    print_memory_usage("line85")
    
    if result.returncode!=0:
        print(f'Something wrong when generate the shorted gene review file for {gene_review_file}.')

    return short_gene_review_file


def short_gencode_file(gencode_file,site_bed_file,tmp_dir,rerun=True):
    tmp=str(uuid.uuid4())

    short_gencode_file_name=os.path.join(tmp_dir,os.path.basename(gencode_file)+tmp+".tmp.short.gene.bed")
    if os.path.exists(short_gencode_file_name) and rerun==False:
        return short_gencode_file_name
    
    if gencode_file.split(".")[-1]=="gz":
        run_shell="zcat %s |grep -v \"#\"|awk '$3==\"gene\"{print $1, $4,$5,$12,$14}' OFS=\"\\t\" |bedtools intersect -a - -b %s -u | grep -v \"MT-\"|sed 's/;//g' > %s" \
                "" % (gencode_file,site_bed_file,short_gencode_file_name)
    else:
        run_shell="cat %s |grep -v \"#\"|awk '$3==\"gene\"{print $1, $4,$5,$12,$14}' OFS=\"\\t\" |bedtools intersect -a - -b %s -u | grep -v \"MT-\" |sed 's/;//g' > %s" \
                "" % (gencode_file,site_bed_file,short_gencode_file_name)
    

    result=subprocess.run(run_shell,shell=True)
    print_memory_usage("line109")

    if result.returncode!=0:
        print(f'Something wrong when generate the shorted gene review file for {short_gencode_file_name}.')
        sys.exit()
    return short_gencode_file_name


def handle_genereview_file(review_file,chrom_colum=1,pos_start=2,pos_end=3,gene_name_colum=4):
    gene_bed_list=list()
    if review_file.split(".")[-1]=="gz":
        f = gzip.open(review_file, 'rb')
        for line in f.readlines():
            s = line.decode().strip().split()
            chr=s[chrom_colum-1]; pos1=int(s[pos_start-1]); pos2=int(s[pos_end-1]); gene_name=int(s[gene_name_colum-1])
            gene_bed_list.append((chr,pos1,pos2,gene_name))
    else:
        f =open(review_file,"r")
        for line in f.readlines():
            s = line.strip().split()
            chr=s[chrom_colum-1]; pos1=int(s[pos_start-1]); pos2=int(s[pos_end-1]); gene_name=int(s[gene_name_colum-1])
            gene_bed_list.append((chr,pos1,pos2,gene_name))
    #print(pos_list)
    return gene_bed_list


def short_mosaic_sites_and_combine_with_germ(ind_geno_file, germ_bed_file, tmp_dir):
    """
    ind_geno_file:
        #chrom	site	ID	germline	mutant	cluster	spot_number	consensus_read_count	genotype	p_mosaic	Gi	vaf
    """
    tmp=str(uuid.uuid4())

    short_ind_file=os.path.join(tmp_dir,os.path.basename(ind_geno_file)+tmp+".tmp.short.mosaic.bed")
    combine_ind_germ_file=os.path.join(tmp_dir,os.path.basename(ind_geno_file)+tmp+".tmp.ind.mosaic.bed")
    # result1=subprocess.run("awk '$13=\"mosaic\"{print $1, $2, $2,$13, $16,$17, $18, \"NA\"}' OFS=\"\t\" %s > %s; cat %s %s > %s" \
    #                       "" % (ind_geno_file,short_ind_file,short_ind_file, germ_bed_file,combine_ind_germ_file),shell=True)
    result1=subprocess.run("awk '$9=\"mosaic\"{print $1, $2, $2,$9, $4,$5, $12, \"NA\"}' OFS=\"\t\" %s > %s; cat %s %s > %s" \
                          "" % (ind_geno_file,short_ind_file,short_ind_file, germ_bed_file,combine_ind_germ_file),shell=True)
    print_memory_usage("line147")

    if result1.returncode != 0:
        print(result1.stderr)

    return short_ind_file, combine_ind_germ_file


### This four function are used to get right quality and base from pysam 
##  not used!!!!!!!!!!!
def get_geno_from_alignment(pos,item,CBtag,UBtag):
    def handle_cigar(ciagr_symbol):
        '''
        ## handel cigar
        # [(0, 76), (2, 1), (0, 33), (3, 139241), (0, 11)]
        # '76M1D33M139241N11M'
        # the 1st is symbol; and the 2nd is count
        # 0: Match; 1: Insertion; 2: deletion; 3: N; 4: S; 5: H; 6: P; 7: =; 8: X
        '''
        seq_length_before = 0
        pos_length_before = 0

        seq_cut_start = None; seq_cut_end = None
        pos_cut=[]
        for cigars, i in zip(ciagr_symbol,range(1,len(ciagr_symbol)+1)):
            symbol = cigars[0]
            count = cigars[1]
            if symbol in [5,6,7,8]:
                # an api for handeling mapping issues "HP=X"
                print(ciagr_symbol)  ## LOG
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
                    elif i ==len(ciagr_symbol):
                        seq_cut_end = seq_length_before
                    else:
                        print(ciagr_symbol) ## LOG
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
                #print(cut_pos_matrix)
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

    Name,geno,quality='','',''
    
    try:
        CB=item.get_tag(CBtag).strip()
        UB=item.get_tag(UBtag).strip()
        Name=CB+"_"+UB
        pos_index=pos-1

        seq_cut, pos_cut = handle_cigar(item.cigar)
        cut_seq=handle_seq(item.seq, seq_cut)
        cut_pos=handle_pos(item.get_reference_positions(), pos_cut)
        raw_index = handle_quality_matrix(cut_pos.index(pos_index),item.seq,cut_seq)
        quality=item.get_forward_qualities()[raw_index]

        if pos_index in cut_pos:
            # effective_DP += 1
            geno = cut_seq[cut_pos.index(pos_index)]            
    except:
        pass

    return Name,geno,quality
        
# util
def check_dir(dir):
    if os.path.exists(dir):
        pass
    else:
        os.mkdir(dir)

# util
def handle_posname(pos_name):
    sitem=pos_name.split("_")
    chrom=str(sitem[0])
    pos=int(sitem[1])
    ref=str(sitem[2])
    alt=str(sitem[3])

    return chrom, pos, ref, alt


def handle_pos_bed(bed_file):
    '''
    input:
    chr1    14623   C       A,<*>
    chr1    14653   C       T,<*>
    
    output:
    list: [("chrX", "119811135", "T", "C,<*>"),...]
    '''
    mutation_identifier_list=[]
    f=open(bed_file,"r")
    for line in f.readlines():
        if line[0]!="#":
            s=line.strip().split()
            chrom=s[0];pos=s[1];ref=s[2];alt=s[3]
            mutation="_".join([str(chrom),str(pos),str(ref),str(alt)])
            mutation_identifier_list.append(mutation)
    
    return mutation_identifier_list


# def read_ind_genotype(result_file):
#     '''
#     input:
#     chr1    14623   C       A,<*>
#     chr1    14653   C       T,<*>
    
#     output:
#     list: [("chrX", "119811135", "T", "C,<*>"),...]
#     '''
#     mutation_identifier_list=[]
#     f=open(result_file,"r")
#     for line in f.readlines():
#         if line[0]!="#":
#             s=line.strip().split()
#             chrom=s[0];pos=s[1];germline=s[15];mutation=s[16]
#             mutation="_".join([str(chrom),str(pos),str(germline),str(mutation)])
#             mutation_identifier_list.append(mutation)
    
#     return mutation_identifier_list


def filter_geno_dict(count_result, scale_ratio=5):
    '''
    example of count_result: {'C,C': 25548, 'G,T': 22144, 'G,C': 178, 'C,T': 55}
    '''
    hSNP_dict=defaultdict(int)
    for k,v in count_result.items():
        new_k=k[-1]
        hSNP_dict[new_k]+=v

    hSNP_rank_list=sorted(hSNP_dict.items(), key = lambda item:item[1], reverse=True)
    # print(dict(hSNP_rank_list))
    h1,h2,h3,h4=0,0,0,0
    try:
        h1=int(hSNP_rank_list[0][1]); h1_g=hSNP_rank_list[0][0]
    except:
        h1_g="NONE"

    try:    
        h2=int(hSNP_rank_list[1][1]); h2_g=hSNP_rank_list[1][0]
    except:
        h2_g="H2none"

    try:
        h3=int(hSNP_rank_list[2][1]); h3_g=hSNP_rank_list[2][0]
    except:
        h3_g=""

    try: 
        h4=int(hSNP_rank_list[3][1]); h4_g=hSNP_rank_list[3][0]
    except:
        h4_g=""
    
    dp=h1+h2+h3+h4
    if dp==0:
        return h1_g,h2_g,{}

    binom_test_p_h1=scipy.stats.binomtest(ceil(100*h1/dp), 100,p=0.5,alternative="two-sided").pvalue
    binom_test_p_h2=scipy.stats.binomtest(ceil(100*h2/dp), 100,p=0.5,alternative="two-sided").pvalue
    # binom_test_result_h3_h2=scipy.stats.binom_test(100*(h1+)/(dp), 100*h2/(h1+h2),p=0.5,alternative="two-sided") 
    # if scale_ratio*(h2+h3+h4) > h1 and h1 > 1/scale_ratio*(h2+h3+h4) and h2 > scale_ratio*(h3+h4): # and sorted([h1_g,h2_g])==sorted([ref,alt])):
    if binom_test_p_h1 >0.05 and  binom_test_p_h2 >0.05:
        ## we hope the count of h1, h2 is similar to 0.5
        short_count_result={k:v for k,v in count_result.items() if k[-1] in [h1_g, h2_g]}
        geno_rank_list=sorted(short_count_result.items(), key = lambda item:item[1], reverse=True)
        geno_count_dict=dict(geno_rank_list)

    else:
        geno_count_dict={}
        # print(f"The hSNP is not like a germline site.\nh1:{h1_g}, h2:{h2_g}\nref:{ref},alt:{alt}")
        pass
    
    # print(h1_g,h2_g,geno_count_dict)
    return h1_g,h2_g,geno_count_dict


def calculate_phased_haplo(geno_count_dict, germline, mutant, h1,h2):
    '''
    Input:
    germline is the ind genotype of pos_initial;
    ref, alt is genotype of the phased candidate allele;
    mutant is the mutated allele

    Output:
    haplo: the kinds of phased haplotype;
    annotated_type: based on the haplotype, the annotation mutation
    '''
    germline=germline.split(",")
    germline_list=list(set(germline))
    haplo=""
    # phased=""
    annotated_type=""
    scale_ratio=5
    # geno_count_dict={'C,C': 25548, 'G,T': 22144, 'G,C': 178, 'C,T': 55}
    # if  len(ref)!=1 or len(alt)!=1 or len(mutant.split(","))!=1:
    #     print(f"The length of candiadate phased sites is not equal to 1, pease check it! Input: phased_ref:{ref}; phased_alt:{alt}; mutatant:{mutant}")
    #     sys.exit()

    if len(germline_list)==1:
        ref_h1,ref_h2,alt_h1,alt_h2=0,0,0,0
        try:
            ref_h1=geno_count_dict[germline_list[0]+","+h1]
        except:
            pass
        try:
            ref_h2=geno_count_dict[germline_list[0]+","+h2]
        except:
            pass
        try:
            alt_h1=geno_count_dict[mutant+","+h1]
        except:
            pass
        try:
            alt_h2=geno_count_dict[mutant+","+h2]
        except:
            pass
        # print(germline_list[0],mutant,ref,alt)
        # print(ref_h1,ref_h2,alt_h1,alt_h2)
        # print(geno_count_dict)
        dp=ref_h1+ref_h2+alt_h1+alt_h2
        mut_allele=alt_h1+alt_h2
        ref_alt_phased=(ref_h1 > scale_ratio*ref_h2 and alt_h2 > scale_ratio*alt_h1) or (ref_h2 > scale_ratio*ref_h1 and alt_h1 > scale_ratio*alt_h2)
        ref_differ=(ref_h1 > scale_ratio*ref_h2 or ref_h2 > scale_ratio*ref_h1)
        alt_differ=(alt_h2 > scale_ratio*alt_h1 or alt_h1> scale_ratio*alt_h2)
        # detail_count=germline_list[0]+","+ref+":"+str(ref_h1)+";"+\
        #              germline_list[0]+","+alt+":"+str(ref_h2)+";"+\
        #              mutant+","+ref+":"+str(alt_h1)+";"+\
        #              mutant+","+alt+":"+str(alt_h2)
        detail_count_list=[germline_list[0]+","+h1+":"+str(ref_h1),
                    germline_list[0]+","+h2+":"+str(ref_h2),
                    mutant+","+h1+":"+str(alt_h1),
                    mutant+","+h2+":"+str(alt_h2)]
        # if (ref_h1<alt_h2*scale_ratio and ref_h1>alt_h2 *1/scale_ratio) and ((ref_h1>ref_h2*5) and (alt_h2>alt_h1*5) or ((ref_h2>scale_ratio*ref_h1) and (alt_h1>scale_ratio*alt_h2))):
        if ref_alt_phased and dp>4:   
            haplo="haplo=2"
            # phased="yes"
            annotated_type="heterzygous"
        # elif (((ref_h1>ref_h2 * 1/scale_ratio) and (ref_h1<ref_h2 * scale_ratio)) and (alt_h1>alt_h2 *scale_ratio or alt_h2>alt_h1*scale_ratio) and (alt_h1+alt_h2>1)):
        elif (not ref_differ) and alt_differ and alt_h1+alt_h2>=2 and ((ref_h1>ref_h2 and alt_h1<alt_h2) or (ref_h2>ref_h1 and alt_h2<alt_h1)):
            haplo="haplo=3"
            # phased="yes"
            annotated_type="mut_mosaic"
        # elif not ((ref_h1<alt_h2*scale_ratio and ref_h1>alt_h2 *1/scale_ratio)) and ((ref_h1>ref_h2*5) and (alt_h2>alt_h1*5) or ((ref_h2>scale_ratio*ref_h1) and (alt_h1>scale_ratio*alt_h2))) or ((ref_h1<ref_h2/5) and (alt_h2<alt_h1/5) and (ref_h2<alt_h1*5 and ref_h2>alt_h1/5)) or (((ref_h1>ref_h2/5) and (ref_h1<ref_h2*5)) and (alt_h1<=alt_h2/5 or alt_h2<=alt_h1/5))  ):
        # elif (not ref_differ) and (not alt_differ) and ((ref_h1+ref_h2)<scale_ratio*(alt_h1+alt_h2) and ((ref_h1+ref_h2)>1/scale_ratio*(alt_h1+alt_h2))):
            # phase[name]['hap>3']=phase[name].get("hap>3",0)+1
        elif (not ref_alt_phased) and ((not ref_differ) and alt_differ):
            haplo="haplo>3"
            # phased="not"
            annotated_type="artifacts"

        else:
            haplo="haplo>3"
            # phased="not"
            annotated_type="artifacts"
            # print("Unknown:",ref_h1,ref_h2,alt_h1,alt_h2)

    elif len(germline_list)==2:
        germline_list=germline
        ref_h1,ref_h2,alt_h1,alt_h2,mut_h1,mut_h2=0,0,0,0,0,0
        try:
            ref_h1=geno_count_dict[germline_list[0]+","+h1]
        except:
            ref_h1=0
        try:
            ref_h2=geno_count_dict[germline_list[0]+","+h2]
        except:
            ref_h2=0
        try:
            alt_h1=geno_count_dict[germline_list[1]+","+h1]
        except:
            alt_h1=0
        try:
            alt_h2=geno_count_dict[germline_list[1]+","+h2]
        except:
            alt_h2=0
        try:
            mut_h1=geno_count_dict[mutant+","+h1]
        except:
            mut_h1=0
        try:
            mut_h2=geno_count_dict[mutant+","+h2]
            # dp=ref_h1+ref_h2+alt_h1+alt_h2+mut_h1+mut_h2
        except:
            mut_h2=0

        # detail_count=germline_list[0]+","+ref+":"+str(ref_h1)+";"+\
        #         germline_list[0]+","+alt+":"+str(ref_h2)+";"+\
        #         germline_list[1]+","+ref+":"+str(alt_h1)+";"+\
        #         germline_list[1]+","+alt+":"+str(alt_h2)+";"+\
        #         mutant+","+ref+":"+str(mut_h1)+";"+\
        #         mutant+","+alt+":"+str(mut_h2)
        detail_count_list=[germline_list[0]+","+h1+":"+str(ref_h1), 
                    germline_list[0]+","+h2+":"+str(ref_h2),
                    germline_list[1]+","+h1+":"+str(alt_h1),
                    germline_list[1]+","+h2+":"+str(alt_h2),
                    mutant+","+h1+":"+str(mut_h1),
                    mutant+","+h2+":"+str(mut_h2)]
        dp=ref_h1+ref_h2+alt_h1+alt_h2+mut_h1+mut_h2
        mut_allele=mut_h1+mut_h2
        ref_alt_phased=(ref_h1 > scale_ratio*ref_h2 and alt_h2 > scale_ratio*alt_h1) or (ref_h2 > scale_ratio*ref_h1 and alt_h1 > scale_ratio*alt_h2)
        ref_alt_differ=((ref_h1+ref_h2) > scale_ratio * (alt_h1+alt_h2))
        ref_differ=(ref_h1 > scale_ratio*ref_h2 or ref_h2 > scale_ratio*ref_h1)
        alt_differ=(alt_h2 > scale_ratio*alt_h1 or alt_h1> scale_ratio*alt_h2)
        mut_differ=(mut_h1 > scale_ratio*mut_h2 or mut_h2> scale_ratio*mut_h1)

        if ref_alt_phased and mut_h1 + mut_h2<2:
            # if ref alr are phased, but the amount of mutant allele is too small, this site will be seemed as a heterzygous, and mutant allele is artifacts
            #haplo="haplo=2" # this will help to tidy the heterzygous sites, but ignore in this version
            haplo="haplo>3"
            # phased="yes"
            annotated_type="artifacts"
        
        elif ref_alt_phased and mut_differ:
            haplo="haplo=3"
            phased="yes"
            annotated_type="mut_mosaic"

        elif ref_alt_phased and (not mut_differ):
            #haplo="haplo=2" # this will help to tidy the heterzygous sites, but ignore in this version
            haplo="haplo>3"
            # haplo=4
            # phased="yes"
            annotated_type="artifacts"

        elif (not ref_differ) and alt_differ and mut_differ:
            # haplo=4
            # phased="yes"
            # if background_ratio*(dp)>(mut_h2+mut_h2) and (not (ref_h2<alt_h1*5 and ref_h2>alt_h1/5)):
            #     annotated_type="alt_mosaic"
            # else:
            # annotated_type="mut_mosaic"
            haplo="haplo>3"
            annotated_type="artifacts"

        elif (not ref_differ) and alt_differ and (not mut_differ):
        #raw: elif ref_alt_phased and (not ref_differ) and alt_differ and (not mut_differ):
        # elif ref_alt_phased and ref_differ and alt_differ and (not mut_differ) and ref_alt_differ:
            # haplo="haplo=3"
            # # phased="yes"
            # annotated_type="alt_mosaic"
            haplo="haplo>3"
            annotated_type="artifacts"   
        
        elif (not ref_alt_phased) and mut_differ:
            # haplo=5
            # phased="no"
            haplo="haplo>3"
            annotated_type="artifacts"

        elif (not ref_alt_phased) and (not mut_differ):
            # haplo=6
            # phased="no"
            haplo="haplo>3"
            annotated_type="artifacts"
        
        else:
            haplo="haplo>3"
            annotated_type="artifacts"
            # print(ref_h1,ref_h2,alt_h1,alt_h2,mut_h1,mut_h2)

    return dp,mut_allele,haplo,annotated_type,detail_count_list


def judge_origin(germline, annotated_type,detail_count_list):
    def get_count(detail_count_list,index):
        count=int(detail_count_list[index].split(":")[-1])
        return count

    mut_origin="NA"
    if len(germline.split(","))==2 and annotated_type=="alt_mosaic" and len(detail_count_list)==6:
        ref=germline
        mut_origin=ref

    elif len(germline.split(","))==2 and annotated_type=="mut_mosaic" and len(detail_count_list)==6:
        ref=germline.split(",")[0]
        alt=germline.split(",")[1]
        mut_h1,mut_h2=get_count(detail_count_list,4),get_count(detail_count_list,5)
        if mut_h1>mut_h2:
            mut_origin=ref if get_count(detail_count_list,0)>get_count(detail_count_list,2) else alt
        elif mut_h1<mut_h2:
            mut_origin=ref if get_count(detail_count_list,1)>get_count(detail_count_list,3) else alt
        else:
            pass
    elif len(germline.split(","))==1 and len(detail_count_list)==4:
        mut_origin=germline

    return mut_origin

def get_candidate_germline(ref_fasta,short_ind_germ_file, in_bam_name,flanking,min_prior,gender,gene_name_index,run_type,bins,line):
    '''
    line_example:
    chr1    601435  724550  RP5-857K21.4    -       ENSG00000230021.3       lincRNA
    
    short_ind_germ_file:
    chr12	52487212	52487212    het C   T   0.1,0.5
    chr12   52487212	52487212    mosaic C   T   0.1,0.5
    '''
    print("====",line)

    out_list=[]
    sline=line.strip().split("\t")
    chr=str(sline[0])
    if gender=="male":
        ignore_chr_list=["chrM", "MT","chrY","chrX"]
    else:
        ignore_chr_list=["chrM", "MT","chrY"]
    
    if chr in ignore_chr_list:
        return []
    
    pos_s=int(sline[1])
    pos_e=int(sline[2])
    gene_name=sline[gene_name_index].replace("\"","")
    chrom_size=int(get_chr_size(ref_fasta+".fai")[chr])
    add_flanking_line="\t".join([chr,str(max(0,pos_s-flanking)),str(min(pos_e+flanking,chrom_size))])
    gene_bed=BedTool(add_flanking_line,from_string=True)
    germ_mosaic_bed=BedTool(short_ind_germ_file)
    with warnings.catch_warnings():
        pos_intersect = germ_mosaic_bed.intersect(gene_bed,u=True)
        warnings.simplefilter("ignore")
 
    
    pos_candidate_dict={"mosaic_pos":[],"informative_SNP":[]}
    count_allele=[]
    for line in pos_intersect:
        new_sline=line.fields
        # print(new_sline)
        mut_name="_".join([new_sline[0],new_sline[1],new_sline[4],new_sline[5]])
        if new_sline[3]=="het":
            prior=new_sline[6]
            if float(prior.split(",")[1])<min_prior:
                continue
            pos_candidate_dict["informative_SNP"].append(mut_name)
            count_allele.append((new_sline[1],new_sline[4][0]))
        elif new_sline[3]=="mosaic":
            pos_candidate_dict["mosaic_pos"].append(mut_name)
            count_allele.append((new_sline[1],new_sline[4][0]))
        else:
            # print(sline)
            print(f"something wrong? please check your input: {short_ind_germ_file}")
    # print(pos_candidate_dict)
    if pos_candidate_dict["mosaic_pos"]==[] or pos_candidate_dict["informative_SNP"]==[]:
        print("no mosaic or no germ")
        return []
    
    print(gene_name,"start")
    # candidate_allele_info=pos_candidate_dict["mosaic_pos"]+pos_candidate_dict["informative_SNP"]
    candidate_allele_info=list(set(count_allele))
    print(len(candidate_allele_info))
    ## offer a chrom pos initial and its candidate allele info
    per_read_dict=defaultdict(dict)
    allele_total_count=defaultdict(int)
    used_barcode=defaultdict(list)

    start=max(0,pos_s-flanking)
    end=min(pos_e+flanking,chrom_size)
    # chunk_size = 1000
    # print(range(start, end, chunk_size))
    # i=0
    # for chunk_start in range(start, end, chunk_size):
    try:
        # in_bam_read=pysam.AlignmentFile(in_bam_name, "rb", reference_filename=ref_fasta)
        with pysam.AlignmentFile(in_bam_name, "rb", reference_filename=ref_fasta) as in_bam_read:
            for reads in in_bam_read.fetch(chr, start, end):
                try:
                    if run_type=="visium":
                        try:
                            CB=reads.get_tag("CB").strip()
                            UB=reads.get_tag("UB").strip()

                        except:
                            continue
                    elif run_type=="stereo":
                        try:
                            Cx_raw=int(reads.get_tag("Cx"))
                            Cy_raw=int(reads.get_tag("Cy"))
                            if bins !=1:
                                Cx=round_to_nearest_bin(Cx_raw,bins)
                                Cy=round_to_nearest_bin(Cy_raw,bins)
                            else:
                                Cx=Cx_raw
                                Cy=Cy_raw
                            UR=reads.get_tag("UR").strip()

                            CB=str(Cx)+"_"+str(Cy)
                            UB=str(UR)
                            # UMI_name=str(Cx_raw)+"_"+str(Cy_raw)+"_"+str(UR)
                        except:
                            continue

                    elif run_type=="ST":
                        try:
                            CB=str(reads.get_tag("B0"))
                            UB=str(reads.get_tag("B3"))
                            # barcode_name=str(CB)
                            # UMI_name=str(UB)
                        except:
                            continue
                    else:
                        # print("type",run_type)
                        continue

                    seq_cut, pos_cut = handle_cigar(reads.cigar)
                    cut_seq=handle_seq(reads.seq, seq_cut)
                    cut_pos=handle_pos(reads.get_reference_positions(), pos_cut)

                    # UMI_name = str(UB)
                    barcode_name = str(CB)+"_"+str(UB)
                except:
                    continue
                # print("read is done: ",barcode_name)   
                for candidate_allele_tuple in candidate_allele_info:
                    candidate_allele=candidate_allele_tuple[0]
                    ref=candidate_allele_tuple[1]

                        # print("Old")

                    pos_index=int(candidate_allele)-1
                    geno=""
                    if pos_index in cut_pos:
                        # effective_DP += 1
                        raw_index = handle_quality_matrix(cut_pos.index(pos_index),reads.seq,cut_seq)
                        quality=reads.get_forward_qualities()[raw_index]
                        geno = cut_seq[cut_pos.index(pos_index)]

                        if geno not in "ATCG":
                            # per_read_dict[barcode_name][candidate_allele]["count"]={}
                            # per_read_dict[barcode_name][candidate_allele]["quality"]={}
                            continue
                        # print(barcode_name,candidate_allele)
                        if candidate_allele not in per_read_dict[barcode_name].keys():
                            # print("New")
                            per_read_dict[barcode_name][candidate_allele]={}
                            per_read_dict[barcode_name][candidate_allele]["count"]=defaultdict(int)
                            per_read_dict[barcode_name][candidate_allele]["quality"]={"A":defaultdict(int), \
                                                                                    "T":defaultdict(int), \
                                                                                    "C":defaultdict(int), \
                                                                                    "G":defaultdict(int)}
                        else:
                            pass
                        per_read_dict[barcode_name][candidate_allele]["count"][geno]+=1
                        per_read_dict[barcode_name][candidate_allele]["quality"][geno][quality]+=1

                        if barcode_name not in used_barcode[candidate_allele]:
                            allele_total_count[candidate_allele]+=1
                            used_barcode[candidate_allele].append(barcode_name)
            # if i==1:
            #     print(per_read_dict)
            del reads
                # print("allele is done: ",barcode_name)   
    except Exception as e:
        print(f"Error in fetch: {e}")
    finally:
        print(gene_name," is done!")  # 确保最终执行        
        # print(i," in ",gene_name," is done!")
            #consensus_read_count, consensus_read_quality=UMI_combination_spot(site_barcode_UMI_dict,chrom,pos,ref)
        # calculater per read
    print("yes1")    
    per_read_genotypes_count=defaultdict(list)
    per_read_genotypes_quality=defaultdict(list)
        # print("stat:",per_read_dict)

    for barcode in per_read_dict.keys():
        # print(barcode)
        # for candidate_allele in per_read_dict[barcode]:
        for candidate_allele,ref in candidate_allele_info:
            if candidate_allele not in per_read_dict[barcode].keys() :
                phred_dict={}
                most_proved_allele,phred=".","."
            else:
                # ref=dict(candidate_allele_info)[candidate_allele][0]
                count_dict=per_read_dict[barcode][candidate_allele]["count"]
                quality_dict=per_read_dict[barcode][candidate_allele]["quality"]
                if count_dict!={}:
                    # print("yes")
                    phred_dict=calculate_UMI_combine_phred(count_dict,quality_dict,weigh=0.5)
                    # print(phred_dict)
                    most_proved_allele,phred=get_most_candidate_allele(phred_dict,ref)
                    # print(most_proved_allele,phred)
                else:
                    # print("no")
                    phred_dict={}
                    most_proved_allele,phred=".","."
            
            per_read_genotypes_count[barcode].append(most_proved_allele)
            per_read_genotypes_quality[barcode].append(phred)
    del per_read_dict
    # print(candidate_allele_info)
    # print(per_read_genotypes_count.values())
    # print(candidate_allele_info[1:])
    # print("###",candidate_allele_info[0])
    # print(pos_candidate_dict)
    # print("yes2")
    candidate_allele_info_only_pos=[s[0] for s in candidate_allele_info]
    for mosaic_pos in pos_candidate_dict["mosaic_pos"]:
        chrom, pos, germline, mutant = handle_posname(mosaic_pos)
    
        for info_SNP in pos_candidate_dict["informative_SNP"]:
            # print("#",mut)
            chrom_germ, pos_germ, _, _=handle_posname(info_SNP)
            if pos_germ == pos and chrom_germ==chrom:
                continue
            
            short_list=[]
            # print("1",candidate_allele_info_only_pos)
            index_mosaic=candidate_allele_info_only_pos.index(str(pos))
            index_germ=candidate_allele_info_only_pos.index(str(pos_germ))
            # print("2",index_mosaic,index_germ)
            # print("3",pos,pos_germ)

            # print(pos,index_mosaic,pos_germ,index_germ)
            # print(mosaic_pos,index_mosaic,info_SNP,index_germ)
            # print(per_read_genotypes_count.values())
            short_dict={}
            # print("====",gene_name,mosaic_pos,info_SNP)
            for barcode in per_read_genotypes_count.keys():
                values=per_read_genotypes_count[barcode]
                # print()
                bases=",".join([values[index_mosaic],values[index_germ]])
                short_list.append(bases)
                # if "." not in bases and bases[0]==mutant:
                if bases[0]==mutant:
                    short_dict[barcode]=bases
                    # print("-",barcode,bases,per_read_genotypes_count[barcode])
                    # print("-",barcode,bases,per_read_genotypes_quality[barcode])

            # for values in per_read_genotypes_count.values():
            #     # print(index_mut,pos,ref,alt,values[index_mut])
            #     bases=",".join([values[index_mosaic],values[index_germ]])
            #     short_list.append(bases)
            #     # print(bases)

            count_result=dict(Counter([bases for bases in short_list if "." not in bases]))
            h1,h2,geno_count_dict=filter_geno_dict(count_result)
            # print(mosaic_pos,info_SNP,count_result)
            # print("filter1",mosaic_pos,info_SNP,h1,h2)
            if h2=="H2none":
                continue
            
            # print(count_result)
            # print("filter2",mosaic_pos,info_SNP,h1,h2)   
            if geno_count_dict!={}:
                new_mut_name="_".join([str(chrom_germ), str(pos_germ),h1,h2])
                total_count=allele_total_count[str(pos)]
                dp,mut_allele,haplo,annotated_type,detail_count_list=calculate_phased_haplo(geno_count_dict, germline, mutant, h1, h2)
                if mut_allele==0 or mut_allele == "0":
                    continue

                if annotated_type != "artifacts" and annotated_type !="heterzygous":
                    mut_origin=judge_origin(germline, annotated_type,detail_count_list)
                else:
                    mut_origin="NA"
                detail_count=";".join(detail_count_list)
                # print("filter3",mosaic_pos,info_SNP,h1,h2)
                out_list.append([chrom, str(pos), germline, mutant,mut_origin, new_mut_name, gene_name,total_count, dp,mut_allele,haplo,annotated_type,detail_count])
                # print([chrom, str(pos), germline, mutant,new_mut_name, gene_name, dp,mut_allele,haplo,annotated_type,detail_count])
                # print("==========================")
                # print("\n",geno_count_dict)
                # print(candidate_allele_info[0],mut,haplo,phased,annotated_type,"\n", geno_count_dict)
            else:
                pass
    del pos_candidate_dict
    print(gene_name,"finished")
    return out_list
                

def main():
    print_memory_usage("程序启动")
    ref_fasta=args.fasta
    tmp_dir=os.path.join(str(args.outdir),"tmp3")
    check_dir(tmp_dir)

    print_memory_usage("程序启动")
    germ_site_bed_file=handel_candidate_informative_SNP_site(args.species,args.germline, tmp_dir)
    print_memory_usage("生成germ_site_bed_file后")
    short_ind_file, combine_ind_germ_file=short_mosaic_sites_and_combine_with_germ(args.indgeno, germ_site_bed_file, tmp_dir)
    print_memory_usage("生成short_ind_file和combine_ind_germ_file后")

    if args.gtexGene:
        #short_gene_review_file=short_bed_file(args.genebed, short_ind_file, tmp_dir)
        short_gencode_file_name=short_bed_file(args.gtexGene, short_ind_file, tmp_dir)
        gene_name_index=3
    if args.gencode:
        short_gencode_file_name=short_gencode_file(args.gencode,short_ind_file,tmp_dir,rerun=True)
        gene_name_index=4


    print("finish preparing work")
    
    # out_file=open(os.path.join(args.outdir,args.sample+"."+args.outsuffix),"w")
    if args.gender in ["M","male"]:
        gender="male"
    else:
        gender="female"

    partial_func=partial(get_candidate_germline,ref_fasta,combine_ind_germ_file, args.bam,args.flanking,args.minprior,gender,gene_name_index,args.seq_type,args.bin)
    thread=int(args.thread)

    short_gene_review=open(short_gencode_file_name,"r")
    
    print_memory_usage("读取gene review文件后")
    lines=short_gene_review.readlines()
    # # print(lines)
    # num_chunks=100
    # chunk_size = len(lines) // num_chunks
    # for i in range(num_chunks):
    #     chunk_start = i * chunk_size
    #     chunk_end = chunk_start + chunk_size if i < num_chunks - 1 else len(lines)
    #     chunk = lines[chunk_start:chunk_end]
    #     # print(lines)

    # with multiprocessing.Pool(thread) as pool:
    #     # results=pool.map(partial_func,chunksize=100) # type: ignore # per result is [[site_info]]
    #     results=pool.map(partial_func, lines,chunksize=1) # type: ignore # per result is [[site_info]]


    # out_file=open(os.path.join(args.outdir,args.sample+"."+args.outsuffix),"w")
    # # print(f"finish running chunk{i}")
    # for lists in results:
    #     for line in lists:
    #         out_text="\t".join([str(i) for i in line])
    #         out_file.write(f'{out_text}\n')

    # out_file.close()

    out_file=os.path.join(args.outdir,args.sample+"."+args.outsuffix)
    COLUMNS = ['#chr', 'pos', 'germline', 'mutant','mut_origin', 'new_mut_name', 'gene_name','total_count', 'dp','mut_allele','haplo','annotated_type','detail_count']
    # for line in lines:
    #     result=partial_func(line)
    #     print(result)
    with multiprocessing.Pool(thread) as pool, open(out_file, "w") as f:
        print_memory_usage("进入多进程前")
        pd.DataFrame(columns=COLUMNS).to_csv(f, header=True, index=False,sep='\t')
        
        for result in pool.imap(partial_func, lines,chunksize=1):
            print_memory_usage("每个多进程结果返回后")
            if result: 
                for per_result in result: 
                    chunk_df = pd.DataFrame([per_result], columns=COLUMNS)
                    chunk_df.to_csv(f, header=False, mode="a", index=False,sep='\t')
            del result
            gc.collect()
    # 在不再需要临时文件时删除
    if os.path.exists(germ_site_bed_file):
        os.remove(germ_site_bed_file)

    print_memory_usage("多进程结束后")
        
## parameters
parser = argparse.ArgumentParser()
parser.add_argument("--fasta","-f", required=True,help="fasta file")
parser.add_argument("--bam","-b", required=False, help="bam file")
parser.add_argument("--sample","-s",required=True, help="sample name")
parser.add_argument("--germline",required=True, help="germline file")
# parser.add_argument("--CBtag",required=False,default="CB", help="CBtag")
# parser.add_argument("--UBtag",required=False,default="UB", help="UBtag")
parser.add_argument("--seq_type", dest='seq_type',default="visium",choices=["visium","stereo","ST"],type=str, required=False, help="Your input sequence type")
parser.add_argument("--bin", required=False,type=int,default=100,help="bin value")
parser.add_argument("--gtexGene",required=False, help="gene bed file downloaded from uscs goldpath")
parser.add_argument("--gencode",required=False, help="gtf file downloaded from gencode")
parser.add_argument("--flanking",required=False,type=int,default=120, help="read length or the flanking length")
parser.add_argument("--minprior",required=False,type=float,default=0.01, help="the min prior for filtering informative SNP")
parser.add_argument("--indgeno",required=True,help="individual genotype file")
parser.add_argument("--gender",required=False,default="female",choices=["F","M","female","male"],help="The gender of input sample, if male, chrX would not be considered")
parser.add_argument("--outdir","-o",required=True,help="output dir")
parser.add_argument("--outsuffix",required=False,default="after_phasing_candidate_list.txt",type=str,help="output suffix")
parser.add_argument("--thread","-t",required=False,default=2, type=int,help="thread")
parser.add_argument("--species",required=False,default="human", type=str,help="human species would filter hSNP, while other not")
args = parser.parse_args()
    

if __name__ == '__main__':
    main()

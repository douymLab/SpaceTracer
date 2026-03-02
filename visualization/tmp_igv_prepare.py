import os
import pandas as pd
import argparse
import pysam
import numpy as np

def read_mutation_list(in_file):
    pos_list=list()
    f =open(in_file,"r")
    for line in f.readlines():
        s = line.strip().split()
        chr=s[0]; pos1=str(s[1]); ref=str(s[2]); alt=str(s[3])
        identifier="_".join([chr,pos1,ref,alt])
        if identifier not in pos_list:
            pos_list.append(identifier)
    return pos_list


def read_mutation_barcode(file):
    barcode_mutinfo_df_colnames=["barcode_name", "cluster","in_tissue","pos_x","pos_y","depth","vaf","mutation_prob","mosaic_likelihood", \
                                    "mutant_allele_num","has_mutant_allele","high_prob","high_likelihood","mutated"]
    barcode_mutinfo_df=pd.read_csv(file, sep='\t', header=None, names=barcode_mutinfo_df_colnames, comment = "#")
    return barcode_mutinfo_df


def check_dir(dir_name):
    if not os.path.exists(dir_name):
        os.mkdir(dir_name)
    else:
        #print("dir exist")
        pass

def read_barcode_list(file):
    # barcode_mutinfo_df_colnames=["barcode_name", "mutant_allele_num"]
    barcode_mutinfo_df_colnames=["barcode_name", "raw_allele_num","raw_mutant_allele_num","allele_num","mutant_allele_num"]
    barcode_mutinfo_df=pd.read_csv(file, sep='\t', header=None, names=barcode_mutinfo_df_colnames, comment = "#")
    return barcode_mutinfo_df


def read_mutation_list2(mutation_list_file):
    file=open(mutation_list_file,"r")
    pos_list=list()
    for line in file.readlines():
        s = line.strip().split()[0]
        pos_list.append(s)
    return pos_list


def substract_bams(identifier,outlist,bam,output_bam_file,flanking=50):
    bam_handle=pysam.AlignmentFile(bam,"r", ignore_truncation=True)
    output_bam=pysam.AlignmentFile(output_bam_file, "wb", header=bam_handle.header)

    chrom,pos,ref,alt=identifier.split("_")
    reads=bam_handle.fetch(chrom,int(pos)-1-flanking,int(pos)+flanking)
    # pos_index = int(pos)-1

    # def round_to_nearest_bins(x):
    #     return int(np.ceil(x / float(bins)) * float(bins))
    
    for item in reads:
        try:
            # Cx=str(round_to_nearest_bins(item.get_tag("Cx")))
            # Cy=str(round_to_nearest_bins(item.get_tag("Cy")))
            barcode_name=item.get_tag("CB").strip()
        except:
            continue
    
        if barcode_name in outlist:
            output_bam.write(item)
    output_bam.close()        
    pysam.index(output_bam_file)
    # pysam.IndexedReads(output_bam)


def detail_bam(identifier,mut_df,cell_info,input_bam,outdir,mut_type):
    out_barcode_name=os.path.join(outdir,"barcode_info", args.sample + "_" + identifier + str(mut_type) + ".barcode.list")
    out_barcode_file=open(out_barcode_name,"w")
    if mut_type==0:
        out_list=list(mut_df[mut_df["mutant_allele_num"]==0]["barcode_name"])
    elif mut_type==1:
        out_list=list(mut_df[mut_df["mutant_allele_num"]>=1]["barcode_name"])
    else:
        return False

    anno_list = list(set([cell_info.loc[bam]["anno"] for bam in out_list if bam in cell_info.index]))
    # print(anno_list)
    for anno in anno_list:
        anno_barcodes = cell_info[cell_info["anno"] == anno].index.tolist()
        current_out = [bam for bam in out_list if bam in anno_barcodes]
        # print(current_out)
        if not current_out:
            continue

        output_bam_file =os.path.join(outdir,"short_bam",f"{anno}_{identifier}_{str(mut_type)}.bam")
        substract_bams(
            identifier=identifier,
            outlist=current_out,
            bam=input_bam,
            output_bam_file=output_bam_file
        )
        
        out_barcode_file.write(f'{os.path.abspath(output_bam_file)}\n')
    out_barcode_file.close()

    return True

def main():
    # prob_df=read_prob_file(args.prob)
    # threshold=float(args.threshold)
    check_dir(args.outdir)
    print("dir check!")

    # if args.mutation:
    #     mutation_list=read_mutation_list(args.mutation)
    # elif args.mutation_list:
    #     mutation_list=read_mutation_list2(args.mutation_list)
    # else:
    #     mutation_list=[]
    mutation=args.mutation
    # bin_info=pd.read_csv(args.cellinfo,names=["binx","biny","anno"],sep="\t")
    # bin_info=pd.read_csv(args.cellinfo,names=["binx","biny","anno"],sep="\t")
    cell_info=pd.read_csv(args.cellinfo,names=["barcode_name","anno"],sep="\t")
    # bin_info["barcode_name"]= bin_info['binx'].astype(str) + "_" + bin_info['biny'].astype(str)
    cell_info.set_index("barcode_name", inplace=True)
    # print(mutation)
    smu=mutation.split("_")
    chrom=smu[0]; pos=smu[1]
    # try:
    if os.path.exists(os.path.join(args.barcodedir,args.sample+"."+mutation+".barcode.mutinfo.txt")):
        mut_file=os.path.join(args.barcodedir,args.sample+"."+mutation+".barcode.mutinfo.txt")
        mut_df=read_mutation_barcode(mut_file)
    elif os.path.exists(os.path.join(args.barcodedir,mutation+".mut.spots.txt")):
        mut_file=os.path.join(args.barcodedir,mutation+".mut.spots.txt")
        mut_df=read_barcode_list(mut_file)
    else:
        print("Opps! cannot find the mutation file!")
        pass
    print(mutation,"done")

    check_dir(os.path.join(args.outdir,"barcode_info"))
    check_dir(os.path.join(args.outdir,"short_bam"))

    if not detail_bam(mutation,mut_df,cell_info,args.bam,args.outdir,0):
        print("not finished 0!")
    if not detail_bam(mutation,mut_df,cell_info,args.bam,args.outdir,1):
        print("not finished 1!")


    check_dir(os.path.join(args.outdir,"bed_info"))     
    bed_name=os.path.join(args.outdir,"bed_info",args.sample + "_" + mutation +".bed")
    bed_file=open(bed_name,"w")
    bed_file.write(f'{chrom}\t{pos}\t{pos}')
    bed_file.close()


## parameters
parser = argparse.ArgumentParser()
parser.add_argument("--barcodedir", required=True,help="barcodedir")
parser.add_argument("--cellinfo", required=True,help="bin annotation file")
parser.add_argument("--sample","-s",required=True, help="sample name")
parser.add_argument("--mutation","-m",required=False, help="one mutation")
#parser.add_argument("--mutation_list",required=False, help="identifier list: chr_pos_ref_alt")
parser.add_argument("--outdir","-o",required=False, default="IGV_Plot", help="output dir")
parser.add_argument("--bam","-b",required=True, type=str, help="The bam file  to grep info")

args = parser.parse_args()


if __name__ == '__main__':
    main()

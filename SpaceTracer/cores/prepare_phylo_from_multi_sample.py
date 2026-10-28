import pandas as pd
import os
from collections import defaultdict
import argparse
import pysam

from SpaceTracer.utils.handle_UMI_combine import (
    calculate_UMI_combine_phred,
    get_most_candidate_allele,
    handle_cigar,
    handle_pos,
    handle_quality_matrix,
    handle_seq,
    handle_seq_type,
)

from SpaceTracer.utils.read_files import handle_barcode
from SpaceTracer.utils.utils import barcode_cell_mapping


def handle_reads_per_pos_allele_count(
    reads,
    pos,
    seq_type="visium",
    bins=50,
    barcode_list=None,
    cell_dict={}
):
    pos_index = pos - 1
    site_barcode_UMI_dict = defaultdict(dict)

    use_filter = barcode_list is not None and len(barcode_list) > 0
    barcode_set = set(barcode_list) if use_filter else set()
    
    for item in reads:
        barcode_name, UMI_name = handle_seq_type(item, seq_type, bins,cell_dict)
        if barcode_name is None or UMI_name is None:
            continue

        if use_filter and barcode_name not in barcode_set:
            continue

        ref_positions = item.get_reference_positions()
        if pos_index not in ref_positions:
            continue

        seq_cut, pos_cut = handle_cigar(item.cigar)
        cut_seq = handle_seq(item.seq, seq_cut)
        cut_pos = handle_pos(ref_positions, pos_cut)

        if pos_index not in cut_pos:
            continue

        cut_idx = cut_pos.index(pos_index)
        geno = cut_seq[cut_idx]
        if geno not in "ATCG":
            continue

        raw_index = handle_quality_matrix(cut_idx, item.seq, cut_seq)
        try:
            quality = item.get_forward_qualities()[raw_index]
        except Exception:
            continue

        if UMI_name not in site_barcode_UMI_dict[barcode_name]:
            site_barcode_UMI_dict[barcode_name][UMI_name] = {
                "count": [0, 0, 0, 0],
                "quality": {
                    "A": defaultdict(int),
                    "T": defaultdict(int),
                    "C": defaultdict(int),
                    "G": defaultdict(int)
                }
            }

        umi_entry = site_barcode_UMI_dict[barcode_name][UMI_name]
        geno_index = "ATCG".index(geno)
        umi_entry["count"][geno_index] += 1
        umi_entry["quality"][geno][quality] += 1
    return site_barcode_UMI_dict


def summarize_barcode_support(site_barcode_UMI_dict, ref, alt, run_type="UMI"):
    allele_index_dict = {"A": 0, "T": 1, "C": 2, "G": 3}
    count_info_dict = {}
    alt_info_dict = defaultdict(int)
    total_count = 0
    alt_count = 0

    for barcode_name, umi_dict in site_barcode_UMI_dict.items():
        per_barcode_list = [0, 0, 0, 0]

        if run_type == "UMI":
            for UMI_name, umi_info in umi_dict.items():
                count_list = umi_info["count"]
                quality_dict = umi_info["quality"]

                phred_dict = calculate_UMI_combine_phred(count_list, quality_dict, weigh=0.5)
                candidate_allele, phred = get_most_candidate_allele(phred_dict, ref)
                per_barcode_list[allele_index_dict[candidate_allele]] += 1

        elif run_type == "read":
            for UMI_name, umi_info in umi_dict.items():
                count_list = umi_info["count"]
                for i in range(4):
                    per_barcode_list[i] += count_list[i]

        else:
            raise ValueError(f"Unsupported run_type: {run_type}")

        alt_support = per_barcode_list[allele_index_dict[alt]]
        total_support = sum(per_barcode_list)

        count_info_dict[barcode_name] = f"{alt_support}/{total_support}"
        alt_info_dict[barcode_name] = alt_support
        alt_count += alt_support
        total_count += total_support

    return count_info_dict, alt_info_dict, alt_count, total_count


def get_identifier_info_from_bam(
    identifier,
    bam_file,
    barcode_list,
    run_type="UMI",
    seq_type="visium",
    bins=50,
    cell_dict={}
):
    chrom, pos, ref, alt = identifier.split("_")
    pos = int(pos)
    with pysam.AlignmentFile(bam_file, "r") as bam_handle:
        reads = bam_handle.fetch(chrom, pos - 1, pos)

        site_barcode_UMI_dict = handle_reads_per_pos_allele_count(
            reads=reads,
            pos=pos,
            seq_type=seq_type,
            bins=bins,
            barcode_list=barcode_list,
            cell_dict=cell_dict
        )

    count_info_dict, alt_info_dict, alt_count, total_count = summarize_barcode_support(
        site_barcode_UMI_dict=site_barcode_UMI_dict,
        ref=ref,
        alt=alt,
        run_type=run_type
    )

    count_values = [
        count_info_dict[key] if key in count_info_dict else "0/0"
        for key in barcode_list
    ]

    posteria_values = [
        "1" if key in alt_info_dict and alt_info_dict[key] > 0
        else "0" if key in alt_info_dict and alt_info_dict[key] == 0
        else "NA"
        for key in barcode_list
    ]

    mut_likelihood_values = [
        "1" if key in alt_info_dict and alt_info_dict[key] > 0
        else "0" if key in alt_info_dict and alt_info_dict[key] == 0
        else "NA"
        for key in barcode_list
    ]

    nomut_likelihood_values = [
        "0" if key in alt_info_dict and alt_info_dict[key] > 0
        else "1" if key in alt_info_dict and alt_info_dict[key] == 0
        else "NA"
        for key in barcode_list
    ]

    return (
        [chrom, str(pos), f'"{ref}"', alt, "1"],
        [alt_count, total_count],
        count_values,
        posteria_values,
        mut_likelihood_values,
        nomut_likelihood_values
    )


def handle_posname(pos_name):
    sitem = pos_name.split("_")
    chrom = str(sitem[0])
    pos = int(sitem[1])
    ref = str(sitem[2])
    alt = str(sitem[3])
    return chrom, pos, ref, alt


def read_mutation_file(mutation_file):
    mutation_identifier_list = []
    with open(mutation_file, "r") as f:
        for line in f:
            sline = line.strip().split()
            if not sline:
                continue
            if len(sline) == 1:
                mutation_identifier_list.append(sline[0])
            else:
                chrom, pos, ref, alt = sline[0], sline[1], sline[2], sline[3]
                mutation_identifier_list.append("_".join([chrom, pos, ref, alt]))
    return mutation_identifier_list


def tidy_func(sample, barcode_list, bam_file, run_type, seq_type, bins,cell_dict, mutation_identifier):
    ind_info, counts, count_values, posteria_values, mut_likelihood_values, nomut_likelihood_values = get_identifier_info_from_bam(
        mutation_identifier, bam_file, barcode_list, run_type, seq_type, bins,cell_dict
    )
    return sample, ind_info, posteria_values, mut_likelihood_values, nomut_likelihood_values, counts, count_values


class MutationExtractor:
    def __init__(
        self,
        samples,
        bams,
        mutlist,
        outprefix,
        barcode_files=None,
        target_barcodes=None,
        bins=100,
        seq_type="visium",
        run_type="UMI",
        cell_dict={},
        min_spot_number=0
    ):
        # list or str (seperated by ,)
        self.samples = samples.split(",") if isinstance(samples, str) else list(samples)
        self.bams = bams.split(",") if isinstance(bams, str) else list(bams)
        self.mutlist = mutlist
        self.outprefix = outprefix
        self.bins = int(bins) if isinstance(bins,int) else None
        self.seq_type = seq_type
        self.run_type = run_type
        self.cell_dict = cell_dict # for only one sample nowdays
        self.min_spot_number = int(min_spot_number)

        if isinstance(barcode_files, str):
            self.barcode_files = barcode_files.split(",")
        else:
            self.barcode_files = list(barcode_files) if barcode_files else None
            
        self.target_barcodes = target_barcodes
        
        self.barcode_lists = []
        self.spot_number = 0
        
        self._load_barcodes()

    def _load_barcodes(self):
        if self.barcode_files:
            self.barcode_lists = []
            self.spot_number = 0
            for i in self.barcode_files:
                barcode_dict = handle_barcode(i)
                self.spot_number += len(list(barcode_dict.keys()))
                self.barcode_lists.append(list(barcode_dict.keys()))
        elif self.target_barcodes:
            self.barcode_lists = []
            barcode_df = pd.read_csv(self.target_barcodes, sep="\t", header=None, names=["sample", "barcode"])
            self.spot_number = barcode_df.shape[0]
            for sample in self.samples:
                self.barcode_lists.append(barcode_df[barcode_df["sample"] == sample]["barcode"].tolist())
        elif self.cell_dict:
            self.barcode_lists = []
            self.barcode_lists.append(list(set(self.cell_dict.values())))
            self.spot_number += len(list(set(self.cell_dict.values())))

        else:
            raise FileNotFoundError(f'The barcode files are not found or not provided! If you run this command for one sample, cell_info file is also suported! Please check!')

    def run(self):
        mutation_identifier_list = read_mutation_file(self.mutlist)
        out_list = []

        out_dir = os.path.dirname(self.outprefix)
        if out_dir and not os.path.exists(out_dir):
            os.makedirs(out_dir, exist_ok=True)

        out_name = f"{self.outprefix}_spot_c_{self.spot_number}.csv"
        
        for mutation_identifier in mutation_identifier_list:
            posteria_values_list = []
            mut_likelihood_values_list = []
            nomut_likelihood_values_list = []
            alt_counts = 0
            total_counts = 0
            count_values_list = []
            ind_info = None
            sample = None
            
            for sample, bam_file, barcode_list in zip(self.samples, self.bams, self.barcode_lists):
                _, ind_info, posteria_values, \
                mut_likelihood_values, nomut_likelihood_values, \
                counts, count_values = tidy_func(
                    sample, barcode_list, bam_file, self.run_type, self.seq_type, self.bins,self.cell_dict, mutation_identifier
                )
                
                posteria_values_list += posteria_values
                mut_likelihood_values_list += mut_likelihood_values
                nomut_likelihood_values_list += nomut_likelihood_values
                alt_counts += counts[0]
                total_counts += counts[1]
                count_values_list += count_values

            # spots with coverage at this identifier (non-NA genotype)
            covered_spot_number = sum(1 for v in posteria_values_list if v != "NA")
            if covered_spot_number <= self.min_spot_number:
                continue
            
            out_list.append(
                [sample] +
                ind_info +
                posteria_values_list +
                mut_likelihood_values_list +
                nomut_likelihood_values_list +
                [f"{alt_counts}/{total_counts}"] +
                count_values_list
            )
        
        # 1. save CSV
        with open(out_name, "w") as out_file:
            for line in out_list:
                write_info = ",".join([str(k) for k in line if k != ""])
                out_file.write(f'{write_info}\n')
        
        # 2. save barcode file
        outsuffix = os.path.basename(self.outprefix)
        out_barcode_path = os.path.join(out_dir, f"{outsuffix}_scid_barcode.txt")
        
        with open(out_barcode_path, "w") as out_barcode_file:
            out_barcode_file.write("scid_basedTree\n")
            for sample, barcode_list in zip(self.samples, self.barcode_lists):
                for barcode in barcode_list:
                    out_barcode_file.write(f"{sample}_{barcode}\n")
        
        return out_name, out_barcode_path, self.spot_number


def main():
    parser = argparse.ArgumentParser(description="Extract mutation information per cell/spot from BAM files.")
    group = parser.add_mutually_exclusive_group(required=False)
    group.add_argument("--barcode_files", help="Path to barcode files, split by comma for multiple samples")
    group.add_argument("--target_barcodes", help="One file containing target barcodes, 1st col: sample, 2nd col: barcode")
    parser.add_argument("--bins", required=False, default=100, type=int, help="The combine bin level in the cluster file")
    parser.add_argument("--seq_type", dest='seq_type', default="visium", choices=["visium","stereo","ST","visium-HD"], type=str, help="Your input sequence type")
    parser.add_argument("--bams", required=True, help="BAM files, split by comma for multiple samples")
    parser.add_argument("--mutlist", required=True, help="The file storing mutation list")
    parser.add_argument("--outprefix", required=True, help="The output file prefix")
    parser.add_argument("--cell_info", required=False, default="", help="The 2-colums file record barcode and cell info,sep by Tab")
    parser.add_argument("--samples", required=True, type=str, help="Sample name, split by comma for multiple samples")
    parser.add_argument("--type", required=False, dest="run_type", default="UMI", choices=["UMI","read"], type=str, help="Do you want to use UMI or read")
    parser.add_argument(
        "--min_spot_number",
        required=False,
        default=0,
        type=int,
        help="Filter out an identifier if its covered spot number is <= this value"
    )
    
    args = parser.parse_args()
    if args.cell_info:
        cell_dict = barcode_cell_mapping(args.cell_info)
    else:
        cell_dict={}

    extractor = MutationExtractor(
        samples=args.samples,
        bams=args.bams,
        mutlist=args.mutlist,
        outprefix=args.outprefix,
        barcode_files=args.barcode_files,
        target_barcodes=args.target_barcodes,
        bins=args.bins,
        seq_type=args.seq_type,
        run_type=args.run_type,
        cell_dict=cell_dict,
        min_spot_number=args.min_spot_number
    )
    
    out_name, out_barcode_file, spot_num = extractor.run()
    out_spot_file=args.outprefix+".spot_num.txt"
    with open(out_spot_file,"w") as f:
        f.write(f'{spot_num}')
    f.close()

    print(f"Extraction completed successfully. Total spots processed: {spot_num}")


if __name__ == '__main__':
    main()

from collections import defaultdict
import gc
from itertools import groupby
from operator import itemgetter
import os
import subprocess
import tempfile

from SpaceTracer.utils.handle_UMI_combine import (
    BASE2IDX,
    IDX2BASE,
    KeyInterner,
    calculate_UMI_combine_phred,
    calculate_UMI_combine_phred_list,
    extract_barcode_and_umi_raw,
    format_barcode_key,
    get_most_candidate_allele,
    handel_barcode_name,
    handle_cigar,
    handle_pos,
    handle_quality_matrix,
    handle_seq,
    handle_seq_type,
)

from SpaceTracer.utils.logger import get_logger

model_name = __name__
logger = get_logger("<core function>: " + model_name)


def UMI_combination_spot_ind_and_judge_error(check_mosaic, check_error, site_barcode_UMI_dict, ref, threshold):
    all_genos = ["A", "T", "C", "G"]
    consensus_read_count = {}
    consensus_read_quality = {}
    raw_count_list = [0, 0, 0, 0]

    lysis_errors = 0
    lysis_alts = []
    error_allele = None

    for barcode in site_barcode_UMI_dict.keys():
        consensus_read_count[barcode] = [0, 0, 0, 0]
        consensus_read_quality[barcode] = {}

        for UMI in site_barcode_UMI_dict[barcode]:
            count_list = site_barcode_UMI_dict[barcode][UMI]["count"]

            if check_mosaic:
                quality_dict = site_barcode_UMI_dict[barcode][UMI]["quality"]
                phred_dict = calculate_UMI_combine_phred(count_list, quality_dict, weigh=0.5)
                candidate_allele, phred = get_most_candidate_allele(phred_dict, ref)
                index = "ATCG".index(candidate_allele)

                consensus_read_count[barcode][index] += 1

                if candidate_allele not in consensus_read_quality[barcode]:
                    consensus_read_quality[barcode][candidate_allele] = {}

                if str(phred) not in consensus_read_quality[barcode][candidate_allele]:
                    consensus_read_quality[barcode][candidate_allele][str(phred)] = 0

                consensus_read_quality[barcode][candidate_allele][str(phred)] += 1

                for index, dp in enumerate(count_list):
                    raw_count_list[index] += dp

            if check_error:
                _, lysis_error, _, lysis_alt = check_errors(count_list, ref, threshold)
                lysis_errors += lysis_error
                lysis_alts += lysis_alt

    if check_error:
        if lysis_errors > 1 or lysis_errors == 0:
            error_allele = None
        else:
            error_allele = lysis_alts[0]

    return consensus_read_count, consensus_read_quality, error_allele


def check_errors(count_list, ref, threshold=3):
    ref_index = "ATCG".index(ref)
    count_above_threshold = sum(1 for index, dp in enumerate(count_list) if index != ref_index and dp >= threshold)

    count_nonzero = sum(1 for dp in count_list if dp > 0)
    pcr_error = 1 if count_above_threshold >= 1 and count_nonzero >= 2 else 0
    pcr_alt = []
    if pcr_error:
        pcr_alt = ["ATCG"[index] for index, dp in enumerate(count_list) if index != ref_index and dp >= threshold]

    lysis_error = 0
    lysis_alt = []
    for index, dp in enumerate(count_list):
        if index != ref_index and dp > threshold:
            if all(count_list[other_idx] == 0 for other_idx in [0, 1, 2, 3] if other_idx != index):
                lysis_error = 1
                lysis_alt.append("ATCG"[index])
                break

    return pcr_error, lysis_error, pcr_alt, lysis_alt


# =========================================================
# 旧版兼容函数
# =========================================================

def handle_reads_per_pos_read_count_and_strand(reads, pos, run_type):
    pos_index = pos - 1
    site_barcode_UMI_dict = defaultdict(dict)
    reverse_dp = 0
    forward_dp = 0

    for item in reads:
        barcode_name, UMI_name = handle_seq_type(item, run_type, 1)
        if barcode_name is None or UMI_name is None:
            continue

        try:
            item.get_reference_positions().index(pos_index)
        except Exception:
            continue

        if item.is_reverse in [True, "TRUE", "true", "True"]:
            reverse_dp += 1
        else:
            forward_dp += 1

        seq_cut, pos_cut = handle_cigar(item.cigar)
        cut_seq = handle_seq(item.seq, seq_cut)
        cut_pos = handle_pos(item.get_reference_positions(), pos_cut)

        if pos_index in cut_pos:
            geno = cut_seq[cut_pos.index(pos_index)]
            if geno not in "ATCG":
                continue

            raw_index = handle_quality_matrix(cut_pos.index(pos_index), item.seq, cut_seq)
            try:
                qualities = item.get_forward_qualities()
                quality = qualities[raw_index]
            except Exception:
                continue

            if UMI_name not in site_barcode_UMI_dict[barcode_name]:
                site_barcode_UMI_dict[barcode_name][UMI_name] = {
                    "count": [0, 0, 0, 0],
                    "quality": {}
                }

            umi_entry = site_barcode_UMI_dict[barcode_name][UMI_name]
            umi_entry["count"]["ATCG".index(geno)] += 1

            if geno not in umi_entry["quality"]:
                umi_entry["quality"][geno] = {}

            umi_entry["quality"][geno][quality] = umi_entry["quality"][geno].get(quality, 0) + 1

    if reverse_dp >= forward_dp:
        major_read_strand = "-"
    elif reverse_dp < forward_dp:
        major_read_strand = "+"
    else:
        major_read_strand = "unknown"

    return site_barcode_UMI_dict, major_read_strand


def scan_region_reads_once_for_targets(reads, sites, run_type, bin_size):
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
        barcode_name, UMI_name = handle_seq_type(item, run_type, bin_size)
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

            if UMI_name not in site_barcode_UMI_dict[barcode_name]:
                site_barcode_UMI_dict[barcode_name][UMI_name] = {
                    "count": [0, 0, 0, 0],
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
        major_read_strand = "-" if reverse_dp >= forward_dp else "+"
        result[pos] = (bucket["site_barcode_UMI_dict"], major_read_strand)

    return result


def summarize_UMI_spot_for_both_mosaic_and_error(
    site_barcode_UMI_dict,
    strand,
    check_mosaic,
    check_error,
    identifier,
    threshold=3
):
    chrom = identifier[0]
    pos = int(identifier[1])
    ref = identifier[2]

    new_list = []
    all_genos = ["A", "T", "C", "G"]

    consensus_read_count, consensus_read_quality, error_allele = UMI_combination_spot_ind_and_judge_error(
        check_mosaic, check_error, site_barcode_UMI_dict, ref, threshold
    )

    if check_mosaic:
        for barcode in consensus_read_count.keys():
            alt = ",".join([alt_geno for alt_geno, dp in zip("ATCG", consensus_read_count[barcode]) if alt_geno != ref and dp != 0])
            if alt == "":
                alt = "."

            consensus_read_count_list = [str(i) for i in consensus_read_count[barcode]]
            consensus_read_qual_dict = defaultdict(str)

            for geno in all_genos:
                if geno not in consensus_read_quality[barcode]:
                    geno_qual_str = "NA"
                else:
                    geno_qual_list = [
                        str(phred) + ":" + str(consensus_read_quality[barcode][geno][phred])
                        for phred in consensus_read_quality[barcode][geno]
                    ]
                    geno_qual_str = ",".join(geno_qual_list)
                consensus_read_qual_dict[geno] = geno_qual_str

            consensus_read_count_str = ",".join(consensus_read_count_list)

            barcode_list = [
                chrom, pos, strand, ref, alt, barcode,
                consensus_read_count_str,
                consensus_read_qual_dict["A"],
                consensus_read_qual_dict["T"],
                consensus_read_qual_dict["C"],
                consensus_read_qual_dict["G"]
            ]
            new_list.append(barcode_list)

    if error_allele:
        error_list = [chrom, pos, ref, error_allele, strand]
    else:
        error_list = []

    return new_list, error_list


def process_single_region_for_umi_combine(bam_handle, region_info, seq_type, bin_size):
    chrom, start, end, sites = region_info

    try:
        reads = bam_handle.fetch(chrom, start, end)
        pos_data_map = scan_region_reads_once_for_targets(reads, sites, seq_type, bin_size)
    except Exception as exc:
        logger.warning(f"[UmiCombine] region failed {chrom}:{start}-{end}: {exc}")
        return []

    results = []
    for pos, ref, alt, check_mosaic, check_error in sites:
        site_barcode_UMI_dict, strand = pos_data_map.get(pos, ({}, "unknown"))
        mosaic_spot_list, error_list = summarize_UMI_spot_for_both_mosaic_and_error(
            site_barcode_UMI_dict=site_barcode_UMI_dict,
            strand=strand,
            check_mosaic=check_mosaic,
            check_error=check_error,
            identifier=(chrom, pos, ref),
        )
        results.append((mosaic_spot_list, error_list))

    del pos_data_map
    gc.collect()
    return results


# =========================================================
# 新版自动切换事件流 + debug
# =========================================================

def generate_region_events_once(
    reads,
    sites,
    run_type,
    bin_size,
    barcode_intern,
    umi_intern,
    cell_dict={}
):
    target_pos0_set = {pos - 1 for pos, *_ in sites}

    for item in reads:
        barcode_key, raw_umi = extract_barcode_and_umi_raw(item, run_type, bin_size)
        if barcode_key is None or raw_umi is None:
            continue

        barcode_name = handel_barcode_name(cell_dict, barcode_key)
        umi_key = (barcode_name, raw_umi)

        try:
            seq_cut, pos_cut = handle_cigar(item.cigar)
            cut_seq = handle_seq(item.seq, seq_cut)
            cut_pos = handle_pos(item.get_reference_positions(), pos_cut)
            qualities = item.get_forward_qualities()
        except Exception:
            continue

        if not cut_pos or not cut_seq or qualities is None:
            continue

        barcode_id = barcode_intern.encode(barcode_name)
        umi_id = umi_intern.encode(umi_key)
        strand_bit = 0 if item.is_reverse in [True, "TRUE", "true", "True"] else 1

        for idx, ref_pos0 in enumerate(cut_pos):
            if ref_pos0 not in target_pos0_set:
                continue

            geno = cut_seq[idx]
            if geno not in BASE2IDX:
                continue

            try:
                raw_index = handle_quality_matrix(idx, item.seq, cut_seq)
                quality = qualities[raw_index]
            except Exception:
                continue

            if quality is None or quality < 0 or quality > 40:
                continue

            pos1 = ref_pos0 + 1
            base_idx = BASE2IDX[geno]
            yield (pos1, barcode_id, umi_id, base_idx, quality, strand_bit)


def write_event_line(fh, event):
    pos1, barcode_id, umi_id, base_idx, qual, strand_bit = event
    fh.write(f"{pos1}\t{barcode_id}\t{umi_id}\t{base_idx}\t{qual}\t{strand_bit}\n")


def sort_events_tsv_int(in_tsv_path, out_sorted_tsv_path, tmp_dir=None):
    cmd = [
        "sort",
        "-t", "\t",
        "-k1,1n",
        "-k2,2n",
        "-k3,3n",
    ]
    if tmp_dir is not None:
        cmd.extend(["-T", tmp_dir])
    cmd.append(in_tsv_path)

    with open(out_sorted_tsv_path, "w") as fout:
        subprocess.run(cmd, stdout=fout, check=True)


def read_sorted_events_tsv_int(sorted_tsv_path):
    with open(sorted_tsv_path, "r") as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            pos, barcode_id, umi_id, base_idx, qual, strand_bit = line.split("\t")
            yield (
                int(pos),
                int(barcode_id),
                int(umi_id),
                int(base_idx),
                int(qual),
                int(strand_bit),
            )


def aggregate_one_umi_group_int(group_iter):
    count_list = [0, 0, 0, 0]
    quality_list = [[0] * 41 for _ in range(4)]
    forward_dp = 0
    reverse_dp = 0

    for _, _, _, base_idx, qual, strand_bit in group_iter:
        count_list[base_idx] += 1
        quality_list[base_idx][qual] += 1
        if strand_bit == 1:
            forward_dp += 1
        else:
            reverse_dp += 1

    return count_list, quality_list, forward_dp, reverse_dp


def summarize_one_umi_group_int(
    pos,
    barcode_id,
    umi_id,
    ref,
    alt,
    count_list,
    quality_list,
    threshold=3,
    weigh=0.5
):
    phred_dict = calculate_UMI_combine_phred_list(
        count_list=count_list,
        quality_list=quality_list,
        weigh=weigh
    )

    candidate_allele, phred = get_most_candidate_allele(phred_dict, ref)
    _, lysis_error, _, lysis_alt = check_errors(count_list, ref, threshold)

    return {
        "pos": pos,
        "barcode_id": barcode_id,
        "umi_id": umi_id,
        "candidate_allele": candidate_allele,
        "phred": phred,
        "lysis_error": lysis_error,
        "lysis_alt": lysis_alt[0] if lysis_alt else None,
        "ref": ref,
        "alt": alt,
    }


def stream_umi_summaries_from_sorted_events_int(sorted_events, site_info_dict, threshold=3, weigh=0.5):
    keyfunc = itemgetter(0, 1, 2)

    for (pos, barcode_id, umi_id), group in groupby(sorted_events, key=keyfunc):
        count_list, quality_list, forward_dp, reverse_dp = aggregate_one_umi_group_int(group)

        site_info = site_info_dict[pos]
        umi_summary = summarize_one_umi_group_int(
            pos=pos,
            barcode_id=barcode_id,
            umi_id=umi_id,
            ref=site_info["ref"],
            alt=site_info["alt"],
            count_list=count_list,
            quality_list=quality_list,
            threshold=threshold,
            weigh=weigh
        )

        umi_summary["strand"] = "-" if reverse_dp >= forward_dp else "+"
        yield umi_summary

        
def stream_barcode_summaries_from_sorted_umi_summaries_int(sorted_umi_summaries, barcode_intern):
    def keyfunc_pos(x):
        return x["pos"]

    for pos, pos_group in groupby(sorted_umi_summaries, key=keyfunc_pos):
        pos_items = list(pos_group)

        # 1) 先在位点级别决定唯一 strand
        strand_plus = 0
        strand_minus = 0
        for item in pos_items:
            if item["strand"] == "+":
                strand_plus += 1
            else:
                strand_minus += 1

        site_major_strand = "-" if strand_minus >= strand_plus else "+"

        # 2) 再按 barcode 聚合，但不丢弃任何 UMI
        pos_items.sort(key=lambda x: x["barcode_id"])

        for barcode_id, barcode_group in groupby(pos_items, key=lambda x: x["barcode_id"]):
            consensus_read_count = [0, 0, 0, 0]
            consensus_read_quality = {base: defaultdict(int) for base in "ATCG"}

            ref = None
            lysis_alts = []

            for item in barcode_group:
                ref = item["ref"]
                allele = item["candidate_allele"]
                phred = item["phred"]

                if allele in BASE2IDX:
                    base_idx = BASE2IDX[allele]
                    consensus_read_count[base_idx] += 1
                    consensus_read_quality[allele][str(phred)] += 1

                if item["lysis_error"] and item["lysis_alt"] is not None:
                    lysis_alts.append(item["lysis_alt"])

            if len(lysis_alts) == 1:
                error_allele = lysis_alts[0]
            else:
                error_allele = None

            barcode_key = barcode_intern.decode(barcode_id)
            barcode_name = format_barcode_key(barcode_key)

            yield {
                "pos": pos,
                "barcode_id": barcode_id,
                "barcode_name": barcode_name,
                "ref": ref,
                "major_strand": site_major_strand,  # 统一为位点级 strand
                "consensus_read_count": consensus_read_count,
                "consensus_read_quality": consensus_read_quality,
                "error_allele": error_allele,
            }


def format_barcode_summary_to_row(chrom, barcode_summary):
    pos = barcode_summary["pos"]
    barcode_name = barcode_summary["barcode_name"]
    ref = barcode_summary["ref"]
    major_strand = barcode_summary["major_strand"]
    counts = barcode_summary["consensus_read_count"]
    qual = barcode_summary["consensus_read_quality"]
    error_allele = barcode_summary["error_allele"]

    alt = ",".join([base for base, dp in zip("ATCG", counts) if base != ref and dp != 0])
    if alt == "":
        alt = "."

    def qual_to_text(base):
        if not qual[base]:
            return "NA"
        return ",".join(
            f"{k}:{qual[base][k]}"
            for k in sorted(qual[base], key=int)
        )

    barcode_row = [
        chrom,
        pos,
        major_strand,
        ref,
        alt,
        barcode_name,
        ",".join(map(str, counts)),
        qual_to_text("A"),
        qual_to_text("T"),
        qual_to_text("C"),
        qual_to_text("G"),
    ]

    if error_allele:
        error_row = [chrom, pos, ref, error_allele, major_strand]
    else:
        error_row = []

    return barcode_row, error_row


def collect_events_with_auto_spill(
    reads,
    sites,
    run_type,
    bin_size,
    barcode_intern,
    umi_intern,
    events_threshold=300000,
    tmp_root=None,
    cell_dict={}
):
    """
    返回:
        mode: "memory" or "disk"
        payload:
            if memory -> list(events)
            if disk   -> sorted_tsv_path
        tmp_ctx:
            TemporaryDirectory对象或None
        debug_info:
            {
                "mode": ...,
                "spill_happened": bool,
                "n_events": int,
                "n_barcodes": int,
                "n_umis": int,
                "events_threshold": int,
            }
    """
    events = []
    tmp_ctx = None
    tmpdir = None
    raw_tsv = None
    fh = None
    spilled = False
    event_count = 0

    event_iter = generate_region_events_once(
        reads=reads,
        sites=sites,
        run_type=run_type,
        bin_size=bin_size,
        barcode_intern=barcode_intern,
        umi_intern=umi_intern,
        cell_dict=cell_dict
    )

    for event in event_iter:
        event_count += 1

        if not spilled:
            events.append(event)
            if len(events) > events_threshold:
                tmp_ctx = tempfile.TemporaryDirectory(dir=tmp_root)
                tmpdir = tmp_ctx.name
                raw_tsv = os.path.join(tmpdir, "events.tsv")
                fh = open(raw_tsv, "w")

                for old_event in events:
                    write_event_line(fh, old_event)

                events = None
                spilled = True
        else:
            write_event_line(fh, event)

    if not spilled:
        debug_info = {
            "mode": "memory",
            "spill_happened": False,
            "n_events": event_count,
            "n_barcodes": len(barcode_intern),
            "n_umis": len(umi_intern),
            "events_threshold": events_threshold,
        }
        return "memory", events, None, debug_info

    fh.close()

    sorted_tsv = os.path.join(tmpdir, "events.sorted.tsv")
    sort_events_tsv_int(raw_tsv, sorted_tsv, tmp_dir=tmpdir)

    debug_info = {
        "mode": "disk",
        "spill_happened": True,
        "n_events": event_count,
        "n_barcodes": len(barcode_intern),
        "n_umis": len(umi_intern),
        "events_threshold": events_threshold,
    }
    return "disk", sorted_tsv, tmp_ctx, debug_info



def process_region_with_auto_events(
    bam_handle,
    region_info,
    seq_type,
    bin_size,
    tmp_root=None,
    threshold=3,
    weigh=0.5,
    events_threshold=300000,
    cell_dict={},
    debug_log=False
):
    chrom, start, end, sites = region_info

    site_info_dict = {
        pos: {
            "ref": ref,
            "alt": alt,
            "check_mosaic": check_mosaic,
            "check_error": check_error,
        }
        for pos, ref, alt, check_mosaic, check_error in sites
    }

    barcode_intern = KeyInterner()
    umi_intern = KeyInterner()

    reads = bam_handle.fetch(chrom, start, end)

    mode = None
    payload = None
    tmp_ctx = None
    debug_info = None
    events = None
    sorted_events = None
    umi_summaries_iter = None
    barcode_summaries_iter = None
    results = []

    mode, payload, tmp_ctx, debug_info = collect_events_with_auto_spill(
        reads=reads,
        sites=sites,
        run_type=seq_type,
        bin_size=bin_size,
        barcode_intern=barcode_intern,
        umi_intern=umi_intern,
        events_threshold=events_threshold,
        tmp_root=tmp_root,
        cell_dict=cell_dict
    )

    try:
        if mode == "memory":
            events = payload
            payload = None

            if debug_log:
                logger.info(
                    "[umi_combine:auto] "
                    f"region={chrom}:{start}-{end} "
                    f"n_sites={len(sites)} "
                    f"mode={debug_info['mode']} "
                    f"spill_happened={debug_info['spill_happened']} "
                    f"n_events={debug_info['n_events']} "
                    f"n_barcodes={debug_info['n_barcodes']} "
                    f"n_umis={debug_info['n_umis']} "
                    f"events_threshold={debug_info['events_threshold']}"
                )

            if not events:
                debug_info["n_result_rows"] = 0
                return [], debug_info

            events.sort(key=itemgetter(0, 1, 2))
            sorted_events = iter(events)

        else:
            sorted_tsv = payload
            payload = None

            if debug_log:
                logger.info(
                    "[umi_combine:auto] "
                    f"region={chrom}:{start}-{end} "
                    f"n_sites={len(sites)} "
                    f"mode={debug_info['mode']} "
                    f"spill_happened={debug_info['spill_happened']} "
                    f"n_events={debug_info['n_events']} "
                    f"n_barcodes={debug_info['n_barcodes']} "
                    f"n_umis={debug_info['n_umis']} "
                    f"events_threshold={debug_info['events_threshold']} "
                    f"tmp_sorted={sorted_tsv}"
                )

            sorted_events = read_sorted_events_tsv_int(sorted_tsv)

        umi_summaries_iter = stream_umi_summaries_from_sorted_events_int(
            sorted_events=sorted_events,
            site_info_dict=site_info_dict,
            threshold=threshold,
            weigh=weigh
        )

        barcode_summaries_iter = stream_barcode_summaries_from_sorted_umi_summaries_int(
            sorted_umi_summaries=umi_summaries_iter,
            barcode_intern=barcode_intern
        )

        for barcode_summary in barcode_summaries_iter:
            barcode_row, error_row = format_barcode_summary_to_row(chrom, barcode_summary)
            results.append((barcode_row, error_row))

        debug_info["n_result_rows"] = len(results)

        if debug_log:
            logger.info(
                "[umi_combine:auto:done] "
                f"region={chrom}:{start}-{end} "
                f"n_sites={len(sites)} "
                f"mode={debug_info['mode']} "
                f"spill_happened={debug_info['spill_happened']} "
                f"n_events={debug_info['n_events']} "
                f"n_barcodes={debug_info['n_barcodes']} "
                f"n_umis={debug_info['n_umis']} "
                f"n_result_rows={debug_info['n_result_rows']}"
            )

        return results, debug_info

    finally:
        # 主动释放大对象引用，减少 region 间残留
        try:
            del reads
        except Exception:
            pass
        try:
            del payload
        except Exception:
            pass
        try:
            del events
        except Exception:
            pass
        try:
            del sorted_events
        except Exception:
            pass
        try:
            del umi_summaries_iter
        except Exception:
            pass
        try:
            del barcode_summaries_iter
        except Exception:
            pass
        try:
            del barcode_intern
        except Exception:
            pass
        try:
            del umi_intern
        except Exception:
            pass
        try:
            del site_info_dict
        except Exception:
            pass

        gc.collect()

        if tmp_ctx is not None:
            tmp_ctx.cleanup()


def process_single_region_for_umi_combine_auto(
    bam_handle,
    region_info,
    seq_type,
    bin_size,
    tmp_root=None,
    threshold=3,
    weigh=0.5,
    events_threshold=300000,
    cell_dict={},
    debug_log=False
):

    chrom, start, end, sites = region_info

    flat_results = None
    debug_info = None
    pos_to_rows = defaultdict(list)
    pos_to_error = {}
    results = []

    try:
        flat_results, debug_info = process_region_with_auto_events(
            bam_handle=bam_handle,
            region_info=region_info,
            seq_type=seq_type,
            bin_size=bin_size,
            tmp_root=tmp_root,
            threshold=threshold,
            weigh=weigh,
            events_threshold=events_threshold,
            cell_dict=cell_dict,
            debug_log=debug_log
        )

        if not flat_results:
            return [([], []) for _ in sites]

        for barcode_row, error_row in flat_results:
            if barcode_row:
                pos = barcode_row[1]
                pos_to_rows[pos].append(barcode_row)
            if error_row:
                pos = error_row[1]
                pos_to_error[pos] = error_row

        for pos, ref, alt, check_mosaic, check_error in sites:
            mosaic_spot_list = pos_to_rows.get(pos, [])
            error_list = pos_to_error.get(pos, [])
            results.append((mosaic_spot_list, error_list))

        return results

    finally:
        try:
            del flat_results
        except Exception:
            pass
        try:
            del debug_info
        except Exception:
            pass
        try:
            del pos_to_rows
        except Exception:
            pass
        try:
            del pos_to_error
        except Exception:
            pass

        gc.collect()

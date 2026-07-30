from collections import defaultdict
from functools import reduce
from math import ceil, log10

from SpaceTracer.utils.logger import get_logger
from SpaceTracer.utils.utils import round_to_nearest_bin

model_name = __name__
logger = get_logger(model_name)

BASE2IDX = {'A': 0, 'T': 1, 'C': 2, 'G': 3}
IDX2BASE = ['A', 'T', 'C', 'G']


def handle_cigar(cigar_symbol):
    '''
    ## handle cigar
    # [(0, 76), (2, 1), (0, 33), (3, 139241), (0, 11)]
    # '76M1D33M139241N11M'
    # the 1st is symbol; and the 2nd is count
    # 0: Match; 1: Insertion; 2: deletion; 3: N; 4: S; 5: H; 6: P; 7: =; 8: X
    output:
    seq_cut may be a tuple, contain the softclip information;
    pos_cut are the indel information: [(before_length, insert number)]
    '''
    seq_length_before = 0
    pos_length_before = 0

    seq_cut_start = None
    seq_cut_end = None
    pos_cut = []

    for cigars, i in zip(cigar_symbol, range(1, len(cigar_symbol) + 1)):
        symbol = cigars[0]
        count = cigars[1]
        if symbol in [5, 6, 7, 8]:
            logger.debug(cigar_symbol)
        elif symbol in [0, 1, 4]:
            seq_length_before += count
            if symbol == 0:
                pos_length_before += count
            elif symbol == 4:
                if i == 1:
                    seq_cut_start = seq_length_before
                elif i == len(cigar_symbol):
                    seq_cut_end = seq_length_before
                else:
                    logger.debug(cigar_symbol)
            elif symbol == 1:
                pos_cut.append((pos_length_before, count))
        else:
            pass

    seq_cut = (seq_cut_start, seq_cut_end)
    return seq_cut, pos_cut


def handle_seq(seq, seq_cut):
    cut_seq = seq[seq_cut[0]:seq_cut[1]]
    return cut_seq


def handle_pos(pos_matrix, pos_cut):
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
            pos = pos_cut[item][0]
            times = pos_cut[item][1]
            cut_pos_matrix = cut_pos_matrix + pos_matrix[start:pos] + [""] * times
            start = pos
        last_pos = pos_cut[-1][0]
        cut_pos_matrix = cut_pos_matrix + pos_matrix[last_pos:]
    return cut_pos_matrix


def handle_quality_matrix(mutation_in_cutseq_index, seq, cut_seq):
    if len(cut_seq[mutation_in_cutseq_index:]) >= len(cut_seq[:mutation_in_cutseq_index]):
        query_str = cut_seq[mutation_in_cutseq_index:]
        raw_index = seq.index(query_str)
    else:
        query_str = cut_seq[:mutation_in_cutseq_index]
        raw_index = seq.index(query_str) + len(query_str)
    return raw_index


def phred_2_q(phred):
    try:
        phred = int(phred)
        q = 1 - (10 ** -(phred / 10))
    except ValueError as e:
        q = 0
        print(f"wrong phred format: {phred}")
        print("Error:", e)
    return q


def q_2_phred(q):
    try:
        q = float(q)
        if q == 1.0:
            phred = 40
        else:
            phred = ceil(-log10(1 - q) * 10)
    except ValueError as e:
        print(f"wrong q format: {q}")
        print("Error:", e)
        phred = 0
    return phred


def trans(qual_list):
    if qual_list == []:
        qual_str = "NA"
    else:
        qual_str = ",".join(qual_list)
    return qual_str


def format_barcode_key(barcode_key):
    if isinstance(barcode_key, tuple):
        return "_".join(map(str, barcode_key))
    return str(barcode_key)


def handel_barcode_name(cell_dict, barcode_key):
    barcode_name = format_barcode_key(barcode_key)
    if cell_dict != {}:
        return cell_dict.get(barcode_name, None)
    else:
        return barcode_name


def extract_barcode_and_umi_raw(read, run_type, bins):
    if run_type == "visium" or run_type =="visium-HD":
        try:
            CB = read.get_tag("CB").strip()
            UB = read.get_tag("UB").strip()
            barcode_key = str(CB)
            umi_key = str(UB)
        except Exception:
            return None, None

    elif run_type == "stereo":
        try:
            Cx_raw = int(read.get_tag("Cx"))
            Cy_raw = int(read.get_tag("Cy"))

            if bins != 1:
                Cx = round_to_nearest_bin(Cx_raw, bins)
                Cy = round_to_nearest_bin(Cy_raw, bins)
            else:
                Cx = Cx_raw
                Cy = Cy_raw

            UR = read.get_tag("UR").strip()
            barcode_key = (Cx, Cy)
            umi_key = str(UR)
        except Exception:
            return None, None

    elif run_type == "ST":
        try:
            CB = str(read.get_tag("B0"))
            UB = str(read.get_tag("B3"))
            barcode_key = str(CB)
            umi_key = str(UB)
        except Exception:
            return None, None

    else:
        raise ValueError(f"The input type is not recognized {run_type}")

    return barcode_key, umi_key


def handle_seq_type(read, run_type, bins, cell_dict={}):
    """
    return:
        barcode_name, UMI_name
    """
    barcode_key, umi_key = extract_barcode_and_umi_raw(read, run_type, bins)
    if barcode_key is None or umi_key is None:
        return None, None

    barcode_name = handel_barcode_name(cell_dict, barcode_key)
    if barcode_name is None:
        return None,None

    UMI_name = barcode_name + "_" + str(umi_key)
    return barcode_name, UMI_name


class KeyInterner:
    def __init__(self):
        self.key2id = {}
        self.id2key = []

    def encode(self, key):
        if key in self.key2id:
            return self.key2id[key]
        idx = len(self.id2key)
        self.key2id[key] = idx
        self.id2key.append(key)
        return idx

    def decode(self, idx):
        return self.id2key[idx]

    def __len__(self):
        return len(self.id2key)


def calculate_UMI_combine_phred_count_dict(count_dict, quality_dict, weigh=0.5):
    all_genos = ["A", "T", "C", "G"]
    pcr_error = 1e-6
    no_pcr_error = (1.0 - pcr_error) ** 100
    rightP = 1.0
    sumP = 0.0
    dp = sum(count_dict.values())

    proP_dict = defaultdict(lambda: 1.0)
    pcrP_dict = defaultdict(float)
    likelihood_dict = defaultdict(float)
    phred_dict = defaultdict(float)

    for geno in count_dict.keys():
        qual_geno_list = [phred_2_q(key) ** int(quality_dict[geno][key]) for key in quality_dict[geno].keys()]
        qual_geno = reduce(lambda x, y: x * y, qual_geno_list)
        proP_dict[geno] *= qual_geno

        for other_geno in quality_dict.keys() - set([geno]):
            other_qual_geno_list = [(1 - phred_2_q(key)) ** int(quality_dict[other_geno][key]) for key in quality_dict[other_geno].keys()]
            if other_qual_geno_list == []:
                continue
            other_qual_geno = reduce(lambda x, y: x * y, other_qual_geno_list)
            proP_dict[geno] *= other_qual_geno

        rightP = rightP * qual_geno

    for geno in all_genos:
        count_geno = 0 if geno not in count_dict.keys() else count_dict[geno]
        ratio = (count_geno + 0.5) / (dp + 0.5 * 4) if dp > 0 else 0.25
        pcrP = 10.0 ** (-6.0 * ratio)
        pcrP_dict[geno] = pcrP

    for geno in all_genos:
        if geno in count_dict.keys():
            base_calling_error = proP_dict[geno]
            no_base_calling_error = rightP
            pcr_error = min([pcrP_dict[char] for char in pcrP_dict.keys() if char != geno])
            likelihood_value = weigh * no_pcr_error * base_calling_error + (1 - weigh) * no_base_calling_error * pcr_error
        else:
            likelihood_value = rightP
            for char in set(all_genos) - set([geno]):
                likelihood_value *= pcrP_dict[char]

        likelihood_dict[geno] = likelihood_value
        sumP += likelihood_value

    for geno in likelihood_dict.keys():
        phred_dict[geno] = 0 if sumP <= 0 else q_2_phred(likelihood_dict[geno] / sumP)

    return phred_dict


def calculate_UMI_combine_phred(count_list, quality_dict, weigh=0.5):
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

        qual_geno_list = [
            phred_2_q(key) ** int(quality_dict[geno][key])
            for key in quality_dict[geno].keys()
        ]
        qual_geno = reduce(lambda x, y: x * y, qual_geno_list)
        proP_dict[geno] *= qual_geno

        for other_geno in quality_dict.keys() - set([geno]):
            other_qual_geno_list = [
                (1 - phred_2_q(key)) ** int(quality_dict[other_geno][key])
                for key in quality_dict[other_geno].keys()
            ]
            if other_qual_geno_list:
                other_qual_geno = reduce(lambda x, y: x * y, other_qual_geno_list)
                proP_dict[geno] *= other_qual_geno

        rightP = rightP * qual_geno

    for index, geno in enumerate("ATCG"):
        count_geno = count_list[index]
        ratio = (count_geno + 0.5) / (dp + 0.5 * 4) if dp > 0 else 0.25
        pcrP = 10.0 ** (-6.0 * ratio)
        pcrP_dict[geno] = pcrP

    for index, geno in enumerate("ATCG"):
        if count_list[index] != 0:
            base_calling_error = proP_dict[geno]
            no_base_calling_error = rightP

            other_pcrP = [pcrP_dict[char] for char in pcrP_dict.keys() if char != geno]
            pcr_error_value = min(other_pcrP) if other_pcrP else 1.0

            likelihood_value = (
                weigh * no_pcr_error * base_calling_error
                + (1 - weigh) * no_base_calling_error * pcr_error_value
            )
        else:
            likelihood_value = rightP
            for char in set(all_genos) - set([geno]):
                likelihood_value *= pcrP_dict[char]

        likelihood_dict[geno] = likelihood_value
        sumP += likelihood_value

    for geno in likelihood_dict.keys():
        phred_dict[geno] = 0 if sumP <= 0 else q_2_phred(likelihood_dict[geno] / sumP)

    return phred_dict


def calculate_UMI_combine_phred_list(count_list, quality_list, weigh=0.5):
    all_genos = ["A", "T", "C", "G"]
    geno2idx = {"A": 0, "T": 1, "C": 2, "G": 3}

    pcr_error = 1e-6
    no_pcr_error = (1.0 - pcr_error) ** 100

    dp = sum(count_list)

    proP_dict = defaultdict(lambda: 1.0)
    pcrP_dict = defaultdict(float)
    likelihood_dict = defaultdict(float)
    phred_dict = defaultdict(float)

    rightP = 1.0

    for geno in all_genos:
        geno_idx = geno2idx[geno]
        dp_count = count_list[geno_idx]
        if dp_count == 0:
            continue

        qual_geno = 1.0
        for q, cnt in enumerate(quality_list[geno_idx]):
            if cnt == 0:
                continue
            qual_geno *= phred_2_q(q) ** int(cnt)

        proP_dict[geno] *= qual_geno

        for other_geno in all_genos:
            if other_geno == geno:
                continue

            other_idx = geno2idx[other_geno]
            other_qual_geno = 1.0
            has_other = False

            for q, cnt in enumerate(quality_list[other_idx]):
                if cnt == 0:
                    continue
                other_qual_geno *= (1 - phred_2_q(q)) ** int(cnt)
                has_other = True

            if has_other:
                proP_dict[geno] *= other_qual_geno

        rightP *= qual_geno

    for geno in all_genos:
        idx = geno2idx[geno]
        count_geno = count_list[idx]
        ratio = (count_geno + 0.5) / (dp + 0.5 * 4) if dp > 0 else 0.25
        pcrP = 10.0 ** (-6.0 * ratio)
        pcrP_dict[geno] = pcrP

    sumP = 0.0

    for geno in all_genos:
        idx = geno2idx[geno]
        if count_list[idx] != 0:
            base_calling_error = proP_dict[geno]
            no_base_calling_error = rightP

            other_pcrP = [pcrP_dict[char] for char in all_genos if char != geno]
            pcr_error_value = min(other_pcrP) if other_pcrP else 1.0

            likelihood_value = (
                weigh * no_pcr_error * base_calling_error
                + (1 - weigh) * no_base_calling_error * pcr_error_value
            )
        else:
            likelihood_value = rightP
            for char in all_genos:
                if char == geno:
                    continue
                likelihood_value *= pcrP_dict[char]

        likelihood_dict[geno] = likelihood_value
        sumP += likelihood_value

    for geno in all_genos:
        phred_dict[geno] = 0 if sumP <= 0 else q_2_phred(likelihood_dict[geno] / sumP)

    return phred_dict


def get_most_candidate_allele(phred_dict, ref_allele):
    rank_list = sorted(phred_dict.items(), key=lambda item: item[1], reverse=True)
    major_allele = rank_list[0][0]
    major_allele_phred = rank_list[0][1]
    ref_allele_phred = phred_dict[ref_allele]

    if major_allele != ref_allele and ref_allele_phred >= major_allele_phred:
        candidate_allele = ref_allele
        phred = ref_allele_phred
    else:
        candidate_allele = major_allele
        phred = major_allele_phred

    return candidate_allele, phred


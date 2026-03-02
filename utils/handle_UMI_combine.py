

from math import ceil, log10

import numpy as np
from SpaceTracer.utils.logger import get_logger

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


def round_to_nearest_bin(x,bins):
    return int(np.ceil(x / bins) * bins)


def handle_seq_type(read,run_type,bins):
    if run_type=="visium":
        try:
            CB=read.get_tag("CB").strip()
            UB=read.get_tag("UB").strip()

            barcode_name=str(CB)
            UMI_name=str(UB)
        except:
            return None,None

    # elif run_type=="stereo":
    #     Cx=str(read.get_tag("Cx"))
    #     Cy=str(read.get_tag("Cy"))
    #     UR=read.get_tag("UR").strip()

    #     barcode_name=Cx+"_"+Cy
    #     UMI_name=str(UR)
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
            UMI_name=str(Cx_raw)+"_"+str(Cy_raw)+"_"+str(UR)
        except:
            return None,None

    elif run_type=="ST":
        try:
            CB=str(read.get_tag("B0"))
            UB=str(read.get_tag("B3"))

            barcode_name=str(CB)
            UMI_name=str(UB)
        except:
            return None,None
    else:
        raise ValueError(f"The input type is not recognized {run_type}")
    
    return barcode_name, UMI_name
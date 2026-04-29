#!/usr/bin/env python3

import sys
from typing import TextIO

from SpaceTracer.utils.logger import get_logger
model_name=__name__
logger = get_logger(model_name)

def get_sorted_alleles_counter(ref,base_counts):
    from collections import Counter,defaultdict
    """ sort the counts"""
    alleles = ['A', 'T', 'C', 'G']
    allele_dict=defaultdict(int)
    for allele in alleles:
        total = base_counts[allele] + base_counts[allele.lower()]
        allele_dict[allele] = total
    
    allele_dict[ref]=base_counts['REF'] + base_counts['ref'] 
    counter = Counter(allele_dict)
    sorted_items = counter.most_common()
    
    # result = [allele1, allele1_count，allele2, allele2_count，allele3, allele3_count，allele4, allele4_count] which are sorted by count
    result = []
    for allele, count in sorted_items:
        if allele!=ref:
            result+=[allele, count]
        else:
            result=[allele, count]+result
    return result

class PileupHandle:
    """ tidy and filter mpileup results (only minimal read depth is required)"""
    def __init__(self, min_read_depth: int = 30):
        """
        min_read_depth: the minimal read depth used, and this filtration is for informative depth, which not include bases with "N" (default: 30)
        """
        self.min_read_depth = min_read_depth
    
    def filter_pileup(self, input_stream: TextIO, output_stream: TextIO):
        """ the run function  """
        for line in input_stream:
            filtered_line = self._filter_line(line)
            
            if filtered_line:
                output_stream.write(filtered_line + '\n')
    
    def _filter_line(self, line: str) -> str:
        """ Core function, for each line (those not pass filtration will return None) """
        line = line.strip()
        if not line:
            return None
        
        parts = line.split('\t')
        
        if len(parts) < 6:
            return None
        
        chrom = parts[0]
        pos = parts[1]
        ref = parts[2]
        depth = int(parts[3])
        bases = parts[4] 
        qualities = parts[5]
        
        base_counts = self._count_bases(bases, qualities)
        info_dp=sum(base_counts.values())
        if info_dp >=self.min_read_depth:
            result=get_sorted_alleles_counter(ref,base_counts)
            output_parts = [chrom,pos,ref,depth,info_dp] + result
            return '\t'.join(str(x) for x in output_parts)
        else:
            return None

    def _count_bases(self, bases: str, qualities: str) -> dict:
        """
        Only bases are under consideration, qualities are ignored

        For each read covering the position, this column contains (from samtools mpileup):
            If this is the first position covered by the read, a “^” character followed by the alignment's mapping quality encoded as an ASCII character.
            A single character indicating the read base and the strand to which the read has been mapped:
            Forward	Reverse	Meaning
            . dot	, comma	Base matches the reference base
            ACGTN	acgtn	Base is a mismatch to the reference base
            >	<	Reference skip (due to CIGAR “N”)
            *	*/#	Deletion of the reference base (CIGAR “D”)
            Deleted bases are shown as “*” on both strands unless --reverse-del is used, in which case they are shown as “#” on the reverse strand.

            If there is an insertion after this read base, text matching “\+[0-9]+[ACGTNacgtn*#]+”: a “+” character followed by an integer giving the length of the insertion and then the inserted sequence. Pads are shown as “*” unless --reverse-del is used, in which case pads on the reverse strand will be shown as “#”.
            If there is a deletion after this read base, text matching “-[0-9]+[ACGTNacgtn]+”: a “-” character followed by the deleted reference bases represented similarly. (Subsequent pileup lines will contain “*” for this read indicating the deleted bases.)
            If this is the last position covered by the read, a “$” character.
        """
        min_base_qual=0
        counts = {
            'A': 0, 'a': 0,
            'C': 0, 'c': 0,
            'G': 0, 'g': 0,
            'T': 0, 't': 0,
            'REF':0, 'ref':0,
            'insertions': 0,
            'deletions': 0
        }
        i = 0  # index of bases
        qual_i = 0  # index of qualities
        
        if bases == "<"*len(bases):
            return counts
        
        while i < len(bases):
            base = bases[i]
            # 1. mark of read start '^'
            if base == '^':
                i += 2  # skip ^ and the quality string
                continue
            
            # 2. mark of read end '$'
            elif base == '$':
                i += 1 # skip $
                continue
            
            # 3. mark of insertion '+'
            elif base == '+':
                i += 1  # skip '+'
                insert_len_str = ''
                while i < len(bases) and bases[i].isdigit():
                    insert_len_str += bases[i]
                    i += 1 # skip '+'
                insert_len = int(insert_len_str) if insert_len_str else 0
                i += insert_len  # skin the string represents insertion length
                counts['insertions'] += 1
                continue
            
            # 4. mark of deletion '-'
            elif base == '-':
                i += 1  # skip '-'
                del_len_str = ''
                while i < len(bases) and bases[i].isdigit():
                    del_len_str += bases[i]
                    i += 1
                del_len = int(del_len_str) if del_len_str else 0
                i += del_len  # skip the string represents deletion length
                counts['deletions'] += 1
                continue
            
            # 5. mark of reference skip (due to CIGAR “N”) 
            elif base in ['<','>']: 
                i += 1 # skip 
                continue

            # 5. base deletion '*' and '#'
            elif base in ['*','#']:
                i += 1  # deletion of the base (CIGAR “D”)
                continue

            # 6. splice junction also represent as "N" 
            elif base=='N' or base=='n':
                i += 1
                continue

            # 7. bases for count (.,ATCGatcg)
            else:
                if qual_i < len(qualities):
                    qual = ord(qualities[qual_i]) - 33  # if the Phred is not 33 please change this!
                    qual_i += 1  # quality index
                else:
                    qual = 0  # not known
                
                # filter quality (not applied, only a api dor further update)
                if qual >= min_base_qual:
                    if base =='.': # in pileup, "." and "," were used to represent ref allele.
                        counts['REF'] += 1
                    elif base==',':
                        counts['ref'] += 1
                    else:
                        counts[base] +=1
                
                i += 1
        return counts



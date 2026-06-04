
from collections import defaultdict
from pathlib import Path
from typing import Dict
import os
import re

from SpaceTracer.steps.base import BaseStep
from SpaceTracer.utils.logger import get_logger

model_name=__name__
logger = get_logger(model_name)

class PriorCalculator(BaseStep):
    """
    get the prior info from gnomAD
    """
    
    BASES = ['A', 'T', 'C', 'G']
    
    def get_inputs(self, context: Dict) -> Dict[str, str]:
        """input"""
        inputs = {
            'filter_mpileup_file': context.get('filter_mpileup_file')
        }

        return inputs
    
    def get_outputs(self,context: Dict) -> Dict[str, str]:
        """output"""
        genome_details=self.config['genome_details']
        species=genome_details['species']
        if species != "human":
            return {'prior_file': 0}
        else:
            return {
                'prior_file': os.path.join(self.work_dir, 'prior.txt')
            }
    
    def _run(self,context: Dict):
        genome_details=self.config['genome_details']
        species=genome_details['species']
        if species != "human":
            logger.info("No population allele fraction information could be provided, prior would be target as same.")

        else:
            gnomad_path=self.config['gnomad_path']
            query_file=self.get_inputs(context)['filter_mpileup_file']
            gnomAD_dict=gnomAD(gnomad_path).get_chrom_list_from_file()
            auto_chrom=self.genome_details['chromosomes']['autosomes']
            sex_chrom=self.genome_details['chromosomes']['sex_chromosomes']

            prior_dict=query_sites_for_prior(query_file, gnomAD_dict,auto_chrom+sex_chrom)

            out_prior_file=self.get_outputs(context)['prior_file']
            try:
                prior_dict_to_file(prior_dict,out_prior_file)
            except:
                raise


def prior_dict_to_file(prior_dict: dict,output_file: str):
    BASES = ['A', 'T', 'C', 'G']
    
    with open(output_file, 'w') as f_out:
        f_out.write("#chrom\tpos\tref\tprior_A\tprior_T\tprior_C\tprior_G\n")
        
        for pos_key, alt_afs in prior_dict.items():
            chrom, pos, ref = pos_key
            # priors = {base: float(alt_afs.get(base, 0.0)) for base in BASES}
            priors = {base: float(alt_afs.get(base, 0.0)) if alt_afs.get(base, 0.0) != '.' else 0.0 for base in BASES}

            other_sum = sum(priors[base] for base in BASES if base != ref)
            
            if other_sum > 1.0:
                scale = 0.99 / other_sum
                for base in BASES:
                    if base != ref:
                        priors[base] *= scale
                priors[ref] = 1.0 - sum(priors[base] for base in BASES if base != ref)
            else:
                priors[ref] = max(0.0, 1.0 - other_sum)
            
            f_out.write(f"{chrom}\t{pos}\t{ref}\t"
                        f"{priors['A']}\t{priors['T']}\t"
                        f"{priors['C']}\t{priors['G']}\n")
                
def query_sites_for_prior(query_file: str, db_dict: dict, auto_chrom: list):
    queries_by_chrom = {}
    with open(query_file) as f:
        for line in f:
            parts = line.strip().split()
            chrom = parts[0]
            pos_previous="" 
            if chrom[0]!="#" and chrom in auto_chrom:
                pos = int(parts[1])
                ref = parts[3]
                
                if chrom not in queries_by_chrom:
                    queries_by_chrom[chrom] = []
                if pos != pos_previous:
                    queries_by_chrom[chrom].append((pos, ref))
                    
                pos_previous=pos
        
    prior_dict=defaultdict(dict)
    for chrom, queries in queries_by_chrom.items():
        if not queries:
            continue
        
        query_set = set()
        for pos, ref in queries:
            key = (pos, ref)
            query_set.add(key)
            prior_dict[(chrom, pos, ref)] = {}
        
        db_file = db_dict[chrom]
        if not Path(db_file).exists():
            continue
        
        with open(db_file) as f_db:
            for line in f_db:
                parts = line.strip().split()
                if len(parts) < 6:
                    continue
                    
                db_pos = int(parts[2])
                db_ref = parts[3]
                db_alt = parts[4]
                db_af = parts[5]
                
                key = (db_pos, db_ref)
                
                if key in query_set:
                    pos_key = (chrom, db_pos, db_ref)
                    prior_dict[pos_key][db_alt] = db_af
        
        del query_set

    return prior_dict


class gnomAD:
    def __init__(self,input_gnomAD):
        self.input_gnomAD=input_gnomAD

    def get_chrom_list_from_file(self) -> Dict[str, str]:
        '''
        Get the chromosome name from the file path recording absolute paths of gnomAD files
        or from the directory containing the gnomAD files.
        Args:
            file_or_dir: the file path or the directory path
        
        Returns:
            {chromosome: file_path}
        '''
        path = Path(self.input_gnomAD)
        
        if path.is_file():
            suffix = path.suffix.lower()
            if suffix in ['.gz', '.bgz', '.vcf']:
                raise ValueError(f'You provided the wrong gnomAD info!',
                                f'Both file list(record all absolute paths of gnomAD files) and directory(contain gnomAD files) are supported.')
            else:
                return self._parse_file_list()
        
        elif path.is_dir():
            return self._parse_directory()
        
        else:
            raise FileNotFoundError(f"Path does not exist: {self.input_gnomAD}")

    def _parse_file_list(self) -> Dict[str, str]:
        """parse the file list containing gnomAD files"""
        gnomad_dict = {}
        
        with open(self.input_gnomAD) as f:
            for line in f:
                file_path = line.strip()
                if not file_path or file_path.startswith('#'):
                    continue
                
                chrom = self._extract_chromosome_from_filename(os.path.basename(file_path))
                if chrom:
                    gnomad_dict[chrom] = file_path
        
        return gnomad_dict

    def _parse_directory(self) -> Dict[str, str]:
        """find all gnomAD files under the directory"""
        gnomad_dict = {}
        path=Path(self.input_gnomAD)
        patterns = [
            "*gnomad*.vcf.gz",
            "*gnomad*.vcf.bgz",
            "*gnomad*.vcf",
            "gnomad.genomes.*.sites.*.vcf.bgz",
            "gnomad.genomes.*.sites.*.vcf.gz",
            "*gnomad*_genome_only_af_*.txt"
        ]
        
        for pattern in patterns:
            for file_path in path.glob(pattern):
                if file_path.is_file():
                    chrom = self._extract_chromosome_from_filename(file_path.name)
                    if chrom and chrom not in gnomad_dict:
                        gnomad_dict[chrom] = str(file_path)
        
        # if the system still cannot find the file, try again 
        if not gnomad_dict:
            for file_path in path.glob("*"):
                if file_path.is_file() and file_path.suffix in ['.gz', '.bgz', '.vcf']:
                    chrom = self._extract_chromosome_from_filename(file_path.name)
                    if chrom and chrom not in gnomad_dict:
                        gnomad_dict[chrom] = str(file_path)
        
        if not gnomad_dict:
            raise ValueError(f"No gnomAD files found in directory: {self.input_gnomAD}")
        
        return gnomad_dict


    def _extract_chromosome_from_filename(self,filename: str) -> str:
        """get the chromosome name from file name"""
        
        # the normal 
        chrom_patterns = [
            r'chr(\d+|X|Y|M|MT)',           # chr1, chrX, chrM
            r'chromosome_(\d+|X|Y|M|MT)',    # chromosome_1
            r'\.(\d+|X|Y|M|MT)\.',           # .1.
            r'_(\d+|X|Y|M|MT)_',             # _1_
        ]
        
        for pattern in chrom_patterns:
            match = re.search(pattern, filename, re.IGNORECASE)
            if match:
                chrom = match.group(0)
                return chrom
        
        return None



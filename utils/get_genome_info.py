import os
import subprocess
from pathlib import Path

from SpaceTracer.utils.logger import get_logger

model_name=__name__
logger = get_logger(model_name)

GENOME_CONFIGS = {
    # Human references
    'hg38': {  
        'species': 'human',
        'description': 'Human (hg38)',
        'chromosomes': {
            'autosomes': [f'chr{i}' for i in range(1, 23)],
            'sex_chromosomes': ['chrX', 'chrY'],
            'mitochondrial': ['chrM'],
            'contigs' : [],
            'length':{}
        }
    },
    
    'hg19': {
        'species': 'human',
        'description': 'Human (hg19)',
        'chromosomes': {
            'autosomes': [str(i) for i in range(1, 23)],  # 注意：无chr前缀
            'sex_chromosomes': ['X', 'Y'],
            'mitochondrial': ['MT', 'M'],
            'contigs' : [],
            'length':{}

        }
    },
    
    # Mouse references
    'mm10': {
        'species': 'mouse',
        'description': 'Mouse ()',
        'chromosomes': {
            'autosomes': [f'chr{i}' for i in range(1, 20)],
            'sex_chromosomes': ['chrX', 'chrY'],
            'mitochondrial': ['chrM'],
            'contigs' : [],
            'length':{}

        }
    },
    
    'mm39': {
        'species': 'mouse',
        'description': 'Mouse (mm39)',
        'chromosomes': {
            'autosomes': [f'chr{i}' for i in range(1, 20)],
            'sex_chromosomes': ['chrX', 'chrY'],
            'mitochondrial': ['chrM'],
            'contigs' : [],
            'length':{}

        }
    },

    # Rat reference
    'rn7': {
        'species': 'rat',
        'description': 'Rat (mRatBN7.2)',
        'chromosomes': {
            'autosomes': [f'chr{i}' for i in range(1, 21)], 
            'sex_chromosomes': ['chrX', 'chrY'],
            'mitochondrial': ['chrM'],
            'contigs' : [],
            'length':{}

        }
    }
}

class GenomeDetails:
    def __init__(self,genome, genome_fasta):
        self.genome=genome
        self.genome_fasta=genome_fasta
    
    # def _find_fai_file(self):
    #     """ find raw fai file """
    #     fasta_path=self.genome_fasta
    #     dir_name = os.path.dirname(fasta_path)
    #     file_name = os.path.basename(fasta_path)
    #     # 1. most common fai file     
    #     path1 = str(fasta_path) + '.fai'
    #     if os.path.exists(path1):
    #         logger.debug(f"Found fai: {path1}")
    #         return path1
        
    #     # if not, try other path
    #     suffixes_to_try = ['.fa', '.fasta', '.fna', '.fas']
    #     # 2. add fai
    #     for suffix in suffixes_to_try:
    #         if fasta_path.endswith(suffix):
    #             path2 = fasta_path+'.fai'
    #             if path2.exists():
    #                 logger.debug(f"Found fai: {path2}")
    #                 return path2
        
    #     # 3. remove .fa .fasta; and add .fai
    #     if fasta_path.endswith(('.fa', '.fasta')):
    #         path3 = os.path.join(dir_name, file_name + '.fai') 
    #         if path3.exists():
    #             logger.debug(f"Found fai: {path3}")
    #             return path3
        
    #     logger.error(f"No fai file found for {fasta_path}")
    #     return None
    
    def _find_fai_file(self):
        """find raw fai file"""
        fasta_path = Path(self.genome_fasta)
        
        path1 = Path(str(fasta_path) + '.fai')
        if path1.exists():
            logger.debug(f"Found fai: {path1}")
            return str(path1)
        
        suffixes_to_try = ['.fa', '.fasta', '.fna', '.fas']
        for suffix in suffixes_to_try:
            if fasta_path.suffix == suffix:
                path2 = Path(str(fasta_path) + '.fai')
                if path2.exists():
                    logger.debug(f"Found fai: {path2}")
                    return str(path2)
        
        if fasta_path.suffix in ('.fa', '.fasta'):
            path3 = fasta_path.with_suffix('.fai')
            if path3.exists():
                logger.debug(f"Found fai: {path3}")
                return str(path3)
        
        logger.error(f"No fai file found for {fasta_path}")
        return None
        
    def _standardize_chrom_name(self, chrom_name: str) -> str:
        """standardize chromosome name"""
        
        # remove other info
        chrom_name = chrom_name.split()[0]  # the first word
        chrom_name = chrom_name.split('|')[0]  # remove | and further word
        
        # standardize
        if chrom_name.upper() in ['M', 'MT']:
            return 'chrM'
        elif chrom_name.upper() == 'X':
            return 'chrX'
        elif chrom_name.upper() == 'Y':
            return 'chrY'
        elif chrom_name.startswith('chr'):
            return chrom_name  
        elif chrom_name.isdigit():
            return f'chr{chrom_name}'  
        
        return chrom_name  # if others, keep the raw description

    def _classify_chromosome_by_name(self,chrom_name: str) -> str:
        """classify chrom as autosome, mitochondria, sex_chromosome and other"""
        # standard name
        std_name = self._standardize_chrom_name(chrom_name)
        
        # judge the chromosome type 
        if std_name == 'chrM':
            return 'mitochondrial'
        elif std_name in ['chrX', 'chrY']:
            return 'sex_chromosome'
        elif std_name.startswith('chr') and std_name[3:].isdigit():
            return 'autosome'
        else:
            return 'other'

    def _get_chromosomes_from_fasta_idx(self):
        """ use samtools faidx to get chromosome info """
        fasta_file=self.genome_fasta
        fai_file = self._find_fai_file()
        if not fai_file:
            logger.info(f"Creating faidx index for {fasta_file}")
            cmd = f"samtools faidx {fasta_file}"
            try:
                subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True)
                fai_file = f"{fasta_file}.fai"
            except subprocess.CalledProcessError as e:
                logger.error(f"samtools faidx failed: {e.stderr}")
                raise RuntimeError(f"samtools faidx failed: {e.stderr.strip()}")

        chromosomes = {}
        with open(fai_file, 'r') as f:
            for line in f:
                # fai：chr_name, length, offset, linebases, linewidth
                parts = line.strip().split('\t')
                if len(parts) >= 2:
                    chrom_name = parts[0]
                    chrom_length = int(parts[1])
                    
                    # use standard_chrom_name to judge the chromosome type
                    standard_chrom_name = self._standardize_chrom_name(chrom_name)
                    chromosomes[chrom_name] = {
                        'length': chrom_length,
                        'type': self._classify_chromosome_by_name(standard_chrom_name)
                    }
        
        return chromosomes
    
    def _get_genome_details(self):
        # 1. from known info
        if self.genome in GENOME_CONFIGS.keys():
            genome_info=GENOME_CONFIGS[self.genome]
            chromosomes=self._get_chromosomes_from_fasta_idx()
            for chrom in chromosomes:
                if chromosomes[chrom]['type']=="other":
                    genome_info['chromosomes']['contigs'].append(chrom)

                genome_info['chromosomes']['length'][chrom]=chromosomes[chrom]['length']

        # 2. only from fai file
        else:
            other_chrom=[]
            genome_info={'species': self.genome,
                'description': 'Other species',
                'chromosomes': {
                    'autosomes': 'auto',  
                    'sex_chromosomes': [],  
                    'mitochondrial': [],
                    'contigs': [],
                    'length': {}
                }
            }
            chromosomes=self._get_chromosomes_from_fasta_idx()
            for chrom in chromosomes:
                genome_info['chromosomes']['length'][chrom]=chromosomes[chrom]['length']

                if chromosomes[chrom]['type']=='autosome':
                    genome_info['chromosomes']['autosomes'].append(chrom)

                elif chromosomes[chrom]['type']=='sex_chromosome':
                    genome_info['chromosomes']['sex_chromosomes'].append(chrom)

                elif chromosomes[chrom]['type']=='mitochondrial':
                    genome_info['chromosomes']['mitochondrial'].append(chrom)
                    
                else:
                    genome_info['chromosomes']['contigs'].append(chrom)
                    
                    other_chrom.append(chrom)

            if genome_info['chromosomes']['autosomes']==[]:
                genome_info['chromosomes']['autosomes']=other_chrom
        return genome_info


def get_chr_size(fai_file):
    chr_sizes=dict()
    for line in open(fai_file,"r"):
        line=line.rstrip()
        fields=line.split('\t')
        chr_sizes[fields[0]]=chr_sizes.get(fields[0],fields[1])

    return chr_sizes

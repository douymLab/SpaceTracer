#!/usr/bin/env python3
"""
Split Mpileup Step - Split mpileup file when the file is so large 
(this will not work when you use the one-step "run" function, but when you run separately)
"""

import os
from pathlib import Path
from typing import Dict
import json

from collections import defaultdict
from SpaceTracer.utils.logger import get_logger
from SpaceTracer.utils.utils import check_dir
model_name=__name__
logger = get_logger(model_name)

class SplitMpileupStep():
    """
    split filtered_mpileup file
    
    1. check filtered_mpileup file size
    2. if file length >= threshold, split file
    3. split by chrom[chrM will be specific handle]
    
    Parameters:
        input_file: the input mpileup file
        output_dir: the output dir
        split_manifest: the manifest file recoding split information
        genome_details: the genome information
        chrom_chunk_size: the chunk size for auto chromosome (default:5000)
        chrM_chunk_size: the chunk size for mitochondrial (default:100)
    """
    def __init__(self, input_file: str, output_dir: str, split_manifest:str,genome_details: Dict, chrom_chunk_size:int, chrM_chunk_size:int, config: dict):
        self.input_file = input_file
        self.output_dir =output_dir
        self.manifest_file=split_manifest
        self.genome_details=genome_details
        self.config=config
        self.chrom_chunk_size=chrom_chunk_size
        self.chrM_chunk_size=chrM_chunk_size

    def get_step_config(self) -> Dict:
        return self.config.get('steps', {}).get('mpileup', {})
    
    def _run(self):
        """do split"""
        try:
            split_files, manifest = self._split_by_chromosome_and_chunk(
                    chrom_chunk_size=self.chrom_chunk_size,
                    chrM_chunk_size=self.chrM_chunk_size
                )
                
            logger.info(f"Split into {len(split_files)} files")
                
            # 3. save manifest
            manifest_file = self.manifest_file
            save_manifest(manifest, manifest_file)
            return True,None
            
        except Exception as e:
            error_msg = f"Error in filter mpileup (steps.mpileup mp_handle.filter_pileup): \n{e}"
            return False, error_msg

    def _split_by_chromosome_and_chunk(self,chrom_chunk_size,chrM_chunk_size):
        mpileup_file=self.input_file
        genome_info=self.genome_details
        split_dir = os.path.join(self.output_dir,"split_for_candidates")
        check_dir(split_dir)
        chrom_config = genome_info['chromosomes']
        autosomes = set(chrom_config['autosomes'])
        sex_chromosomes = set(chrom_config['sex_chromosomes'])
        mitochondrial = set(chrom_config['mitochondrial'])
        
        handlers = {} 
        split_files = []
        manifest = {
            'original_file': str(mpileup_file),
            'genome_info': genome_info,
            'chromosome_groups': defaultdict(dict)
        }
        
        def get_handler(chrom: str):
            """ build a handler for split """
            if chrom not in handlers:
                # check the chunk size for each chrom
                if chrom in mitochondrial:
                    chrom_type = 'mitochondrial'
                    chunk_size = int(chrM_chunk_size)
                elif chrom in sex_chromosomes:
                    chrom_type = 'sex_chromosome'
                    chunk_size = int(chrom_chunk_size)
                elif chrom in autosomes:
                    chrom_type = 'autosome'
                    chunk_size = int(chrom_chunk_size)
                else:
                    logger.error(f'There may be something wrong in display genome_info: {self.genome_info}' )
                
                handlers[chrom] = {
                    'display_name': chrom,
                    'chrom_type': chrom_type,
                    'chunk_size': chunk_size,
                    'count': 0,
                    'chunk_idx': 0,
                    'files': [],
                    'file_obj': None
                }
                
                # open the first file for each chrom
                file_path = os.path.join(split_dir ,f"{chrom}_chunk0000.mpileup")
                handlers[chrom]['file_obj'] = open(file_path, 'w')
                handlers[chrom]['files'].append(str(file_path))
                split_files.append(file_path)
            
            return handlers[chrom]
        
        # open the file
        with open(mpileup_file, 'r') as f:
            for line_num, line in enumerate(f, 1):
                if line.startswith('#'):
                    continue
                
                fields = line.strip().split('\t')
                if len(fields) < 3:
                    continue
                
                chrom = fields[0]
                
                handler = get_handler(chrom)
                handler['count'] += 1
                
                # do we need new file for each chrom
                if handler['count'] > handler['chunk_size']:
                    handler['count'] = 1
                    handler['chunk_idx'] += 1
                    
                    # close old split file
                    handler['file_obj'].close()
                    
                    # open new split file
                    file_path = os.path.join(split_dir, f"{handler['display_name']}_chunk{handler['chunk_idx']:04d}.mpileup")
                    handler['file_obj'] = open(file_path, 'w')
                    handler['files'].append(str(file_path))
                    split_files.append(file_path)
                
                # write file
                handler['file_obj'].write(line)
        
        # close all file and renew the manifest
        for chrom, handler in handlers.items():
            if handler['file_obj']:
                handler['file_obj'].close()
            
            manifest['chromosome_groups'][chrom] = {
                'display_name': handler['display_name'],
                'chrom_type': handler['chrom_type'],
                'total_records': handler['count'] + handler['chunk_idx'] * handler['chunk_size'],
                'num_chunks': handler['chunk_idx'] + 1,
                'chunk_size': handler['chunk_size'],
                'files': handler['files']
            }
        
        return split_files, manifest


def save_manifest(manifest: Dict, output_manifest: Path):
    """save manifest into JSON file"""
    with open(output_manifest, 'w') as f:
        json.dump(manifest, f, indent=2)
    
    logger.info(f"Manifest saved to: {output_manifest}")
    

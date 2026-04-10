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
    def __init__(self, input_file: str, output_dir: str, split_manifest:str,genome_details: Dict, config: dict):
        self.input_file = input_file
        self.output_dir =output_dir
        self.manifest_file=split_manifest
        self.genome_details=genome_details
        self.config=config


    def get_step_config(self) -> Dict:
        return self.config.get('steps', {}).get('mpileup', {})
    
    def _run(self,chrom_chunk_size,chrM_chunk_size,read_len,max_cost):
        """do split"""
        try:
            split_files, manifest = self._split_by_chromosome_and_chunk(
                    chrom_chunk_size=chrom_chunk_size,
                    chrM_chunk_size=chrM_chunk_size,
                    read_len=read_len,
                    max_cost=max_cost
                )
                
            logger.info(f"Split into {len(split_files)} files")
                
            # 3. save manifest
            manifest_file = self.manifest_file
            save_manifest(manifest, manifest_file)
            return True,None
            
        except:
            raise 

    def _split_by_chromosome_and_chunk(self,chrom_chunk_size,chrM_chunk_size,read_len,max_cost=1000):
        mpileup_file=self.input_file
        genome_info=self.genome_details
        split_dir = os.path.join(self.output_dir,"split_for_candidates")
        check_dir(split_dir)
        chrom_config = genome_info['chromosomes']

        autosomes = set(chrom_config['autosomes'])
        sex_chromosomes = set(chrom_config['sex_chromosomes'])
        mitochondrial = set(chrom_config['mitochondrial'])
        contigs = set(chrom_config['contigs'])
        
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
                elif chrom in contigs:
                    chrom_type = 'contig'
                    chunk_size = int(chrom_chunk_size)

                else:
                    logger.error(f'There may be something wrong in find {chrom} in genome_info: {genome_info}' )
                    # raise 
                
                # handlers[chrom] = {
                #     'display_name': chrom,
                #     'chrom_type': chrom_type,
                #     'chunk_size': chunk_size,
                #     'count': 0,
                #     'chunk_idx': 0,
                #     'files': [],
                #     'file_obj': None,
                #     'chunk_start_pos': None,
                #     'last_pos': None,
                #     'max_depth': 0
                # }
                handlers[chrom] = {
                    'display_name': chrom,
                    'chrom_type': chrom_type,
                    'chunk_size': chunk_size,
                    'chunk_idx': 0,
                    'files': [],
                    'chunks': [],
                    'file_obj': None,

                    # 当前 chunk 统计
                    'count': 0,
                    'chunk_start_pos': None,
                    'last_pos': None,
                    'max_depth': 0,
                    'depth_sum': 0,
                }

                
                # open the first file for each chrom
                file_path = os.path.join(split_dir ,f"{chrom}_chunk0000.mpileup")
                handlers[chrom]['file_obj'] = open(file_path, 'w')
                handlers[chrom]['files'].append(str(file_path))
                split_files.append(file_path)
            
            return handlers[chrom]
        
        # open the file
        # with open(mpileup_file, 'r') as f:
        #     for line_num, line in enumerate(f, 1):
        #         if line.startswith('#'):
        #             continue

        #         fields = line.strip().split('\t')
        #         if len(fields) < 9:
        #             continue

        #         chrom = fields[0]
        #         pos = int(fields[1])
        #         depth = int(fields[5])

        #         # 1st: split by chrom
        #         handler = get_handler(chrom)

        #         # 如果这是该 chrom 当前 chunk 的第一条记录
        #         if handler['chunk_start_pos'] is None:
        #             next_count = 1
        #             next_span = 1
        #             next_max_depth = depth
        #         else:
        #             next_count = handler['count'] + 1
        #             next_span = pos - handler['chunk_start_pos'] + 1
        #             next_max_depth = max(handler['max_depth'], depth)

        #         # 成本估计
        #         next_cost = next_max_depth * (next_span / read_len)

        #         # 2nd and 3rd: split by variants number and cost
        #         # 注意：只有当前 chunk 已经有内容时，才需要先切再写当前行
        #         if handler['count'] > 0 and (
        #             next_count > handler['chunk_size'] or next_cost > max_cost
        #         ):
        #             finalize_chunk(handler,read_len)
        #             handler['file_obj'].close()
        #             handler['chunk_idx'] += 1

        #             file_path = os.path.join(
        #                 split_dir,
        #                 f"{handler['display_name']}_chunk{handler['chunk_idx']:04d}.mpileup"
        #             )
        #             handler['file_obj'] = open(file_path, 'w')
        #             handler['files'].append(str(file_path))
        #             split_files.append(file_path)

        #             # 新 chunk 从当前行开始
        #             handler['count'] = 0
        #             handler['chunk_start_pos'] = None
        #             handler['last_pos'] = None
        #             handler['max_depth'] = 0

        #             # 重新计算当前行作为新 chunk 第一条记录时的状态
        #             next_count = 1
        #             next_span = 1
        #             next_max_depth = depth

        #         # 写当前行
        #         handler['file_obj'].write(line)

        #         # 更新当前 chunk 状态
        #         if handler['chunk_start_pos'] is None:
        #             handler['chunk_start_pos'] = pos

        #         handler['count'] = next_count
        #         handler['last_pos'] = pos
        #         handler['max_depth'] = next_max_depth
                
        # # close all file and renew the manifest
        # for chrom, handler in handlers.items():
        #     if handler['file_obj']:
        #         handler['file_obj'].close()
            
        #     manifest['chromosome_groups'][chrom] = {
        #         'display_name': handler['display_name'],
        #         'chrom_type': handler['chrom_type'],
        #         'total_records': sum(x['records'] for x in handler['chunks']),
        #         'num_chunks': len(handler['chunks']),
        #         'chunk_size': handler['chunk_size'],
        #         'files': handler['files'],
        #         'chunks': handler['chunks'],
        #         'total_cost': sum(x['cost'] for x in handler['chunks']),
        #     }
        
        # return split_files, manifest
        with open(mpileup_file, 'r') as f:
            for line_num, line in enumerate(f, 1):
                if line.startswith('#'):
                    continue

                fields = line.strip().split('\t')
                if len(fields) < 9:
                    continue

                chrom = fields[0]
                pos = int(fields[1])
                depth = int(fields[5])

                # 1st: split by chrom
                handler = get_handler(chrom)

                # 如果这是该 chrom 当前 chunk 的第一条记录
                if handler['chunk_start_pos'] is None:
                    next_count = 1
                    next_span = 1
                    next_max_depth = depth
                else:
                    next_count = handler['count'] + 1
                    next_span = pos - handler['chunk_start_pos'] + 1
                    next_max_depth = max(handler['max_depth'], depth)

                # 成本估计
                next_cost = next_max_depth * (next_span / read_len)
                # 或者更稳：
                # next_cost = next_max_depth * max(1.0, next_span / read_len)

                # 2nd and 3rd: split by variants number and cost
                if handler['count'] > 0 and (
                    next_count > handler['chunk_size'] or next_cost > max_cost
                ):
                    finalize_chunk(handler, read_len)
                    handler['file_obj'].close()
                    handler['chunk_idx'] += 1

                    file_path = os.path.join(
                        split_dir,
                        f"{handler['display_name']}_chunk{handler['chunk_idx']:04d}.mpileup"
                    )
                    handler['file_obj'] = open(file_path, 'w')
                    handler['files'].append(str(file_path))
                    split_files.append(file_path)

                    # 新 chunk 从当前行开始
                    handler['count'] = 0
                    handler['chunk_start_pos'] = None
                    handler['last_pos'] = None
                    handler['max_depth'] = 0

                    # 重新计算当前行作为新 chunk 第一条记录时的状态
                    next_count = 1
                    next_span = 1
                    next_max_depth = depth

                # 写当前行
                handler['file_obj'].write(line)

                # 更新当前 chunk 状态
                if handler['chunk_start_pos'] is None:
                    handler['chunk_start_pos'] = pos

                handler['count'] = next_count
                handler['last_pos'] = pos
                handler['max_depth'] = next_max_depth

        # close all file and renew the manifest
        for chrom, handler in handlers.items():
            if handler['count'] > 0:
                finalize_chunk(handler, read_len)

            if handler['file_obj']:
                handler['file_obj'].close()

            manifest['chromosome_groups'][chrom] = {
                'display_name': handler['display_name'],
                'chrom_type': handler['chrom_type'],
                'total_records': sum(x['records'] for x in handler['chunks']),
                'num_chunks': len(handler['chunks']),
                'chunk_size': handler['chunk_size'],
                'files': handler['files'],
                'chunks': handler['chunks'],
                'total_cost': sum(x['cost'] for x in handler['chunks']),
            }

        return split_files, manifest


def save_manifest(manifest: Dict, output_manifest: Path):
    """save manifest into JSON file"""
    with open(output_manifest, 'w') as f:
        json.dump(manifest, f, indent=2)
    
    logger.info(f"Manifest saved to: {output_manifest}")
    

def finalize_chunk(handler,read_len):
    if handler['count'] == 0 or handler['chunk_start_pos'] is None or handler['last_pos'] is None:
        return

    span_bp = handler['last_pos'] - handler['chunk_start_pos'] + 1
    mean_depth = handler['depth_sum'] / handler['count'] if handler['count'] > 0 else 0.0
    cost = mean_depth * max(1.0, span_bp / read_len)

    chunk_meta = {
        'chunk_idx': handler['chunk_idx'],
        'file': handler['files'][-1],
        'records': handler['count'],
        'start_pos': handler['chunk_start_pos'],
        'end_pos': handler['last_pos'],
        'span_bp': span_bp,
        'max_depth': handler['max_depth'],
        'mean_depth': mean_depth,
        'cost': cost,
    }
    handler['chunks'].append(chunk_meta)

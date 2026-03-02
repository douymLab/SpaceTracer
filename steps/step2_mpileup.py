#!/usr/bin/env python3
"""
Mpileup Step - This step will help you to get candidate sites by mpileup
"""

from typing import Dict
import os
import subprocess
import pandas as pd

from SpaceTracer.steps.base import BaseStep
from SpaceTracer.utils.logger import get_logger

from SpaceTracer.cores.mpileup_01_handle import PileupHandle
from SpaceTracer.cores.mpileup_02_filter import FilterCandidatesStep
from SpaceTracer.cores.mpileup_03_split import SplitMpileupStep, save_manifest

model_name=__name__
logger = get_logger("[mpileup]: "+model_name)

class MpileupStep(BaseStep):
    """
    mpileup (The definition of name, context, output_dir, work_dir, config, step_dir can be found in base.py.)
    
    input:
        - in_filter_bam: the bam file after filtration
        - reference: the reference fasta file
        - regions_file (optional): the reference callable region, if not, all genome will be scanned

    output:
        - mpileup_file: mpileup out put file
    
    parameters:
        - min_depth: minimal depth (default: 30)
        - max_depth: maximum depth (default: 200000) 
        - min_mapq: minimal mapping quality (default: 0)
        - min_baseq: minimal base quality (default: 0)
        - excl-flags: exclude reads with specific flags (default: 0) # Note: We provide this parameter for flexibility, but 0 is strongly recommended
    
    """
    
    def get_inputs(self, context: Dict) -> Dict[str, str]:
        """input"""
        inputs = {
            'in_filter_bam': context.get('bam_file'),
            'reference': self.config.get('genome_fasta')
        }
        if context.get('regions_file'):
            inputs['regions'] = context.get('regions_file')
        return inputs
    
    def get_outputs(self,context: Dict) -> Dict[str, str]:
        """output"""
        return {
            'mpileup_file': os.path.join(self.step_dir, 'raw_mpileup.txt'),
            'filter_mpileup_file': os.path.join(self.step_dir, 'filter_mpileup.txt'),
            'manifest_file': os.path.join(self.work_dir, 'split_manifest.json')

        }
    
    def get_step_config(self) -> Dict:
        return self.config.get('steps', {}).get('mpileup', {})
        
    def _run(self, context: Dict) -> Dict:
        """
        run mpileup
        use samtools mpileup and python file to handle the mpileup result
        """
        # parameters:
        inputs=self.get_inputs(context)
        processed_bam = inputs['in_filter_bam']
        reference = inputs['reference']
        regions_file = context.get('regions_file')
        raw_output_file = self.get_outputs(context)['mpileup_file']
        filter_output_file=self.get_outputs(context)['filter_mpileup_file']

        step_config = self.get_step_config()
        enable_split = step_config.get('enable_split', True)
        split_threshold = step_config.get('split_threshold', 100000)
        self.chrom_chunk_size = step_config.get('chrom_chunk_size', 5000)
        self.chrM_chunk_size = step_config.get('chrM_chunk_size', 100)
        print("step_config_chunk",self.chrom_chunk_size,self.chrM_chunk_size)

        logger.debug(f"parameters: inputs-> {inputs}")
        try:
            min_depth=int(step_config.get('mpileup').get('min_depth'))
        except:
            min_depth=30

        samtools_cmd=self._build_samtools_mpileup(processed_bam, reference, regions_file, step_config)
        handle_finish=self._run_and_handle_mpileup_results(samtools_cmd, min_depth ,raw_output_file)
        
        filter_finish,filter_log_info=self._filter_mpileup_results(context)
        if not filter_finish:
            logger.error(filter_log_info, exc_info=True)
            # raise # RuntimeError(f"Processing failed: {filter_log_info}")
                
        # 1. get line number and decide whether to split files
        total_lines = self._count_lines(filter_output_file)
        should_split = enable_split and total_lines > split_threshold
        logger.info(f"Mpileup file has {total_lines:,} lines")

        should_split = enable_split and total_lines > split_threshold
        if should_split:
            split_finish,split_log_info=self._split_mpileup_results(context)
            if not split_finish:
                logger.error(split_log_info, exc_info=True)
                # raise # "Error in mpileup split"

        else:
            manifest={}
            manifest['chromosome_groups']={}
            manifest['chromosome_groups']['all']={}
            manifest['chromosome_groups']['all']['files']=filter_output_file
            manifest_file = os.path.join(self.work_dir , 'split_manifest.json')
            save_manifest(manifest,manifest_file)
        
    def _build_samtools_mpileup(self, bam_file, reference, regions_file, config):
        print("###############",config)
        max_depth = config.get('max_depth', 200000)
        min_mapq = config.get('min_mapq', 0)
        min_baseq = config.get('min_baseq', 0)
        exclude_flag =  config.get('exclude_flag', 0)
        
        cmd_parts = [
            'samtools', 'mpileup', '-s -B '
            f'-f {reference}',
            f'--max-depth {max_depth}',
            f'--min-MQ {min_mapq}',
            f'--min-BQ {min_baseq}',
            f'--excl-flags {exclude_flag}'  
        ]
        
        print(cmd_parts)
        if regions_file:
            cmd_parts.append(f'-l {regions_file}')
        
        cmd_parts.append(bam_file)
        cmd = ' '.join(str(p) for p in cmd_parts)
        return cmd

    def _run_and_handle_mpileup_results(self, samtools_cmd, min_depth ,output_file):
        """ 1st step run samtools mpileup """ 
        try:
            samtools_proc = subprocess.Popen(
                samtools_cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                bufsize=1,
                shell=True
            )
            
            mp_handle = PileupHandle(min_depth)
            with open(output_file, 'w') as output_stream:
                mp_handle.filter_pileup(samtools_proc.stdout, output_stream)
            
            samtools_proc.stdout.close()
            stderr = samtools_proc.stderr.read()
            samtools_proc.wait()
            
            if samtools_proc.returncode != 0:
                raise

            return True, None
        
        except Exception as e:
            raise
    
    def _filter_mpileup_results(self,context):
        mp_filter=FilterCandidatesStep(self.get_outputs(context)['mpileup_file'],self.get_outputs(context)['filter_mpileup_file'], self.config)
        finish,info=mp_filter._run()
        return finish,info

    def _split_mpileup_results(self,context):
        input_file=self.get_outputs(context)['filter_mpileup_file']
        genome_details=self.genome_details
        manifest_file = os.path.join(self.work_dir , 'split_manifest.json')
        mp_split=SplitMpileupStep(input_file,self.work_dir,manifest_file,genome_details,self.chrom_chunk_size,self.chrM_chunk_size, self.config)
        finish,info=mp_split._run()
        return finish,info
        
    #### other functions
    def _load_regions(self, regions_file):
        """load regions of the bed file"""
        df = pd.read_csv(regions_file, sep='\t', header=None,
                        names=['chrom', 'start', 'end'])
        
        return [(row['chrom'], row['start'], row['end']) 
                for _, row in df.iterrows()]
    
    
    def _count_lines(self, file_path: str) -> int:
        """count file line number"""
        result = subprocess.run(
            ['wc', '-l', str(file_path)],
            capture_output=True,
            text=True,
            check=True
        )
        line_count = int(result.stdout.strip().split()[0])
        return line_count

#!/usr/bin/env python3
"""
Filter Candidates Step - filter and define the candidate errors and mosaics
"""

import subprocess
from typing import Dict

from SpaceTracer.utils.logger import get_logger
model_name=__name__
logger = get_logger(model_name)

class FilterCandidatesStep:
    """
    Filter mpileup result:
    
    The format of raw_mpileup file:
        chrom pos ref mp_depth info_depth ref_allele ref_count 
        alt1 alt1_count alt2 alt2_count alt3 alt3_count
    
    The format of filtered_mpileup file:
        chrom pos type ref alt1 total_depth ref_depth alt1_depth vaf
    
    Rule of filtrationL
        candidate_error: 
            - info_depth >= 100
            - ref_count / info_depth > 0.95
            - alt1_count >= 5
        
        candidate_somatic:
            - alt1_count >= 5
            - info_depth >= 30
            - 0.001 <= alt1_count / info_depth <= 0.6
    
    Params:
        # candidate_error 
        error_min_depth: 100
        error_min_vaf: 0.95
        error_min_alt_count: 5
        
        # candidate_somatic 
        somatic_min_depth: 30
        somatic_min_alt_count: 5
        somatic_min_vaf: 0.001
        somatic_max_vaf: 0.6
    """
    def __init__(self, input_file: str, output_file: str,config: dict):
        self.input_file = input_file
        self.output_file = output_file
        self.config=config

    def get_step_config(self) -> Dict:
        return self.config.get('steps', {}).get('mpileup', {})
    
    def _run(self):
        """run filter"""
        raw_mpileup = self.input_file
        step_config = self.get_step_config()
        
        # output file
        filtered_file = self.output_file
        
        # parameters
        params = {
            # for candidate_error
            'error_min_depth': step_config.get('error_min_depth', 100),
            'error_min_vaf': step_config.get('error_min_vaf', 0.95),
            'error_min_alt_count': step_config.get('error_min_alt_count', 5),
            
            # for candidate_somatic
            'somatic_min_depth': step_config.get('somatic_min_depth', 30),
            'somatic_min_alt_count': step_config.get('somatic_min_alt_count', 5),
            'somatic_min_vaf': step_config.get('somatic_min_vaf', 0.001),
            'somatic_max_vaf': step_config.get('somatic_max_vaf', 0.6),
        }
        
        logger.debug(f"Filtering candidates with params: {params}")
        
        finish,error_info = self._process_with_awk(raw_mpileup, filtered_file, params)
        return finish,error_info

    
    def _process_with_awk(self, input_file: str, output_file: str, params: dict) -> bool:
        """filter the raw mpileup by awk"""
        error_depth = params['error_min_depth']
        error_vaf = params['error_min_vaf']
        error_count = params['error_min_alt_count']
        
        somatic_depth = params['somatic_min_depth']
        somatic_count = params['somatic_min_alt_count']
        somatic_min_vaf = params['somatic_min_vaf']
        somatic_max_vaf = params['somatic_max_vaf']
        
        awk_script = f'''
        BEGIN {{
            OFS = "\\t";
            error_depth = {error_depth};
            error_vaf = {error_vaf};
            error_count = {error_count};
            somatic_depth = {somatic_depth};
            somatic_count = {somatic_count};
            somatic_min_vaf = {somatic_min_vaf};
            somatic_max_vaf = {somatic_max_vaf};
            
            print "#chrom", "pos", "type", "ref", "alt1", 
                "total_depth", "ref_depth", "alt1_depth", "vaf";
        }}
        
        !/^#/ {{
            chrom = $1;
            pos = $2;
            ref = $3;
            total_depth = $5 ;   # info_depth
            ref_count = $7 ;
            alt1 = $8;
            alt1_count = $9 ;
            
            # calculate vaf
            if (total_depth > 0) {{
                vaf = alt1_count / total_depth;
                ref_vaf = ref_count / total_depth;

                # check candidate error
                is_error = (ref!="N" && 
                        total_depth >= error_depth && 
                        ref_vaf >= error_vaf && 
                        alt1_count >= error_count);
                
                # check candidate somatic
                is_somatic = (ref!="N" && 
                            alt1_count >= somatic_count && 
                            total_depth >= somatic_depth && 
                            vaf >= somatic_min_vaf && 
                            vaf <= somatic_max_vaf);
                
                if (is_error && is_somatic) {{
                    type = "candidate_error,candidate_somatic";
                }} else if (is_error) {{
                    type = "candidate_error";
                }} else if (is_somatic) {{
                    type = "candidate_somatic";
                }} else {{
                    next;  # skip not output
                }}
                
                # output
                print chrom, pos, type, ref, alt1, 
                    total_depth, ref_count, alt1_count, vaf;
            }}
        }}
        '''
        
        try:
            # run awk 
            cmd = ['awk', awk_script, input_file]
            # logger.debug(f"awk command {cmd}")
            
            with open(output_file, 'w') as f_out:
                result = subprocess.run(
                    cmd,
                    stdout=f_out,
                    stderr=subprocess.PIPE,
                    text=True,
                    encoding='utf-8'
                )
            
            # chech output
            if result.returncode != 0:
                error_msg = f"{result.returncode}"
                if result.stderr:
                    error_msg += f"{result.stderr[:500]}"
                return False, error_msg
            
            # debug out
            if result.stderr:
                logger.debug(f"awk output: {result.stderr.strip()}")
            
            return True, None
            
        except FileNotFoundError:
            error_info="where is awk???"
            return False,error_info
        except Exception as e:
            return False, e


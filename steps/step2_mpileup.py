#!/usr/bin/env python3
"""
Mpileup Step - This step will help you to get candidate sites by mpileup
"""

from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
import shutil
from typing import Dict
import os
import subprocess
import pandas as pd

from SpaceTracer.steps.base import BaseStep
from SpaceTracer.utils.logger import get_logger
from SpaceTracer.utils.get_read_level_feature import detect_read_length

from SpaceTracer.cores.mpileup_01_handle import PileupHandle
from SpaceTracer.cores.mpileup_02_filter import FilterCandidatesStep
from SpaceTracer.cores.mpileup_03_split import SplitMpileupStep

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
            'db_path': os.path.join(self.work_dir, 'split_chunk.db')

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
        regions_file = context.get('regions_file','')
        raw_output_file = self.get_outputs(context)['mpileup_file']
        filter_output_file=self.get_outputs(context)['filter_mpileup_file']

        step_config = self.get_step_config()
        enable_split = step_config.get('enable_split', True)
        split_threshold = step_config.get('split_threshold', 100000)
        self.chrom_chunk_size = step_config.get('chrom_chunk_size', 5000)
        self.chrM_chunk_size = step_config.get('chrM_chunk_size', 100)
        print("step_config_chunk",self.chrom_chunk_size,self.chrM_chunk_size)

        self.readLen=detect_read_length(processed_bam)
        self.max_cost= 150 * 20000
        logger.debug(f"parameters: inputs-> {inputs}")
        try:
            min_depth=int(step_config.get('mpileup').get('min_depth'))
        except:
            min_depth=30

        ## samtools_cmd=self._build_samtools_mpileup(processed_bam, reference, regions_file, step_config)
        ## handle_finish=self._run_and_handle_mpileup_results(samtools_cmd, min_depth ,raw_output_file)
        handle_finish=self._run_and_merge_mpileup_results(processed_bam, reference, regions_file, raw_output_file, min_depth, step_config)
        
        filter_finish,filter_log_info=self._filter_mpileup_results(context)

        total_lines = self._count_lines(filter_output_file)
        should_split = enable_split and total_lines > split_threshold
        print("++++++++++++++++++++++++++++++++++++++++",should_split)
        logger.info(f"Mpileup file has {total_lines:,} lines")

        if should_split:
            split_finish=self._split_mpileup_results(context)

        else:
            split_finish = self._write_single_chunk_db(context, filter_output_file)
            # manifest={}
            # manifest['chromosome_groups']={}
            # manifest['chromosome_groups']['all']={}
            # manifest['chromosome_groups']['all']['files'] = [filter_output_file]
            # manifest_file = os.path.join(self.work_dir , 'split_manifest.json')
            # save_manifest(manifest,manifest_file)

            
    def _write_single_chunk_db(self, context, input_file):
        db_path = self.get_outputs(context)['db_path']

        mp_split = SplitMpileupStep(
            input_file=input_file,
            output_dir=self.work_dir,
            db_path=db_path,
            genome_details=self.genome_details,
            config=self.config
        )

        conn = mp_split._init_chunk_index_db(db_path)
        cur = conn.cursor()

        cur.execute("""
        INSERT INTO chunks (
            chunk_id, chrom, chrom_type, chunk_idx, chunk_file,
            start_pos, end_pos, records, span_bp, max_depth, mean_depth, cost
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """, (
            "all_chunk0000",
            "all",
            "whole_file",
            0,
            str(input_file),
            None,
            None,
            None,
            None,
            None,
            None,
            None,
        ))

        cur.execute("""
        INSERT INTO chroms (
            chrom, chrom_type, total_records, num_chunks, chunk_size, total_cost
        ) VALUES (?, ?, ?, ?, ?, ?)
        """, (
            "all",
            "whole_file",
            None,
            1,
            1,
            None,
        ))

        conn.commit()
        conn.close()

        return True

    def _split_bed_by_total_bases(self, regions_file, result_dir, n_parts):
        lines = []
        with open(regions_file, "r") as f:
            for line in f:
                if line.strip() == "" or line.startswith("#"):
                    continue
                parts = line.strip().split()
                if len(parts) < 3:
                    continue
                lines.append(line)

        if not lines:
            raise ValueError("region_file has no valid intervals")

        total = len(lines)
        n_parts = min(n_parts, total)  # 防止 n_parts 比行数还多
        chunk_size = (total + n_parts - 1) // n_parts  # 向上取整，保证覆盖所有行

        region_files = []
        for part_idx in range(n_parts):
            chunk = lines[part_idx * chunk_size : (part_idx + 1) * chunk_size]
            if not chunk:
                break
            bed_path = os.path.join(result_dir, f"part_{part_idx + 1}.bed")
            with open(bed_path, "w") as out:
                out.writelines(chunk)
            region_files.append(bed_path)

        return region_files

    def _make_bed_splits(self, chrom_length_dict,regions_file, window_size=100000000):
        split_regions_dir = os.path.join(self.step_dir, "split_regions")
        os.makedirs(split_regions_dir, exist_ok=True)

        if regions_file:
            # if provide the bed file
            n_parts = max(1, self.threads * 2)
            return self._split_bed_by_total_bases(regions_file, split_regions_dir,n_parts)
        
        else:
            region_files = []
            for chrom, length in chrom_length_dict.items():
                length = int(length)
                for start in range(0, length, window_size):
                    end = min(start + window_size, length)
                    bed_path = os.path.join(split_regions_dir, f"{chrom}_{start}_{end}.bed")
                    with open(bed_path, "w") as f:
                        f.write(f"{chrom}\t{start}\t{end}\n")
                    region_files.append(bed_path)

            return region_files

    def _prepare_regions_and_outputs(self, regions_file, config):
        split_results_dir = os.path.join(self.step_dir, "split_results")
        os.makedirs(split_results_dir, exist_ok=True)

        chrom_length_dict=self.genome_details['chromosomes']['length']
        # window_size = config.get("window_size", 100000000)
        region_files = self._make_bed_splits(chrom_length_dict,regions_file, window_size=100000000)

        out_files = []
        for bed in region_files:
            name = os.path.splitext(os.path.basename(bed))[0]
            out_files.append(os.path.join(split_results_dir, f"{name}.mpileup"))

        return region_files, out_files


    def _build_samtools_mpileup(self, bam_file, reference, regions_file, config):
        max_depth = config.get('max_depth', 200000)
        min_mapq = config.get('min_mapq', 0)
        min_baseq = config.get('min_baseq', 0)
        exclude_flag = config.get('exclude_flag', 0)

        cmd_parts = [
            'samtools', 'mpileup',
            '-s', '-B',
            '-f', reference,
            '--max-depth', str(max_depth),
            '--min-MQ', str(min_mapq),
            '--min-BQ', str(min_baseq),
            '--excl-flags', str(exclude_flag)
        ]

        if regions_file:
            cmd_parts += ['-l', regions_file]

        cmd_parts.append(bam_file)
        return cmd_parts


    def _run_mpileup_one(self, cmd, min_depth, out_file):
        samtools_proc = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            bufsize=-1,
            shell=True
        )
        mp_handle = PileupHandle(min_depth)
        with open(out_file, 'w') as output_stream:
            mp_handle.filter_pileup(samtools_proc.stdout, output_stream)
        samtools_proc.stdout.close()
        stderr = samtools_proc.stderr.read()
        samtools_proc.wait()
        if samtools_proc.returncode != 0:
            raise RuntimeError(stderr)
        return out_file


    def _run_mpileup_parallel(self, bam, reference, region_files, min_depth, out_files, config):
        with ProcessPoolExecutor(max_workers=self.threads) as exe:
            futures = []
            for reg, out in zip(region_files, out_files):
                cmd = self._build_samtools_mpileup(bam, reference, reg, config)
                futures.append(exe.submit(self._run_mpileup_one, cmd, min_depth, out))

            for f in as_completed(futures):
                f.result()

    def _merge_mpileup_files(self, out_files, merged_file):
        os.makedirs(os.path.dirname(merged_file), exist_ok=True)
        with open(merged_file, "w") as w:
            for f in out_files:
                with open(f, "r") as r:
                    shutil.copyfileobj(r, w)
        return merged_file

    def _run_and_merge_mpileup_results(self, bam, reference, regions_file, output_file, min_depth, config):
        try:
            region_files, out_files = self._prepare_regions_and_outputs(regions_file, config)
            self._run_mpileup_parallel(bam, reference, region_files, min_depth, out_files, config)
            merged = output_file
            self._merge_mpileup_files(out_files, merged)
            return True

        except:
            raise

        
    def _build_samtools_mpileup(self, bam_file, reference, regions_file, config):
        # print("###############",config)
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
        
        if regions_file:
            cmd_parts.append(f'-l {regions_file}')
        # print(cmd_parts)
        
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
        db_path=self.get_outputs(context)['db_path']
        genome_details=self.genome_details
        # db_path = os.path.join(self.work_dir , 'split_chunk.db')
        mp_split=SplitMpileupStep(input_file,self.work_dir,db_path,genome_details, self.config)
        mp_split._run(self.chrom_chunk_size,self.chrM_chunk_size,self.readLen,self.max_cost)

        return True
        
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

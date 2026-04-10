#!/usr/bin/env python3
"""
bam process: here we need to filter the bam
"""

import os
from pathlib import Path
import glob
import subprocess
from typing import Dict
from SpaceTracer.steps.base import BaseStep
from SpaceTracer.utils.logger import get_logger
from SpaceTracer.utils.read_files import handle_barcode
from SpaceTracer.utils.utils import check_dir

model_name=__name__
logger = get_logger(model_name)

class BamProcessingStep(BaseStep):
    """
    bam_processing (The definition of name, context, output_dir, work_dir, config, step_dir can be found in base.py.)
    
    input:
        - raw_bam: the raw bam file input
        - barcode_file (optional): the barcode you want to keep in the bam file. 
    
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
        """ input """
        inputs = {
            'raw_bam': context.get('bam_file')
        }
        return inputs
    

    def optional_parameters(self, context: Dict) -> Dict[str, str]:
        """ That's optional parameters """
        parameters={}
        if self.config.get('tissue_position'):
            parameters['tissue_position'] = self.config.get('tissue_position')
        else:
            parameters['tissue_position'] = None
        return parameters
            

    def get_outputs(self,context: Dict) -> Dict[str, str]:
        """output"""
        return {
            'in_bam': os.path.join(self.step_dir, 'IN.bam'),
            'in_filter_bam': os.path.join(self.step_dir, 'IN_filter.bam')
        }
    
    def get_step_config(self) -> Dict:
        return self.config.get('steps', {}).get('bam_processing', {})
    

    def _run(self, context: Dict):
        """
        run bam processing
        use sinto and samtools 
        """
        # parameters:
        inputs=self.get_inputs(context)
        parameters=self.optional_parameters(context)
        raw_bam = inputs['raw_bam']
        
        barcode_file = parameters['tissue_position']
        barcode_key=self.config.get('barcode_key')
        threads=self.threads

        step_config=self.get_step_config()
        nm_threshold=step_config["nm_threshold"]
        mapq_threshold=step_config["mapq_threshold"]

        outputs=self.get_outputs(context)
        in_bam=outputs['in_bam']
        in_filter_bam=outputs['in_filter_bam']
        print(raw_bam,self.step_dir)

        if barcode_file:
            in_barcode_file=os.path.join(self.step_dir,"in_barcode.txt")
            barcode_dict=handle_barcode(barcode_file)
            with open(in_barcode_file,"w") as f1:
                for barcode in barcode_dict.keys():
                    f1.write(f'{barcode}\tIN\n')
            
            finish=self._run_sinto_filterbarcodes(raw_bam,in_barcode_file,barcode_key,self.step_dir)

        else:
            link_finish=self._link_raw_bam(raw_bam,self.step_dir,threads)
        
        filter_finish=self._filter_bam(in_bam, in_filter_bam, threads,nm_threshold, mapq_threshold)
    
            
    def _run_sinto_filterbarcodes(self, bam_file: str, barcode_file: str, barcode_key:str,
                                output_dir: str) -> bool:
        """
        run sinto filterbarcodes commands 

        bam_file: input bam file 
        barcode_file: the in_barcode_file
        barcode_key: for visium data, the key word record barcode info is "CB"
        output_dir: save dir 
        """
        threads = self.threads
        output_dir = str(output_dir)

        if os.path.exists(output_dir):
            for file_path in glob.glob(os.path.join(output_dir, 'IN*')):
                logger.warning(f"Removing existing file: {file_path}")
                os.remove(file_path)
                
        cmd = [
            'sinto', 'filterbarcodes',
            '-b', str(bam_file),
            '-c', str(barcode_file),
            '--barcodetag', str(barcode_key),
            '--outdir', output_dir,
            '-p', str(threads)
        ]
        logger.debug(f"run sinto: {' '.join(cmd)}")

        result = subprocess.run(
            cmd,  
            capture_output=True,
            shell=False,  
            encoding='utf-8',
            check=False  
        )

        cmd_str=' '.join(cmd)
        if result.returncode != 0:
            raise RuntimeError(
                    f"samtools mpileup failed with command: {cmd_str}. Please try to rm -rf {output_dir}, and try again.",)

        return True

    def _link_raw_bam(self, raw_bam: str,output_dir: str, threads: int):
        """
        If the barcode file not provided, all reads will be used as in_tissue reads, so link the raw file as in_bam_file
        raw_bam: raw bam file   
        output_dir: save dir
        """
        logger.debug(f"link for raw bam: {raw_bam}")
        output_dir=Path(output_dir)
        raw_bam=Path(raw_bam)
        in_bam = output_dir / "IN.bam"

        raw_bam_abs = raw_bam.resolve()
        
        if in_bam.exists() or in_bam.is_symlink():
            in_bam.unlink()
        os.symlink(raw_bam_abs, in_bam)
        logger.debug(f"link bam file: {raw_bam_abs} -> {in_bam}")
        
        possible_index_extensions = ['.bai', '.csi', raw_bam.suffix + '.bai']
        
        for ext in possible_index_extensions:
            index_file = raw_bam.with_suffix(ext)
            if index_file.exists():
                linked_index = output_dir / "IN.bam.bai"
                if linked_index.exists() or linked_index.is_symlink():
                    linked_index.unlink()
                os.symlink(index_file.resolve(), linked_index)
                logger.debug(f"link bai file: {index_file} -> {linked_index}")
        
        bam_index=self._check_index_bam(raw_bam,threads)
        if not bam_index:
            raise RuntimeError(
                "BAM index file cannot be copied or reindexed. "
            )

        if not in_bam.exists():
            raise FileNotFoundError(
                f"Linked file does not exist: {in_bam}. "
            )
        
        if in_bam.is_symlink():
            target = in_bam.resolve()
            if target != raw_bam_abs:
                raise OSError(f"Ops! Why the link bam is not equal with raw bam?")

        return True

    def _copy_file(self, raw_bam: Path,  target_path: Path) -> Path:
        import shutil
        try:
            shutil.copy2(raw_bam, target_path)
        except Exception as e:
            raise

    def _check_index_bam(self, bam_file: str, threads: int):
        """
        Index BAM file using samtools 
        bam_file: input bam file to index
        """
        bam_file = Path(bam_file)

        possible_index = [
            bam_file.with_suffix('.bam.bai'),
            bam_file.with_suffix('.bai')
        ]
        
        for index_file in possible_index:
            if index_file.exists():
                logger.debug(f"Index file created: {index_file}")
                return True

        cmd = [
            'samtools', 'index',
            '-@', str(threads),
            str(bam_file)
        ]
        
        logger.debug(f"run samtools index: {' '.join(cmd)}")
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            encoding='utf-8',
            check=False
        )
        
        if result.returncode != 0:
            raise RuntimeError(f'Samtools index is wrong with command: {cmd}')

        return True

    def _filter_bam(self,bam_file: str, output_bam: str, threads: int, nm_threshold: int, mapq_threshold: int ):
        """
        Run samtools view to filter BAM file by NM tag and MAPQ score
        
        bam_file: input bam file
        step_dir: output directory for filtered bam
        nm_threshold: maximum allowed NM (number of mismatches) value, default 5
        mapq_threshold: minimum MAPQ score, default 255
        context: context dict containing threads, etc.
        """
        
        cmd = [
            'samtools', 'view',
            '-e', f'"[nM] <= {nm_threshold}"',
            '-q', str(mapq_threshold),
            '-o', str(output_bam),
            '-@', str(threads),
            str(bam_file)
        ]
        
        logger.debug(f"run samtools view filter: {' '.join(cmd)}")
        
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            encoding='utf-8',
            check=False
        )
        
        if result.returncode != 0:
            raise RuntimeError(
                f"samtools view failed with command: {cmd}")
        
        if not Path(output_bam).exists():
            raise FileNotFoundError(
                f"Output BAM file was not created: {output_bam}")
        
        if Path(output_bam).stat().st_size == 0:
            logger.warning(f"Filtered BAM file is empty: {output_bam}")

        index_finish=self._check_index_bam(str(output_bam),threads)

        return True


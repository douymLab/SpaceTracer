#!/usr/bin/env python3
"""
Base Step Class - The Class object for each step 
(in this step, we define the raw structure of each step)
"""

import logging
from abc import ABC, abstractmethod
from pathlib import Path
import sys
from typing import Dict
import os

from SpaceTracer.utils.get_genome_info import GenomeDetails
from SpaceTracer.utils.utils import check_dir

logger = logging.getLogger(__name__)

class BaseStep(ABC):
    """
    Each step should has:
    1. _run() - core
    2. get_inputs() - input definition
    3. get_outputs() - output definition
    4. validate_inputs() - make sure input is validated (optional)
    5. validate_outputs() - make sure output is validated (optional)
    """
    
    def __init__(self, name: str, context: Dict, checkpoint_manager):
        """
        Initialization steps
        
        Args:
            name: The step name
            context: the context will be updated across pipeline
            checkpoint_manager: checkpoint manager
        """
        self.name = name
        self.context = context
        self.checkpoint = checkpoint_manager
        
        # get the basic context from contig
        self.config = context['config']
        self.work_dir = self.config.get('output_dir')
        self.threads=int(self.config.get('run').get('threads'))
        self.memory=self.config.get('run').get('memory')
        
        self.genome_details=self.config['genome_details']
        
        # 创建步骤专用目录
        self.step_dir = os.path.join(self.work_dir, self.name)
        check_dir(self.step_dir)
        
        logger.debug(f"Initialized step: {self.name}")
    
    def execute(self, context: Dict, skip_validation: bool = True) -> Dict:
        # validate input
        if not self.validate_inputs(context):
            raise ValueError(f"Step {self.name}: Input validation failed")
        
        # run func
        self._run(context)
        result_context = self.get_outputs(context)
        
        # validate output
        if not skip_validation and not self.validate_outputs(result_context):
            raise ValueError(f"Something wrong in {self.name}")

        # update context
        context.update(result_context)
        return context

    @abstractmethod
    def _run(self, context: Dict) -> Dict:
        pass
    
    def get_inputs(self, context: Dict) -> Dict[str, str]:
        return {}
    
    def get_outputs(self, context: Dict) -> Dict[str, str]:
        return {}
    
    def validate_inputs(self, context: Dict) -> bool:
        inputs = self.get_inputs(context)
        # logger.debug(f'###The inputs is :{inputs}')
        for input_name, input_path in inputs.items():
            if isinstance(input_path,int): # the int output was used to pass some step 
                return True
            
            elif isinstance(input_path,str) and not Path(input_path).exists():
                logger.error(f"Input file not found: {input_name} -> {input_path}")
                return False
        
        return True
    
    def validate_outputs(self, context: Dict) -> bool:
        outputs = self.get_outputs(context)
        
        for output_name, output_path in outputs.items():
            if isinstance(output_path,int): # the int output was used to pass some step 
                continue
            
            path = Path(output_path)
            if not os.path.exists(output_path):
                logger.error(f"Output file not created: {output_name} -> {output_path}")
                return False
            
            if path.stat().st_size == 0:
                logger.error(f"Output file is empty: {output_name} -> {output_path}")
                return False
        
        return True
    
    def load_outputs(self, context: Dict) -> Dict:
        outputs = self.get_outputs(context)
        
        for output_name, output_path in outputs.items():
            context[output_name] = output_path
        
        return context
    
    def run_command(self, cmd: str, description: str = None) -> int:
        import subprocess
        
        if description:
            logger.debug(f"Running: {description}")
        logger.debug(f"Command: {cmd}")
        
        result = subprocess.run(
            cmd,
            shell=True,
            capture_output=True,
            text=True
        )
        
        if result.returncode != 0:
            logger.error(f"Command failed with code {result.returncode}")
            logger.error(f"STDERR: {result.stderr}")
            raise RuntimeError(f"Command failed: {cmd}")
        
        return result.returncode
    
    def get_step_config(self) -> Dict:
        """get the step configurations from config file"""
        return self.config.get('steps', {}).get(self.name, {})

    def get_executor(self) -> str:
        step_cfg = self.get_step_config()
        return step_cfg.get("executor", self.config.get("executor", "internal"))



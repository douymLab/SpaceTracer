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
        
        # 从context获取基本信息
        self.config = context['config']
        self.work_dir = self.config.get('output_dir')
        self.threads=int(self.config.get('run').get('threads'))

        
        # 创建步骤专用目录
        self.step_dir = os.path.join(self.work_dir, self.name)
        check_dir(self.step_dir)
        
        logger.debug(f"Initialized step: {self.name}")
    
    def execute(self, context: Dict) -> Dict:
        """
        执行步骤
        
        Args:
            context: 输入上下文
        
        Returns:
            更新后的上下文
        """
        logger.info(f"[{self.name}] Executing...")
        # 验证输入
        if not self.validate_inputs(context):
            raise ValueError(f"Step {self.name}: Input validation failed")
        
        # 运行核心逻辑
        self._run(context)
        result_context = self.get_outputs(context)
        
        # 验证输出
        if not self.validate_outputs(result_context):
            logger.error(f"Something wrong in {self.name}")
            sys.exit()
        
        # 更新context
        context.update(result_context)
        
        return context
    
    @abstractmethod
    def _run(self, context: Dict) -> Dict:
        """
        核心处理逻辑（子类必须实现）
        
        Args:
            context: 输入上下文
        
        Returns:
            包含输出的字典
        """
        pass
    
    def get_inputs(self, context: Dict) -> Dict[str, str]:
        """
        定义输入文件
        
        Returns:
            {input_name: file_path} 字典
        """
        return {}
    
    def get_outputs(self, context: Dict) -> Dict[str, str]:
        """
        定义输出文件
        
        Returns:
            {output_name: file_path} 字典
        """
        return {}
    
    def validate_inputs(self, context: Dict) -> bool:
        """
        验证输入文件
        
        默认实现：检查所有输入文件是否存在
        """
        inputs = self.get_inputs(context)
        logger.debug(f'###The inputs is :{inputs}')
        for input_name, input_path in inputs.items():
            if not Path(input_path).exists():
                logger.error(f"Input file not found: {input_name} -> {input_path}")
                return False
        
        return True
    
    def validate_outputs(self, context: Dict) -> bool:
        """
        验证输出文件
        
        默认实现：检查所有输出文件是否存在且非空
        """
        outputs = self.get_outputs(context)
        
        for output_name, output_path in outputs.items():
            if isinstance(output_path,int):
                continue
            
            path = Path(output_path)
            # print(output_path)
            if not os.path.exists(output_path):
                logger.error(f"Output file not created: {output_name} -> {output_path}")
                return False
            
            if path.stat().st_size == 0:
                logger.error(f"Output file is empty: {output_name} -> {output_path}")
                return False
        
        return True
    
    def load_outputs(self, context: Dict) -> Dict:
        """
        从之前的运行中加载输出
        
        用于checkpoint恢复
        """
        outputs = self.get_outputs(context)
        
        # 将输出路径添加到context
        for output_name, output_path in outputs.items():
            context[output_name] = output_path
        
        return context
    
    def run_command(self, cmd: str, description: str = None) -> int:
        """
        运行外部命令
        
        Args:
            cmd: 命令字符串
            description: 命令描述
        
        Returns:
            返回码
        """
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
        """获取该步骤的配置"""
        return self.config.get('steps', {}).get(self.name, {})

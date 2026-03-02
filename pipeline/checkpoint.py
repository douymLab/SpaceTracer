#!/usr/bin/env python3
"""
Checkpoint Manager
"""

import json
import logging
from pathlib import Path
from datetime import datetime
from typing import Dict, Any, Optional, List

logger = logging.getLogger(__name__)

class CheckpointManager:
    """
    检查点管理器
    
    功能：
    1. 记录每个步骤的完成/失败状态
    2. 记录输入/输出文件信息
    3. 支持断点续跑
    4. 检测输入文件是否变化
    """
    
    def __init__(self, output_dir: Path, disabled: bool = False):
        """
        初始化检查点管理器
        
        Args:
            output_dir: 输出目录
            disabled: 是否禁用检查点（force模式）
        """
        self.output_dir = Path(output_dir)
        self.checkpoint_file = self.output_dir / '.pipeline_checkpoints.json'
        self.disabled = disabled
        
        self.checkpoints = self._load_checkpoints()
    
    def _load_checkpoints(self) -> Dict:
        """从文件加载检查点"""
        if self.disabled:
            return {}
        
        if self.checkpoint_file.exists():
            try:
                with open(self.checkpoint_file, 'r') as f:
                    checkpoints = json.load(f)
                logger.info(f"Loaded {len(checkpoints)} checkpoints from {self.checkpoint_file}")
                return checkpoints
            except Exception as e:
                logger.warning(f"Failed to load checkpoints: {e}")
                return {}
        
        return {}
    
    def _save_checkpoints(self):
        """保存检查点到文件"""
        if self.disabled:
            return
        
        try:
            # 确保目录存在
            self.checkpoint_file.parent.mkdir(parents=True, exist_ok=True)
            
            with open(self.checkpoint_file, 'w') as f:
                json.dump(self.checkpoints, f, indent=2)
            logger.debug(f"Saved checkpoints to {self.checkpoint_file}")
        except Exception as e:
            logger.error(f"Failed to save checkpoints: {e}")
    
    def _get_file_state(self, file_path: str) -> Optional[Dict]:
        """获取文件的状态信息"""
        path = Path(file_path)
        if not path.exists():
            return None
        
        stat = path.stat()
        return {
            'path': str(path.resolve()),
            'mtime': stat.st_mtime,
            'size': stat.st_size,
        }
    
    def mark_complete(self, step_name: str, 
                    # inputs: Dict[str, str], 
                    # outputs: Dict[str, str],
                    context: Dict = None):
        """
        标记步骤完成
        
        Args:
            step_name: 步骤名称
            inputs: 输入文件路径字典
            outputs: 输出文件路径字典
            context: 其他元数据（参数、版本等）
        """
        if self.disabled:
            return
        
        logger.info(f"Marking step '{step_name}' as complete")
        
        # # 记录输入文件的状态
        # input_states = {}
        # for name, path in inputs.items():
        #     if path:
        #         state = self._get_file_state(path)
        #         if state:
        #             input_states[name] = state
        
        # # 记录输出文件的状态
        # output_states = {}
        # for name, path in outputs.items():
        #     if path:
        #         state = self._get_file_state(path)
        #         if state:
        #             output_states[name] = state
        
        self.checkpoints[step_name] = {
            'status': 'complete',
            'timestamp': datetime.now().isoformat(),
            # 'inputs': input_states,
            # 'outputs': output_states,
            'context': context or {},
        }
        
        self._save_checkpoints()
        logger.debug(f"Marked {step_name} as complete")
    
    def mark_failed(self, step_name: str, error_message: str, 
                    inputs: Dict[str, str] = None):
        """
        标记步骤失败
        
        Args:
            step_name: 步骤名称
            error_message: 错误信息
            inputs: 当时的输入文件（可选）
        """
        if self.disabled:
            return
        
        logger.error(f"Marking step '{step_name}' as failed: {error_message}")
        
        checkpoint = {
            'status': 'failed',
            'timestamp': datetime.now().isoformat(),
            'error': error_message,
        }
        
        # 如果有输入文件，也记录下来
        if inputs:
            input_states = {}
            for name, path in inputs.items():
                if path:
                    state = self._get_file_state(path)
                    if state:
                        input_states[name] = state
            checkpoint['inputs'] = input_states
        
        self.checkpoints[step_name] = checkpoint
        self._save_checkpoints()
    
    def is_complete(self, step_name: str, 
                    current_inputs: Dict[str, str] = None) -> bool:
        """
        检查步骤是否已完成且输入未变化
        
        Args:
            step_name: 步骤名称
            current_inputs: 当前输入文件（用于检查变化）
        
        Returns:
            True 如果已完成且输入未变
        """
        if self.disabled:
            return False
        
        # 1. 检查是否有记录
        if step_name not in self.checkpoints:
            return False
        
        checkpoint = self.checkpoints[step_name]
        
        # 2. 检查状态
        if checkpoint.get('status') != 'complete':
            return False
        
        # 3. 检查输出文件是否都存在
        outputs = checkpoint.get('outputs', {})
        for name, state in outputs.items():
            if not Path(state['path']).exists():
                logger.warning(f"Output file missing for {step_name}: {state['path']}")
                return False
        
        # 4. 如果提供了当前输入，检查输入是否有变化
        if current_inputs:
            saved_inputs = checkpoint.get('inputs', {})
            
            # 检查输入文件数量是否一致
            current_input_keys = {k for k, v in current_inputs.items() if v}
            saved_input_keys = set(saved_inputs.keys())
            
            if current_input_keys != saved_input_keys:
                logger.info(f"Input files changed for {step_name}: "
                            f"current={current_input_keys}, saved={saved_input_keys}")
                return False
            
            # 检查每个输入文件是否有变化
            for name, current_path in current_inputs.items():
                if not current_path:  # 跳过None值
                    continue
                    
                saved_state = saved_inputs.get(name)
                if not saved_state:
                    logger.info(f"New input file for {step_name}: {name}")
                    return False
                
                # 检查文件是否存在
                current_path_obj = Path(current_path)
                if not current_path_obj.exists():
                    logger.warning(f"Current input file missing: {current_path}")
                    return False
                
                # 检查路径是否相同（考虑软链接）
                if str(current_path_obj.resolve()) != str(Path(saved_state['path']).resolve()):
                    logger.info(f"Input file path changed for {step_name}: {name}")
                    return False
                
                # 检查文件修改时间
                current_stat = current_path_obj.stat()
                if abs(current_stat.st_mtime - saved_state['mtime']) > 1:  # 允许1秒误差
                    logger.info(f"Input file modified for {step_name}: {name}")
                    return False
                
                # 检查文件大小
                if current_stat.st_size != saved_state['size']:
                    logger.info(f"Input file size changed for {step_name}: {name}")
                    return False
        
        logger.debug(f"Step '{step_name}' is complete and inputs unchanged")
        return True
    
    def is_failed(self, step_name: str) -> bool:
        """检查步骤是否失败"""
        if self.disabled:
            return False
        
        checkpoint = self.checkpoints.get(step_name, {})
        return checkpoint.get('status') == 'failed'
    
    def get_step_status(self, step_name: str) -> str:
        """获取步骤状态"""
        if self.disabled:
            return 'disabled'
        
        checkpoint = self.checkpoints.get(step_name, {})
        return checkpoint.get('status', 'not_started')
    
    def get_outputs(self, step_name: str) -> Dict[str, str]:
        """获取步骤的输出文件路径"""
        checkpoint = self.checkpoints.get(step_name, {})
        if checkpoint.get('status') != 'complete':
            return {}
        
        outputs = checkpoint.get('outputs', {})
        return {name: state['path'] for name, state in outputs.items()}
    
    def get_error(self, step_name: str) -> Optional[str]:
        """获取步骤的错误信息"""
        checkpoint = self.checkpoints.get(step_name, {})
        return checkpoint.get('error')
    
    def clear(self):
        """清除所有检查点"""
        self.checkpoints = {}
        if self.checkpoint_file.exists():
            self.checkpoint_file.unlink()
        logger.info("Cleared all checkpoints")
    
    def clear_step(self, step_name: str):
        """清除特定步骤的检查点"""
        if step_name in self.checkpoints:
            del self.checkpoints[step_name]
            self._save_checkpoints()
            logger.info(f"Cleared checkpoint for {step_name}")
    
    def get_all_checkpoints(self) -> Dict:
        """获取所有检查点"""
        return self.checkpoints.copy()
    
    def get_completed_steps(self) -> List[str]:
        """获取所有已完成的步骤"""
        return [name for name, cp in self.checkpoints.items() 
                if cp.get('status') == 'complete']
    
    def get_failed_steps(self) -> List[str]:
        """获取所有失败的步骤"""
        return [name for name, cp in self.checkpoints.items() 
                if cp.get('status') == 'failed']
    
    def print_summary(self):
        """打印检查点摘要"""
        logger.info("=" * 50)
        logger.info("Checkpoint Summary")
        logger.info("=" * 50)
        
        completed = self.get_completed_steps()
        failed = self.get_failed_steps()
        
        logger.info(f"Total checkpoints: {len(self.checkpoints)}")
        logger.info(f"Completed: {len(completed)}")
        logger.info(f"Failed: {len(failed)}")
        
        if completed:
            logger.info("\nCompleted steps:")
            for step in completed:
                timestamp = self.checkpoints[step].get('timestamp', 'unknown')
                logger.info(f"  - {step} ({timestamp})")
        
        if failed:
            logger.info("\nFailed steps:")
            for step in failed:
                error = self.checkpoints[step].get('error', 'unknown error')
                logger.info(f"  - {step}: {error}")
        
        logger.info("=" * 50)
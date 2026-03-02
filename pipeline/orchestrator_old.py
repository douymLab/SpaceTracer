#!/usr/bin/env python3
"""
Pipeline Orchestrator - 管理整个流程的执行
"""

import logging
import time
from pathlib import Path
from typing import Dict, List, Optional
from datetime import datetime
import os

from SpaceTracer.utils.get_genome_info import GenomeDetails
from SpaceTracer.pipeline.checkpoint import CheckpointManager
from SpaceTracer.pipeline.validator import Validator
# from SpaceTracer.steps import (
#     MappingStep,
#     BamProcessingStep,
#     Mpileup,
#     ReadsProcessingStep,
#     GenotypingStep,
#     FilteringStep,
#     FeatureExtractionStep,
#     PredictionStep,
#     FinalFilterStep
# )

# from SpaceTracer.steps.bam_processing import BamProcessingStep
from SpaceTracer.steps.step2_mpileup import MpileupStep

logger = logging.getLogger(__name__)

class PipelineOrchestrator:
    """
    流程调度器
    
    职责：
    1. 初始化所有步骤
    2. 管理步骤执行顺序
    3. 检查点管理
    4. 错误处理和恢复
    5. 生成运行报告
    """
    
    # 定义所有步骤及其顺序
    # STEPS = [
    #     ('mapping', MappingStep),
    #     ('bam_processing', BamProcessingStep),
    #     ('mpileup', MpileupStep),
    #     ('reads_processing', ReadsProcessingStep),
    #     ('genotyping', GenotypingStep),
    #     ('filtering', FilteringStep),
    #     ('feature_extraction', FeatureExtractionStep),
    #     ('prediction', PredictionStep),
    #     ('final_filter', FinalFilterStep)
    # ]
    STEPS = [
        # ('bam_processing',BamProcessingStep),
        ('mpileup', MpileupStep)
    ]

    def __init__(self, bam_file: str, regions_file: Optional[str],
                    output_dir: Path, config: Dict, 
                    resume: bool = False, force: bool = True):
        """
        初始化Pipeline
        
        Args:
            bam_file: 输入BAM文件路径
            regions_file: 区域BED文件路径（可选）
            output_dir: 输出目录
            config: 配置字典
            resume: 是否从检查点恢复
            force: 是否强制重新运行所有步骤
        """
        ### the initial step to load files and configs
        self.bam_file = bam_file
        self.regions_file = regions_file if regions_file else None
        self.output_dir = output_dir
        self.config = config
        self.resume = resume
        self.force = force
        
        # check output dir
        self.work_dir =  output_dir #os.path.join(self.output_dir ,'work')
        os.makedirs(self.work_dir, exist_ok=True)
        
        # 初始化检查点管理器
        self.checkpoint = CheckpointManager(self.output_dir, disabled=force)
        # 初始化验证器
        self.validator = Validator(self.config.get('skip_validation', False))
        
        # 初始化步骤
        self.steps = self._initialize_steps()
        
        # 运行统计
        self.stats = {
            'start_time': None,
            'end_time': None,
            'step_times': {},
            'step_status': {}
        }
    
    def _initialize_steps(self) -> Dict:
        """初始化所有步骤实例"""
        steps = {}
        
        Genome=GenomeDetails(self.config.get('genome'),self.config.get('genome_fasta'))
        genome_details=Genome._get_genome_details()
        context = {
            'bam_file': self.bam_file,
            'regions_file': self.regions_file,
            'output_dir': self.output_dir,
            'work_dir': self.work_dir,
            'genome_fasta': self.config.get('genome_fasta'),
            'genome_details': genome_details,
            'config': self.config
        }
        
        for step_name, step_class in self.STEPS:
            steps[step_name] = step_class(step_name, context, self.checkpoint)
        
        return steps
    
    def run(self, start_from: Optional[str] = None,
            stop_at: Optional[str] = None,
            only_steps: Optional[List[str]] = None) -> Dict:
        """
        运行pipeline
        
        Args:
            start_from: 从哪个步骤开始
            stop_at: 在哪个步骤停止
            only_steps: 只运行指定的步骤列表
        
        Returns:
            运行结果字典
        """
        logger.info("=" * 70)
        logger.info("Starting SpaceTracer Pipeline")
        logger.info("=" * 70)
        
        logger.debug(f"{self.config}")
        self.stats['start_time'] = time.time()
        
        # 确定要运行的步骤
        steps_to_run = self._plan_execution(start_from, stop_at, only_steps)
        
        logger.info(f"Steps to run: {', '.join(steps_to_run)}")
        logger.info("=" * 70)
        
        # 运行每个步骤
        
        for step_name in steps_to_run:
            step = self.steps[step_name]
            context=step.context
            # 运行步骤
            logger.info(f"[{step_name}] Starting...")
            step_start = time.time()
            
            try:
                context = step.execute(context)
                print(context)
                # 检查是否已完成
                if self.resume and not self.force:
                    if self.checkpoint.is_complete(step_name):
                        logger.info(f"[{step_name}] Already completed (skipping)")
                        # 加载该步骤的输出到context
                        context = step.load_outputs(context)
                        self.stats['step_status'][step_name] = 'skipped'
                        continue
                    
                step_elapsed = time.time() - step_start
                self.stats['step_times'][step_name] = step_elapsed
                self.stats['step_status'][step_name] = 'success'
                
                logger.info(f"[{step_name}] ✓ Completed in {step_elapsed:.2f}s")
                
                # 标记完成
                self.checkpoint.mark_complete(step_name, context)
                
                # 验证输出
                if not self.config.get('skip_validation', False):
                    if not self.validator.validate_step_output(step_name, context):
                        raise ValueError(f"Step {step_name} output validation failed")
                
            except Exception as e:
                logger.error(f"[{step_name}] ✗ Failed: {e}")
                self.stats['step_status'][step_name] = 'failed'
                self.checkpoint.mark_failed(step_name, str(e))
                raise 
            logger.debug(f"{context}")
            
        self.stats['end_time'] = time.time()
        
        # 最终结果
        results = {
            'final_vcf': context.get('final_vcf'),
            'total_variants': context.get('total_variants', 0),
            'elapsed_time': self.stats['end_time'] - self.stats['start_time'],
            'stats': self.stats
        }
        
        logger.info("=" * 70)
        logger.info("Pipeline completed successfully!")
        logger.info("=" * 70)
        
        return results
    
    def _plan_execution(self, start_from: Optional[str],
                    stop_at: Optional[str],
                    only_steps: Optional[List[str]]) -> List[str]:
        """
        计划要执行的步骤
        
        Returns:
            步骤名称列表
        """
        all_step_names = [name for name, _ in self.STEPS]
        
        # 如果指定了only_steps，只运行这些步骤
        if only_steps:
            return [s for s in only_steps if s in all_step_names]
        
        # 确定开始和结束位置
        start_idx = 0
        if start_from:
            if start_from not in all_step_names:
                raise ValueError(f"Unknown step: {start_from}")
            start_idx = all_step_names.index(start_from)
        
        end_idx = len(all_step_names)
        if stop_at:
            if stop_at not in all_step_names:
                raise ValueError(f"Unknown step: {stop_at}")
            end_idx = all_step_names.index(stop_at) + 1
        
        return all_step_names[start_idx:end_idx]
    
    def generate_report(self, output_path: Path):
        """生成HTML报告"""
        logger.info(f"Generating report: {output_path}")
        
        html = f"""
<!DOCTYPE html>
<html>
<head>
    <title>Mutation Caller Report</title>
    <style>
        body {{ font-family: Arial, sans-serif; margin: 40px; }}
        h1 {{ color: #333; }}
        table {{ border-collapse: collapse; width: 100%; margin: 20px 0; }}
        th, td {{ border: 1px solid #ddd; padding: 12px; text-align: left; }}
        th {{ background-color: #4CAF50; color: white; }}
        .success {{ color: green; }}
        .failed {{ color: red; }}
        .skipped {{ color: orange; }}
    </style>
</head>
<body>
    <h1>Mutation Caller Pipeline Report</h1>
    
    <h2>Summary</h2>
    <ul>
        <li><strong>Input BAM:</strong> {self.bam_file}</li>
        <li><strong>Output Directory:</strong> {self.output_dir}</li>
        <li><strong>Total Time:</strong> {self.stats['end_time'] - self.stats['start_time']:.2f}s</li>
        <li><strong>Completion Time:</strong> {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</li>
    </ul>
    
    <h2>Step Details</h2>
    <table>
        <tr>
            <th>Step</th>
            <th>Status</th>
            <th>Time (s)</th>
        </tr>
"""
        
        for step_name, _ in self.STEPS:
            status = self.stats['step_status'].get(step_name, 'not_run')
            elapsed = self.stats['step_times'].get(step_name, 0)
            
            status_class = status
            html += f"""
        <tr>
            <td>{step_name}</td>
            <td class="{status_class}">{status}</td>
            <td>{elapsed:.2f}</td>
        </tr>
"""
        
        html += """
    </table>
</body>
</html>
"""
        
        with open(output_path, 'w') as f:
            f.write(html)



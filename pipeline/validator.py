#!/usr/bin/env python3
"""
Validator - 验证Pipeline各个步骤的输入和输出
"""

import logging
from pathlib import Path
from typing import Dict, Any, Optional, List
import subprocess

logger = logging.getLogger(__name__)

class Validator:
    """
    验证器类
    
    功能：
    1. 验证文件存在性
    2. 验证文件格式
    3. 验证文件内容完整性
    4. 验证数据质量
    """
    
    def __init__(self, skip_validation: bool = False):
        """
        初始化验证器
        
        Args:
            skip_validation: 是否跳过验证（用于快速测试）
        """
        self.skip_validation = skip_validation
        
        if skip_validation:
            logger.warning("Validation is DISABLED - use with caution!")
    
    def validate_step_output(self, step_name: str, context: Dict) -> bool:
        """
        验证步骤输出
        
        Args:
            step_name: 步骤名称
            context: 上下文字典
        
        Returns:
            True if validation passed
        """
        if self.skip_validation:
            return True
        
        logger.debug(f"Validating output for step: {step_name}")
        
        # 根据步骤名称选择验证方法
        validators = {
            'mapping': self._validate_mapping,
            'bam_processing': self._validate_bam_processing,
            'mpileup': self._validate_mpileup,
            'reads_processing': self._validate_reads_processing,
            'genotyping': self._validate_genotyping,
            'filtering': self._validate_filtering,
            'feature_extraction': self._validate_feature_extraction,
            'prediction': self._validate_prediction,
            'final_filter': self._validate_final_filter,
        }
        
        validator_func = validators.get(step_name)
        
        if validator_func:
            try:
                return validator_func(context)
            except Exception as e:
                logger.error(f"Validation failed for {step_name}: {e}")
                return False
        else:
            logger.warning(f"No specific validator for step: {step_name}, using basic validation")
            return self._validate_basic(context)
    
    # ==================== 基础验证方法 ====================
    
    def _validate_basic(self, context: Dict) -> bool:
        """基础验证：检查关键输出文件是否存在"""
        
        # 检查所有路径类型的值
        for key, value in context.items():
            if isinstance(value, (str, Path)):
                path = Path(value)
                
                # 如果看起来像文件路径
                if any(ext in str(path) for ext in ['.bam', '.vcf', '.txt', '.tsv', '.bed']):
                    if not path.exists():
                        logger.error(f"Output file missing: {key} -> {path}")
                        return False
                    
                    if path.stat().st_size == 0:
                        logger.error(f"Output file is empty: {key} -> {path}")
                        return False
        
        return True
    
    def validate_file_exists(self, filepath: Path, min_size: int = 0) -> bool:
        """
        验证文件存在且非空
        
        Args:
            filepath: 文件路径
            min_size: 最小文件大小（字节）
        
        Returns:
            True if valid
        """
        filepath = Path(filepath)
        
        if not filepath.exists():
            logger.error(f"File does not exist: {filepath}")
            return False
        
        size = filepath.stat().st_size
        
        if size < min_size:
            logger.error(f"File too small ({size} bytes, min: {min_size}): {filepath}")
            return False
        
        return True
    
    def validate_bam_file(self, bam_file: Path, check_index: bool = True) -> bool:
        """
        validate bam file
        bam_file: the path of bam file
        check_index: check bam index or not?
        """
        bam_file = Path(bam_file)
        
        if not self.validate_file_exists(bam_file, min_size=1000):
            return False
        
        if check_index:
            bai_file = Path(str(bam_file) + '.bai')
            if not bai_file.exists():
                logger.warning(f"BAM index not found: {bai_file}")
        
        try:
            result = subprocess.run(
                ['samtools', 'quickcheck', str(bam_file)],
                capture_output=True,
                text=True,
                timeout=60
            )
            
            if result.returncode != 0:
                logger.error(f"BAM file is corrupted: {bam_file}")
                logger.error(f"samtools output: {result.stderr}")
                return False
            
            logger.debug(f"✓ BAM file valid: {bam_file}")
            return True
            
        except FileNotFoundError:
            logger.warning("samtools not found, skipping BAM format validation")
            return True
        except subprocess.TimeoutExpired:
            logger.error(f"BAM validation timed out: {bam_file}")
            return False
        except Exception as e:
            logger.error(f"Failed to validate BAM: {e}")
            return False
    
    def validate_vcf_file(self, vcf_file: Path) -> bool:
        """
        验证VCF文件
        
        Args:
            vcf_file: VCF文件路径
        
        Returns:
            True if valid
        """
        vcf_file = Path(vcf_file)
        
        # 检查文件存在
        if not self.validate_file_exists(vcf_file, min_size=50):
            return False
        
        # 检查文件格式
        try:
            with open(vcf_file) as f:
                first_line = f.readline()
                
                # VCF文件应该以##开头（header）或直接是数据行
                if not (first_line.startswith('#') or '\t' in first_line):
                    logger.error(f"Invalid VCF format: {vcf_file}")
                    return False
                
                # 统计行数
                line_count = sum(1 for line in f if not line.startswith('#'))
                
                if line_count == 0:
                    logger.warning(f"VCF file has no data lines: {vcf_file}")
                    # 不算错误，可能是没有variant
                
                logger.debug(f"✓ VCF file valid: {vcf_file} ({line_count} variants)")
                return True
                
        except Exception as e:
            logger.error(f"Failed to validate VCF: {e}")
            return False
    
    def validate_tsv_file(self, tsv_file: Path, required_columns: Optional[List[str]] = None,
                         min_rows: int = 0) -> bool:
        """
        验证TSV文件
        
        Args:
            tsv_file: TSV文件路径
            required_columns: 必需的列名列表
            min_rows: 最小行数
        
        Returns:
            True if valid
        """
        tsv_file = Path(tsv_file)
        
        # 检查文件存在
        if not self.validate_file_exists(tsv_file, min_size=10):
            return False
        
        try:
            with open(tsv_file) as f:
                lines = f.readlines()
                
                if len(lines) == 0:
                    logger.error(f"TSV file is empty: {tsv_file}")
                    return False
                
                # 检查列名（如果提供）
                if required_columns:
                    header = lines[0].strip().split('\t')
                    
                    for col in required_columns:
                        if col not in header:
                            logger.error(f"Missing required column '{col}' in {tsv_file}")
                            return False
                
                # 检查最小行数
                data_lines = [line for line in lines if not line.startswith('#')]
                
                if len(data_lines) - 1 < min_rows:  # -1 for header
                    logger.error(f"TSV has only {len(data_lines)-1} rows, minimum is {min_rows}")
                    return False
                
                logger.debug(f"✓ TSV file valid: {tsv_file} ({len(data_lines)-1} rows)")
                return True
                
        except Exception as e:
            logger.error(f"Failed to validate TSV: {e}")
            return False
    
    # ==================== 步骤特定验证方法 ====================
    
    def _validate_mapping(self, context: Dict) -> bool:
        """验证Mapping步骤输出"""
        
        mapped_bam = context.get('mapped_bam')
        
        if not mapped_bam:
            logger.error("Mapping step did not produce mapped_bam")
            return False
        
        # 验证BAM文件
        if not self.validate_bam_file(mapped_bam, check_index=True):
            return False
        
        # 检查BAM是否有reads
        try:
            result = subprocess.run(
                ['samtools', 'view', '-c', str(mapped_bam)],
                capture_output=True,
                text=True,
                timeout=60
            )
            
            read_count = int(result.stdout.strip())
            
            if read_count == 0:
                logger.error("Mapped BAM has no reads")
                return False
            
            logger.info(f"✓ Mapped BAM contains {read_count:,} reads")
            return True
            
        except Exception as e:
            logger.warning(f"Could not count reads: {e}")
            return True  # 不算致命错误
    
    def _validate_bam_processing(self, context: Dict) -> bool:
        """验证BAM Processing步骤输出"""
        
        processed_bam = context.get('in_filter_bam')
        
        if not processed_bam:
            logger.error("BAM processing did not produce processed_bam")
            return False
        
        return self.validate_bam_file(processed_bam, check_index=True)
    
    def _validate_mpileup(self, context: Dict) -> bool:
        """验证Mpileup步骤输出"""
        
        # 检查mpileup文件
        mpileup_file = context.get('mpileup_file')
        vcf_file = context.get('vcf_file')
        
        if mpileup_file:
            if not self.validate_file_exists(mpileup_file, min_size=100):
                return False
            
            # 检查格式
            try:
                with open(mpileup_file) as f:
                    first_line = f.readline().strip()
                    
                    # mpileup格式应该有至少6列
                    fields = first_line.split('\t')
                    
                    if len(fields) < 6:
                        logger.error(f"Invalid mpileup format (expected ≥6 columns, got {len(fields)})")
                        return False
                    
                    line_count = sum(1 for _ in f) + 1
                    logger.info(f"✓ Mpileup file has {line_count:,} positions")
                    
            except Exception as e:
                logger.error(f"Failed to validate mpileup file: {e}")
                return False
        
        # 检查VCF文件
        if vcf_file:
            if not self.validate_vcf_file(vcf_file):
                return False
        
        return True
    
    def _validate_reads_processing(self, context: Dict) -> bool:
        """验证Reads Processing步骤输出"""
        
        reads_data = context.get('reads_data')
        
        if not reads_data:
            logger.error("Reads processing did not produce reads_data")
            return False
        
        # 验证为TSV文件
        return self.validate_tsv_file(reads_data, min_rows=1)
    
    def _validate_genotyping(self, context: Dict) -> bool:
        """验证Genotyping步骤输出"""
        
        genotypes = context.get('genotypes')
        
        if not genotypes:
            logger.error("Genotyping did not produce genotypes file")
            return False
        
        # 验证TSV文件
        required_cols = ['chrom', 'pos', 'ref', 'alt', 'genotype']
        
        if not self.validate_tsv_file(genotypes, required_columns=required_cols, min_rows=0):
            return False
        
        # 检查统计信息
        stats = context.get('genotype_stats', {})
        call_rate = stats.get('call_rate', 0)
        
        if call_rate < 0.01:  # 少于1%的call rate可能有问题
            logger.warning(f"Low genotype call rate: {call_rate:.2%}")
        
        return True
    
    def _validate_filtering(self, context: Dict) -> bool:
        """验证Filtering步骤输出"""
        
        filtered_variants = context.get('filtered_variants')
        
        if not filtered_variants:
            logger.error("Filtering did not produce filtered_variants")
            return False
        
        return self.validate_tsv_file(filtered_variants, min_rows=0)
    
    def _validate_feature_extraction(self, context: Dict) -> bool:
        """验证Feature Extraction步骤输出"""
        
        features_file = context.get('features_file')
        
        if not features_file:
            logger.error("Feature extraction did not produce features_file")
            return False
        
        if not self.validate_tsv_file(features_file, min_rows=0):
            return False
        
        # 检查是否有足够的features
        try:
            with open(features_file) as f:
                header = f.readline().strip().split('\t')
                
                # 应该有多个feature列
                if len(header) < 5:
                    logger.warning(f"Only {len(header)} features extracted (might be too few)")
                
                logger.info(f"✓ Extracted {len(header)} features")
                
        except Exception as e:
            logger.warning(f"Could not check feature count: {e}")
        
        return True
    
    def _validate_prediction(self, context: Dict) -> bool:
        """验证Prediction步骤输出"""
        
        predictions = context.get('predictions')
        
        if not predictions:
            logger.error("Prediction did not produce predictions file")
            return False
        
        # 验证预测文件
        if not self.validate_tsv_file(predictions, min_rows=0):
            return False
        
        # 检查预测分数范围
        try:
            with open(predictions) as f:
                header = f.readline()
                
                for i, line in enumerate(f):
                    if i >= 100:  # 只检查前100行
                        break
                    
                    fields = line.strip().split('\t')
                    
                    # 假设最后一列是预测分数
                    try:
                        score = float(fields[-1])
                        
                        if score < 0 or score > 1:
                            logger.warning(f"Prediction score out of range [0,1]: {score}")
                    except (ValueError, IndexError):
                        pass
            
            logger.debug("✓ Prediction scores validated")
            
        except Exception as e:
            logger.warning(f"Could not validate prediction scores: {e}")
        
        return True
    
    def _validate_final_filter(self, context: Dict) -> bool:
        """验证Final Filter步骤输出"""
        
        final_vcf = context.get('final_vcf')
        
        if not final_vcf:
            logger.error("Final filter did not produce final_vcf")
            return False
        
        # 验证VCF文件
        if not self.validate_vcf_file(final_vcf):
            return False
        
        # 统计最终variant数量
        try:
            with open(final_vcf) as f:
                variant_count = sum(1 for line in f if not line.startswith('#'))
                
                logger.info(f"✓ Final VCF contains {variant_count:,} variants")
                
                if variant_count == 0:
                    logger.warning("No variants in final output (might be expected)")
        
        except Exception as e:
            logger.warning(f"Could not count variants: {e}")
        
        return True
    
    # ==================== 辅助验证方法 ====================
    
    def validate_reference_genome(self, fasta_file: Path) -> bool:
        """
        验证参考基因组FASTA文件
        
        Args:
            fasta_file: FASTA文件路径
        
        Returns:
            True if valid
        """
        fasta_file = Path(fasta_file)
        
        # 检查文件存在
        if not self.validate_file_exists(fasta_file, min_size=1000000):
            return False
        
        # 检查fai索引
        fai_file = Path(str(fasta_file) + '.fai')
        
        if not fai_file.exists():
            logger.warning(f"FASTA index not found: {fai_file}")
            logger.info("You may need to run: samtools faidx {fasta_file}")
        
        # 检查FASTA格式
        try:
            with open(fasta_file) as f:
                first_line = f.readline()
                
                if not first_line.startswith('>'):
                    logger.error(f"Invalid FASTA format: {fasta_file}")
                    return False
            
            logger.debug(f"✓ Reference genome valid: {fasta_file}")
            return True
            
        except Exception as e:
            logger.error(f"Failed to validate reference genome: {e}")
            return False
    
    def validate_bed_file(self, bed_file: Path, min_regions: int = 0) -> bool:
        """
        验证BED文件
        
        Args:
            bed_file: BED文件路径
            min_regions: 最小区域数量
        
        Returns:
            True if valid
        """
        bed_file = Path(bed_file)
        
        # 检查文件存在
        if not self.validate_file_exists(bed_file, min_size=10):
            return False
        
        try:
            with open(bed_file) as f:
                region_count = 0
                
                for line in f:
                    if line.startswith('#') or line.strip() == '':
                        continue
                    
                    fields = line.strip().split('\t')
                    
                    # BED文件至少应该有3列：chr, start, end
                    if len(fields) < 3:
                        logger.error(f"Invalid BED format (expected ≥3 columns): {line.strip()}")
                        return False
                    
                    # 验证start和end是数字
                    try:
                        start = int(fields[1])
                        end = int(fields[2])
                        
                        if start >= end:
                            logger.warning(f"Invalid region (start >= end): {line.strip()}")
                    except ValueError:
                        logger.error(f"Invalid coordinates in BED file: {line.strip()}")
                        return False
                    
                    region_count += 1
                
                if region_count < min_regions:
                    logger.error(f"BED file has only {region_count} regions, minimum is {min_regions}")
                    return False
                
                logger.debug(f"✓ BED file valid: {bed_file} ({region_count} regions)")
                return True
                
        except Exception as e:
            logger.error(f"Failed to validate BED file: {e}")
            return False
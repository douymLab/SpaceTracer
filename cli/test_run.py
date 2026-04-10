# test_pipeline.py

from SpaceTracer.pipeline.orchestrator import PipelineOrchestrator
from SpaceTracer.config.config_loader import LoadConfig


from SpaceTracer.utils.logger import setup_logger
setup_logger('INFO') #DEBUG

def test_two_steps():
    
    # 加载配置
    config = LoadConfig().load_config(
        # genome='hg38',
        # genome_fasta='/storage/douyanmeiLab/yangzhirui/Reference/Cellranger/refdata-gex-GRCh38-2020-A/fasta/genome.fa',
        # custom_config='/storage/douyanmeiLab/yangzhirui/SpaceTracer_new/test/test_config.yaml'  # 你的测试配置
        # custom_config='/storage/douyanmeiLab/yangzhirui/SpaceTracer_new/BCC-1_test/config.yaml'
        custom_config='/storage/douyanmeiLab/yangzhirui/SpaceTracer_new/BCC-1_test_cmd/config.yaml'
    )
    
    # 创建pipeline（只注册两个步骤）
    pipeline = PipelineOrchestrator(
        config=config,
        force=True
    )
    
    # 🔥 临时修改：只运行这两个步骤
    # pipeline.STEP_CLASSES = {
    #     # 'bam_processing':BamProcessingStep,
    #     'mpileup': MpileupStep
    # }
    
    # 运行
        # 'bam_processing':BamProcessingStep,
        # 'mpileup': MpileupStep,
        # 'umi_combine': UMICombineStep,
        # 'prior': PriorCalculator
        
        
    results = pipeline.run(only_steps=["umi_combine"])
    
    print(f"✅ Pipeline completed!")
    # print(f"   Variants called: {results['stats']['genotyping']['n_variants']}")
    
if __name__ == '__main__':
    test_two_steps()

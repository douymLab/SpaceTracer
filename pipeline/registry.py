from SpaceTracer.steps.step0_cluster import ClusterStep
from SpaceTracer.steps.step1_bam_processing import BamProcessingStep
from SpaceTracer.steps.step2_mpileup import MpileupStep
from SpaceTracer.steps.step3_UMI_combine import UMICombineStep
from SpaceTracer.steps.step3_cell_number import CellNumStep
from SpaceTracer.steps.step3_get_prior import PriorCalculator
from SpaceTracer.steps.step4_genotyping import GenotypingStep
from SpaceTracer.steps.step5_RNA_level_feature import RNAFeatureStep
from SpaceTracer.steps.step5_read_feature import ReadFeatureStep
from SpaceTracer.steps.step5_spatial_feature import SpatialFeatureStep
from SpaceTracer.steps.step5_mappability_feature import MappabilityFeatureStep
from SpaceTracer.steps.step6_merge_all_features import MergeFeatureStep
from SpaceTracer.steps.step7_mutation_prediction import MutationPredictionStep


STEP_CLASSES = {
    "cluster": ClusterStep,
    "bam_processing": BamProcessingStep,
    "mpileup": MpileupStep,
    "umi_combine": UMICombineStep,
    "cell_num": CellNumStep,
    "prior": PriorCalculator,
    "genotyping": GenotypingStep,
    "spatial_feature": SpatialFeatureStep,
    "mappability_feature": MappabilityFeatureStep,
    "read_feature": ReadFeatureStep,
    "RNA_feature": RNAFeatureStep,
    "merge_feature": MergeFeatureStep,
    # "mutation_prediction": MutationPredictionStep
}
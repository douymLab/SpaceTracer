# SpaceTracer/pipeline/registry.py
from __future__ import annotations

import importlib
from typing import Dict, Type, Any

STEP_CLASS_PATHS: Dict[str, str] = {
    "cluster": "SpaceTracer.steps.step0_cluster:ClusterStep",
    "bam_processing": "SpaceTracer.steps.step1_bam_processing:BamProcessingStep",
    "mpileup": "SpaceTracer.steps.step2_mpileup:MpileupStep",
    "umi_combine": "SpaceTracer.steps.step3_UMI_combine:UMICombineStep",
    "cell_num": "SpaceTracer.steps.step3_cell_number:CellNumStep",
    "prior": "SpaceTracer.steps.step3_get_prior:PriorCalculator",
    "genotyping": "SpaceTracer.steps.step4_genotyping:GenotypingStep",
    "spatial_feature": "SpaceTracer.steps.step5_spatial_feature:SpatialFeatureStep",
    "mappability_feature": "SpaceTracer.steps.step5_mappability_feature:MappabilityFeatureStep",
    "read_feature": "SpaceTracer.steps.step5_read_feature:ReadFeatureStep",
    "RNA_feature": "SpaceTracer.steps.step5_RNA_level_feature:RNAFeatureStep",
    "phasing": "SpaceTracer.steps.step5_phasing:PhasingCandidateStep",
    "merge_feature": "SpaceTracer.steps.step6_merge_all_features:MergeFeatureStep",
    "mutation_prediction": "SpaceTracer.steps.step7_mutation_prediction:MutationPredictionStep",
    "phylogeny": "SpaceTracer.steps.step8_phylosolid:PhylogenyBuildStep"
}

def _load_symbol(path: str):
    mod_path, symbol = path.split(":")
    mod = importlib.import_module(mod_path)
    return getattr(mod, symbol)

def get_step_class(step_name: str):
    try:
        path = STEP_CLASS_PATHS[step_name]
    except KeyError:
        raise KeyError(f"Unknown step: {step_name}. Available: {list(STEP_CLASS_PATHS)}")
    return _load_symbol(path)

def get_step_classes() -> Dict[str, Any]:
    return {name: get_step_class(name) for name in STEP_CLASS_PATHS}

# Algorithm Walkthrough

This page explains how SpaceTracer converts spatial transcriptomics alignments into lineage-informative somatic SNV evidence.

## Pipeline design

SpaceTracer is organized as a directed acyclic graph (DAG) of processing steps. Each step produces context outputs consumed by downstream steps.

The current core execution flow (aligned to `SpaceTracer/steps`) is:

1. `step0_cluster.py` -> `cluster`
2. `step1_bam_processing.py` -> `bam_processing`
3. `step2_mpileup.py` -> `mpileup`
4. `step3_UMI_combine.py` -> `umi_combine`
5. `step3_cell_number.py` -> `cell_num`
6. `step3_get_prior.py` -> `prior`
7. `step4_genotyping.py` -> `genotyping`
8. `step5_spatial_feature.py`, `step5_mappability_feature.py`, `step5_read_feature.py`, `step5_RNA_level_feature.py` -> feature branches
9. `step5_phasing.py` -> `phasing`
10. `step6_merge_all_features.py` -> `merge_feature`
11. `step7_mutation_prediction.py` -> `mutation_prediction`

For a detailed per-step reference (inputs, parameters, and outputs), see [Step Reference Overview](steps/overview.md). For practical rerun/debug patterns, see [Single-Step Debug Cookbook](steps/debug-cookbook.md).

## Step-by-step meaning

### 1) `cluster`

Builds or loads spot/domain grouping information used in downstream genotype and spatial inference.

### 2) `bam_processing`

Prepares BAM-level data for robust pileup and candidate detection (sorting/filtering/index-friendly preprocessing).

### 3) `mpileup`

Generates base-level evidence from aligned reads across the genome/chunks.

### 4) `umi_combine`

Aggregates read evidence at UMI level to reduce read-level technical noise and improve confidence.

### 5) `cell_num`

Estimates spot/cell-level support statistics needed for later probabilistic genotyping.

### 6) `prior`

Builds prior information for mutation likelihood estimation.

### 7) `genotyping`

Combines evidence and priors to infer genotype-level mutation signals.

### 8) Feature branches

From genotyping outputs, SpaceTracer computes multiple complementary feature families:

- `spatial_feature`: neighborhood/tissue-structure signal
- `mappability_feature`: regional mappability/confounder signal
- `read_feature`: read-level quality/bias signal
- `RNA_feature`: RNA-level context (including expression-related cues)

### 9) `phasing`

Refines candidate evidence with RNA-informed phasing information and cluster-level event summaries.

### 10) `merge_feature`

Merges all feature families (plus phasing outputs) into an integrated feature representation for downstream prioritization.

## Step detail index

- [cluster](steps/cluster.md)
- [bam_processing](steps/bam-processing.md)
- [mpileup](steps/mpileup.md)
- [umi_combine](steps/umi-combine.md)
- [cell_num](steps/cell-num.md)
- [prior](steps/prior.md)
- [genotyping](steps/genotyping.md)
- [spatial_feature](steps/spatial-feature.md)
- [mappability_feature](steps/mappability-feature.md)
- [read_feature](steps/read-feature.md)
- [RNA_feature](steps/rna-feature.md)
- [phasing](steps/phasing.md)
- [merge_feature](steps/merge-feature.md)
- [mutation_prediction](steps/mutation-prediction.md)

## Why this structure works

SpaceTracer combines orthogonal information to suppress false positives:

- read/UMI evidence reduces sequencing artifacts
- prior modeling stabilizes genotype inference
- mappability and RNA-level features handle context-specific noise
- spatial features preserve in situ biological structure

This combination improves confidence in mosaic SNV discovery for lineage analysis.

## Parallel execution and checkpoints

SpaceTracer supports:

- parallel execution for independent feature steps
- checkpoint-aware resume (skip completed steps)
- partial execution with `--start-from` and `--stop-at`
- explicit subset execution with `--only-steps` (no automatic dependency completion outside the listed subset)

These capabilities make iterative analysis and parameter tuning practical on real datasets.

## Practical interpretation

When reading results, think in three layers:

1. **Evidence layer**: pileup + UMI + genotype calls
2. **Feature layer**: spatial/read/mappability/RNA features
3. **Integration layer**: merged features for final candidate prioritization

This helps diagnose whether a candidate variant is supported by strong multi-modal evidence or likely a technical artifact.

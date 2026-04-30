# Step-by-step guide

This page is aligned to the code layout in `SpaceTracer/steps` so users can map tutorials directly to implementation files.

For exact input/parameter interpretation, open each step page under this section.

## Code-aligned execution map

| Step tutorial | Runtime step name | Purpose | Code step file |
| --- | --- | --- | --- |
| [Step 0: Cluster Preparation](step0-cluster.md) | `cluster` | Load or compute spot/domain clusters. | `step0_cluster.py` |
| [Step 1: BAM Processing](step1-bam-processing.md) | `bam_processing` | Build in-tissue filtered BAM for downstream analysis. | `step1_bam_processing.py` |
| [Step 2: Pileup Candidate Sites](step2-mpileup.md) | `mpileup` | Generate candidate loci and chunk metadata. | `step2_mpileup.py` |
| [Step 3: Count and Prior Construction](step3-count-and-prior-construction.md) | `umi_combine`<br>`cell_num`<br>`prior` | Aggregate UMI-level count evidence, build cell-number support values, and construct prior information for candidate-site genotyping. | `step3_UMI_combine.py`<br>`step3_cell_number.py`<br>`step3_get_prior.py` |
| [Step 4: Genotyping](step4-genotyping.md) | `genotyping` | Infer genotype evidence at multiple levels. | `step4_genotyping.py` |
| [Step 5: Feature Extraction and Phasing](step5-feature-extraction.md) | `spatial_feature`<br>`mappability_feature`<br>`read_feature`<br>`RNA_feature`<br>`phasing` | Extract spatial, mappability, read-level, and RNA-level features, then refine candidate sites with phasing evidence. | `step5_spatial_feature.py`<br>`step5_mappability_feature.py`<br>`step5_read_feature.py`<br>`step5_RNA_level_feature.py`<br>`step5_phasing.py` |
| [Step 6: Feature Merge and Filtration](step6-merge-filtration.md) | `merge_feature` | Merge features and apply filtration tags/switches. | `step6_merge_all_features.py` |
| [Step 7: Mutation Prediction](step7-mutation-prediction.md) | `mutation_prediction` | Run model inference and export VCF results. | `step7_mutation_prediction.py` |

Post-processing outside `SpaceTracer/steps`:

- [Post-step: Phylogeny (optional)](post-phylogeny.md)

## How to run

### Recommended: one-command workflow

```bash
SpaceTracer run --config config.yaml
```

### Advanced: inspect by code step order

Use these pages in code order:

- [Step 0: Cluster Preparation](step0-cluster.md)
- [Step 1: BAM Processing](step1-bam-processing.md)
- [Step 2: Pileup Candidate Sites](step2-mpileup.md)
- [Step 3: Count and Prior Construction](step3-count-and-prior-construction.md)
- [Step 4: Genotyping](step4-genotyping.md)
- [Step 5: Feature Extraction and Phasing](step5-feature-extraction.md)
- [Step 6: Feature Merge and Filtration](step6-merge-filtration.md)
- [Step 7: Mutation Prediction](step7-mutation-prediction.md)
- [Post-step: Phylogeny (optional)](post-phylogeny.md)

## Important outputs by stage

- Code step 4: genotype evidence files (for example `ind_geno_filter_file`, `spot_geno_file`, `cluster_vaf_file`)
- Code step 6: merged/filtered feature tables (`all_feature.txt` and parquet mirror)
- Code step 7: predicted VCF outputs from mutation prediction
- Optional post-step phylogeny: lineage analysis from high-confidence loci

For practical rerun recipes (for example, rerun only RNA or spatial branches), see [Single-Step Debug Cookbook](debug-cookbook.md).

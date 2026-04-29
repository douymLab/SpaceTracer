# Step-by-step guide

This page combines the algorithm-level view and step-level execution map, so you can understand both **why** and **how** the full workflow runs.

For exact input/parameter interpretation, open each step page under this section.

## Main algorithm logic

1. **Candidate generation**  
   Generate candidate variant sites from mpileup and basic depth/VAF filters.
2. **Allele counting and quality control**  
   Build spot, cluster, and individual-level count tables after quality checks.
3. **Bayesian/statistical genotyping**  
   Infer genotypes using allele counts, prior frequencies, and confidence thresholds.
4. **Feature extraction**  
   Extract spatial, read-level, mappability, and transcriptomic features.
5. **Mutation prediction**  
   Use random-forest models to classify true somatic mutations.
6. **Artifact removal**  
   Remove recurrent artifacts using dbSNP, editing databases, imprinted regions, and PON resources.

## Full workflow steps

| Step | Main SpaceTracer modules | Purpose |
| --- | --- | --- |
| Step0 | `cluster`, `bam_processing`, `mpileup` | Data pre-processing: cluster metadata preparation, BAM filtering, and mpileup candidate generation. |
| Step1 | `umi_combine`, `cell_num`, `prior` | Quality-control-oriented evidence aggregation and prior preparation for robust genotyping. |
| Step2 | `genotyping` | Spot/cluster/individual-level genotype inference. |
| Step3 | `spatial_feature`, `mappability_feature`, `read_feature`, `RNA_feature`, `phasing`, `merge_feature` | Multi-source feature extraction plus phasing refinement, followed by integrated feature matrix generation. |
| Step4 | `mutation_prediction` | Model-based mutation prediction and VCF export. |
| Step5 | external filtering stage | Remove recurrent artifacts using public/PoN-style resources. |
| Step6 | external phylogeny tools | Optional lineage tree reconstruction from high-confidence mutation calls. |

## How to run

### Recommended: one-command workflow

```bash
SpaceTracer run --config config.yaml
```

### Advanced: script-by-script workflow

Use the per-step pages directly:

- [Step 0: Data pre-processing](step0-preprocessing.md)
- [Step 1: Quality control](step1-quality-control.md)
- [Step 2: Genotyping](step2-genotyping.md)
- [Step 3: Feature extraction](step3-feature-extraction.md)
- [Step 4: Mutation prediction](step4-mutation-prediction.md)
- [Step 5: Remove recurrent artifacts](step5-remove-recurrent-artifacts.md)
- [Step 6: Phylogeny (optional)](step6-phylogeny.md)

## Important outputs by stage

- Step1/2: genotype evidence files (for example `ind_geno_filter_file`, `spot_geno_file`, `cluster_vaf_file`)
- Step3: merged feature table (`all_feature.txt` and parquet mirror)
- Step4: predicted VCF outputs from mutation prediction
- Step5/6: curated high-confidence loci used for downstream lineage/phylogeny analyses

For practical rerun recipes (for example, rerun only RNA or spatial branches), see [Single-Step Debug Cookbook](debug-cookbook.md).

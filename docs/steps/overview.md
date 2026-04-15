# Step-by-step guide

This page combines the algorithm-level view and step-level execution map, so you can understand both **why** and **how** the full workflow runs.

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

| Step | Main script | Purpose |
| --- | --- | --- |
| Step0 | `0_other_func.py`, `others/get_umiCount_cellNum.py` | Data preparation: BAM filtering, mpileup, prior generation, cluster/cell-number utilities. |
| Step1 | `1_run_data_process.py` | Allele counting with quality/background filtering. |
| Step2 | `2_run_genotyper.py` | Spot/cluster/individual-level genotype inference. |
| Step3 | `3_run_get_features.py` | Multi-source feature extraction (spatial, read, mappability, RNA). |
| Step4 | `4_run_model_predict.py` | Random-forest based mutation prediction. |
| Step5 | `5_remove_recurrent_mutations.sh` | Remove recurrent artifacts using public/PoN resources. |
| Step6 | external `PhyloSOLID` | Optional phylogeny reconstruction from called mutations. |

## How to run

### Recommended: one-command workflow (Snakemake)

```bash
snakemake --configfile config_example.txt -s run_snakemake_for_conda
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

- Step1: `*.spot_count.out`, `*.cluster.count.out`, `*.ind_filter.count.out`
- Step2: `*.ind_genotype_filter.out`, `*.spot_genotype.out`, `*.cluster_vaf.out`
- Step3: `*.features.add_hFDR.txt`
- Step4: `predict/results/*pred_truesites*.txt`
- Step5: `$model/pred.FINAL.txt`

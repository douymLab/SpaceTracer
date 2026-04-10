# `mutation_prediction`

## Purpose

Runs model-based mutation classification from merged features and produces VCF outputs.

## Current status

- Implemented in `steps/step7_mutation_prediction.py`
- Currently **not enabled in the active default DAG** in `pipeline/dag.py` / orchestrator registration

## Upstream dependency (when enabled)

- `merge_feature`

## Required inputs

- `combine_feature` (`all_feature.txt` or equivalent)

## Key parameters

From `steps.mutation_prediction`:

- `model_dir`
- `model_name`

The current implementation also uses internal defaults for many model options (training mode, plotting, SHAP/PCA export, RF hyperparameters, and class balancing strategy).

## Outputs

- `raw_pred_vcf`: `<step_dir>/Sample_total_pred_truesites.vcf`
- `final_vcf`: `<step_dir>/Sample_total_pred_truesites_PASS.vcf`

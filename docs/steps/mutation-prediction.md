# `mutation_prediction`

## Purpose

Runs model-based mutation classification from merged features and produces VCF outputs.

## Upstream

- `merge_feature`

## Required inputs

- `combine_feature_parquet` (`all_feature.parquet` from `merge_feature`)

### Input interpretation

| Input key | Source step | Required | Interpretation |
| --- | --- | --- | --- |
| `combine_feature_parquet` | `merge_feature` output | Yes | Integrated feature matrix used for model inference; schema must match model expectations. |
| `model_dir` | top-level config | Yes | Directory containing trained model artifacts for prediction. |
| `model_name` | top-level config | Yes | Selected model artifact name/version under `model_dir`. |

## Parameters

### Effective runtime keys (current implementation)

| Parameter | Location | Type | Interpretation |
| --- | --- | --- | --- |
| `model_dir` | top-level config | path string | Directory containing trained model artifacts used for inference. |
| `model_name` | top-level config | string | Model identifier/name loaded from `model_dir`. |
| `random_seed` | `steps.mutation_prediction` | integer | Random seed for reproducible prediction-related operations. |
| `plot` | `steps.mutation_prediction` | boolean | Whether to generate mutation-prediction plots. |

Common pretrained model names shipped in `SpaceTracer/models`:

- `spatial_free_model`
- `spatial_preserved_model`

## Tuning notes

- Ensure model artifacts match the feature schema in `combine_feature`.
- Keep model versioning explicit (`model_dir` + `model_name`) for reproducibility.

## Outputs

- `raw_pred_vcf`: `<step_dir>/results/<sample>_<model_name>_total_pred_truesites.vcf`
- `final_vcf`: `<step_dir>/results/<sample>_<model_name>_total_pred_truesites_PASS.vcf`
- `final_mutation_list`: `<step_dir>/results/<sample>_<model_name>_total_pred_truesites_PASS_mutation_list.txt`
- mitochondrial mutation list: `<step_dir>/results/<sample>_<model_name>_total_pred_truesites_MITO_mutation_list.txt`

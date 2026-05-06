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
| `model_dir` | config (`steps.mutation_prediction`) | Yes | Directory containing trained model artifacts for prediction. |
| `model_name` | config (`steps.mutation_prediction`) | Yes | Selected model artifact name/version under `model_dir`. |

## Parameters

### Effective runtime keys (current implementation)

| Parameter | Location | Type | Interpretation |
| --- | --- | --- | --- |
| `model_dir` | `steps.mutation_prediction` | path string | Directory containing trained model artifacts used for inference. |
| `model_name` | `steps.mutation_prediction` | string | Model identifier/name loaded from `model_dir`. |

Common pretrained model names shipped in `SpaceTracer_new_github/models`:

- `spatial_free_model`
- `spatial_feature_preserved_model`

### Present in template but not fully wired in current step runner

| Parameter | Location | Note |
| --- | --- | --- |
| `random_seed` | `steps.mutation_prediction` | Template key exists; current step uses internal constant seed. |
| `plot` | `steps.mutation_prediction` | Template key exists; current step uses internal plotting setting. |

## Tuning notes

- Ensure model artifacts match the feature schema in `combine_feature`.
- Keep model versioning explicit (`model_dir` + `model_name`) for reproducibility.
- If extending this step, expose internal hardcoded options as config keys progressively.

## Outputs

- `raw_pred_vcf`: `<step_dir>/results/Sample_total_pred_truesites.vcf`
- `final_vcf`: `<step_dir>/results/Sample_total_pred_truesites_PASS.vcf`

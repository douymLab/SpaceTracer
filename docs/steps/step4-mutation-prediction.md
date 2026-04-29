# Step 4: Mutation prediction

This stage runs model-based classification on merged features and exports predicted mutation loci.

## Current status

`mutation_prediction` is implemented but not enabled in the active default DAG unless explicitly integrated in the run path.

## Inputs

- `combine_feature_parquet` (for example `all_feature.parquet` from `merge_feature`)

### Input interpretation

| Input | Interpretation |
| --- | --- |
| `combine_feature_parquet` | Integrated feature matrix used by mutation classifier inference. |
| `model_dir` / `model_name` | Model artifacts selected for prediction. |

## Outputs

- raw prediction VCF
- filtered/pass VCF

## Detailed reference

- [mutation_prediction](mutation-prediction.md)

# Step 4: Mutation prediction

This stage runs model-based classification on merged features and exports predicted mutation loci.

## Current status

`mutation_prediction` is implemented but not enabled in the active default DAG unless explicitly integrated in the run path.

## Inputs

- `combine_feature` (for example `all_feature.txt`)

## Outputs

- raw prediction VCF
- filtered/pass VCF

## Detailed reference

- [mutation_prediction](mutation-prediction.md)

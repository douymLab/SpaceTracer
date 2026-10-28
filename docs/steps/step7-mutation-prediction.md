# Step 7: Mutation Prediction

## Code file

- `SpaceTracer/steps/step7_mutation_prediction.py`

## Runtime step name

- `mutation_prediction`

## What this step does

- loads merged features from Step 6
- runs each model named in `model_name` (`model_dir` + model name)
- exports one VCF and mutation-list set per model

## Detailed reference page

- [mutation_prediction](mutation-prediction.md)

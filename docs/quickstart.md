# Quick Start

This page gives the fastest path to run SpaceTracer on your data.

## 1) Prepare environment

Complete [Installation](installation.md) first.

## 2) Prepare inputs and references

At minimum you need:

- aligned BAM file from spatial transcriptomics
- tissue position file (for Visium)
- reference genome FASTA and required resources from [Resources](resources.md)
- writable output directory

## 3) Edit config file

Create `config.yaml` from the template in [Configuration](configuration.md), then update all sample-specific paths and run parameters.

Pretrained mutation-prediction models are already provided under `SpaceTracer_new_github/models`:

- `spatial_free_model`
- `spatial_feature_preserved_model`

Set these via `model_dir` and `model_name` in your config.
For full parameter/input meaning, use [Config Reference](config-reference.md).

!!! warning
    SpaceTracer validates input/resource paths on startup. Every configured file path must exist before running.

## 4) Run full workflow

```bash
SpaceTracer run --config config.yaml
```

Fallback command:

```bash
python -m SpaceTracer.cli.run --config config.yaml
```

## 5) Resume or rerun when needed

```bash
SpaceTracer run --config config.yaml --start-from genotyping
SpaceTracer run --config config.yaml --stop-at merge_feature
SpaceTracer run --config config.yaml --only-steps "RNA_feature,phasing,merge_feature" --force
SpaceTracer run --config config.yaml --force
```

Notes:

- `--start-from <step>` skips earlier completed steps and resumes from the specified step.
- `--stop-at <step>` runs up to that step and then exits.
- `--force` reruns requested steps even if checkpoint metadata says they are complete.
- `--only-steps` runs exactly the listed steps (topologically ordered within that subset) and does **not** auto-include external dependencies.
- If required upstream outputs are missing, `--only-steps` will fail with missing-key or missing-file errors.

## 6) Check outputs and step details

Use these pages for interpretation and debugging:

- [Outputs](outputs.md)
- [Step-by-step guide](steps/overview.md) (code-aligned with `SpaceTracer/steps`)
- [Single-Step Debug Cookbook](steps/debug-cookbook.md)

Common output files from downstream feature and prediction stages include:

- `output_dir/all_feature.txt`
- `output_dir/all_feature.parquet`
- `output_dir/mutation_prediction/Sample_total_pred_truesites.vcf`
- `output_dir/mutation_prediction/Sample_total_pred_truesites_PASS.vcf`

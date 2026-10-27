# Quick Start

This page gives the fastest path to run SpaceTracer on your data.

## 1) Prepare environment

Complete [Installation](installation.md) first.

Make sure you have installed the package into your environment (for CLI usage):

```bash
pip install .
spacetracer --help
```

## 2) Prepare inputs and references

At minimum you need:

- aligned BAM file from spatial transcriptomics
- tissue position file (for Visium)
- reference genome FASTA and required resources from [Resources](resources.md)
- writable output directory

## 3) Edit config file

Create `config.yaml` from the templates described in [Configuration](configuration.md), then update all sample-specific paths and run parameters.

For a standard Visium run, start from `SpaceTracer/config/config.example.yaml`. For Visium-HD, manual input paths, or external cluster/cell annotations, start from `SpaceTracer/config/config.advanced.example.yaml`.

Pretrained mutation-prediction models are provided under `SpaceTracer/models`:

- `spatial_free_model`
- `spatial_preserved_model`

Set these via `model_dir` and `model_name` in your config.
For step-level default parameters, use [Default Step Config](default-step-config.md).

!!! warning
    SpaceTracer validates input/resource paths on startup. Every configured file path must exist before running.

## 4) Run full workflow

```bash
spacetracer run --config config.yaml
```

Fallback command:

```bash
python -m SpaceTracer.cli.run --config config.yaml
```

## 5) Resume or rerun when needed

```bash
spacetracer run --config config.yaml --start-from genotyping
spacetracer run --config config.yaml --stop-at merge_feature
spacetracer run --config config.yaml --stop-at mutation_prediction
spacetracer run --config config.yaml --only-steps "RNA_feature,phasing,merge_feature" --force
spacetracer run --config config.yaml --force
```

Notes:

- `--start-from <step>` skips earlier completed steps and resumes from the specified step.
- `--stop-at <step>` runs up to that step and then exits.
- `--force` reruns requested steps even if checkpoint metadata says they are complete.
- `--only-steps` runs exactly the listed steps (topologically ordered within that subset) and does **not** auto-include external dependencies.
- If required upstream outputs are missing, `--only-steps` will fail with missing-key or missing-file errors.
- A full run continues through `phylogeny`. Use `--stop-at mutation_prediction` if you only want mutation VCF output and want to skip tree building.

Available step names:

```text
cluster, bam_processing, mpileup, umi_combine, cell_num, prior,
genotyping, spatial_feature, mappability_feature, read_feature,
RNA_feature, phasing, merge_feature, mutation_prediction, phylogeny
```

## 6) Check outputs and step details

Use these pages for interpretation and debugging:

- [Outputs](outputs.md)
- [Step-by-step guide](steps/overview.md) (code-aligned with `SpaceTracer/steps`)
- [Single-Step Debug Cookbook](steps/debug-cookbook.md)

Common output files from downstream feature and prediction stages include:

- `output_dir/all_feature.txt`
- `output_dir/all_feature.parquet`
- `output_dir/mutation_prediction/results/<sample>_<model_name>_total_pred_truesites.vcf`
- `output_dir/mutation_prediction/results/<sample>_<model_name>_total_pred_truesites_PASS.vcf`
- `output_dir/mutation_prediction/results/<sample>_<model_name>_total_pred_truesites_PASS_mutation_list.txt`
- `output_dir/phylogeny/tree/mutation_integrator/phylo/final_cleaned_M_full_basedPivots.filtered_sites_inferred.tree_scphylo.pdf`

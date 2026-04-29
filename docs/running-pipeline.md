# Running SpaceTracer End-to-End

This tutorial shows a practical way to run SpaceTracer with your own sample and inspect progress.

## What you need first

Complete the setup in:

- [Installation](installation.md)
- [Resources](resources.md)

You also need:

- one input BAM file from your spatial transcriptomics run
- a tissue position file (for Visium)
- a writable output directory
- required reference resources (genome, dbSNP, mappability, etc.)

## 1) Prepare a working folder

```bash
mkdir -p tutorial_run/{input,resources,output}
cd tutorial_run
```

Place your files in:

- `input/`: BAM and tissue position file
- `resources/`: reference files

## 2) Create a sample configuration file

Create `config.yaml` in your run directory.

You can start from the template on the [Configuration](configuration.md) page and update all paths.

!!! warning
    SpaceTracer validates input/resource paths on startup. Every configured file path must exist before running.

## 3) Run the full pipeline

Run with your config:

```bash
SpaceTracer run --config config.yaml
```

If your environment does not expose the `SpaceTracer` executable yet, use:

```bash
python -m SpaceTracer.cli.run --config config.yaml
```

## 4) Resume, rerun, or run partial workflow

SpaceTracer tracks completed steps and can skip them automatically.

### Force rerun all requested steps

```bash
SpaceTracer run --config config.yaml --force
```

### Start from a specific step

```bash
SpaceTracer run --config config.yaml --start-from genotyping
```

### Stop at a specific step

```bash
SpaceTracer run --config config.yaml --stop-at merge_feature
```

The valid step names are documented in [Algorithm Walkthrough](algorithm.md).
Detailed per-step inputs/outputs are available in [Step Reference Overview](steps/overview.md).
Debug rerun recipes are available in [Single-Step Debug Cookbook](steps/debug-cookbook.md).

## 5) Check run outputs

The pipeline writes intermediate and merged results under your configured `output_dir`.

Typical artifacts include:

- step-level intermediate files (for pileup, UMI combine, genotyping, and features)
- checkpoint metadata for completed steps
- merged feature matrix files for downstream mutation prediction and interpretation

See [Outputs & Troubleshooting](outputs.md) for output interpretation and common failure modes.

## 6) Recommended first QC checks

After a successful run:

- confirm expected sample/barcode counts in intermediate summaries
- inspect candidate variant counts after `genotyping`
- inspect merged feature file dimensions after `merge_feature`
- verify no required step was skipped unexpectedly

## 7) Example development workflow

For parameter tuning, iterate with partial runs:

1. run once with full pipeline
2. adjust a small set of parameters in `config.yaml`
3. rerun from an affected step using `--start-from`
4. compare candidate/feature statistics across runs

This gives much faster iteration than rerunning from scratch each time.

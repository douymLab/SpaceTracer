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
SpaceTracer run --config config.yaml --force
```

## 6) Check outputs and step details

Use these pages for interpretation and debugging:

- [Outputs](outputs.md)
- [Step-by-step guide](steps/overview.md)
- [Single-Step Debug Cookbook](steps/debug-cookbook.md)

# Quick Start

This page gives the fastest path to your first SpaceTracer run.

## 1) Install and prepare resources

Complete:

- [Installation](installation.md)
- [Resources](resources.md)

## 2) Create a config file

Create `config.yaml` using the minimal template from [Configuration](configuration.md).

At minimum, provide valid paths for:

- input BAM and tissue position file (`input_details`)
- all required references (`resource_details`)
- `output_dir`

## 3) Run SpaceTracer

```bash
SpaceTracer run --config config.yaml
```

Fallback command:

```bash
python -m SpaceTracer.cli.run --config config.yaml
```

## 4) Resume or rerun if needed

```bash
SpaceTracer run --config config.yaml --start-from genotyping
SpaceTracer run --config config.yaml --force
```

## 5) Understand outputs

For result interpretation and debugging, see:

- [Running Tutorial](running-pipeline.md)
- [Algorithm Walkthrough](algorithm.md)
- [Outputs & Troubleshooting](outputs.md)

# Quick Start

This page gives the fastest path to run SpaceTracer on demo or real data.

## 1) Prepare environment

Complete [Installation](installation.md) first.

## 2) Prepare inputs and references

At minimum you need:

- filtered or raw aligned BAM (from spatial transcriptomics)
- tissue position file (Visium)
- reference genome FASTA
- annotation/resources listed in [Resources](resources.md)

## 3) Edit config file

Start from:

- `config_example.txt`

Update all paths and sample-specific settings described in [Configuration](configuration.md).

## 4) Run full workflow

```bash
snakemake --configfile config_example.txt -s run_snakemake_for_conda
```

For Docker workflows:

```bash
snakemake --configfile config_example.txt -s run_snakemake_for_docker
```

## 5) Check final output

Main final result:

`$savePATH/$sample/$model/pred.FINAL.txt`

See [Outputs](outputs.md) for important intermediate files, predicted-site format examples, and visualization/validation examples (lineage tree, spatial scatter, and IGV views).

# SpaceTracer

SpaceTracer is a computational framework for detecting somatic single-nucleotide variants (SNVs) directly from spatial transcriptomics data and reconstructing cellular lineages within native tissue architecture.

![flowchart](figures/flowchart.png)

## Overview

Spatial transcriptomics reveals tissue organization, but lineage-tracing methods applicable to human tissues remain limited. SpaceTracer addresses this challenge by leveraging naturally occurring somatic SNVs to enable perturbation-free lineage tracing directly from spatial transcriptomics data.

SpaceTracer supports the analysis of both **nuclear SNVs** and **mitochondrial SNVs**, and enables the study of lineage spread, migration, and lineage-associated transcriptional changes in situ.


## Release notes

* **2026-02-25 · Version 1.1.0**
  Updated genotype calculation, improved random forest features, and added additional filtering steps to improve accuracy.

* **2025-03-31 · Version 1.0.0**
  Initial release of SpaceTracer.

## Get started

To start using SpaceTracer:

1. Follow [Installation](installation.md)
2. Read [Quick Start](quickstart.md)
3. Run your first dataset with [Running Tutorial](running-pipeline.md)
4. Understand method details in [Algorithm Walkthrough](algorithm.md)
5. Explore each stage in [Step Reference Overview](steps/overview.md)
6. Use practical rerun recipes in [Single-Step Debug Cookbook](steps/debug-cookbook.md)

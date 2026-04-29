# Step 1: Quality control

This stage converts raw pileup evidence into quality-controlled counting and prior inputs used by genotyping.

---

## Step 1.1 - UMI-aware counting (`umi_combine`)

**What it does:** aggregates read-level evidence at UMI/locus level to reduce technical noise.

**Detailed reference:** [umi_combine](umi-combine.md)

---

## Step 1.2 - Cell-number estimation (`cell_num`)

**What it does:** builds cell-count support values used by downstream genotype modeling.

**Detailed reference:** [cell_num](cell-num.md)

---

## Step 1.3 - Population prior construction (`prior`)

**What it does:** builds prior information for candidate loci from resource datasets.

**Detailed reference:** [prior](prior.md)

## Typical outputs

- `spot_count_file`, `error_count_file`
- `cell_num`
- `prior_file`

## Input/parameter interpretation entry points

- `umi_combine` inputs/parameters: [umi_combine](umi-combine.md)
- `cell_num` inputs/parameters: [cell_num](cell-num.md)
- `prior` inputs/parameters: [prior](prior.md)

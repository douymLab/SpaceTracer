# `cell_num`

## Purpose

Computes or forwards spot-level cell number information used in spot genotype inference.

## Upstream

- `cluster`
- `umi_combine`

## Required inputs

- `in_filter_bam`
- `cluster` information from config/context
- `sequence_type` and optional `bins`

### Input interpretation

| Input key | Source | Required | Interpretation |
| --- | --- | --- | --- |
| `in_filter_bam` | `bam_processing` output | Yes | Filtered BAM used to estimate or validate spot-level cell support. |
| `cluster` | `cluster` step output/file | Yes | Spot grouping used for cell-number modeling context. |
| `sequence_type` | top-level config | Yes | Controls mode-specific behavior (Visium vs non-Visium). |
| `bins` (`bin_size`) | top-level config | Conditional | Used in non-Visium/binned workflows when applicable. |

## Parameters

This step is mostly controlled by:

- `steps.cell_number` (integer or file path)
  - `0`: compute from BAM/UMI profiles
  - positive integer: use constant
  - file path: use provided table

### Parameter interpretation

| Parameter | Type | Typical/default | Interpretation |
| --- | --- | --- | --- |
| `steps.cell_number` | int/path/null | unset or `0` | Main mode switch: compute, fixed value, or passthrough external table. |

## Outputs

- `cell_num` context key (consumed by `genotyping`)
  - generated file: `output_dir/refined_cell_num.txt` (columns: `barcode`, `cluster`, `UMI_counts`, `refined_cell_num`), or
  - passthrough integer/path from config when `steps.cell_number` is a fixed value or external table

## Tuning notes

- If cell number is already known, pass it directly to avoid recomputation.

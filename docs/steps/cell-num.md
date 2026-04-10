# `cell_num`

## Purpose

Computes or forwards spot-level cell number information used in spot genotype inference.

## Upstream dependencies

- `cluster`
- `umi_combine`

## Required inputs

- `in_filter_bam`
- `cluster` information from config/context
- `sequence_type` and optional `bins`

## Key parameters

This step is mostly controlled by:

- `steps.cell_number` (integer or file path)
  - `0`: compute from BAM/UMI profiles
  - positive integer: use constant
  - file path: use provided table

## Outputs

- `cell_num` key
  - generated file: `output_dir/refined_umi_read_cellNum.txt`, or
  - passthrough integer/path from config

## Notes

- If cell number is already known, pass it directly to avoid recomputation.

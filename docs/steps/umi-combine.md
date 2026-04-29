# `umi_combine`

## Purpose

Aggregates read-level evidence into UMI-level counts across candidate loci.

## Upstream dependencies

- `mpileup`

## Required inputs

- `in_filter_bam`
- `filter_mpileup_file`
- `db_path` (chunk database from `mpileup`)
- `sequence_type`

## Key parameters

Most controls are currently internal defaults in the step implementation (for example chunk-level parallel workers and buffer sizes), while threads are inherited from `run.threads`.

## Outputs

- `spot_count_file`: `output_dir/umi_combine/spot.count.parquet`
- `error_count_file`: `output_dir/umi_combine/error.count.parquet`

## Notes

- This step is parallelized by chunk files loaded from `split_chunk.db`.
- Outputs are parquet and used by both `genotyping` and `RNA_feature`.

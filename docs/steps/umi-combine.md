# `umi_combine`

## Purpose

Aggregates read-level evidence into UMI-level counts across candidate loci.

## Upstream

- `mpileup`

## Required inputs

- `in_filter_bam`
- `filter_mpileup_file`
- `db_path` (chunk database from `mpileup`)
- `sequence_type`

### Input interpretation

| Input key | Source | Required | Interpretation |
| --- | --- | --- | --- |
| `in_filter_bam` | `bam_processing` output | Yes | Filtered BAM used to aggregate read evidence into UMI-level counts. |
| `filter_mpileup_file` | `mpileup` output | Yes | Candidate loci table that defines UMI aggregation targets. |
| `db_path` / chunk manifest | `mpileup` output | Yes | Chunk metadata controlling parallel chunk processing. |
| `sequence_type` | top-level config | Yes | Controls mode-specific handling during evidence aggregation. |

## Parameters

Most controls are currently internal defaults in the step implementation (for example chunk-level parallel workers and buffer sizes), while threads are inherited from `run.threads`.

### Parameter interpretation

| Parameter area | Interpretation |
| --- | --- |
| `run.threads` / runtime parallel settings | Controls chunk-level parallel processing throughput. |
| internal filtering flags | Duplicate/secondary/QC-fail/supplementary filtering behavior is defined in implementation defaults unless exposed in your version config. |

## Outputs

- `spot_count_file`: `output_dir/umi_combine/spot.count.parquet`
- `error_count_file`: `output_dir/umi_combine/error.count.parquet`

## Tuning notes

- This step is parallelized by chunk files loaded from `split_chunk.db`.
- Outputs are parquet and used by both `genotyping` and `RNA_feature`.

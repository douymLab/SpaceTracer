# `bam_processing`

## Purpose

Filters and prepares BAM input for pileup generation.

## Upstream

None (DAG root step).

## Required inputs

- `bam_file` (from config/input resolution)
- optional `tissue_position` for in-tissue barcode extraction
- `barcode_key` (for example `CB`)

### Input interpretation

| Input key | Source | Required | Interpretation |
| --- | --- | --- | --- |
| `bam_file` | `input_details.bam_file` or `spaceranger_dir` | Yes | Primary aligned BAM. Missing file aborts step execution. |
| `tissue_position` | `input_details.tissue_position` or `spaceranger_dir/spatial/*` | Conditional | Required for Visium in-tissue barcode filtering; otherwise this stage can pass BAM through. |
| `barcode_key` | `input_details.barcode_key` | Yes | BAM tag containing barcode identity used by barcode-based filtering logic. |

## Parameters

From `steps.bam_processing`:

- `nm_threshold`: max mismatch tag (`nM`) allowed
- `mapq_threshold`: minimum MAPQ

### Parameter interpretation

| Parameter | Type | Typical/default | Interpretation |
| --- | --- | --- | --- |
| `nm_threshold` | int | `5` (template) | Maximum allowed mismatch count per read; lower values are stricter. |
| `mapq_threshold` | int | `255` (template) | Minimum mapping quality threshold; lower values are more permissive. |

## Outputs

- `in_bam`: `output_dir/bam_processing/IN.bam`
- `in_filter_bam`: `output_dir/bam_processing/IN_filter.bam`

## Tuning notes

- If tissue positions are provided, `sinto filterbarcodes` is used to keep in-tissue barcodes.
- Otherwise, raw BAM is symlinked as `IN.bam`.
- Filtered BAM is indexed for downstream steps.

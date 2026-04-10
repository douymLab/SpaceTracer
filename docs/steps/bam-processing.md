# `bam_processing`

## Purpose

Filters and prepares BAM input for pileup generation.

## Upstream dependencies

None (DAG root step).

## Required inputs

- `bam_file` (from config/input resolution)
- optional `tissue_position` for in-tissue barcode extraction
- `barcode_key` (for example `CB`)

## Key parameters

From `steps.bam_processing`:

- `nm_threshold`: max mismatch tag (`nM`) allowed
- `mapq_threshold`: minimum MAPQ

## Outputs

- `in_bam`: `output_dir/bam_processing/IN.bam`
- `in_filter_bam`: `output_dir/bam_processing/IN_filter.bam`

## Notes

- If tissue positions are provided, `sinto filterbarcodes` is used to keep in-tissue barcodes.
- Otherwise, raw BAM is symlinked as `IN.bam`.
- Filtered BAM is indexed for downstream steps.

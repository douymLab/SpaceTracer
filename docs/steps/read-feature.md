# `read_feature`

## Purpose

Computes read-level quality and bias features for each candidate mutation.

## Upstream dependencies

- `genotyping`

## Required inputs

- `bam_file`
- `ind_geno_filter_mutation_list`

## Key parameters

From `steps.read_feature`:

- `cell_info` (optional barcode-to-cell mapping)
- `downsample`
- `downsample_target_depth`
- `max_region_size`
- `max_variants_per_region`
- `seed`

## Outputs

- `read_feature`: `output_dir/read_feature/read_feature.txt`
- parquet mirror: `read_feature.parquet`

## Notes

- Regions are batched and processed in multiprocessing workers.
- Feature extraction uses mutation-local read information and produces one row per candidate.

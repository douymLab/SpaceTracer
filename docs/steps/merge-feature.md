# `merge_feature`

## Purpose

Merges RNA, spatial, read, and mappability feature tables into one integrated feature matrix and computes filtration tags.

## Upstream dependencies

- `spatial_feature`
- `mappability_feature`
- `read_feature`
- `RNA_feature`

## Required inputs

- `RNA_feature`
- `spatial_feature`
- `read_feature`
- `mappability_feature`

## Key parameters

From `steps.feature_filtration` (group-based enablement):

- enabled filtration groups or defaults
- implied thresholds in step logic (for example hFDR, mismatch bias, mappability, ALT count, population AF)

## Outputs

- `combine_feature`: `output_dir/all_feature.txt`
- parquet mirror: `output_dir/all_feature.parquet`

## Notes

- Multi-index keys are `#chrom`, `pos`, `ref`, `alt`.
- Filtration summary is consolidated into the `Filtration` column (`PASS` or semicolon-separated tags).

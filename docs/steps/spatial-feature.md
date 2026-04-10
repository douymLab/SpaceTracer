# `spatial_feature`

## Purpose

Extracts mutation-level spatial consistency and neighborhood-distribution features from spot genotypes.

## Upstream dependencies

- `genotyping`

## Required inputs

- `ind_geno_filter_mutation_list`
- `spot_geno_file`
- `tissue_position`

## Key parameters

From `steps.spatial_feature`:

- `alpha`
- `thr_r2`
- `thr_prob`
- `thr_likelihood`
- `thr_vaf`
- `plot_supp`
- `fig_size`
- `method`
- `num_directions`

## Outputs

- `spatial_feature`: `output_dir/spatial_feature/spatial_feature.txt`
- parquet mirror: `spatial_feature.parquet`

## Notes

- This step is parallelized across mutation identifiers.
- Spatial metrics include distribution tests and VAF consistency summaries.

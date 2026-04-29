# `spatial_feature`

## Purpose

Extracts mutation-level spatial consistency and neighborhood-distribution features from spot genotypes.

## Upstream

- `genotyping`

## Required inputs

- `ind_geno_filter_mutation_list`
- `spot_geno_file`
- `tissue_position`

### Input interpretation

| Input key | Source | Required | Interpretation |
| --- | --- | --- | --- |
| `ind_geno_filter_mutation_list` | `genotyping` output | Yes | Candidate mutation list used to iterate spatial-feature extraction targets. |
| `spot_geno_file` | `genotyping` output | Yes | Spot-level genotype/probability table for spatial consistency tests. |
| `tissue_position` | input details | Yes (Visium) | Spot coordinate matrix required for neighborhood and directional spatial metrics. |

## Parameters

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

### Parameter interpretation highlights

| Parameter | Interpretation |
| --- | --- |
| `alpha` | Significance threshold used in statistical spatial tests. |
| `thr_r2`, `thr_prob`, `thr_likelihood`, `thr_vaf` | Core cutoffs that define mutation-spatial consistency criteria. |
| `method` | Spatial modeling method (for example LDA-style directional aggregation). |
| `num_directions` | Direction granularity for directional spatial summaries. |
| `plot_supp`, `fig_size` | Controls supplementary plotting behavior and figure size. |

## Outputs

- `spatial_feature`: `output_dir/spatial_feature/spatial_feature.txt`
- parquet mirror: `spatial_feature.parquet`

## Tuning notes

- This step is parallelized across mutation identifiers.
- Spatial metrics include distribution tests and VAF consistency summaries.

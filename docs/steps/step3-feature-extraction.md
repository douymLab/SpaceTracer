# Step 3: Feature extraction

This stage computes multi-source features for genotype-supported loci and merges them into one feature matrix.

## Feature branches

- `spatial_feature`: neighborhood and tissue-structure evidence
- `mappability_feature`: genomic mappability/confounder evidence
- `read_feature`: read-level quality and bias evidence
- `RNA_feature`: transcript-level contextual evidence

## Integration step

- `phasing`: RNA-informed phasing refinement before final merge
- `merge_feature`: merges all branch outputs and emits `all_feature.txt`

## Inputs

- Step 2 genotype outputs (`spot_geno_file`, `ind_geno_filter_file`, related context)
- branch-specific resources (for example mappability/reference annotations)

### Sub-function input/parameter map

| Sub-function | Critical inputs | Key parameters |
| --- | --- | --- |
| `spatial_feature` | `ind_geno_filter_mutation_list`, `spot_geno_file`, `tissue_position` | `alpha`, `thr_r2`, `thr_prob`, `thr_likelihood`, `thr_vaf`, `method`, `num_directions` |
| `mappability_feature` | `ind_geno_filter_mutation_list`, `mappability_path` | mostly resource-driven (few exposed thresholds) |
| `read_feature` | `bam_file`, `ind_geno_filter_mutation_list` | `cell_info`, `downsample`, `downsample_target_depth`, `max_region_size`, `max_variants_per_region`, `seed` |
| `RNA_feature` | `germline_file`, `ind_geno_filter_file`, `error_count_file`, RNA/filter resources | `min_count_for_germline`, `min_prior_for_germline`, `default_range_of_gene`, `p_threshold`, `previous_base` |
| `phasing` | `in_filter_bam`, `merged_germline_file`, `merged_ind_geno_filter_file` | `minprior`, `min_dp`, `min_total_dp`, `alpha`, `phasing_pad`, `merge_gap`, `max_target`, `seed` |
| `merge_feature` | branch feature outputs + `phasing_result` + `cluster_event_result` | `steps.feature_filtration.*` group switches |

See each branch page for exact semantics:

## Outputs

- branch feature tables
- `output_dir/all_feature.txt`
- `output_dir/all_feature.parquet`

## Detailed references

- [spatial_feature](spatial-feature.md)
- [mappability_feature](mappability-feature.md)
- [read_feature](read-feature.md)
- [RNA_feature](rna-feature.md)
- [phasing](phasing.md)
- [merge_feature](merge-feature.md)

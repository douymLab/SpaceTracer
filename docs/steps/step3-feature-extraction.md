# Step 3: Feature extraction

This stage computes multi-source features for genotype-supported loci and merges them into one feature matrix.

## Feature branches

- `spatial_feature`: neighborhood and tissue-structure evidence
- `mappability_feature`: genomic mappability/confounder evidence
- `read_feature`: read-level quality and bias evidence
- `RNA_feature`: transcript-level contextual evidence

## Integration step

- `merge_feature`: merges all branch outputs and emits `all_feature.txt`

## Inputs

- Step 2 genotype outputs (`spot_geno_file`, `ind_geno_filter_file`, related context)
- branch-specific resources (for example mappability/reference annotations)

## Outputs

- branch feature tables
- `output_dir/all_feature.txt`
- `output_dir/all_feature.parquet`

## Detailed references

- [spatial_feature](spatial-feature.md)
- [mappability_feature](mappability-feature.md)
- [read_feature](read-feature.md)
- [RNA_feature](rna-feature.md)
- [merge_feature](merge-feature.md)

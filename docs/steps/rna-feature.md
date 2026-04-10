# `RNA_feature`

## Purpose

Builds RNA-context features and annotations (ASE, hFDR, editing/imprinted/PON tags, sequence-context metrics).

## Upstream dependencies

- `genotyping`

## Required inputs

- `germline_file`
- `ind_geno_filter_file`
- `error_count_file`
- resources: `gene_bed`, `dbsnp_vcf_file`, `imprinted_bed`, `editing_bed`, `PON_file`, `genome_fasta`, `reference_error_profile`

## Key parameters

From `steps.RNA_feature`:

- `min_count_for_germline`
- `min_prior_for_germline`
- `default_range_of_gene`
- `p_threshold`
- `previous_base`

## Outputs

- `RNA_feature`: `output_dir/RNA_feature/RNA_feature.txt`
- parquet mirror: `RNA_feature.parquet`

## Notes

- This step merges several RNA-oriented annotation and statistical modules.
- Output fields are heavily used in final merged filtration rules.

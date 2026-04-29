# `RNA_feature`

## Purpose

Builds RNA-context features and annotations (ASE, hFDR, editing/imprinted/PON tags, sequence-context metrics).

## Upstream

- `genotyping`

## Required inputs

- `germline_file`
- `ind_geno_filter_file`
- `error_count_file`
- resources: `gene_bed`, `dbsnp_vcf_file`, `imprinted_bed`, `editing_bed`, `PON_file`, `genome_fasta`, `reference_error_profile`

### Input interpretation

| Input key | Source | Required | Interpretation |
| --- | --- | --- | --- |
| `germline_file` | `genotyping` output | Yes | Germline context used for RNA-level contrast and filtering logic. |
| `ind_geno_filter_file` | `genotyping` output | Yes | Candidate loci table entering RNA annotation and tests. |
| `error_count_file` | `umi_combine` output | Yes | Error-profile/count evidence used in RNA-context filtering features. |
| RNA/filter resources | `resource_details.*` | Yes | Annotation resources for dbSNP/editing/imprinted/PON/reference-error-derived features. |

## Parameters

From `steps.RNA_feature`:

- `min_count_for_germline`
- `min_prior_for_germline`
- `default_range_of_gene`
- `p_threshold`
- `previous_base`

### Parameter interpretation highlights

| Parameter | Interpretation |
| --- | --- |
| `min_count_for_germline`, `min_prior_for_germline` | Controls germline evidence sufficiency for RNA-context decisions. |
| `default_range_of_gene` | Window size for gene-context assignment when annotation ambiguity exists. |
| `p_threshold` | Main significance threshold for RNA-derived statistical tests. |
| `previous_base` | Sequence-context offset used in motif/context feature extraction. |

## Outputs

- `RNA_feature`: `output_dir/RNA_feature/RNA_feature.txt`
- parquet mirror: `RNA_feature.parquet`

## Tuning notes

- This step merges several RNA-oriented annotation and statistical modules.
- Output fields are heavily used in final merged filtration rules.

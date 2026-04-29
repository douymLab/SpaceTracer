# Step 2: Genotyping

This stage infers mutation evidence at individual, cluster, and spot levels from quality-controlled counts and priors.

Script/module stage: `genotyping`

## What it does

- combines allele counts with prior information
- applies confidence and filtering thresholds
- emits genotype outputs used by all feature branches

## Inputs

- count outputs from Step 1
- `cluster` context
- `prior_file`
- `cell_num`

## Outputs

- `ind_geno_filter_file`
- `ind_geno_filter_mutation_list`
- `germline_file`
- `cluster_vaf_file`
- `spot_geno_file`

## Detailed reference

- [genotyping](genotyping.md)

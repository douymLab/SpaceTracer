# `mappability_feature`

## Purpose

Annotates candidate loci with mappability information to flag low-uniqueness regions.

## Upstream dependencies

- `genotyping`

## Required inputs

- `ind_geno_filter_mutation_list`
- `mappability_path`

## Key parameters

No major per-step thresholds are currently exposed; performance mainly scales with `run.threads`.

## Outputs

- `mappability_feature`: `output_dir/mappability_feature/mappability_feature.txt`
- parquet mirror: `mappability_feature.parquet`

## Notes

- Internally processed by chromosome groups.
- Output includes `mappabilityScore` used in final filtration logic.

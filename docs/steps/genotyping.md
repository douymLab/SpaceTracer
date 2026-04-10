# `genotyping`

## Purpose

Performs cluster-level and spot-level genotype inference from UMI counts and priors.

## Upstream dependencies

- `cluster`
- `prior`
- `cell_num`

## Required inputs

- `spot_count_file`
- `prior_file`
- `cluster` information
- `cell_num`

## Key parameters

From `steps.genotyping`:

- `alpha`
- `epsQ`
- `epsAF`
- `mu`
- `thr_dp`
- `pop_vaf`
- `filter_oneallele`

## Outputs

Main outputs:

- `ind_geno_filter_file`
- `ind_geno_filter_mutation_list`
- `germline_file`
- `cluster_vaf_file`
- `spot_geno_file`

When `run.keep_intermediates` is true, extra intermediate count/genotype files are also emitted.

## Notes

- This step is the key bridge from evidence aggregation to multi-feature extraction.
- Downstream feature steps all depend directly on these genotype outputs.

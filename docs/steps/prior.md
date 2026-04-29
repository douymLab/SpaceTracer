# `prior`

## Purpose

Computes base priors per candidate locus using gnomAD reference data.

## Upstream dependencies

- `umi_combine`

## Required inputs

- `gnomad_path`
- `filter_mpileup_file`
- genome autosome metadata from `genome_details`

## Key parameters

No major step-level tunables are exposed in config; this step uses query-driven prior lookup and normalizes per-base priors (`A/T/C/G`) per site.

## Outputs

- `prior_file`: `output_dir/prior.txt`

## Notes

- gnomAD input can be either a directory or a file list.
- Missing chromosome resources are tolerated, but reduce prior completeness.

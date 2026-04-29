# `prior`

## Purpose

Computes base priors per candidate locus using gnomAD reference data.

## Upstream

- `umi_combine`

## Required inputs

- `gnomad_path`
- `filter_mpileup_file`
- genome autosome metadata from `genome_details`

### Input interpretation

| Input key | Source | Required | Interpretation |
| --- | --- | --- | --- |
| `gnomad_path` | `resource_details.gnomad_path` | Yes | Population AF resource directory/list used for prior lookup. |
| `filter_mpileup_file` | `mpileup` output | Yes | Candidate-locus file used as prior-query target list. |
| `genome_details` | config loader | Yes | Chromosome metadata used for per-chromosome resource querying. |

## Parameters

No major step-level tunables are exposed in config; this step uses query-driven prior lookup and normalizes per-base priors (`A/T/C/G`) per site.

### Parameter interpretation

| Parameter area | Interpretation |
| --- | --- |
| `gnomad_path` content/coverage | Primary determinant of prior completeness and quality. |
| chromosome metadata (`genome_details`) | Controls which chromosomes are queried and merged into prior output. |
| `run.threads` (indirect) | Affects overall pipeline throughput; this step itself is mostly I/O/resource-lookup driven. |

## Outputs

- `prior_file`: `output_dir/prior.txt`

## Tuning notes

- gnomAD input can be either a directory or a file list.
- Missing chromosome resources are tolerated, but reduce prior completeness.

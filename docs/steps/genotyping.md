# `genotyping`

## Purpose

Performs cluster-level and spot-level genotype inference from UMI counts and priors.

## Upstream

- `cluster`
- `prior`
- `cell_num`

## Required inputs

- `spot_count_file`
- `prior_file`
- `cluster` information
- `cell_num`

### Input interpretation

| Input key | Source step/config | Required | Interpretation |
| --- | --- | --- | --- |
| `spot_count_file` | `umi_combine` manifest output | Yes | Chunk manifest pointing to spot-level count parquet files used for genotype inference. |
| `prior_file` | `prior` output (or fixed/empty mode) | Yes | Prior-frequency table for genotype calculations; behavior differs if fixed/empty prior is used. |
| `cluster` | `cluster` output (`cluster_file`) | Yes | Spot-to-cluster mapping used for cluster-level aggregation before individual calls. |
| `cell_num` | `cell_num` step output in pipeline context | Yes | Per-spot or constant cell-number support used in spot genotype inference (`refined_cell_num.txt` when `steps.cell_number` is unset or `0`, otherwise the configured integer or file path). |

## Parameters (`steps.genotyping`)

| Parameter | Type | Typical/default | Interpretation |
| --- | --- | --- | --- |
| `alpha` | float | `0.05` | Statistical significance threshold used in allele-level filtering logic. |
| `epsQ` | int | `20` | Quality-to-error conversion scale for UMI/read evidence aggregation. |
| `epsAF` | float | `0.003` | Allele-frequency error floor used during cluster allele filtering. |
| `mu` | float | `1e-5` | Prior mutation-rate term used in individual genotype inference. |
| `thr_dp` | int | `1000` | Depth threshold for robust genotype calling/retention. |
| `pop_vaf` | float | `1e-5` | Population-AF threshold used in genotype filtering logic. |
| `filter_oneallele` | bool | `true` | If true, applies one-allele style filtering for stricter genotype selection. |

## Tuning notes

- `alpha`, `epsAF`, and `mu` jointly control strictness of candidate retention.
- Raise `thr_dp` for more conservative calls on noisy/high-depth data.
- Lower `pop_vaf` for stricter rare-variant emphasis.
- Keep `filter_oneallele=true` unless you explicitly want a more permissive candidate set.

## Outputs

Main outputs:

- `ind_geno_filter_file`
- `ind_geno_filter_mutation_list`
- `germline_file`
- `cluster_vaf_file`
- `spot_geno_file`

When `run.keep_intermediates` is true, extra intermediate count/genotype files are also emitted.

## Tuning notes

- This step is the key bridge from Step 3 count/prior construction to Step 5 multi-feature extraction.
- Downstream feature steps all depend directly on these genotype outputs.

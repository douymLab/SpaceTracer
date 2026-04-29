# `mappability_feature`

## Purpose

Annotates candidate loci with mappability information to flag low-uniqueness regions.

## Upstream

- `genotyping`

## Required inputs

- `ind_geno_filter_mutation_list`
- `mappability_path`

### Input interpretation

| Input key | Source | Required | Interpretation |
| --- | --- | --- | --- |
| `ind_geno_filter_mutation_list` | `genotyping` output | Yes | Candidate mutation list used for mappability lookup. |
| `mappability_path` | `resource_details.mappability_path` | Yes | Mappability track/resource used to annotate locus uniqueness. |

## Parameters

No major per-step thresholds are currently exposed; performance mainly scales with `run.threads`.

### Parameter interpretation

| Parameter area | Interpretation |
| --- | --- |
| `run.threads` / runtime parallel controls | Main performance knob for chromosome-group processing throughput. |
| `mappability_path` resolution | Most behavior is resource-driven; choice/version of mappability resource strongly affects output scores. |

## Outputs

- `mappability_feature`: `output_dir/mappability_feature/mappability_feature.txt`
- parquet mirror: `mappability_feature.parquet`

## Tuning notes

- Internally processed by chromosome groups.
- Output includes `mappabilityScore` used in final filtration logic.

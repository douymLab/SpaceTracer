# `read_feature`

## Purpose

Computes read-level quality and bias features for each candidate mutation.

## Upstream

- `genotyping`

## Required inputs

- `bam_file`
- `ind_geno_filter_mutation_list`

### Input interpretation

| Input key | Source | Required | Interpretation |
| --- | --- | --- | --- |
| `bam_file` | `input_details.bam_file` | Yes | Read-level evidence source for mismatch/mapQ/baseQ/query-position features. |
| `ind_geno_filter_mutation_list` | `genotyping` output | Yes | Candidate locus list used to define per-region read feature extraction targets. |

## Parameters

From `steps.read_feature`:

- `cell_info` (optional barcode-to-cell mapping)
- `downsample`
- `downsample_target_depth`
- `max_region_size`
- `max_variants_per_region`
- `seed`

### Parameter interpretation highlights

| Parameter | Interpretation |
| --- | --- |
| `cell_info` | Optional barcode-to-cell map for cell-aware feature derivation. |
| `downsample`, `downsample_target_depth` | Controls depth normalization for robust comparisons across loci. |
| `max_region_size`, `max_variants_per_region` | Region partition controls for balancing runtime and memory. |
| `seed` | Reproducibility for stochastic operations (for example downsampling). |

## Outputs

- `read_feature`: `output_dir/read_feature/read_feature.txt`
- parquet mirror: `read_feature.parquet`

## Tuning notes

- Regions are batched and processed in multiprocessing workers.
- Feature extraction uses mutation-local read information and produces one row per candidate.

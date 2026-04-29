# `merge_feature`

## Purpose

Merges RNA, spatial, read, and mappability feature tables into one integrated feature matrix and computes filtration tags.

## Upstream

- `spatial_feature`
- `mappability_feature`
- `read_feature`
- `RNA_feature`

## Required inputs

- `RNA_feature`
- `spatial_feature`
- `read_feature`
- `mappability_feature`

### Input interpretation

| Input key | Source step | Required | Interpretation |
| --- | --- | --- | --- |
| `RNA_feature` | `RNA_feature` | Yes | RNA-level annotation feature table. |
| `spatial_feature` / `spatial_feature_results` | `spatial_feature` | Yes | Spatial feature table or chunk manifest to be merged. |
| `read_feature` / `read_feature_results` | `read_feature` | Yes | Read-level feature table or chunk manifest to be merged. |
| `mappability_feature` | `mappability_feature` | Yes | Mappability/confounder feature table. |
| `phasing_result` | `phasing` | Recommended | Phasing summary that enriches merged features before filtration. |
| `cluster_event_result` | `phasing` | Recommended | Cluster-event flags used by filtration grouping (`CLUSTER_EVENTS`). |

## Parameters (`steps.feature_filtration`)

This step uses group-based filtration switches. If the section is omitted, the default enabled-group set in code is applied.

| Group key | Interpretation |
| --- | --- |
| `ASE` | Filter loci flagged by allele-specific expression evidence. |
| `hFDR` | Filter loci with high hFDR-style artifact risk. |
| `imprinted` | Filter loci in imprinted-gene contexts. |
| `homopolymer` | Filter loci affected by homopolymer-related artifacts. |
| `PON` | Filter loci present in panel-of-normals evidence. |
| `RNA_editing` | Filter loci overlapping known RNA-editing contexts. |
| `ABNORMAL_MISMATCHES` | Filter loci with high/bias mismatch patterns. |
| `LOW_READ_DIVERSITY` | Filter loci with poor alternative-read diversity. |
| `HIGH_MULTIPLE_MAPPIN` | Filter loci with high multiple-mapping signal. |
| `WIDE_DISTRIBUTION` | Filter loci with broad suspicious mutant-probability distribution. |
| `NEAR_READ_END` | Filter loci dominated by near-read-end signal. |
| `CLUSTER_EVENTS` | Filter loci associated with cluster-like event patterns. |
| `LOW_MAPQ` | Filter low mapping-quality loci. |
| `LOW_BASEQ` | Filter low base-quality loci. |

## Tuning notes

- Start with defaults for production runs.
- Disable only specific groups when investigating why candidates were filtered.
- Document any filtration-group changes to preserve run comparability.

## Outputs

- `combine_feature`: `output_dir/all_feature.txt`
- parquet mirror: `output_dir/all_feature.parquet`

## Tuning notes

- Multi-index keys are `#chrom`, `pos`, `ref`, `alt`.
- Filtration summary is consolidated into the `Filtration` column (`PASS` or semicolon-separated tags).

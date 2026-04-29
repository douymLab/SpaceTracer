# `mpileup`

## Purpose

Builds pileup evidence from filtered BAM, filters candidate loci, and creates chunk metadata for later parallel processing.

## Upstream

- `bam_processing`

## Required inputs

- `in_filter_bam` (BAM after `bam_processing`)
- `genome_fasta`
- optional `regions_file`

### Input interpretation

| Input key | Source step/config | Required | Interpretation |
| --- | --- | --- | --- |
| `in_filter_bam` | `bam_processing` output | Yes | Filtered BAM used as pileup evidence source. Missing file will stop this step. |
| `genome_fasta` | `resource_details.genome_fasta` | Yes | Reference FASTA required by pileup and downstream filtering logic. |
| `regions_file` | top-level `regions_file` | No | If provided, limits analysis to target regions (useful for focused reruns/debugging). |

## Parameters (`steps.mpileup`)

| Parameter | Type | Typical/default | Interpretation |
| --- | --- | --- | --- |
| `min_depth` | int | `30` (template) | Minimum depth required for candidate retention after pileup filtering. Higher value reduces low-support candidates. |
| `max_depth` | int | `200000` (template) | Upper depth cap to avoid extreme-depth artifacts and unstable loci. |
| `min_mapq` | int | `0` | Minimum read mapping quality in pileup collection/filtering. |
| `min_baseq` | int | `0` | Minimum base quality used for pileup evidence. |
| `exclude_flag` | int | `0` | SAM flag mask to exclude reads. Use cautiously because aggressive masks may remove true evidence. |
| `enable_split` | bool | `true` | If enabled, large filtered pileup outputs are split into chunks for downstream parallel steps. |
| `split_threshold` | int | `10` (template) | Line-count threshold to trigger splitting. |
| `chrom_chunk_size` | int | `10` (template) | Chunk-size control for autosome splitting. |
| `chrM_chunk_size` | int | `10` (template) | Chunk-size control for mitochondrial chromosome splitting. |

## Tuning notes

- Increase `min_depth` when false positives are driven by sparse coverage.
- Keep `max_depth` high enough to retain informative deep loci, but not so high that PCR hotspots dominate.
- Use `regions_file` for targeted reruns/debugging to reduce runtime.
- Use `enable_split` for large genomes/runs to improve downstream parallel throughput.

## Outputs

- `mpileup_file`: `output_dir/mpileup/raw_mpileup.txt`
- `filter_mpileup_file`: `output_dir/mpileup/filter_mpileup.txt`
- `manifest_path`: `output_dir/mpileup/umi_combine_chunk_manifest.tsv`

## Tuning notes

- This step can split large genome tasks and emit chunk manifests for downstream processing.
- Chunk metadata (`umi_combine_chunk_manifest.tsv`) is consumed directly by `umi_combine`.

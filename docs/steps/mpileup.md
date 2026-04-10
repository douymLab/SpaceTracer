# `mpileup`

## Purpose

Builds pileup evidence from filtered BAM, filters candidate loci, and creates chunk metadata for later parallel processing.

## Upstream dependencies

- `bam_processing`

## Required inputs

- `in_filter_bam` (BAM after `bam_processing`)
- `genome_fasta`
- optional `regions_file`

## Key parameters

From `steps.mpileup`:

- `min_depth`
- `max_depth`
- `min_mapq`
- `min_baseq`
- `exclude_flag`
- `enable_split`
- `split_threshold`
- `chrom_chunk_size`
- `chrM_chunk_size`

## Outputs

- `mpileup_file`: `output_dir/mpileup/raw_mpileup.txt`
- `filter_mpileup_file`: `output_dir/mpileup/filter_mpileup.txt`
- `db_path`: `output_dir/split_chunk.db`

## Notes

- This step can split large genome tasks and create chunk-level database indexes.
- Chunk metadata (`split_chunk.db`) is consumed directly by `umi_combine`.

# Configuration Guide

This page explains the main fields in the SpaceTracer config YAML and provides a minimal template you can copy.

## Minimal runnable template

Save as `config.yaml` and replace all placeholder paths.

```yaml
genome: "hg38"
sequence_type: "visium"

input_details:
  bam_file: "/absolute/path/to/possorted_genome_bam.bam"
  tissue_position: "/absolute/path/to/tissue_positions_list.csv"
  barcode_key: "CB"

resource_details:
  genome_fasta: "/absolute/path/to/genome.fa"
  gnomad_path: "/absolute/path/to/gnomad"
  mappability_path: "/absolute/path/to/mappability"
  gene_bed: "/absolute/path/to/gene_region.bed"
  dbsnp_vcf_file: "/absolute/path/to/dbSNP.vcf"
  imprinted_bed: "/absolute/path/to/imprinted_gene_region.bed"
  editing_bed: "/absolute/path/to/editing.bed"
  PON_file: "/absolute/path/to/PON.txt"
  reference_error_profile: "/absolute/path/to/reference_error_profile.txt"

run:
  threads: 8
  memory: "32G"
  chunk_size: 1000000
  keep_intermediates: false
  skip_validation: false

output_dir: "/absolute/path/to/output"

steps:
  cluster:
    cluster_file: null
  cell_number: 0
```

## Key sections

### `input_details`

- `bam_file`: aligned BAM file
- `tissue_position`: required for Visium mode
- `barcode_key`: BAM tag used as barcode key (usually `CB` for Visium)

### `resource_details`

All resource paths are validated at startup. Missing files will stop the run immediately.

### `run`

- `threads`: CPU threads used by the pipeline
- `memory`: memory string for runtime tools
- `chunk_size`: genomic chunk size used in splitting/parallelization
- `keep_intermediates`: keep or clean temporary files
- `skip_validation`: disable output validation checks (use carefully)

### `steps`

- `cluster.cluster_file`: if provided and exists, SpaceTracer uses this cluster file directly
- `cell_number`: can be a fixed integer or a file path

If `cluster_file` is not provided and `sequence_type` is `visium`, SpaceTracer can compute clusters internally.

For full step-level parameter details, see:

- [Step Reference Overview](steps/overview.md)
- [Single-Step Debug Cookbook](steps/debug-cookbook.md)
- [cluster](steps/cluster.md)
- [genotyping](steps/genotyping.md)
- [RNA_feature](steps/rna-feature.md)
- [merge_feature](steps/merge-feature.md)

## Override behavior and path shortcuts

SpaceTracer supports directory-level shortcuts:

- `spaceranger_dir`: auto-resolves BAM and Visium tissue positions
- `resource_dir`: auto-resolves expected resource file names in one directory

You can still override any specific file in `input_details` or `resource_details`.

## Common config mistakes

- relative paths that point to the wrong working directory
- using `sequence_type: visium` without a valid tissue position file
- missing one of the required files under `resource_details`
- setting too many threads for your machine

## Suggested practice

- keep one validated baseline config per dataset
- copy and edit it for parameter experiments
- record changed parameters per run (for reproducibility)

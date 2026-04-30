# Configuration Guide

This page explains the main fields in the SpaceTracer config YAML and provides a minimal template you can copy.

For full parameter-by-parameter interpretation, see [Config Reference](config-reference.md).

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
- `memory`: memory limit string parsed as `<integer>G` (for example `32G`)
- `skip_validation`: disable output validation checks (use carefully)

### `runtime`

- `max_parallel`: controls parallel execution width for independent DAG layers.

### `steps`

- `cluster.cluster_file`: if provided and exists, SpaceTracer uses this cluster file directly
- `cell_number`: can be a fixed integer or a file path

If `cluster_file` is not provided and `sequence_type` is `visium`, SpaceTracer can compute clusters internally.

Model settings for mutation prediction:

- `model_dir`: directory containing trained models (for example `SpaceTracer_new_github/models`)
- `model_name`: model identifier (for example `spatial_free_model` or `spatial_feature_preserved_model`)

For full step-level parameter details, see:

- [Step Reference Overview](steps/overview.md)
- [Single-Step Debug Cookbook](steps/debug-cookbook.md)
- [cluster](steps/cluster.md)
- [bam_processing](steps/bam-processing.md)
- [mpileup](steps/mpileup.md)
- [umi_combine](steps/umi-combine.md)
- [cell_num](steps/cell-num.md)
- [prior](steps/prior.md)
- [genotyping](steps/genotyping.md)
- [spatial_feature](steps/spatial-feature.md)
- [mappability_feature](steps/mappability-feature.md)
- [read_feature](steps/read-feature.md)
- [RNA_feature](steps/rna-feature.md)
- [feature_filtration / merge_feature](steps/merge-feature.md)
- [phasing](steps/phasing.md)
- [mutation_prediction](steps/mutation-prediction.md)

## Override behavior and path shortcuts

SpaceTracer supports directory-level shortcuts:

- `spaceranger_dir`: auto-resolves BAM and Visium tissue positions
- `resource_dir`: auto-resolves expected resource file names in one directory

You can still override any specific file in `input_details` or `resource_details`.

## Active step names for CLI control

Use these names with `--start-from`, `--stop-at`, and `--only-steps`:

`cluster`, `bam_processing`, `mpileup`, `umi_combine`, `cell_num`, `prior`, `genotyping`, `spatial_feature`, `mappability_feature`, `read_feature`, `RNA_feature`, `phasing`, `merge_feature`, `mutation_prediction`

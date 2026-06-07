# Configuration Guide

This page explains how to choose and edit the SpaceTracer input config YAML used by `spacetracer run --config config.yaml`. For step-level default parameters, see [Default Step Config](default-step-config.md).

## Which template should I use?

SpaceTracer v2.1.0 keeps the main user config small and moves most step-specific parameters into a separate default file.

| Template | Best for | What to edit |
| --- | --- | --- |
| `SpaceTracer/config/config.example.yaml` | Standard Visium runs with a Space Ranger `outs` directory and a prepared resource directory. | Sample name, `spaceranger_dir`, `output_dir`, `resource_dir`, `model_name`, and compute resources. |
| `SpaceTracer/config/config.advanced.example.yaml` | Visium-HD, manual BAM/tissue-position/resource paths, external cluster files, or custom model locations. | Optional top-level shortcuts such as `cluster_method`, `cell_number`, `cell_info`, `model_dir`, `input_details`, and `resource_details`. |
| `SpaceTracer/config/default_step_config.yaml` | Advanced step tuning. | Default step parameter values. This is not the run input config; see [Default Step Config](default-step-config.md). |

## Minimal Visium config

Copy `SpaceTracer/config/config.example.yaml` to your project directory as `config.yaml` and replace the placeholder paths:

```yaml
genome: hg38
sequence_type: visium
spaceranger_dir: /path/to/sample/outs
output_dir: /path/to/output
sample: Sample
resource_dir: /path/to/resource_dir
model_name: spatial_free_model

run:
  threads: 60
  memory: "400G"
```

The two shipped model names are:

- `spatial_free_model`
- `spatial_preserved_model`

## Advanced config

Use `SpaceTracer/config/config.advanced.example.yaml` when the shortcut directories are not enough:

```yaml
genome: hg38
sequence_type: visium
spaceranger_dir: /path/to/sample/outs
output_dir: /path/to/output
sample: BCC-1
resource_dir: /path/to/resource_dir

# Optional advanced settings
model_name: spatial_preserved_model
model_dir: /path/to/SpaceTracer/models
cluster_file: /path/to/cluster_info/result/sample_cluster.txt
cluster_method: SpaGCN # choices: default, SpaGCN,
cell_number: /path/to/cellNum.txt # both int and path are supported
cell_info: /path/to/cell_info.txt # colum1 is barcode, colum2 is cell

run:
  threads: 60
  memory: "400G"
  max_parallel: 100
  keep_intermediates: true
  skip_validation: false

logging:
  level: "INFO"
  log_file: "pipeline.log"

# If spaceranger_dir is not available, use input_details instead
input_details:
  bam_file: /path/to/possorted_genome_bam.bam
  tissue_position: /path/to/tissue_positions_list.csv
  barcode_key: CB
  barcode_mapping: /path/to/barcode_mappings.parquet

# If resource_dir is not available, use resource_details instead
resource_details:
  genome_fasta: /path/to/genome.fa
  gnomad_path: /path/to/gnomad_af
  mappability_path: /path/to/k24.umap.bedgraph
  gene_bed: /path/to/gene_region.bed
  dbsnp_vcf_file: /path/to/Homo_sapiens_assembly38.dbsnp138.vcf.gz
  imprinted_bed: /path/to/imprinted_genes.region.bed
  editing_bed: /path/to/known_editing.bed
  PON_file: /path/to/PON.txt
  reference_error_profile: /path/to/Artifacts.Sigprofile.txt
```

If `spaceranger_dir` or `resource_dir` is not available, fill in the corresponding `input_details` or `resource_details` blocks from the advanced template.

## Config Fields

### Top-level fields

- `sample`: sample label used in output filenames.
- `genome`: genome build label, such as `hg38` or `mm10`.
- `sequence_type`: spatial platform mode, such as `visium` or `visium-HD`.
- `spaceranger_dir`: Space Ranger `outs` directory. SpaceTracer can resolve BAM, tissue positions, and related Visium files from this directory.
- `resource_dir`: directory containing the prepared reference resources from [Resources](resources.md).
- `output_dir`: root directory for all pipeline outputs.
- `model_name`: mutation prediction model name, usually `spatial_free_model` or `spatial_preserved_model`.
- `model_dir`: directory containing the trained model folders. Use this when the models are not in the default location.
- `cluster_file`: existing cluster/domain assignment file. If this is provided, SpaceTracer reuses it instead of generating a new cluster file.
- `cluster_method`: clustering strategy used when SpaceTracer needs to generate clusters. Supported values are `default`, `SpaGCN`, and `GraphST`.
- `cell_number`: fixed cell number or path to a cell-number file.
- `cell_info`: optional barcode-to-cell annotation file used by read-level feature extraction. The first column is barcode and the second column is cell.

### `run`

- `threads`: CPU threads used by the pipeline.
- `memory`: memory limit string parsed as `<integer>G`, for example `400G`.
- `max_parallel`: maximum parallel execution width for independent DAG layers.
- `keep_intermediates`: keep intermediate files when supported by a step.
- `skip_validation`: disable validation checks when needed for debugging.

### `logging`

- `level`: logging level label.
- `log_file`: optional pipeline log file name.

### `input_details`

Use this section when you do not want SpaceTracer to infer paths from `spaceranger_dir`.

- `bam_file`: aligned BAM file.
- `tissue_position`: Visium tissue-position table. SpaceTracer supports both `tissue_positions_list.csv` and `tissue_positions.csv`.
- `barcode_key`: BAM tag used as barcode key, usually `CB`.
- `barcode_mapping`: Visium-HD barcode mapping parquet file when required.

### `resource_details`

Use this section when you do not want SpaceTracer to infer paths from `resource_dir`. All configured resource paths are validated at startup.

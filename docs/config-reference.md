# Config Reference

This page is the full reference companion to [Configuration](configuration.md).
Use `configuration.md` for a minimal runnable template, and use this page for complete parameter interpretation.

## Complete template (verbatim from `config/default_config.yaml`)

```yaml
sample_name:
genome: hg38
sequence_type: visium
spaceranger_dir:
resource_dir: SpaceTracer/resources/hg38
model_dir: SpaceTracer/models/
model_name: "spatial_free_model"
bin_size:
regions_file:
output_dir:

run:
  threads: 16
  memory: "96G"
  skip_validation: true

input_details:
  bam_file:
  tissue_position:
  barcode_key: CB

resource_details:
  genome_fasta:
  gnomad_path:
  mappability_path:
  gene_bed:
  dbsnp_vcf_file:
  imprinted_bed:
  editing_bed:
  PON_file:
  reference_error_profile:

steps:
  cluster:
    cluster_file:
    ncluster: 8
    plot: true
    method: "SpaGCN"
    init_method: "louvain"
    data_type: "Visium"
    h5_file_name: "filtered_feature_bc_matrix.h5"
    histology: true
    spot_area: 49
    weight_histology: 1
    distance_threshold: 2
    min_samples: 1
    num_threshold: 30
    percentage: 0.5
    seed: 100
    tol: 5e-3
    lr: 0.05
    max_epochs: 200
    graphst_tool: "louvain"
    radius: 6
    refinement: true

  bam_processing:
    nm_threshold: 5
    mapq_threshold: 255

  mpileup:
    min_depth: 30
    max_depth: 200000
    min_mapq: 0
    min_baseq: 0
    exclude_flag: 0
    enable_split: true
    split_threshold: 10
    chrom_chunk_size: 10
    chrM_chunk_size: 10

  cell_number:

  UMI_combine:
    filter_duplicates: true
    filter_secondary: true
    filter_qcfail: true
    filter_supplementary: true
    min_read_quality: 20

  genotyping:
    alpha: 0.05
    epsQ: 20
    epsAF: 0.003
    mu: 1e-5
    thr_dp: 1000
    pop_vaf: 1e-5
    filter_oneallele: true

  spatial_feature:
    alpha: 0.05
    thr_r2: 0.3
    thr_prob: 0.9
    thr_likelihood: 0.9
    thr_vaf: 0
    plot_supp: false
    fig_size: 5
    method: LDA
    num_directions: 8

  read_feature:
    cell_info: None
    downsample: true
    downsample_target_depth: 2000
    max_region_size: 20000
    max_variants_per_region: 100
    seed: 42

  RNA_feature:
    min_count_for_germline: 50
    min_prior_for_germline: 0.0001
    default_range_of_gene: 150
    p_threshold: 0.05
    previous_base: 5

  phasing:
    minprior: 0.01
    min_dp: 20
    min_total_dp: 50
    alpha: 0.05
    phasing_pad: 1000
    merge_gap: 200
    max_target: 200000
    seed: 42

  feature_filtration:
    ASE: true
    hFDR: true
    imprinted: true
    homopolymer: true
    PON: true
    RNA_editing: true
    ABNORMAL_MISMATCHES: true
    LOW_READ_DIVERSITY: true
    HIGH_MULTIPLE_MAPPIN: true
    WIDE_DISTRIBUTION: true
    NEAR_READ_END: true
    CLUSTER_EVENTS: true
    LOW_MAPQ: true
    LOW_BASEQ: true

  mutation_prediction:
    random_seed: 42
    plot: true
```

## Top-level fields

### `sample_name`

- **Type**: string
- **Purpose**: Sample label used in downstream naming and outputs.

### `genome`

- **Type**: string
- **Purpose**: Genome build label used across steps and resources.
- **Example**: `"hg38"`

### `sequence_type`

- **Type**: string
- **Purpose**: Selects input mode behavior (for example Visium-specific handling).
- **Example**: `"visium"`

### `spaceranger_dir`

- **Type**: path string
- **Purpose**: Shortcut input root (`<prefix>/outs`) for auto-resolving BAM and Visium tissue-position files.

### `resource_dir`

- **Type**: path string
- **Purpose**: Shortcut directory for auto-resolving common resource files.

### `model_dir`, `model_name`

- **Type**: string/path
- **Purpose**: Model location and model identifier used by mutation prediction.

### `bin_size`

- **Type**: integer or null
- **Purpose**: Bin-size setting used for non-Visium workflows (for example stereo-seq).

### `regions_file`

- **Type**: path string or null
- **Purpose**: Restrict analysis to a target-region file when provided.

### `output_dir`

- **Type**: path string
- **Purpose**: Root output directory for all step outputs and checkpoints.

## `input_details`

- `bam_file`: aligned BAM input path.
- `tissue_position`: Visium tissue-position table path.
- `barcode_key`: BAM tag used as barcode key (commonly `CB`).

## `resource_details`

- `genome_fasta`: reference FASTA.
- `gnomad_path`: population-frequency resource path.
- `mappability_path`: mappability resource path.
- `gene_bed`: gene annotation BED.
- `dbsnp_vcf_file`: dbSNP VCF.
- `imprinted_bed`: imprinted-region BED.
- `editing_bed`: RNA-editing BED/resource.
- `PON_file`: panel-of-normals file.
- `reference_error_profile`: reference error profile file.

## `run`

- `threads`: total CPU threads for execution.
- `memory`: memory limit string (`<integer>G` format).
- `skip_validation`: disables output validation checks when `true`.

## `steps`

This namespace contains step-specific parameters.

### `steps.cluster`

- `cluster_file`, `ncluster`, `plot`, `method`, `init_method`, `data_type`, `h5_file_name`, `histology`, `spot_area`, `weight_histology`, `distance_threshold`, `min_samples`, `num_threshold`, `percentage`, `seed`, `tol`, `lr`, `max_epochs`, `graphst_tool`, `radius`, `refinement`.

### `steps.bam_processing`

- `nm_threshold`, `mapq_threshold`.

### `steps.cell_number` / top-level `steps.cell_number`

- fixed integer or file-backed setting depending on workflow mode.

### `steps.UMI_combine`

- `filter_duplicates`, `filter_secondary`, `filter_qcfail`, `filter_supplementary`, `min_read_quality`.

### `steps.genotyping`

- `alpha`, `epsQ`, `epsAF`, `mu`, `thr_dp`, `pop_vaf`, `filter_oneallele`.

### `steps.mpileup`

Common keys:

- `min_depth`
- `max_depth`
- `min_mapq`
- `min_baseq`
- `exclude_flag`
- `enable_split`
- `split_threshold`
- `chrom_chunk_size`
- `chrM_chunk_size`

See [mpileup step](steps/mpileup.md) for details.

### `steps.spatial_feature`

- `alpha`, `thr_r2`, `thr_prob`, `thr_likelihood`, `thr_vaf`, `plot_supp`, `fig_size`, `method`, `num_directions`.

### `steps.read_feature`

- `cell_info`, `downsample`, `downsample_target_depth`, `max_region_size`, `max_variants_per_region`, `seed`.

### `steps.RNA_feature`

- `min_count_for_germline`, `min_prior_for_germline`, `default_range_of_gene`, `p_threshold`, `previous_base`.

### `steps.phasing`

- `minprior`, `min_dp`, `min_total_dp`, `alpha`, `phasing_pad`, `merge_gap`, `max_target`, `seed`.

### `steps.feature_filtration`

- `ASE`, `hFDR`, `imprinted`, `homopolymer`, `PON`, `RNA_editing`, `ABNORMAL_MISMATCHES`, `LOW_READ_DIVERSITY`, `HIGH_MULTIPLE_MAPPIN`, `WIDE_DISTRIBUTION`, `NEAR_READ_END`, `CLUSTER_EVENTS`, `LOW_MAPQ`, `LOW_BASEQ`.

### `steps.merge_feature`

- behavior is tied to merged feature generation and filtration tags; see [merge_feature step](steps/merge-feature.md).

### `steps.mutation_prediction`

- `random_seed`, `plot` (plus model settings from top-level `model_dir`/`model_name`).

See [mutation_prediction step](steps/mutation-prediction.md).

## CLI step names (for `--start-from` / `--stop-at` / `--only-steps`)

`cluster`, `bam_processing`, `mpileup`, `umi_combine`, `cell_num`, `prior`, `genotyping`, `spatial_feature`, `mappability_feature`, `read_feature`, `RNA_feature`, `phasing`, `merge_feature`, `mutation_prediction`

## Important quick guide

Use this as a fast checklist for both parameter tuning and step-input handoff. For full field lists, see [Configuration](configuration.md) and the `steps.*` sections above.

- Runtime and input wiring: `run.*`, `sequence_type`, `spaceranger_dir`, `input_details.*`, `resource_details.*`
- Candidate detection and chunking: [`steps.mpileup`](#stepsmpileup) -> [mpileup step](steps/mpileup.md)
- Genotype confidence and inputs: [`steps.genotyping`](#stepsgenotyping) -> [genotyping step](steps/genotyping.md)
- Phasing behavior and inputs: [`steps.phasing`](#stepsphasing) -> [phasing step](steps/phasing.md)
- Filtration and feature merge handoff: [`steps.feature_filtration`](#stepsfeature_filtration), [`steps.merge_feature`](#stepsmerge_feature) -> [merge_feature step](steps/merge-feature.md)
- Final prediction and model artifacts: [`steps.mutation_prediction`](#stepsmutation_prediction), `model_dir`, `model_name` -> [mutation_prediction step](steps/mutation-prediction.md)

!!! note
    Keep one baseline config per dataset and change only a few important parameters per experiment to preserve comparability.

## Recommended usage pattern

1. Start with the minimal template in [Configuration](configuration.md).
2. Add only the step parameters you need to override.
3. Keep one validated baseline config per dataset.
4. Track parameter changes per run for reproducibility.

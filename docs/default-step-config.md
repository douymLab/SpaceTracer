# Default Step Config

`SpaceTracer/config/default_step_config.yaml` stores the default parameter values used by SpaceTracer steps.

This file is for advanced step tuning. It does **not** replace `config.yaml`: `config.yaml` is still the input file passed to `spacetracer run --config config.yaml`. In the current version, if advanced users want to modify step parameter values, they should edit `SpaceTracer/config/default_step_config.yaml` directly.

For workflow order and detailed step behavior, see:

- [Step-by-step guide](steps/overview.md)
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
- [phasing](steps/phasing.md)
- [feature_filtration / merge_feature](steps/merge-feature.md)
- [mutation_prediction](steps/mutation-prediction.md)
- [phylogeny](steps/phylogeny.md)

## Complete default file

```yaml
steps:
  cluster:
    ncluster: 8
    plot: true
    method: "SpaGCN"
    init_method: "louvain"
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
    split_threshold: 100000
    chrom_chunk_size: 20000
    chrM_chunk_size: 5000
    max_cost: 200000

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
    LOW_VAF: True
    HIGH_VAF: True
    LOW_SPOT_NUM: True
    LOW_ALT_SPOT_NUM: True
    ASE: True
    hFDR: True
    imprinted: True
    homopolymer: True
    PON: True
    RNA_editing: True
    ABNORMAL_MISMATCHES: True
    LOW_READ_DIVERSITY: True
    HIGH_MULTIPLE_MAPPIN: True
    WIDE_DISTRIBUTION: True
    NEAR_READ_END: True
    CLUSTER_EVENTS: True
    LOW_MAPQ: True
    LOW_BASEQ: True
    INDEL_PROPORTION: True
    ALT_ALLELE_COUNT: True
    POPULATION_AF: True
    CONTIG: True
    MITOCHONDRIA: True
    CLUSTERED_NOISE(RNA_editing): True
    CLUSTERED_NOISE: True

  mutation_prediction:
    random_seed: 42
    plot: true
```

## `cluster`

Controls automatic spot/domain clustering when no `cluster_file` is provided in `config.yaml`.

- `ncluster`: target number of clusters.
- `plot`: whether to generate clustering plots.
- `method`: clustering backend. Supported values include `default`, `SpaGCN`, and `GraphST`.
- `init_method`: initialization method for clustering, such as `louvain`.
- `h5_file_name`: Space Ranger feature matrix file name.
- `histology`: whether to use histology image information.
- `spot_area`: expected spot area used by spatial clustering.
- `weight_histology`: weight for histology information.
- `distance_threshold`: spatial distance threshold for domain refinement.
- `min_samples`: minimum samples used in neighborhood/domain filtering.
- `num_threshold`: minimum spot count threshold for domain filtering.
- `percentage`: spatial/histology weighting parameter used by the clustering backend.
- `seed`: random seed for reproducible clustering.
- `tol`: optimization tolerance.
- `lr`: learning rate.
- `max_epochs`: maximum optimization epochs.
- `graphst_tool`: community detection tool used by GraphST.
- `radius`: neighborhood radius used by GraphST.
- `refinement`: whether to refine cluster assignments after initial clustering.

## `bam_processing`

Controls BAM-level filtering before pileup and read-feature extraction.

- `nm_threshold`: maximum allowed edit distance/mismatch tag value.
- `mapq_threshold`: minimum mapping quality threshold.

## `mpileup`

Controls candidate-site pileup generation and chunking.

- `min_depth`: minimum depth required for a site.
- `max_depth`: maximum depth allowed for a site.
- `min_mapq`: minimum read mapping quality used by pileup.
- `min_baseq`: minimum base quality used by pileup.
- `exclude_flag`: SAM flag mask for reads excluded from pileup.
- `enable_split`: whether to split large pileup work into chunks.
- `split_threshold`: threshold for enabling chunk splitting.
- `chrom_chunk_size`: chunk size for standard chromosomes.
- `chrM_chunk_size`: chunk size for mitochondrial chromosome processing.
- `max_cost`: maximum chunk cost used to control workload size.

## `UMI_combine`

Controls read and UMI filtering before combined UMI evidence is produced.

- `filter_duplicates`: remove duplicate reads.
- `filter_secondary`: remove secondary alignments.
- `filter_qcfail`: remove reads marked as QC-failed.
- `filter_supplementary`: remove supplementary alignments.
- `min_read_quality`: minimum read quality for UMI evidence.

## `genotyping`

Controls genotype inference thresholds.

- `alpha`: statistical significance threshold.
- `epsQ`: sequencing error quality parameter.
- `epsAF`: allele-frequency error parameter.
- `mu`: prior mutation rate.
- `thr_dp`: depth threshold.
- `pop_vaf`: population allele-frequency threshold.
- `filter_oneallele`: whether to filter one-allele-only evidence.

## `spatial_feature`

Controls spatial feature extraction from genotype evidence.

- `alpha`: significance threshold.
- `thr_r2`: spatial correlation threshold.
- `thr_prob`: probability threshold.
- `thr_likelihood`: likelihood threshold.
- `thr_vaf`: variant allele frequency threshold.
- `plot_supp`: whether to generate supplementary plots.
- `fig_size`: output figure size.
- `method`: spatial feature method, such as `LDA`.
- `num_directions`: number of spatial directions used for directional features.

## `read_feature`

Controls read-level feature extraction.

- `downsample`: whether to downsample high-depth regions.
- `downsample_target_depth`: target depth after downsampling.
- `max_region_size`: maximum region size processed for read features.
- `max_variants_per_region`: maximum number of variants processed per region.
- `seed`: random seed for reproducible downsampling.

## `RNA_feature`

Controls RNA-level feature extraction.

- `min_count_for_germline`: minimum count used when evaluating germline-like evidence.
- `min_prior_for_germline`: minimum prior used for germline filtering.
- `default_range_of_gene`: default gene-flanking range.
- `p_threshold`: p-value threshold.
- `previous_base`: number of previous bases considered for sequence-context features.

## `phasing`

Controls RNA-informed phasing and cluster-event refinement.

- `minprior`: minimum prior required for phasing.
- `min_dp`: minimum depth for a phased site.
- `min_total_dp`: minimum total depth.
- `alpha`: significance threshold.
- `phasing_pad`: padding around target sites for phasing.
- `merge_gap`: maximum gap for merging nearby phased regions.
- `max_target`: maximum number of target sites.
- `seed`: random seed.

## `feature_filtration`

Controls which filtration tags are applied during feature merge and filtration.

- `LOW_VAF`: filter low variant allele frequency calls.
- `HIGH_VAF`: filter unusually high variant allele frequency calls.
- `LOW_SPOT_NUM`: filter calls supported by too few spots.
- `LOW_ALT_SPOT_NUM`: filter calls supported by too few mutant spots.
- `ASE`: filter allele-specific expression artifacts.
- `hFDR`: filter calls failing hFDR-related criteria.
- `imprinted`: filter imprinted-region artifacts.
- `homopolymer`: filter homopolymer-context artifacts.
- `PON`: filter panel-of-normals artifacts.
- `RNA_editing`: filter known RNA-editing sites.
- `ABNORMAL_MISMATCHES`: filter abnormal mismatch-pattern artifacts.
- `LOW_READ_DIVERSITY`: filter calls with low read diversity.
- `HIGH_MULTIPLE_MAPPIN`: filter high multiple-mapping artifacts.
- `WIDE_DISTRIBUTION`: filter calls with overly broad spatial distribution.
- `NEAR_READ_END`: filter calls close to read ends.
- `CLUSTER_EVENTS`: filter cluster-level event artifacts.
- `LOW_MAPQ`: filter low mapping-quality evidence.
- `LOW_BASEQ`: filter low base-quality evidence.
- `INDEL_PROPORTION`: filter calls with high nearby indel proportion.
- `ALT_ALLELE_COUNT`: filter calls with insufficient alternate allele count.
- `POPULATION_AF`: filter common population variants.
- `CONTIG`: filter calls on configured contig chromosomes.
- `MITOCHONDRIA`: filter mitochondrial calls.
- `CLUSTERED_NOISE(RNA_editing)`: filter clustered A-to-G/T-to-C noise patterns consistent with RNA-editing-like artifacts.
- `CLUSTERED_NOISE`: filter other clustered technical-noise patterns.

## `mutation_prediction`

Controls final model-based mutation prediction behavior.

- `random_seed`: random seed for reproducible prediction-related operations.
- `plot`: whether to generate mutation prediction plots.

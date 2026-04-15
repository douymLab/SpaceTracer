# Step 3: Feature extraction

Script: `3_run_get_features.py`

**What it does:** extract spatial-, read-, and site-level features for each candidate site.

## How to run (template)

```bash
python 3_run_get_features.py \
  --fasta <FASTAFILE> \
  --raw_bam <RAWBAMFILE> \
  --filter_bam <FILTERBAMFILE> \
  --gender <GENDER> \
  --outdir <OUTPUTDIRECTORY> \
  --outprefix <OUTPREFIX> \
  --thread <NUM_THREADS> \
  --germline <GERMLINEFILE> \
  --ind_genotype <INDGENOFILE> \
  --spot_genotype <SPOTGENOFILE> \
  --barcodes <BARCODESFILE> \
  --species <SPECIES> \
  --readLen <READLENGTH> \
  --prior <PRIORFILE> \
  --h5ad <H5ADFILE>  \
  --spaceranger_result_dir <SPACERANGERDIRECTORY> \
  --ind_count_file <INDCOUNTFILE> \
  --mappbablity_file <MAPPABILITYFILE> \
  --gff3_file <GFF3FILE> \
  --vaf_cluster_file <CLUSTERVAFFILE> \
  --gtexGene <GTEXTFILE> \
  --artifact_signature <ERRORSIGFILE> \
  --reference_error_profile <LYSISSIGFILE>
```

## Parameters

| Parameter | Description |
| --- | --- |
| `--fasta` | Reference FASTA |
| `--raw_bam` | Unfiltered BAM |
| `--filter_bam` | Filtered BAM |
| `--gender` | `F`, `M`, `female`, or `male` |
| `--outdir` | Output directory |
| `--outprefix` | Optional (default `sample`) |
| `--thread` | Optional (default 2) |
| `--germline` | Germline file |
| `--ind_genotype` | Individual genotype |
| `--spot_genotype` | Spot genotype |
| `--barcodes` | Spatial positions |
| `--species` | Optional (default `human`) |
| `--readLen` | Optional (default 150) |
| `--prior` | Optional. Population AF/prior |
| `--h5ad` | Optional. AnnData `h5ad` |
| `--spaceranger_result_dir` | Optional. Space Ranger `outs` parent |
| `--ind_count_file` | Optional. Individual allele counts |
| `--mappbablity_file` | Optional. Mappability |
| `--gff3_file` | Optional. GFF3 annotation |
| `--vaf_cluster_file` | Optional. Per-cluster VAF |
| `--gtexGene` | Optional. GTEx expression table |
| `--artifact_signature` | Optional. Sample artifact signature |
| `--reference_error_profile` | Optional. Reference lysis error profile |

## Example outputs

- `demo/output/features_dir/demo.spatial_feature.txt`
- `demo/output/features_dir/demo.phase_beforeUMIcombination.txt`
- `demo/output/features_dir/demo.phase_afterUMIcombination.txt`
- `demo/output/features_dir/demo.features.txt`
- `demo/output/features_dir/demo.features.add_hFDR.txt`

## Demo command

```bash
python 3_run_get_features.py \
  --fasta genome.fa \
  --raw_bam demo_output/demo/bam_filter/IN.bam \
  --filter_bam demo_output/demo/bam_filter/IN_filter.bam \
  --gender male \
  --outdir demo_output/demo/features_dir \
  --outprefix demo \
  --thread 2 \
  --germline demo_output/demo/geno_files/demo.germ_genotype.out \
  --ind_genotype demo_output/demo/geno_files/demo.ind_genotype_filter.out \
  --spot_genotype demo_output/demo/geno_files/demo.spot_genotype.out \
  --readLen 120 \
  --prior demo_input/demo.prior.out \
  --h5ad demo_input/demo_results.h5ad  \
  --spaceranger_result_dir demo_input/Spaceranger_result/outs/ \
  --ind_count_file demo_output/demo/counts_files/demo.ind_filter.count.out \
  --mappbablity_file Resources/demo.k24.umap.bedgraph \
  --gff3_file Resources/demo.gencode.v44.annotation.exon.sort.gff3 \
  --vaf_cluster_file demo_output/demo/geno_files/demo.cluster_vaf.out \
  --gtexGene Resources/demo.gtexGene.txt \
  --artifact_signature demo_inputs/lysis_errors.Sigprofile.txt \
  --reference_error_profile Resources/lysis_Signatures_split3.Sigprofile \
  --barcodes demo_input/Spaceranger_result/outs/spatial/tissue_positions.csv \

```

This step may take about 42 seconds on the demo data.

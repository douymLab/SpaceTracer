# Configuration

SpaceTracer Snakemake workflow reads parameters from a text config file (for example `config_example.txt`).

## Full example (`config_example.txt`)

```yaml
sample: "demo"
threads: 1
readLen: 120
spacerangerResult: "demo_input/Spaceranger_result"  #absolute path
srcPath: ${SpaceTracer_path}
GENOME_FA: "Resources/genome.fa"
gff3_file: "Resources/wgEncodeGencodeExonSupportV44.sort.bed"
mappbablity_file: "Resources/k24.umap.bedgraph"
prior_file: "demo_input/demo.prior.out" # use your personal prior file
error_profile: "demo_input/lysis_errors.SigProfile.txt"
gnomad_list: "" # if prior_file is provided, gnomad_list is not required
savePATH: "output/"
SPECIES: "human"
gender: "male"
gtexGene: "Resources/gtexGene.txt"
reference_error_SigProfile: "Resources/split_3_lysis_errors.SigProfile.txt"
cluster: "demo_input/demo_cluster.txt" ## you can provide personal cluster result. If not, we would run cluster automaticly
h5ad_file: "demo_input/demo_results.h5ad" ## you can provide personal cluster result. If not, we would run cluster automaticly
PON_file: "Resources/PON_33samples.txt"
dbSNP: "Resources/Homo_sapiens_assembly38.dbsnp138.vcf"
RNA_editing_path: "Resources/RNA_editing.bed"
imprinted_genes_bed: "imprinted_genes.region.bed"
```

## Core parameters

| Parameter | Description |
| --- | --- |
| `sample` | Sample name used for output prefixes. |
| `savePATH` | Output root directory. |
| `threads` | Number of CPU threads. |
| `readLen` | Sequencing read length. |
| `spacerangerResult` | Directory containing Space Ranger outputs. |
| `srcPath` | Path to SpaceTracer repository. |
| `GENOME_FA` | Reference genome FASTA path. |
| `gff3_file` | Genome annotation GFF3 file. |
| `mappbablity_file` | Mappability BEDGraph file. |
| `prior_file` | Population allele frequency file (optional but recommended). |
| `genomad_file` | gnomAD annotation path. |
| `SPECIES` | Species label. |
| `gender` | `F`, `M`, `female`, or `male`. |
| `gtexGene` | GTEx gene expression reference. |
| `error_SigProfile` | Artifact/error signature profile. |
| `cluster` | Spot cluster file (optional). |
| `h5ad_file` | h5ad file with spatial metadata (optional). |

## Practical tips

- use absolute file paths to avoid path issues
- verify all required files exist before starting Snakemake
- start with low threads for first dry run, then scale up
- keep one validated config per dataset for reproducibility

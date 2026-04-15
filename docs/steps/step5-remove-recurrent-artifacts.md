# Step 5: Remove recurrent artifacts

Script: `5_remove_recurrent_mutations.sh`

**What it does:** filter predicted sites using public resources and Panel of Normals (PON).  
Resource bundle: [Figshare](https://figshare.com/s/c7836f53c4eafb556ee1)

## How to run (template)

```bash
bash 5_remove_recurrent_mutations.sh \
  --src_path <SRC_PATH> \
  --pred_sites <PRED_SITES> \
  --feature_path <FEATURE_PATH> \
  --outpath <OUTPUT_PATH> \
  --germ_path <GERM_PATH> \
  --genome_fa <GENOME_FA> \
  --rna_editing_path <RNA_EDITING_PATH> \
  --gtexGene <GTEX_GENE> \
  --dbSNP <DB_SNP> \
  --imprinted_genes <IMPRINTED_GENES> \
  --PON_file <PON_FILE>
```

## Parameters

| Parameter | Description |
| --- | --- |
| `--src_path` | Working/source path |
| `--pred_sites` | Predicted mutation sites |
| `--feature_path` | Feature file used in classification |
| `--outpath` | Output directory |
| `--germ_path` | Germline information file |
| `--genome_fa` | Reference FASTA |
| `--rna_editing_path` | RNA editing data |
| `--gtexGene` | GTEx gene expression file |
| `--dbSNP` | dbSNP VCF |
| `--imprinted_genes` | Imprinted regions BED |
| `--PON_file` | Panel of Normals file |

## Example output

- `demo/output/$model/pred.FINAL.txt`

## Demo command

```bash
bash 5_remove_recurrent_mutations.sh \
  --src_path ${SpaceTracer_path} \
  --pred_sites demo/output/predict/results/demo_total_pred_truesites.txt \
  --feature_path demo/output/features_dir/demo.features.add_hFDR.txt \
  --outpath demo/output/$model \
  --germ_path $germ_path \
  --genome_fa $genome_fa \
  --rna_editing_path $RNA_editing_path \
  --gtexGene $gtexGene \
  --dbSNP $dbSNP \
  --imprinted_genes $imprinted_genes \
  --PON_file $PON_file
```

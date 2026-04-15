# Step 1: Quality control

Script: `1_run_data_process.py`

**What it does:** for each candidate site, count alleles per spot/cluster/individual after quality, background-error, and allele-frequency checks.

## How to run (template)

```bash
python 1_run_data_process.py \
  --posfile <POSFILE> \
  --bam <BAMFILE> \
  --outdir <OUTPUTDIRECTORY> \
  --outprefix <OUTPREFIX> \
  --barcodes <BARCODESFILE> \
  --cellcluster <CELLCLUSTERFILE> \
  --thread <NUM_THREADS>
```

## Parameters

| Parameter | Description |
| --- | --- |
| `--posfile` | Position file |
| `--bam` | BAM file |
| `--outdir` | Output directory |
| `--outprefix` | Optional (default `sample`). Output prefix |
| `--barcodes` | Spots spatial location file |
| `--cellcluster` | Spots cluster file |
| `--cellpos` | Optional (required for Stereo-seq). Barcode/name map |
| `--platform` | Optional (`visium`, `stereo`, `ST`; default `visium`) |
| `--thread` | Optional (default 2) |
| `--epsQ` | Optional (default 20). Consensus read quality threshold |
| `--alpha` | Optional (default 0.05). Confidence level |
| `--epsAF` | Optional (default 0.003). Background AF test threshold |

## Example outputs

- `demo/output/count_files/demo.spot_count.out`
- `demo/output/count_files/demo.cluster_filter.count.out`
- `demo/output/count_files/demo.cluster.count.out`
- `demo/output/count_files/demo.ind_filter.count.out`

## Demo command

```bash
python 1_run_data_process.py  \
  --posfile demo_output/demo/mpileup.filter.result \
  --outprefix demo \
  --bam demo_output/demo/bam_filter/IN_filter.bam \
  --barcodes  demo_input/Spaceranger_result/outs/spatial/tissue_positions.csv \
  --outdir demo_output/demo/counts_files \
  --thread 2 \
  --cellcluster demo_input/demo_cluster.txt
```

This step may take about 2 seconds on the demo data.

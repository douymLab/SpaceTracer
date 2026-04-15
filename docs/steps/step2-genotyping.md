# Step 2: Genotyping

Script: `2_run_genotyper.py`

**What it does:** infer genotypes at spot, cluster, and individual levels.

## How to run (template)

```bash
python 2_run_genotyper.py \
  --spot_count <SPOTCOUNTFILE> \
  --cluster_count <CLUSTERCOUNTFILE> \
  --ind_count <INTCOUNTFILE> \
  --cluster <CLUSTERFILE> \
  --outdir <OUTPUTDIRECTORY> \
  --outprefix <OUTPREFIX> \
  --prior <PRIORFILE> \
  --cellnum_file <CELLNUMFILE>
```

## Parameters

| Parameter | Description |
| --- | --- |
| `--spot_count` | Spot allele count file |
| `--cluster_count` | Cluster-level allele counts |
| `--ind_count` | Individual-level allele counts |
| `--cluster` | Spot cluster file |
| `--outdir` | Output directory |
| `--outprefix` | Optional (default `sample`) |
| `--prior` | Optional. Population allele frequency file |
| `--cellnum_file` | Optional. Estimated cells per spot |
| `--cell_num` | Optional (default 20). Cells per spot if no `cellnum_file` |
| `--epsQ` | Optional (default 20) |
| `--alpha` | Optional (default 0.05) |
| `--epsAF` | Optional (default 0.003) |
| `--mu` | Optional (default 1e-7). Population mutation rate prior |
| `--max_dp` | Optional (default 1000). Downsample if depth exceeds |
| `--vaf` | Optional (default 1e-5). Avoid zero AF |
| `--min_dp` | Optional (default 30). Min depth |
| `--max_vaf` | Optional (default 0.3). Max VAF to keep |
| `--min_vaf` | Optional (default 1e-5). Min VAF |

## Example outputs

- `demo/output/geno_files/demo.germ_genotype.out`
- `demo/output/geno_files/demo.ind_genotype.out`
- `demo/output/geno_files/demo.ind_genotype_filter.out`
- `demo/output/geno_files/demo.cluster_vaf.out`
- `demo/output/geno_files/demo.spot_genotype.out`

## Demo command

```bash
python 2_run_genotyper.py \
  --spot_count demo_output/demo/counts_files/demo.spot.count.out \
  --cluster_count demo_output/demo/counts_files/demo.cluster.count.out \
  --ind_count demo_output/demo/counts_files/demo.ind_filter.count.out \
  --cluster demo_input/demo_cluster.txt \
  --outprefix demo \
  --prior demo_input/demo.prior.out \
  --cellnum_file demo_input/refined_umi_read_cellNum.txt  \
  --outdir demo_output/demo/geno_files
```

This step may take about 4 seconds on the demo data.

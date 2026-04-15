# Step 0: Data pre-processing

This optional stage prepares raw/public data before the core SpaceTracer steps.

If your inputs are already ready for Step 1, you can skip Step 0.

---

## Step 0.1 - Filter BAM

**What it does:** keep in-tissue reads, then filter by mapping quality and mismatch count.

**Inputs**

| Item | Description |
| --- | --- |
| Space Ranger outputs | e.g. `demo_input/Spaceranger_result/outs/possorted_genome_bam.bam`, `.../spatial/tissue_positions.csv` |

**Outputs**

| Path | Description |
| --- | --- |
| `demo_output/demo/bam_filter/demo_barcode.txt` | IN/OUT barcode list |
| `demo_output/demo/bam_filter/IN.bam` | In-tissue BAM |
| `demo_output/demo/bam_filter/IN_filter.bam` | Filtered BAM (indexed) |

**Demo command**

```bash
# create output directory
mkdir -p demo_output/demo/bam_filter/
# filter bam file (1. reads not in tissue; 2. reads with low Mapping Quality; 3. reads with high mismatches)
sed 's/,/\t/g' demo_input/Spaceranger_result/outs/spatial/tissue_positions.csv | sed '1d' | awk '{if ($2==0) print $1, "OUT"; else print $1, "IN"}' OFS="\t" >demo_output/demo/bam_filter/demo_barcode.txt;sinto filterbarcodes -b demo_input/Spaceranger_result/outs/possorted_genome_bam.bam -c demo_output/demo/bam_filter/demo_barcode.txt --barcodetag "CB" --outdir demo_output/demo/bam_filter -p 4;samtools index demo_output/demo/bam_filter/IN.bam; samtools view -e "[nM] <= 5" -q 255 -o demo_output/demo/bam_filter/IN_filter.bam demo_output/demo/bam_filter/IN.bam; samtools index demo_output/demo/bam_filter/IN_filter.bam
```

---

## Step 0.2 - Mpileup and candidate sites

**What it does:** run `samtools mpileup` + Java `PileupFilter`, then apply candidate-site filters.

**Inputs**

| Item | Description |
| --- | --- |
| `IN_filter.bam` | From Step 0.1 |
| `${genome.fa}` | Reference genome (absolute path recommended) |

**Outputs**

| Path | Description |
| --- | --- |
| `demo_output/demo/mpileup.result` | Raw mpileup table |
| `demo_output/demo/mpileup.filter.result` | Candidate sites |

**Demo command**

```bash
# mpile-up (please offer the absolute path of ${genome.fa} )
samtools mpileup demo_output/demo/bam_filter/IN_filter.bam -s --excl-flags 0 -B -Q 0 -q 0 -d 200000 -f ${genome.fa} | java -classpath others/java_mpileup/ PileupFilter --minbasequal=0 --minmapqual=0 --asciibase=33 --filtered=1| awk '$3!="N"' |cut -f 1-3,8-15 >demo_output/demo/mpileup.result

# filter mpile-up results to get initial candidate sites
awk '$6+$7+$9+$10 >= 5 && $4+$5+$6+$7+$9+$10>=30' demo_output/demo/mpileup.result |awk '($6+$7+$9+$10)/($4+$5+$6+$7+$9+$10)<=0.6 && ($6+$7+$9+$10)/($4+$5+$6+$7+$9+$10)>=0.001' - > demo_output/demo/mpileup.filter.result
```

---

## Step 0.3 - Population allele frequency (prior)

**What it does:** build prior information from gnomAD-style resources.

**Inputs**

| Item | Description |
| --- | --- |
| `--posfile` | Candidate positions from Step 0.2 |
| `${annovar_info.txt}` | Chromosome-wise gnomAD info list |

**Parameters**

| Parameter | Description |
| --- | --- |
| `--outprefix` | Output prefix |
| `--posfile` | Position file |
| `--annovar` | annovar/gnomAD info file |
| `--thread` | Threads |
| `--outdir` | Output directory |

**Outputs**

Files under `--outdir` with prefix `--outprefix` (used later as `--prior`).

**Demo command**

```bash
# the ${annovar_info.txt} was a file containing the gnomad files, which seperated by chromosome. And these files are accessable in https://figshare.com/s/c7836f53c4eafb556ee1, with 
python 0_other_func.py prior \
  --outprefix demo \
  --posfile demo_output/mpileup.filter.result  \
  --annovar ${annovar_info.txt} \
  --thread 2 \
  --outdir demo_output/demo
```

---

## Step 0.4 - Calculate cell number per spot

**What it does:** estimate cells per spot from UMI/read information.

**Inputs**

| Item | Description |
| --- | --- |
| `--bam` | Filtered BAM |
| `--cluster` | Spot cluster file |

**Parameters**

| Parameter | Description |
| --- | --- |
| `--bam` | Input BAM |
| `--run` | Run estimation |
| `--cluster` | Cluster assignments |
| `--outdir` | Output directory |

**Outputs**

Files under `--outdir` (used downstream as `cellnum_file`).

**Demo command**

```bash
python others/get_umiCount_cellNum.py \
  --bam demo_output/demo/bam_filter/IN_filter.bam \
  --run \
  --cluster demo_input/demo_cluster.txt \
  --outdir demo_output/demo
```

---

## Step 0.5 - Cluster spots

**What it does:** compute spot clusters (example method: SpaGCN).

**Parameters**

| Parameter | Description |
| --- | --- |
| `--indir` | Space Ranger `outs` directory |
| `--method` | e.g. `SpaGCN` |
| `--ncluster` | Number of clusters |
| `--sample` | Sample name |
| `--outdir` | Output directory |

**Demo command**

```bash
python 0_other_func.py cluster \
  --indir demo_input/Spaceranger_result/outs \
  --method SpaGCN \
  --ncluster 6 \
  --sample demo \
  --outdir demo_input/Spaceranger_result/outs
```

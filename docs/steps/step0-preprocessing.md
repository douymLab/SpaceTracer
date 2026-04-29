# Step 0: Data pre-processing

This stage prepares raw inputs before genotype inference.

If your project already has prepared candidate-locus and metadata inputs, this stage can be shortened.

---

## Step 0.1 - Cluster preparation (`cluster`)

**What it does:** loads or computes tissue-domain/cluster assignments and cell-number context.

**Detailed reference:** [cluster](cluster.md)

---

## Step 0.2 - BAM filtering (`bam_processing`)

**What it does:** retains in-tissue reads and generates filtered/indexed BAM files for pileup.

**Detailed reference:** [bam_processing](bam-processing.md)

---

## Step 0.3 - Candidate generation (`mpileup`)

**What it does:** runs pileup and candidate filtering, then prepares chunk metadata for downstream processing.

**Detailed reference:** [mpileup](mpileup.md)

## Typical outputs

- `output_dir/bam_processing/IN.bam`
- `output_dir/bam_processing/IN_filter.bam`
- `output_dir/mpileup/raw_mpileup.txt`
- `output_dir/mpileup/filter_mpileup.txt`
- `output_dir/split_chunk.db`

# Public Resources

This page lists reference resources commonly used with SpaceTracer.

## Start here: prebuilt reference packages (recommended)

**Recommended for first-time users** to avoid manual reference setup.

- Source: [Zenodo resource package (SpaceTracer v2.1.0)](https://zenodo.org/records/19896967)
- Main archive: `resources.tar.tar` (~11.5 GB)
- Demo archive: `demo.zip` (~4.4 MB)
- Includes:
  - `mm10_resources.tar.zst` (mouse, mm10 / GRCm38)
  - `hg38_resources.tar.zst` (human, hg38 / GRCh38)

Download from Zenodo with `wget` (append `?download=1` so the server serves the file directly):

```bash
# Full reference package (~11.5 GB) — pick a working directory first
wget -O resources.tar "https://zenodo.org/records/19896967/files/resources.tar.tar?download=1"

# Optional: small demo package from the same record
wget -O demo.zip "https://zenodo.org/records/19896967/files/demo.zip?download=1"
```

If you prefer `curl`, use:

```bash
curl -L -o resources.tar "https://zenodo.org/records/19896967/files/resources.tar.tar?download=1"
```

Extract references:

```bash
tar -xf resources.tar

# mm10
mkdir -p mm10 && zstd -d mm10_resources.tar.zst | tar -xf - -C mm10

# hg38
mkdir -p hg38 && zstd -d hg38_resources.tar.zst | tar -xf - -C hg38
```

Extract the demo package:

```bash
unzip demo.zip
```

The demo package contains `human_visium_demo_input.tar` and `human_visium_demo_output.tar`. See [Demo](demo.md) for extraction, config editing, and run instructions.

If you need custom references, use the resource links below to build your own.

## Human reference genome example

### GRCh37 / hg19

```bash
wget ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz
wget ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz.fai
```

### GRCh38 / hg38

```bash
wget ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa
wget ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa.fai
```

NCBI reference page: [NCBI human datasets](https://www.ncbi.nlm.nih.gov/datasets/taxonomy/9606/).

## Annotation and supporting resources

- gene annotation (Gencode): [GENCODE human](https://www.gencodegenes.org/human/#)
- mappability resources: [UCSC Genome Browser](https://genome.ucsc.edu/)
- GTEx expression: [GTEx Portal](https://gtexportal.org/home/) or [dbGaP website](https://www.ncbi.nlm.nih.gov/gap/)
- dbSNP bundle: [GATK resources](https://gatk.broadinstitute.org/)


# Public Resources

This page lists public datasets and reference files commonly used with SpaceTracer.

## Reference data

### Human reference genome

Human reference genome files can be downloaded from public repositories.

#### GRCh37 / hg19

??? note "Download commands"
    ```bash
    wget ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz
    wget ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz.fai
    ```

#### GRCh38 / hg38

??? note "Download commands"
    ```bash
    wget ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa
    wget ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa.fai
    ```

Reference genome files are also available from the [NCBI human datasets page](https://www.ncbi.nlm.nih.gov/datasets/taxonomy/9606/).

### Human genome annotation

Human genome annotation files can be downloaded directly from the [GENCODE website](https://www.gencodegenes.org/human/#).

### dbSNP

The dbSNP138 VCF file for hg38 (`Homo_sapiens_assembly38.dbsnp138.vcf`) can be obtained from the [Broad Institute GATK Resource Bundle](https://gatk.broadinstitute.org/).

## Supporting data

### Mappability score

Mappability resources can be obtained from the [UCSC Genome Browser](https://genome.ucsc.edu/) or downloaded using the commands below.

#### Umap score (k = 24, GRCh37 / hg19)

??? note "Download commands"
    ```bash
    wget https://bismap.hoffmanlab.org/raw/hg19.umap.tar.gz
    tar -zxvf hg19.umap.tar.gz
    ```

#### Umap score (k = 24, GRCh38 / hg38)

??? note "Download commands"
    ```bash
    wget https://bismap.hoffmanlab.org/raw/hg38.umap.tar.gz
    tar -zxvf hg38.umap.tar.gz
    ```

### GTEx gene expression data

GTEx gene expression data can be accessed through the [GTEx Portal](https://gtexportal.org/home/) or the [dbGaP website](https://www.ncbi.nlm.nih.gov/gap/).

## Resource collection

For convenience, we provide a [Figshare collection](https://figshare.com/s/c7836f53c4eafb556ee1) that includes many of the resources listed on this page.

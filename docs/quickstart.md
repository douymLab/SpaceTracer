# Quick Start

This page provides a minimal workflow for running SpaceTracer.

## Before you begin

Make sure you have:

* cloned the SpaceTracer repository
* prepared the required environment
* downloaded the required reference resources

If you have not completed these steps yet, please see the [Installation](installation.md) and [Resources](resources.md) pages.

## Step 1: Clone the repository

```bash
git clone https://github.com/douymLab/SpaceTracer.git
cd SpaceTracer
```

## Step 2: Prepare the environment

You can either use the Docker image or create the Conda environment.

### Option 1: Docker

```bash
docker pull xiayh17/spacetracer
docker run -it -v $(pwd):/mnt/workflow xiayh17/spacetracer bash
```

### Option 2: Conda

```bash
conda env create -f environment.yaml
conda activate SpaceTracer
```

## Step 3: Prepare input files

Before running SpaceTracer, make sure your input data and reference files are ready.

Typical inputs include:

* spatial transcriptomics data
* reference genome
* genome annotation
* other required resource files described in the [Resources](resources.md) page

## Step 4: Run SpaceTracer

Run SpaceTracer with your input files and output directory.

```bash
python SpaceTracer.py \
    --input <input_data> \
    --reference <reference_genome> \
    --annotation <annotation_file> \
    --output <output_directory>
```

Replace the placeholders above with your actual file paths and parameters.

## Step 5: Check the output

After the run is complete, the output directory will contain the SpaceTracer results for downstream analysis.

Typical outputs may include:

* detected somatic SNVs
* intermediate result files
* summary tables
* lineage-related analysis results

## Next steps

For more details about required inputs, parameters, and output formats, please refer to the full tutorial and documentation pages.

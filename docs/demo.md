# Demo

This page explains how to run the human Visium demo after downloading `demo.zip` from the SpaceTracer Zenodo record.

## 1) Download and extract the demo package

```bash
wget -O demo.zip "https://zenodo.org/records/19896967/files/demo.zip?download=1"
unzip demo.zip
```

The extracted package contains two tar archives:

- `human_visium_demo_input.tar`
- `human_visium_demo_output.tar`

Extract both archives:

```bash
tar -xf human_visium_demo_input.tar
tar -xf human_visium_demo_output.tar
```

After extraction, you should have two folders:

- `human_visium_demo`: input files and demo config for running SpaceTracer.
- `human_visium`: example output files for comparison.

## 2) Prepare reference resources

Download and extract the reference resources described in [Resources](resources.md). For this human demo, prepare the `hg38` resource directory and use that path as `resource_dir` in the demo config.

## 3) Edit the demo config

Open:

```text
human_visium_demo/config.yaml
```

The demo config is already set up for the packaged example, but several paths must match where you extracted the demo and resources. The `output_dir: human_visium` setting writes demo results to the `human_visium` folder unless you change it.

A complete edited config looks like this:

```yaml
genome: hg38
sequence_type: visium
spaceranger_dir:
output_dir: human_visium
sample: demo
resource_dir: /path/to/resource_dir

# Optional advanced settings
model_name: spatial_preserved_model
cluster_file: /path/to/human_visium_demo/cluster.txt
cluster_method: default
cell_number: /path/to/human_visium_demo/cell_num.txt

# If spaceranger_dir is not available, use input_details instead
input_details:
  bam_file: /path/to/human_visium_demo/demo_human.bam
  tissue_position: /path/to/human_visium_demo/tissue_positions_list.csv
  barcode_key:
  barcode_mapping:

run:
  threads: 60
  memory: "400G"
  max_parallel: 100
  skip_validation: false
```

You can use absolute paths, or relative paths if you run SpaceTracer from the extracted `demo` folder.

For example, if your current directory is `demo` and it contains `human_visium_demo`, the config can use:

```yaml
output_dir: human_visium
resource_dir: /path/to/hg38
cluster_file: human_visium_demo/cluster.txt
cell_number: human_visium_demo/cell_num.txt

input_details:
  bam_file: human_visium_demo/demo_human.bam
  tissue_position: human_visium_demo/tissue_positions_list.csv
```

## 4) Run the demo

From the directory where `human_visium_demo` is located, run:

```bash
spacetracer run --config human_visium_demo/config.yaml
```

PhyloSOLID currently requires at least 3 mutations to construct a phylogeny tree. The demo output contains only 1 mutation, so the `phylogeny` step is skipped automatically even when the full workflow is run.

To run only through mutation prediction and skip PhyloSOLID tree building:

```bash
spacetracer run --config human_visium_demo/config.yaml --stop-at mutation_prediction
```

## 5) Check outputs

The demo config writes outputs to `human_visium` unless you change `output_dir`.

Important outputs include:

- `human_visium/all_feature.txt`
- `human_visium/all_feature.parquet`
- `human_visium/mutation_prediction/results/<sample>_<model_name>_total_pred_truesites.vcf`
- `human_visium/mutation_prediction/results/<sample>_<model_name>_total_pred_truesites_PASS.vcf`
- `human_visium/mutation_prediction/results/<sample>_<model_name>_total_pred_truesites_PASS_mutation_list.txt`

You can compare your generated output folder with the extracted `human_visium` example output.

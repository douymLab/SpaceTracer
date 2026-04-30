# SpaceTracer

SpaceTracer is an open-source algorithm capable of accurately detecting mosaic SNVs, including both nuclear SNVs and mitochondria SNVs, directly from spatial transcriptomics data. 

![flowchart](./figures/flowchart.png)

## Release Notes

- 2026/04/30: Version 2.0.0 ([changelog](https://github.com/douymLab/SpaceTracer/releases/tag/v2.0.0))  
This release introduces one-command workflow execution, integrates lysis error calculation into the main pipeline, and significantly improves running speed.

- 2026/02/25: Version 1.1.0  
This release focuses on updating the genotype calculation and enhancing the features used in the random forest model. Additionally, we've added more filter steps to improve the accuracy of the results.
- 2025/03/31: Version 1.0.0  
This is the initial version of SpaceTracer.

## Documentation

**Full tutorial website:** [SpaceTracer Docs](https://douymLab.github.io/SpaceTracer/)  
Includes installation, resources, configuration, quick start, outputs, and step-by-step pages.

## Quick Start

SpaceTracer requires Python `3.9`, `samtools`, `bedtools`, and `vcftools`.

### 1) Install SpaceTracer

```bash
git clone https://github.com/douymLab/SpaceTracer.git
cd SpaceTracer
```

### 2) Prepare environment (recommended: Conda)

```bash
conda env create -f environment.yml
conda activate SpaceTracer
pip install .
spacetracer --help
```

If install failed, you can try to add the pythonpath, and run the command by SpaceTracer/cli/run.py --config config.yaml

### 3) Run the pipeline (most efficient path)

Run with your config:

```bash
spacetracer run --config config.yaml
```

Fallback command:

```bash
python -m SpaceTracer.cli.run --config config.yaml
```

### 4) Configure inputs/resources/model

Start from the docs template and update sample/resource paths.

SpaceTracer provides two pretrained mutation-prediction models under `models/`:

- `spatial_free_model`
- `spatial_feature_preserved_model`

For complete parameter descriptions and resource guidance, see:

- Docs home: [SpaceTracer Docs](https://douymLab.github.io/SpaceTracer/)
- Configuration: [Configuration Guide](https://douymLab.github.io/SpaceTracer/configuration/)
- Config reference: [Config Reference](https://douymLab.github.io/SpaceTracer/config-reference/)
- Resources: [Resources](https://douymLab.github.io/SpaceTracer/resources/)
- Quick start: [Quick Start](https://douymLab.github.io/SpaceTracer/quickstart/)

You can also download the packaged resource files (for mm10 and hg38) and the `demo_input` archive to quickly get started with SpaceTracer from [zendo].

### 5) Main outputs

Common downstream outputs:

- `output_dir/all_feature.txt`
- `output_dir/all_feature.parquet`
- `output_dir/mutation_prediction/Sample_total_pred_truesites.vcf`
- `output_dir/mutation_prediction/Sample_total_pred_truesites_PASS.vcf`

### Optional: advanced step-by-step mode

Advanced users can run/debug by step:

- [Step-by-step guide](https://douymLab.github.io/SpaceTracer/steps/overview/) (code-aligned with `SpaceTracer/steps`)
- [Single-Step Debug Cookbook](https://douymLab.github.io/SpaceTracer/steps/debug-cookbook/)

## Contact:

If you have any questions please contact us:  
Zhirui Yang: [yangzhirui@westlake.edu.cn](mailto:yangzhirui@westlake.edu.cn)  
Mengdie Yao: [yaomengdie@westlake.edu.cn](mailto:yaomengdie@westlake.edu.cn)  
Yanmei Dou: [douyanmei@westlake.edu.cn](mailto:douyanmei@westlake.edu.cn)

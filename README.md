# SpaceTracer

SpaceTracer is an open-source algorithm capable of accurately detecting mosaic SNVs, including both nuclear SNVs and mitochondria SNVs, directly from spatial transcriptomics data. 

> **Active development:** We are actively maintaining and extending SpaceTracer. New releases and documentation updates are posted in [Release Notes](#release-notes) and on [GitHub Releases](https://github.com/douymLab/SpaceTracer/releases).

![flowchart](./figures/flowchart.png)

## Release Notes

*SpaceTracer is actively maintained. Check the Release Notes below and GitHub Releases for the latest updates.*

- 2026/06/08: Version 2.1.1 ([changelog](https://github.com/douymLab/SpaceTracer/releases/tag/v2.1.1))
Fixed a PhyloSOLID input formatting issue so the integrated phylogeny step can run successfully.

- 2026/06/05: Version 2.1.0 ([changelog](https://github.com/douymLab/SpaceTracer/releases/tag/v2.1.0))
This release updates the packaged demo, improves Space Ranger 4.0.1 clustering compatibility, simplifies the config templates, and adds the integrated PhyloSOLID phylogeny step.

- 2026/04/30: Version 2.0.0 ([changelog](https://github.com/douymLab/SpaceTracer/releases/tag/v2.0.0))  
This release introduces one-command workflow execution, integrates lysis error calculation into the main pipeline, and significantly improves running speed.

- 2026/02/25: Version 1.1.0  
This release focuses on updating the genotype calculation and enhancing the features used in the random forest model. Additionally, we've added more filter steps to improve the accuracy of the results.
- 2025/03/31: Version 1.0.0  
This is the initial version of SpaceTracer.

## Documentation

**Full tutorial website:** [SpaceTracer Docs](https://douymLab.github.io/SpaceTracer/)  
Includes installation, resources, configuration, demo, quick start, outputs, and step-by-step pages.

## Quick Start

SpaceTracer requires Python `3.9`, `samtools`, `bedtools`, and `bcftools`.

### 1) Install SpaceTracer

```bash
git clone https://github.com/douymLab/SpaceTracer.git
cd SpaceTracer
```

### 2) Prepare environment (recommended: Conda)

```bash
conda env create -f environment.yaml
conda activate SpaceTracer
pip install .
spacetracer --help
```

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
- `spatial_preserved_model`

For complete parameter descriptions and resource guidance, see:

- Docs home: [SpaceTracer Docs](https://douymLab.github.io/SpaceTracer/)
- Configuration: [Configuration Guide](https://douymLab.github.io/SpaceTracer/configuration/)
- Default step config: [Default Step Config](https://douymLab.github.io/SpaceTracer/default-step-config/)
- Resources: [Resources](https://douymLab.github.io/SpaceTracer/resources/)
- Demo: [Demo](https://douymLab.github.io/SpaceTracer/demo/)
- Quick start: [Quick Start](https://douymLab.github.io/SpaceTracer/quickstart/)

You can also download packaged resources (mm10/hg38) and the human Visium demo from Zenodo:

```bash
wget -O resources.tar "https://zenodo.org/records/19896967/files/resources.tar.tar?download=1"
wget -O demo.zip "https://zenodo.org/records/19896967/files/demo.zip?download=1"
```

### 5) Main outputs

Common downstream outputs:

- `output_dir/all_feature.txt`
- `output_dir/all_feature.parquet`
- `output_dir/mutation_prediction/results/Sample_total_pred_truesites.vcf`
- `output_dir/mutation_prediction/results/Sample_total_pred_truesites_PASS.vcf`
- `output_dir/phylogeny/tree/mutation_integrator/phylo/final_cleaned_M_full_basedPivots.filtered_sites_inferred.tree_scphylo.pdf`

### Optional: advanced step-by-step mode

Advanced users can run/debug by step:

- [Step-by-step guide](https://douymLab.github.io/SpaceTracer/steps/overview/) (code-aligned with `SpaceTracer/steps`)
- [Single-Step Debug Cookbook](https://douymLab.github.io/SpaceTracer/steps/debug-cookbook/)

## Contact:

If you have any questions please contact us:  
- Zhirui Yang: [yangzhirui@westlake.edu.cn](mailto:yangzhirui@westlake.edu.cn)
- Mengdie Yao: [yaomengdie@westlake.edu.cn](mailto:yaomengdie@westlake.edu.cn)
- Qing Yang: [yangqing@westlake.edu.cn](mailto:yangqing@westlake.edu.cn) (PhyloSOLID)
- Yanmei Dou: [douyanmei@westlake.edu.cn](mailto:douyanmei@westlake.edu.cn) (Corresponding author)

# SpaceTracer

SpaceTracer is an open-source algorithm capable of accurately detecting mosaic SNVs, including both nuclear SNVs and mitochondria SNVs, directly from spatial transcriptomics data. 

flowchart

## Release Notes

- 2025/03/31: Version 1.0.0  
This is the initial version of SpaceTracer.
- 2026/02/25: Version 1.1.0  
This release focuses on updating the genotype calculation and enhancing the features used in the random forest model. Additionally, we've added more filter steps to improve the accuracy of the results.

We expect to release the **next version** of SpaceTracer soon. This release will 
feature substantial performance improvements, the **most important update** in this release is significantly improved running speed.

**Key Improvements for Upcoming Release**

- Integrated Lysis Error Calculation   
Lysis error calculation for single samples will be incorporated directly into the full algorithm pipeline, providing more comprehensive and streamlined analysis.
- One-Command Execution Mode   
A simplified execution option will allow users to run SpaceTracer with a single command, removing the Snakemake dependency and making the workflow more accessible.
- Improved Running Speed   
Major performance optimizations will significantly accelerate processing compared with the previous version.

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
```

Alternative (container):  
Docker image is available at [xiayh17/spacetracer](https://hub.docker.com/r/xiayh17/spacetracer).

```bash
docker pull xiayh17/spacetracer
docker run -it -v $(pwd):/mnt/workflow xiayh17/spacetracer bash
```

### 3) Run the pipeline (most efficient path)

SpaceTracer is designed for one-command execution through Snakemake:

```bash
# conda-based run
snakemake --configfile config_example.txt -s run_snakemake_for_conda

# docker-based run
snakemake --configfile config_example.txt -s run_snakemake_for_docker
```

### 4) Configure inputs/resources

Start from `config_example.txt` and update paths for your sample and resources.

For complete parameter descriptions and resource guidance, see:

- Docs home: [SpaceTracer Docs](https://douymLab.github.io/SpaceTracer/)
- Configuration: [Configuration Guide](https://douymLab.github.io/SpaceTracer/configuration/)
- Resources: [Resource Guide](https://douymLab.github.io/SpaceTracer/resources/)

### 5) Final output

Final high-confidence somatic mutation list:

`$savePATH/$sample/$model/pred.FINAL.txt`

### Optional: advanced step-by-step mode

Advanced users can run all the steps manually for debugging/custom workflows; see **Step Reference** in the tutorial website.

## Contact:

If you have any questions please contact us:  
Zhirui Yang: [yangzhirui@westlake.edu.cn](mailto:yangzhirui@westlake.edu.cn)  
Mengdie Yao: [yaomengdie@westlake.edu.cn](mailto:yaomengdie@westlake.edu.cn)  
Yanmei Dou: [douyanmei@westlake.edu.cn](mailto:douyanmei@westlake.edu.cn)
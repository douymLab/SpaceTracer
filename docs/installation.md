## Installation

SpaceTracer requires **Python 3.9**, `samtools`, `bedtools` and `bcftools`.

### Clone the repository

First, clone the SpaceTracer repository from GitHub:

```bash
git clone https://github.com/douymLab/SpaceTracer.git
cd SpaceTracer
```

### Prepare dependencies

Before running SpaceTracer, make sure all required dependencies are available in your environment. We provide two recommended ways to prepare the environment.


#### Option 1: Conda environment

You can create your own Conda environment using the provided `environment.yaml` file:

```bash
conda env create -f environment.yaml
conda activate SpaceTracer
```

Because `environment.yaml` contains many dependencies, we **strongly recommend using `mamba`** to create the environment:

```bash
mamba env create -f environment.yaml
conda activate SpaceTracer
```

!!! note
    During environment creation, the `pot` or `sorted-nearest` packages may fail to install on some systems. If this happens, activate the `SpaceTracer` environment and install them separately:

    ```bash
    mamba install -c conda-forge pot
    pip install sorted-nearest
    ```


Install SpaceTracer into the activated environment (this provides the `spacetracer` command):

```bash
pip install .
spacetracer --help
```

!!! note
    Conda must be installed before running these commands.
    Please refer to the [official Conda documentation](https://docs.conda.io/en/latest/).
    If you plan to use mamba, please install it first from the [mamba documentation](https://mamba.readthedocs.io/en/latest/).


#### Option 2: Docker image

We provide a Docker image, `spacetracer`, with all required dependencies preinstalled. The image is available on [Docker Hub](https://hub.docker.com/r/xiayh17/spacetracer):

```bash
docker pull xiayh17/spacetracer
docker run -it -v $(pwd):/mnt/workflow xiayh17/spacetracer bash
```

You can also download the container image file `spacetracer_latest.sif` from [Figshare](https://figshare.com/s/c7836f53c4eafb556ee1).

This command runs the Docker container in interactive mode (`-it`), mounts the current working directory (`$(pwd)`) to `/mnt/workflow` inside the container, and starts a Bash shell.

!!! note
    Docker must be installed on your machine before using this option.
    If Docker is not installed, please follow the [official Docker installation guide](https://docs.docker.com/get-started/get-docker/).


#### Additional Python packages

If needed, install extras from:

```bash
pip install -r requirements_pip.txt
```

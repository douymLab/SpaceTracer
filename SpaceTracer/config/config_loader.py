import json
import os
from pathlib import Path
import yaml
from typing import Dict, Any
from multiprocessing import cpu_count
from importlib.resources import files

from SpaceTracer.utils.get_genome_info import GenomeDetails


def check_file_exist(path):
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(
            f'Path {path} not exist! Please check your command or config file!'
        )
    return True


def deep_update(base: dict, override: dict) -> dict:
    """
    Recursively merge override into base.
    Only skip values that are None.
    """
    for key, value in override.items():
        if (
            key in base
            and isinstance(base[key], dict)
            and isinstance(value, dict)
        ):
            deep_update(base[key], value)
        else:
            if value is not None:
                base[key] = value
    return base


class LoadConfig:
    def __init__(self):
        self.config = {}

    def _load_yaml_config(self, config_path: str = None, required: bool = False) -> Dict[str, Any]:
        if not config_path:
            if required:
                raise FileNotFoundError("User config file is required but not provided.")
            return {}

        if not os.path.exists(config_path):
            raise FileNotFoundError(f'Cannot find config file: {config_path}')

        with open(config_path, "r", encoding="utf-8") as f:
            config = yaml.safe_load(f)

        if config is None:
            return {}
        if not isinstance(config, dict):
            raise ValueError(f"Config file must contain a yaml mapping: {config_path}")
        return config

    def _load_builtin_step_config(self) -> Dict[str, Any]:
        try:
            config_path = files("SpaceTracer.config").joinpath("default_step_config.yaml")
            with open(config_path, "r", encoding="utf-8") as f:
                config = yaml.safe_load(f)

            if config is None:
                return {}
            if not isinstance(config, dict):
                raise ValueError("Built-in default_step_config.yaml must be a yaml mapping.")
            return config
        except FileNotFoundError:
            return {}

    def _normalize_cli_kwargs(self, kwargs: Dict[str, Any]) -> Dict[str, Any]:
        config = {}
        for key, value in kwargs.items():
            if key == "custom_config":
                continue
            if value is None:
                continue

            if isinstance(value, str) and value.lower() in ["true", "false"]:
                value = json.loads(value.lower())

            config[key] = value
        return config

    def _replace_config(self, config_key: str, replace_dict: Dict[str, Any]):
        if not isinstance(replace_dict, dict):
            return
        if config_key in replace_dict and replace_dict[config_key] is not None:
            self.config[config_key] = replace_dict[config_key]

    def _check_required_values(self, check_list):
        for key in check_list:
            if self.config.get(key) is None or self.config.get(key) == "":
                raise ValueError(f'You have not provided required config: {key}')

    def _validate_basic_inputs(self):
        self._check_required_values([
            "genome",
            "sequence_type",
            "output_dir",
            "sample",
        ])

    def _resolve_input_files(self):
        input_details = self.config.get("input_details", {})
        sequence_type = self.config.get("sequence_type")
        spaceranger_dir = self.config.get("spaceranger_dir")

        if sequence_type == "stereo":
            if self.config.get("bin_size") is None:
                self.config["bin_size"] = 100
            self._replace_config("bam_file", input_details)
            self._replace_config("barcode_mapping", input_details)
            self.config["barcode_key"] = None

            self.config["tissue_position"] = ""

            if self.config.get("cell_number") is None:
                self.config["cell_number"] = 1

        elif sequence_type == "visium-HD":
            if self.config.get("bin_size") is None:
                self.config["bin_size"] = 16

            if spaceranger_dir:
                in_dir = Path(spaceranger_dir)
                self.config["bam_file"] = str(in_dir / "possorted_genome_bam.bam")
                # self.config["tissue_position"] = str(
                #     in_dir / "binned_outputs/square_002um/spatial/tissue_positions_list.csv"
                # )
                spatial_dir = in_dir / "binned_outputs/square_002um/spatial/"
                candidates = [
                    spatial_dir / "tissue_positions_list.csv",
                    spatial_dir / "tissue_positions.csv",
                    spatial_dir / "tissue_positions.parquet"
                ]
                for f in candidates:
                    if f.exists():
                        self.config["tissue_position"] = str(f)
                        break
                else:
                    raise FileNotFoundError(
                        f"Neither 'tissue_positions_list.csv', 'tissue_positions.csv' nor 'tissue_positions.parquet' was found in {spatial_dir}. "
                        f"Consider using the advanced config to provide the tissue position file path."
                    )
                self.config["barcode_mapping"] = str(in_dir / "barcode_mappings.parquet")

            else:
                self._replace_config("bam_file", input_details)
                self._replace_config("tissue_position", input_details)
                self._replace_config("barcode_mapping", input_details)
                
            if input_details.get("barcode_key") is not None:
                self.config["barcode_key"] = input_details["barcode_key"]
            else: 
                self.config["barcode_key"] = "CB"

            if self.config.get("cell_number") is None:
                self.config["cell_number"] = 1
                
        elif sequence_type == "visium":
            if input_details.get("barcode_key") is not None:
                    self.config["barcode_key"] = input_details["barcode_key"]
            else: 
                self.config["barcode_key"] = "CB"

            self.config["bin_size"] = None

            if spaceranger_dir:
                in_dir = Path(spaceranger_dir)
                self.config["bam_file"] = str(in_dir / "possorted_genome_bam.bam")
                spatial_dir = in_dir / "spatial"
                candidates = [
                    spatial_dir / "tissue_positions_list.csv",
                    spatial_dir / "tissue_positions.csv"
                ]
                for f in candidates:
                    if f.exists():
                        self.config["tissue_position"] = str(f)
                        break
                else:
                    raise FileNotFoundError(
                        f"Neither 'tissue_positions_list.csv' nor 'tissue_positions.csv' was found in {spatial_dir}. "
                        f"Consider using the advanced config to provide the tissue position file path."
                    )
            else:
                self._replace_config("bam_file", input_details)
                self._replace_config("tissue_position", input_details)


        elif sequence_type == "ST":
            if input_details.get("barcode_key") is not None:
                self.config["barcode_key"] = input_details["barcode_key"]
            else: 
                self.config["barcode_key"] = "B0"
            self._replace_config("bam_file", input_details)
            self._replace_config("tissue_position", input_details)
            self._replace_config("barcode_mapping", input_details)

        else:
            self._replace_config("bam_file", input_details)
            self._replace_config("tissue_position", input_details)
            self._replace_config("barcode_mapping", input_details)
            self._replace_config("barcode_key", input_details)

    def _resolve_resource_files(self):
        resource_dir = self.config.get("resource_dir")
        resource_details = self.config.get("resource_details", {})

        if resource_dir:
            in_dir = Path(resource_dir)
            self.config["genome_fasta"] = str(in_dir / "genome.fa")
            self.config["gnomad_path"] = str(in_dir / "gnomad_af")
            self.config["mappability_path"] = str(in_dir / "k24.umap.bedgraph")
            self.config["gene_bed"] = str(in_dir / "gene_region.bed")
            self.config["dbsnp_vcf_file"] = str(in_dir / "Homo_sapiens_assembly38.dbsnp138.vcf.gz")
            self.config["imprinted_bed"] = str(in_dir / "imprinted_genes.region.bed")
            self.config["editing_bed"] = str(in_dir / "known_editing.bed")
            self.config["PON_file"] = str(in_dir / "PON.txt")
            self.config["reference_error_profile"] = str(in_dir / "Artifacts.Sigprofile.txt")

        for key in [
            "genome_fasta",
            "gnomad_path",
            "mappability_path",
            "gene_bed",
            "dbsnp_vcf_file",
            "imprinted_bed",
            "editing_bed",
            "PON_file",
            "reference_error_profile",
        ]:
            self._replace_config(key, resource_details)

    def _resolve_model_defaults(self):
        if not self.config.get("model_name"):
            self.config["model_name"]="spatial_preserved_model"
        if not self.config.get("model_dir"):
            self.config["model_dir"] = str(files("SpaceTracer").joinpath("models"))

    def _sync_top_level_to_steps(self):
        steps = self.config.setdefault("steps", {})
        steps.setdefault("cluster", {})
        steps.setdefault("mutation_prediction", {})
        steps.setdefault("read_feature", {})

        if self.config.get("cluster_file") is not None:
            steps["cluster"]["cluster_file"] = self.config["cluster_file"]

        if self.config.get("cluster_method") is not None:
            steps["cluster"]["method"] = self.config["cluster_method"]

        if self.config.get("cell_number") is not None:
            steps["cell_number"] = self.config["cell_number"]

        if self.config.get("cell_info") is not None:
            steps["read_feature"]["cell_info"] = self.config["cell_info"]

        if self.config.get("model_name") is not None:
            steps["mutation_prediction"]["model_name"] = self.config["model_name"]

        if self.config.get("model_dir") is not None:
            steps["mutation_prediction"]["model_dir"] = self.config["model_dir"]

    def _validate_resolved_inputs(self):
        if not self.config.get("bam_file"):
            raise ValueError(
                "bam_file is required. Please provide spaceranger_dir or input_details.bam_file"
            )

        if not self.config.get("tissue_position") and self.config.get("sequence_type") != "stereo":
            raise ValueError(
                "tissue_position is required for this sequence_type. "
                "Please provide spaceranger_dir or input_details.tissue_position"
            )

        if not self.config.get("model_name"):
            raise ValueError(
                "model_name is required. Please provide model_name in config.yaml"
            )

    def _parse_memory_to_bytes(self, memory_value):
        if isinstance(memory_value, int):
            return memory_value

        memory_str = str(memory_value).strip().upper()
        if memory_str.endswith("G"):
            return int(float(memory_str[:-1]) * (1024 ** 3))
        if memory_str.endswith("M"):
            return int(float(memory_str[:-1]) * (1024 ** 2))
        if memory_str.endswith("K"):
            return int(float(memory_str[:-1]) * 1024)
        return int(memory_str)

    def load_config(self, **kwargs) -> Dict[str, Any]:
        DEFAULT_CONFIG = {
            "sample": 'Sample',
            "genome": None,
            "resource_dir": None,
            "sequence_type": None,
            "bin_size": None,
            "spaceranger_dir": None,
            "regions_file": None,
            "output_dir": None,

            "cluster_file": None,
            "cluster_method": None,
            "cell_number": None,
            "cell_info": None,
            "model_name": None,
            "model_dir": None,

            "bam_file": None,
            "tissue_position": None,
            "barcode_key": None,
            "barcode_mapping": None,

            "genome_fasta": None,
            "gnomad_path": None,
            "mappability_path": None,
            "gene_bed": None,
            "dbsnp_vcf_file": None,
            "imprinted_bed": None,
            "editing_bed": None,
            "PON_file": None,
            "reference_error_profile": None,

            "input_details": {
                "bam_file": None,
                "tissue_position": None,
                "barcode_key": None,
                "barcode_mapping": None,
            },

            "resource_details": {
                "genome_fasta": None,
                "gnomad_path": None,
                "mappability_path": None,
                "gene_bed": None,
                "dbsnp_vcf_file": None,
                "imprinted_bed": None,
                "editing_bed": None,
                "PON_file": None,
                "reference_error_profile": None,
            },

            "run": {
                "threads": 4,
                "memory": "32G",
                "max_parallel": 100,
                "keep_intermediates": True,
                "skip_validation": False,
            },

            "logging": {
                "level": "INFO",
                "log_file": "pipeline.log",
            },

            "steps": {},
        }

        builtin_step_config = self._load_builtin_step_config()
        user_config = self._load_yaml_config(kwargs.get("custom_config"), required=True)
        cli_config = self._normalize_cli_kwargs(kwargs)

        config = {}
        deep_update(config, DEFAULT_CONFIG)
        deep_update(config, builtin_step_config)
        deep_update(config, user_config)
        deep_update(config, cli_config)

        self.config = config

        self._validate_basic_inputs()
        self._resolve_input_files()
        self._resolve_resource_files()
        self._resolve_model_defaults()
        self._sync_top_level_to_steps()
        self._validate_resolved_inputs()

        self.config["run"]["threads"] = min(int(self.config["run"]["threads"]), cpu_count())
        self.config["run"]["memory"] = self._parse_memory_to_bytes(self.config["run"]["memory"])

        for key in [
            "regions_file",
            "spaceranger_dir",
            "cluster_file",
            "cell_info",
            "bam_file",
            "barcode_mapping",
            "resource_dir",
            "genome_fasta",
            "gnomad_path",
            "mappability_path",
            "gene_bed",
            "dbsnp_vcf_file",
            "imprinted_bed",
            "editing_bed",
            "PON_file",
            "reference_error_profile",
        ]:
            if self.config.get(key):
                check_file_exist(self.config[key])

        if self.config.get("tissue_position") not in [None, ""]:
            check_file_exist(self.config["tissue_position"])

        Genome = GenomeDetails(
            self.config.get("genome"),
            self.config.get("genome_fasta")
        )
        self.config["genome_details"] = Genome._get_genome_details()
        self.config.pop("input_details", None)
        self.config.pop("resource_details", None)

        return self.config

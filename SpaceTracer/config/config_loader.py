import json
import os
from pathlib import Path
import yaml
from typing import Dict, Any
from multiprocessing import cpu_count

from SpaceTracer.utils.get_genome_info import GenomeDetails

def check_file_exist(path):
    path=Path(path)
    if not path.exists():
        raise FileExistsError(f'Path {path} not exist! Please check your command or config file!')
    else:
        return True

class LoadConfig:
    def _replace_config(self, config_key, replace_dict):
        if config_key in replace_dict.keys():
            if replace_dict[config_key]:
                check_file_exist(replace_dict[config_key])
                self.config[config_key]=replace_dict[config_key]


    def _check_input_details(self):
        """
        If spaceranger dir provided, no matter in cli command or in configure file, 
        the bam file and tissue position file (not required for stereo-seq) will use the spaceranger result.
        """
        
        if "spaceranger_dir" in self.config.keys():
            in_dir=Path(self.config["spaceranger_dir"])
            self.config["bam_file"]=in_dir/"possorted_genome_bam.bam"
            self.config["tissue_position"]=in_dir/"spatial/tissue_positions_list.csv"
        else:
            self._replace_config("bam_file",self.config["input_details"])
            self._replace_config("tissue_position",self.config["input_details"])
        
        check_file_exist(self.config["bam_file"])
        if self.config["sequence_type"]=="visium":
            check_file_exist(self.config["tissue_position"])

        if self.config["input_details"]['barcode_key']:
            self.config['barcode_key']=self.config["input_details"]['barcode_key']
        else:
            if self.config['sequence_type']=='visium':
                self.config['barcode_key']="CB"


    def _check_resource_details(self):
        """
        If any resource_details provided in configure file, 
        the path will be replaced, rather than use the resource_dir.
        """
        if "resource_dir" in self.config.keys():
            in_dir=Path(self.config["resource_dir"])
            self.config["genome_fasta"]=str(in_dir/"genome.fa")
            self.config["gnomad_path"]=str(in_dir/"gnomad_af")
            self.config["mappability_path"]=str(in_dir/"k24.umap.bedgraph")
            self.config["gene_bed"]=str(in_dir/"gene_region.bed")
            self.config["dbsnp_vcf_file"]=str(in_dir/"Homo_sapiens_assembly38.dbsnp138.vcf.gz")
            self.config["imprinted_bed"]=str(in_dir/"imprinted_genes.region.bed")
            self.config["editing_bed"]=str(in_dir/"known_editing.bed")
            self.config["PON_file"]=str(in_dir/"PON.txt")
            self.config["reference_error_profile"]=str(in_dir/"Artifacts.Sigprofile.txt")

        details=self.config["resource_details"]

        resource_list=["genome_fasta","gnomad_path","mappability_path","gene_bed",
                        "dbsnp_vcf_file","imprinted_bed","editing_bed","PON_file","reference_error_profile"]
        for key in resource_list:
            self._replace_config(key, details)

        Genome=GenomeDetails(self.config.get('genome'),self.config.get('genome_fasta'))
        self.config['genome_details']=Genome._get_genome_details()
        
    def _check_required_values(self,check_list):
        for key in check_list:
            if not self.config[key]:
                raise ValueError(f'You have not provided the input {key}! Please check!')
            
    def load_config(self, **kwargs) -> Dict[str, Any]:
        """
        Load configuration with parameter priority.
        
        Priority order (from lowest to highest):
        1. Default parameters (built-in defaults)
        2. User custom parameters (from custom config file)
        3. Command-line parameters (highest priority, override all others)
        
        Args:
            genome: Genome name (e.g., 'hg38', 'mm10')
            **kwargs: Command-line parameters that will have highest priority
        
        Returns:
            Configuration dictionary with merged parameters
        """
        
        # Reference genome settings

        # 1. Default parameters (built-in defaults)
        DEFAULT_CONFIG = {'run':
            {
                'threads': 4,
                'memory': '32G',
                'skip_validation': False,
                'sequence_type': 'visium'
            }
        }

        # priority order
        config_sources = [
            DEFAULT_CONFIG, # 1. default
            self._load_custom_config(kwargs.get('custom_config')),  # 2 custom parameters fom config file
            {k: v for k, v in kwargs.items() if v is not None},  # 3. command-line parameters
        ]
        
        # combine parameters
        config = {}
        for source in config_sources:
            for key, value in source.items():
                if value is not None and value != "None":
                    if isinstance(value, str) and value.lower() in ['true', 'false']:
                        config[key] = json.loads(value.lower())
                    else:
                        config[key] = value
                        

        if 'bin_size' in config.keys() and config['sequence_type']!='visium':
            config['bin_size']=config.get('bin_size',100)
        else:
            config['bin_size']=None

        config['run']['threads']=min(config['run']['threads'],cpu_count())
        config['run']['memory']=int(config['run']['memory'].split("G")[0])*(1024**3)
        config['model_used']=config.get('model_used','spatial_preserved_model' )
        
        self.config=config

        self._check_input_details()
        self._check_resource_details()

        return self.config


    def _load_custom_config(self,custom_config_path: str = None) -> Dict[str, Any]:
        """load config file"""
        if not custom_config_path or not os.path.exists(custom_config_path):
            raise FileNotFoundError(f'Cannot find {custom_config_path}, please check!')
        
        try:
            with open(custom_config_path) as f:
                config = yaml.safe_load(f)
                return config if isinstance(config, dict) else {}
        except Exception:
            raise
            # return {}
        


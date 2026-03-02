import json
import os
from pathlib import Path
import yaml
from typing import Dict, Any

def check_file_exist(path):
    # if path=="":
    #     raise FileExistsError(f'You did not provide {path}, Please check your command or config file!')

    path=Path(path)
    if not path.exists():
        raise FileExistsError(f'Path {path} not exist! Please check your command or config file!')
    else:
        return True

class LoadConfig:
    def _replace_config(self, config_key, replace_dict):
        if config_key in replace_dict.keys():
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
            self.config["genome_fasta"]=in_dir/"genome.fa"
            self.config["gnomad_path"]=in_dir/"gnomad"

        details=self.config["resource_details"]
        self._replace_config("genome_fasta",details)
        self._replace_config("gnomad_path",details)
        
        print(self.config["genome_fasta"])

        check_file_exist(self.config["genome_fasta"])    
        check_file_exist(self.config["gnomad_path"])

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
        DEFAULT_CONFIG = {
            'threads': 4,
            'memory': '32G',
            'chunk_size': 1000000,
            'keep_intermediates': False,
            'skip_validation': False,
            'sequence_type': 'visium' 
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
            
        self.config=config

        self._check_input_details()
        self._check_resource_details()

        return self.config


    def _load_custom_config(self,custom_config_path: str = None) -> Dict[str, Any]:
        """load config file"""
        if not custom_config_path or not os.path.exists(custom_config_path):
            return {}
        
        try:
            with open(custom_config_path) as f:
                config = yaml.safe_load(f)
                return config if isinstance(config, dict) else {}
        except Exception:
            return {}
        

step_default_config={}


# def load_config(genome: str, **kwargs) -> Dict[str, Any]:
#     """
#     Load configuration with parameter priority.
    
#     Priority order (from lowest to highest):
#     1. Default parameters (built-in defaults)
#     2. User custom parameters (from custom config file)
#     3. Command-line parameters (highest priority, override all others)
    
#     Args:
#         genome: Genome name (e.g., 'hg38', 'mm10')
#         **kwargs: Command-line parameters that will have highest priority
    
#     Returns:
#         Configuration dictionary with merged parameters
#     """
    
#     # Reference genome settings
#     genome: null  # 将从命令行指定

#     # 1. Default parameters (built-in defaults)
#     DEFAULT_CONFIG = {
#         'threads': 4,
#         'memory': '32G',
#         'chunk_size': 1000000,
#         'keep_intermediates': False,
#         'skip_validation': False,
#         'sequence_type': 'visium' 
#     }

#     # priority order
#     config_sources = [
#         DEFAULT_CONFIG, # 1. default
#         load_custom_config(kwargs.get('custom_config')),  # 2 custom parameters fom config file
#         {k: v for k, v in kwargs.items() if v is not None},  # 3. command-line parameters
#     ]
    
#     # combine parameters
#     config = {'genome': genome}
#     for source in config_sources:
#         for key, value in source.items():
#             if value is not None:  # skip None value
#                 config[key] = value


#     if 'barcode_key' not in config.keys():
#         if config['sequence_type']=='visium':
#             config['barcode_key']="CB"

#     if 'bin_size' in config.keys() and config['sequence_type']!='visium':
#         config['bin_size']=config.get('bin_size',100)
#     else:
#         config['bin_size']=None
        
#     return config




# def load_custom_config(custom_config_path: str = None) -> Dict[str, Any]:
#     """load config file"""
#     if not custom_config_path or not os.path.exists(custom_config_path):
#         return {}
    
#     try:
#         with open(custom_config_path) as f:
#             config = yaml.safe_load(f)
#             return config if isinstance(config, dict) else {}
#     except Exception:
#         return {}


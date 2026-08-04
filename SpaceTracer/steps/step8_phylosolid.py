import os
from pathlib import Path
import subprocess
from SpaceTracer.utils.utils import barcode_cell_mapping
from SpaceTracer.steps.base import BaseStep
from SpaceTracer.cores.prepare_phylo_from_multi_sample import MutationExtractor
from SpaceTracer.utils.logger import get_logger

model_name=__name__
logger = get_logger(model_name)


class PhylogenyBuildStep(BaseStep):

    def get_inputs(self, context):
        inputs={
            "in_filter_bam": context.get("in_filter_bam"),
            'final_mutation_list':context.get('final_mutation_list')
        }
        return inputs

    def get_outputs(self, context):
        final_mutation_list=self.get_inputs(context)['final_mutation_list']

        with open(final_mutation_list, "r") as f:
            mut_number = len([line for line in f if line.strip()])

        if mut_number<3:
            logger.info(f"The pred mutation number({mut_number}) is smaller than 3. The phylogeny step will be skipped.")

            outputs={
                'final_mutation_list': context.get('final_mutation_list'),
            }
        else:
            outputs={
                'final_mutation_list': context.get('final_mutation_list'),
                'tree_pdf': os.path.join(self.step_dir,"tree/mutation_integrator/phylo/final_cleaned_M_full_basedPivots.filtered_sites_inferred.tree_scphylo.pdf")
            }
        return outputs


    def optional_parameters(self):
        """ That's optional parameters """
        parameters={}
        if self.config.get('tissue_position'):
            parameters['tissue_position'] = self.config.get('tissue_position')
        else:
            parameters['tissue_position'] = None

        parameters['cell_info']=self.config.get("steps", {}).get("read_feature", {}).get('cell_info','')
        parameters['target_barcodes']=self.config.get("steps", {}).get("phylogeny", {}).get('target_barcodes',None)

        return parameters
                

    def build_input_matrix(self, filter_bam, final_mutation_list, cell_dict, 
                        tissue_position_barcode, target_barcodes):
        if target_barcodes is not None and target_barcodes: 
            extractor = MutationExtractor(
                samples=self.sample,
                bams=filter_bam,
                mutlist=final_mutation_list,
                outprefix=str(self.step_dir) + "/" + self.sample,  
                cell_dict=cell_dict,
                target_barcodes=target_barcodes,
                bins=self.bin_size,                                    
                seq_type=self.seq_type,                       
                run_type="UMI"                     
            )
        
        elif tissue_position_barcode is not None and tissue_position_barcode:
            extractor = MutationExtractor(
                samples=self.sample,
                bams=filter_bam,
                mutlist=final_mutation_list,
                outprefix=str(self.step_dir) + "/" + self.sample,  
                cell_dict=cell_dict,
                barcode_files=tissue_position_barcode,  # 注意参数名
                bins=self.bin_size,                                    
                seq_type=self.seq_type,                       
                run_type="UMI"                     
            )
        
        else:
            raise FileNotFoundError(
                f'Neither target_barcodes ({target_barcodes}) nor '
                f'tissue_position_barcode ({tissue_position_barcode}) was found!'
            )
        
        out_name, out_barcode_file, spot_num = extractor.run()
        return out_name, out_barcode_file, spot_num 


    def prepare_data(self, matrix_path, barcode_path, spot_num, out_dir):
        log_file_path = self.log_dir / f"{self.sample}_phylosolid_prepare.log"
        r_script_path = (
            Path(__file__).resolve().parent.parent 
            / "third_party" 
            / "phylosolid" 
            / "rawPreprocess_spatial.extract_identifier_sites.R"
        )

        data_output_path = Path(out_dir) / "data"
        data_output_path.mkdir(parents=True, exist_ok=True)

        cmd = [
            "Rscript",
            str(r_script_path),
            "--cellnum", str(spot_num),
            "--inputfile", str(matrix_path),
            "--outputpath", str(data_output_path),
            "--scid_file", str(barcode_path),
            "-r", "no",
            "-t", "0.5",
            "-s", self.sample
        ]

        try:
            with open(log_file_path, "w") as log_file:
                subprocess.run(
                    cmd,
                    stdout=log_file,
                    stderr=subprocess.STDOUT,
                    check=True
                )
            logger.info("Rscript completed successfully.")
        except subprocess.CalledProcessError as e:
            logger.error(f"ERROR: Rscript failed! Detail logs: {log_file_path}")
            
            raise e
            
        return data_output_path


    def build_tree(self, data_output_path,tree_output_path):
        log_file_path = self.log_dir / f"{self.sample}_phylosolid_build_tree.log"
        
        py_script_path = (
            Path(__file__).resolve().parent.parent 
            / "third_party" 
            / "phylosolid" 
            / "run_phylosilid_fullTree_spacetracer.py"
        )

        cmd = [
            "python",
            str(py_script_path),
            "-s", self.sample,
            "-i", str(data_output_path),
            "-o", str(tree_output_path) + "/",
            "-c", "None",
            "--is_predict_germ", "no",
            "--is_detect_passtree_by_dp", "no",
            "--is_filter_quality", "yes"
        ]

        logger.info(f"Running build_tree command: {' '.join(cmd)}")
        logger.info(f"Tree building log will be saved to: {log_file_path}")

        try:
            with open(log_file_path, "w") as log_file:
                subprocess.run(
                    cmd,
                    stdout=log_file,
                    stderr=subprocess.STDOUT,
                    check=True
                )
            logger.info("Build tree completed successfully.")
        except subprocess.CalledProcessError as e:
            logger.error(f"ERROR: build_tree failed! Detail logs: {log_file_path}")
            raise e
        return tree_output_path


    def _run(self, context):

        self.log_dir = Path(self.step_dir) / "logs"
        self.log_dir.mkdir(parents=True, exist_ok=True)

        inputs=self.get_inputs(context)
        filter_bam=inputs["in_filter_bam"]
        final_mutation_list=inputs["final_mutation_list"]
        with open(final_mutation_list, "r") as f:
            mut_number = len([line for line in f if line.strip()])

        if mut_number>=3:
        
            parameters=self.optional_parameters()
            tissue_position_barcode=parameters["tissue_position"]
            target_barcodes=parameters["target_barcodes"]
            cell_info_file = parameters.get("cell_info","")
            
            if cell_info_file:
                cell_dict = barcode_cell_mapping(cell_info_file)
            else:
                cell_dict = {}
                
            if not tissue_position_barcode and not target_barcodes:
                logger.info(f'The current version of spacetracer not support build the tree automaticlly, please check the update of spacetracer.')

            else:
                tree_output_path = Path(self.step_dir) / "tree"
                tree_output_path.mkdir(parents=True, exist_ok=True)

                matrix_path, barcode_path, spot_num=self.build_input_matrix(filter_bam,final_mutation_list,cell_dict,tissue_position_barcode,target_barcodes)
                data_output_path=self.prepare_data(matrix_path, barcode_path, spot_num,self.step_dir)
                tree_output_path = self.build_tree(data_output_path,tree_output_path)



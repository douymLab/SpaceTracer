import os
from pathlib import Path
import subprocess
import glob
from SpaceTracer.utils.utils import barcode_cell_mapping
from SpaceTracer.steps.base import BaseStep
from SpaceTracer.cores.prepare_phylo_from_multi_sample import MutationExtractor
from SpaceTracer.utils.logger import get_logger

model_name=__name__
logger = get_logger(model_name)


class PhylogenyBuildStep(BaseStep):

    def _normalize_mutation_list_paths(self, mutation_lists):
        if not mutation_lists:
            return []
        if isinstance(mutation_lists, (str, Path)):
            return [str(mutation_lists)]
        return [str(path) for path in mutation_lists if path]

    def _merged_mutation_list_path(self):
        return os.path.join(self.step_dir, f"{self.sample}_merged_mutation_list.txt")

    def _merge_mutation_lists(self, mutation_lists, output_file=None):
        """Merge mutation list files and write sorted unique lines (sort -u)."""
        paths = self._normalize_mutation_list_paths(mutation_lists)
        if output_file is None:
            output_file = self._merged_mutation_list_path()

        out_dir = os.path.dirname(output_file)
        if out_dir:
            os.makedirs(out_dir, exist_ok=True)

        unique_lines = set()
        for path in paths:
            if not path or not os.path.exists(path):
                logger.warning(f"Mutation list not found, skip: {path}")
                continue
            with open(path, "r", encoding="utf-8", errors="ignore") as f:
                for line in f:
                    line = line.strip()
                    if line:
                        unique_lines.add(line)

        with open(output_file, "w", encoding="utf-8") as out:
            for line in sorted(unique_lines):
                out.write(line + "\n")

        logger.info(
            f"Merged {len(paths)} mutation list(s) into {output_file} "
            f"({len(unique_lines)} unique mutations)"
        )
        return output_file

    def _extract_prefix_from_mutation_list(self, mutation_list_path):
        """Extract prefix name from mutation list filename."""
        filename = os.path.basename(mutation_list_path)
        prefix = filename.removesuffix("_total_pred_truesites_PASS_mutation_list.txt")

        return prefix

    def get_inputs(self, context):
        inputs={
            "in_filter_bam": context.get("in_filter_bam"),
            'final_mutation_list':context.get('final_mutation_list')
        }
        return inputs

    def get_outputs(self, context):
        final_mutation_list = self._merged_mutation_list_path()
        if not os.path.exists(final_mutation_list):
            mutation_lists = self.get_inputs(context).get("final_mutation_list")
            if mutation_lists:
                final_mutation_list = self._merge_mutation_lists(mutation_lists)

        with open(final_mutation_list, "r") as f:
            mut_number = len([line for line in f if line.strip()])

        outputs = {
            'final_mutation_list': final_mutation_list,
        }

        if mut_number < 3:
            logger.info(f"The pred mutation number({mut_number}) is smaller than 3. The phylogeny step will be skipped.")
        else:
            tree_pdf_pattern = os.path.join(self.step_dir, "*_tree/04_final_results/phylo/*_inferred.tree_scphylo.pdf")
            tree_pdfs = glob.glob(tree_pdf_pattern)
            
            if tree_pdfs:
                outputs['tree_pdfs'] = tree_pdfs
                outputs['tree_pdf'] = tree_pdfs[0]
            else:
                outputs['tree_pdf'] = ""

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

    def build_input_matrix(self, filter_bam, mutation_list, cell_dict, 
                        tissue_position_barcode, target_barcodes, output_prefix):
        if target_barcodes is not None and target_barcodes: 
            extractor = MutationExtractor(
                samples=self.sample,
                bams=filter_bam,
                mutlist=mutation_list,
                outprefix=output_prefix,
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
                mutlist=mutation_list,
                outprefix=output_prefix,
                cell_dict=cell_dict,
                barcode_files=tissue_position_barcode,
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

    def prepare_data(self, matrix_path, barcode_path, spot_num, out_dir, prefix_name):
        log_file_path = self.log_dir / f"{prefix_name}_phylosolid_prepare.log"
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
            "-s", prefix_name
        ]

        try:
            with open(log_file_path, "w") as log_file:
                subprocess.run(
                    cmd,
                    stdout=log_file,
                    stderr=subprocess.STDOUT,
                    check=True
                )
            logger.info(f"Rscript prepare_data completed successfully for {prefix_name}.")
        except subprocess.CalledProcessError as e:
            logger.error(f"ERROR: Rscript failed for {prefix_name}! Detail logs: {log_file_path}")
            raise e
            
        return data_output_path

    def build_tree(self, data_output_path, tree_output_path, prefix_name):
        log_file_path = self.log_dir / f"{prefix_name}_phylosolid_build_tree.log"
        
        py_script_path = (
            Path(__file__).resolve().parent.parent 
            / "third_party" 
            / "phylosolid" 
            / "run_phylosilid_fullTree_spacetracer.py"
        )

        cmd = [
            "python",
            str(py_script_path),
            "-s", prefix_name,
            "-i", str(data_output_path),
            "-o", str(tree_output_path) + "/",
            "-c", "None",
            "--is_predict_germ", "no",
            "--is_detect_passtree_by_dp", "no",
            "--is_filter_quality", "yes"
        ]

        logger.info(f"Running build_tree command for {prefix_name}: {' '.join(cmd)}")
        logger.info(f"Tree building log will be saved to: {log_file_path}")

        try:
            with open(log_file_path, "w") as log_file:
                subprocess.run(
                    cmd,
                    stdout=log_file,
                    stderr=subprocess.STDOUT,
                    check=True
                )
            logger.info(f"Build tree completed successfully for {prefix_name}.")
        except subprocess.CalledProcessError as e:
            logger.error(f"ERROR: build_tree failed for {prefix_name}! Detail logs: {log_file_path}")
            raise e
        return tree_output_path

    def _run(self, context):
        self.log_dir = Path(self.step_dir) / "logs"
        self.log_dir.mkdir(parents=True, exist_ok=True)

        inputs = self.get_inputs(context)
        filter_bam = inputs["in_filter_bam"]
        mutation_lists = inputs["final_mutation_list"]
        final_mutation_list = self._merge_mutation_lists(mutation_lists)

        mutation_lists_normalized = self._normalize_mutation_list_paths(mutation_lists)
        all_mutation_lists = mutation_lists_normalized + [final_mutation_list]

        parameters = self.optional_parameters()
        tissue_position_barcode = parameters["tissue_position"]
        target_barcodes = parameters["target_barcodes"]
        cell_info_file = parameters.get("cell_info", "")
        
        if cell_info_file:
            cell_dict = barcode_cell_mapping(cell_info_file)
        else:
            cell_dict = {}

        if not tissue_position_barcode and not target_barcodes:
            logger.info(f'The current version of spacetracer does not support building the tree automatically, please check the update of spacetracer.')
            return

        for one_mutation_list in all_mutation_lists:
            if not os.path.exists(one_mutation_list):
                logger.warning(f"Mutation list not found, skip: {one_mutation_list}")
                continue

            prefix_name = self._extract_prefix_from_mutation_list(one_mutation_list)
            
            with open(one_mutation_list, "r") as f:
                mut_number = len([line for line in f if line.strip()])

            logger.info(f"Processing {prefix_name}: {mut_number} mutations")

            if mut_number < 3:
                logger.info(f"The mutation number ({mut_number}) for {prefix_name} is smaller than 3. Skipping tree building.")
                continue

            tree_output_path = Path(self.step_dir) / f"{prefix_name}_tree"
            tree_output_path.mkdir(parents=True, exist_ok=True)

            try:
                output_prefix = str(tree_output_path / prefix_name)
                matrix_path, barcode_path, spot_num = self.build_input_matrix(
                    filter_bam, one_mutation_list, cell_dict, 
                    tissue_position_barcode, target_barcodes, output_prefix
                )
                
                data_output_path = self.prepare_data(
                    matrix_path, barcode_path, spot_num, tree_output_path, prefix_name
                )
                
                self.build_tree(data_output_path, tree_output_path, prefix_name)
                
                logger.info(f"Successfully built tree for {prefix_name}")
            except Exception as e:
                logger.error(f"Failed to build tree for {prefix_name}: {str(e)}")
                continue

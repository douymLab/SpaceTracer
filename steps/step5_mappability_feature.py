import os
import pandas as pd
import gc

from SpaceTracer.steps.base import BaseStep
from SpaceTracer.cores.mappability_features import mappabilityFeatures
from SpaceTracer.utils.logger import get_logger
from SpaceTracer.utils.parallel import parallel_map
from SpaceTracer.utils.utils import load_manifest_tsv

model_name = __name__
logger = get_logger(model_name)


def merge_text_files(file_list, output_file, deduplicate=True):
    """
    合并多个文本文件为一个文件。
    默认去重，并跳过空文件/不存在文件。
    """
    valid_files = [f for f in file_list if f and os.path.exists(f) and os.path.getsize(f) > 0]

    out_dir = os.path.dirname(output_file)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    if not valid_files:
        open(output_file, "w").close()
        return output_file

    seen = set()
    with open(output_file, "w", encoding="utf-8") as out:
        for file in valid_files:
            with open(file, "r", encoding="utf-8", errors="ignore") as f:
                for line in f:
                    line = line.strip()
                    if not line:
                        continue
                    if deduplicate:
                        if line in seen:
                            continue
                        seen.add(line)
                    out.write(line + "\n")

    return output_file


def collect_files_from_manifest(manifest_file: str, column_name: str):
    rows = load_manifest_tsv(manifest_file)
    files = []
    for row in rows:
        if not row:
            continue
        f = row.get(column_name, "")
        if f:
            files.append(f)
    return files


class MappabilityFeatureStep(BaseStep):
    def get_inputs(self, context):
        return {
            "genotype_results": context.get("genotype_results", "")
        }

    def get_outputs(self, context):
        return {
            "mappability_feature": os.path.join(self.step_dir, "mappability_feature.txt"),
            "merged_ind_geno_filter_mutation_list": os.path.join(
                self.step_dir, "merged_ind_geno_filter_mutation_list.txt"
            ),
        }

    def _prepare_mutation_list(self, manifest_file, merged_output_file):
        """
        从 genotype chunk manifest 中提取 ind_geno_filter_mutation_list 并合并。
        """
        if not manifest_file:
            return None

        manifest_file = str(manifest_file)
        if not os.path.exists(manifest_file):
            return None

        file_list = collect_files_from_manifest(manifest_file, "ind_geno_filter_mutation_list")
        merged_file = merge_text_files(file_list, merged_output_file, deduplicate=True)
        return merged_file

    def _run(self, context):
        mappability_path = self.config["mappability_path"]

        # 预初始化一次，触发 chrom parquet 准备逻辑
        mappability_infos = mappabilityFeatures(mappability_path)
        del mappability_infos
        gc.collect()

        inputs = self.get_inputs(context)
        outputs = self.get_outputs(context)

        genotype_results = inputs["genotype_results"]
        if not genotype_results or not os.path.exists(genotype_results):
            raise FileNotFoundError(f"genotype_results manifest not found: {genotype_results}")

        merged_mutation_list = self._prepare_mutation_list(
            genotype_results,
            outputs["merged_ind_geno_filter_mutation_list"]
        )

        if (
            not merged_mutation_list
            or (not os.path.exists(merged_mutation_list))
            or os.path.getsize(merged_mutation_list) == 0
        ):
            logger.warning("No valid mutation list found, writing empty mappability feature output.")
            self._write_empty_outputs(outputs["mappability_feature"])
            return {
                "mappability_feature": outputs["mappability_feature"],
                "merged_ind_geno_filter_mutation_list": outputs["merged_ind_geno_filter_mutation_list"],
            }

        mutation_identifier_list = pd.read_csv(
            merged_mutation_list,
            sep="_",
            header=None,
            names=["chrom", "pos", "ref", "alt"]
        )

        if mutation_identifier_list.empty:
            logger.warning("Mutation list is empty after reading, writing empty mappability feature output.")
            self._write_empty_outputs(outputs["mappability_feature"])
            return {
                "mappability_feature": outputs["mappability_feature"],
                "merged_ind_geno_filter_mutation_list": outputs["merged_ind_geno_filter_mutation_list"],
            }

        chrom_groups = list(mutation_identifier_list.groupby("chrom"))

        results = []
        def process_chrom_wrapper(chrom_group):
            chrom, group = chrom_group
            return process_chromosome_mutations_for_mappable(chrom, group, mappability_path)

        results = parallel_map(
            list(chrom_groups),
            worker_fn=process_chrom_wrapper,
            max_workers=self.threads,
            desc="mappability_feature",
            raise_on_error=True
        )

        del chrom_groups
        del mutation_identifier_list
        gc.collect()

        results = [df for df in results if df is not None and not df.empty]
        gc.collect()

        # with ThreadPoolExecutor(max_workers=self.threads) as executor:
        #     future_to_chrom = {
        #         executor.submit(
        #             process_chromosome_mutations_for_mappable,
        #             chrom,
        #             group,
        #             mappability_path
        #         ): chrom
        #         for chrom, group in chrom_groups
        #     }

            # for future in tqdm(as_completed(future_to_chrom), total=len(future_to_chrom), desc="mappability_feature"):
            #     chrom = future_to_chrom[future]
            #     try:
            #         result_df = future.result()
            #         if result_df is not None and not result_df.empty:
            #             results.append(result_df)
            #     except Exception as e:
            #         logger.error(f"Failed to process {chrom}: {e}")
            #         raise

        if not results:
            logger.warning("No mappability results generated, writing empty outputs.")
            self._write_empty_outputs(outputs["mappability_feature"])
            return {
                "mappability_feature": outputs["mappability_feature"],
                "merged_ind_geno_filter_mutation_list": outputs["merged_ind_geno_filter_mutation_list"],
            }

        mutation_with_mappability = pd.concat(results, ignore_index=True)

        out_mappability_features = outputs["mappability_feature"]
        colnames = ["#chrom", "pos", "ref", "alt", "mappabilityScore"]
        mutation_with_mappability.columns = colnames
        mutation_with_mappability.to_csv(out_mappability_features, index=False, sep="	")

        mutation_with_mappability = mutation_with_mappability.set_index(["#chrom", "pos", "ref", "alt"])
        parquet_file = str(out_mappability_features).replace(".txt", ".parquet")
        mutation_with_mappability.to_parquet(
            parquet_file,
            index=True,
            engine="pyarrow",
            compression="snappy"
        )
        del mutation_with_mappability
        gc.collect()

        return {
            "mappability_feature": out_mappability_features,
            "merged_ind_geno_filter_mutation_list": outputs["merged_ind_geno_filter_mutation_list"],
        }

    def _write_empty_outputs(self, output_txt):
        empty_df = pd.DataFrame(columns=["#chrom", "pos", "ref", "alt", "mappabilityScore"])
        empty_df.to_csv(output_txt, index=False, sep="	")

        parquet_file = str(output_txt).replace(".txt", ".parquet")
        empty_df = empty_df.set_index(["#chrom", "pos", "ref", "alt"])
        empty_df.to_parquet(
            parquet_file,
            index=True,
            engine="pyarrow",
            compression="snappy"
        )


def process_chromosome_mutations_for_mappable(chrom: str, mutations: pd.DataFrame, mappability_path: str):
    """
    get the mappability score for positions in one chromosome
    """
    try:
        mappability_infos = mappabilityFeatures(mappability_path)
        chrom_data = mappability_infos._load_mappability_for_chrom(chrom)
        pos_list = mutations["pos"].tolist()
        mappability_scores = mappability_infos._query_mappability(
            chrom_data,
            pos_list
        )

        result_df = mutations.copy()
        result_df["mappability"] = mappability_scores
        
        del pos_list
        del mappability_scores
        del chrom_data
        del mappability_infos
        gc.collect()

        return result_df

    except Exception as e:
        logger.error(f"Error processing mappability scores for {chrom}: {e}")
        raise


from functools import partial
import multiprocessing
import os
from pathlib import Path
import pandas as pd
from tqdm import tqdm
from typing import Dict

from SpaceTracer.cores.spatial_features import handle_per_line_for_sf
from SpaceTracer.steps.base import BaseStep
from SpaceTracer.utils.read_files import handle_barcode, load_spot_genotypes_data
from SpaceTracer.utils.utils import check_dir, load_manifest_tsv,save_manifest_tsv
from SpaceTracer.utils.parallel import parallel_map

global_spot_geno_df = None
global_barcode_dict = None
global_barcode_dir = None
global_step_dir = None
global_in_name = None
global_alpha = None
global_thr_r2 = None
global_thr_prob = None
global_thr_likelihood = None
global_thr_vaf = None
global_plot_supp = None
global_fig_size = None
global_method = None
global_num_directions = None


def _str2bool(v):
    if isinstance(v, bool):
        return v
    if isinstance(v, str):
        return v.lower() in ("1", "true", "yes", "y", "t")
    return bool(v)


def init_worker(
    barcode_dir,
    step_dir,
    in_name,
    spot_geno_df,
    barcode_dict,
    alpha,
    thr_r2,
    thr_prob,
    thr_likelihood,
    thr_vaf,
    plot_supp,
    fig_size,
    method,
    num_directions
):
    global global_spot_geno_df, global_barcode_dict
    global global_barcode_dir, global_step_dir, global_in_name
    global global_alpha, global_thr_r2, global_thr_prob, global_thr_likelihood, global_thr_vaf
    global global_plot_supp, global_fig_size, global_method, global_num_directions

    global_spot_geno_df = spot_geno_df
    global_barcode_dict = barcode_dict
    global_barcode_dir = barcode_dir
    global_step_dir = step_dir
    global_in_name = in_name
    global_alpha = alpha
    global_thr_r2 = thr_r2
    global_thr_prob = thr_prob
    global_thr_likelihood = thr_likelihood
    global_thr_vaf = thr_vaf
    global_plot_supp = plot_supp
    global_fig_size = fig_size
    global_method = method
    global_num_directions = num_directions


def worker(mutation_id):
    return handle_per_line_for_sf(
        global_barcode_dir,
        global_step_dir,
        global_in_name,
        global_spot_geno_df,
        global_barcode_dict,
        global_alpha,
        global_thr_r2,
        global_thr_prob,
        global_thr_likelihood,
        global_thr_vaf,
        global_plot_supp,
        global_fig_size,
        global_method,
        global_num_directions,
        mutation_id
    )


class SpatialFeatureStep(BaseStep):
    def get_inputs(self, context):
        return {
            "genotype_results": context.get("genotype_results", "")
        }


    def get_outputs(self, context):
        sample = context.get("sample", "sample")
        out_dir = os.path.join(self.step_dir, sample)
        os.makedirs(out_dir, exist_ok=True)
        return {
            "spatial_feature_results": os.path.join(out_dir, "spatial_feature_chunk_manifest.tsv")
        }

    def get_step_config(self):
        return self.config.get("steps", {}).get("spatial_feature", {})

    def _get_chunk_output_paths(self, context: Dict, chunk: str) -> Dict[str, str]:
        sample = context.get("sample", "sample")
        chunk_dir = os.path.join(self.step_dir, sample, chunk)
        os.makedirs(chunk_dir, exist_ok=True)

        return {
            "spatial_feature_txt": os.path.join(chunk_dir, "spatial_feature.txt"),
            "spatial_feature_parquet": os.path.join(chunk_dir, "spatial_feature.parquet"),
            "barcode_dir": os.path.join(chunk_dir, "barcode_dir"),
        }

    def _run_one_chunk(
        self,
        row: Dict[str, str],
        tissue_positions: str,
        alpha: float,
        thr_r2: float,
        thr_prob: float,
        thr_likelihood: float,
        thr_vaf: float,
        plot_supp: bool,
        fig_size: int,
        method: str,
        num_directions: int,
        inner_processes: int,
        context: Dict,
    ) -> Dict[str, str]:

        chunk = row["chunk"]
        spot_genotype_file = row["spot_geno_file"]
        mutation_list_file = row["ind_geno_filter_mutation_list"]

        if not os.path.exists(spot_genotype_file):
            raise FileNotFoundError(f"spot_geno_file not found for chunk={chunk}: {spot_genotype_file}")
        if not os.path.exists(mutation_list_file):
            raise FileNotFoundError(f"mutation_list_file not found for chunk={chunk}: {mutation_list_file}")

        outputs = self._get_chunk_output_paths(context, chunk)
        out_spatial_features = outputs["spatial_feature_txt"]
        parquet_file = outputs["spatial_feature_parquet"]
        barcode_dir = Path(outputs["barcode_dir"])
        check_dir(barcode_dir)

        mutation_identifier_list = pd.read_csv(
            mutation_list_file,
            header=None,
            sep="\t"
        ).iloc[:, 0].tolist()

        # 空 chunk 直接生成空结果
        colnames = [
            '#chrom', 'pos', 'ref', 'alt', 'pass_spatial_test',
            'early_mutation', 'late_mutation',
            'mut_vs_nonmut_spots_KS_s', 'mut_vs_nonmut_spots_KS_p',
            'mut_vs_nonmut_spots_MI_s', 'mut_vs_nonmut_spots_MI_p',
            'mut_spots_prop', 'mut_spots_prop_by_probablity', 'mut_spots_prop_by_likelihood', 'mut_spots_prop_by_vaf',
            'all_spots_vaf_mean', 'all_spots_vaf_max', 'mut_spots_vaf_mean', 'mut_spots_vaf_median',
            'num_spots', 'num_mut_spots',
            'alt_vs_total_dp_r2', 'alt_vs_total_dp_paired_stat', 'alt_vs_total_dp_paired_wilcoxon_p', 'alt_vs_total_dp_paired_wilcoxon_rbc'
        ]

        if len(mutation_identifier_list) == 0:
            with open(out_spatial_features, "w") as f:
                f.write("\t".join(colnames) + "\n")
            empty_df = pd.DataFrame(columns=colnames)
            empty_df.to_parquet(parquet_file, index=False, engine="pyarrow", compression="snappy")
            return {
                "chunk": chunk,
                "ind_geno_filter_mutation_list": mutation_list_file,
                "spot_geno_file": spot_genotype_file,
                "spatial_feature_txt": out_spatial_features,
                "spatial_feature_parquet": parquet_file,
            }

        spot_geno_df = load_spot_genotypes_data(spot_genotype_file)
        barcode_dict = handle_barcode(tissue_positions)

        in_name = None

        # fork 只在 linux/mac 下更自然；如果你环境固定 linux，这样可以
        try:
            multiprocessing.set_start_method("fork", force=True)
        except RuntimeError:
            pass

        with (
            multiprocessing.Pool(
                processes=inner_processes,
                initializer=init_worker,
                initargs=(
                    barcode_dir,
                    self.step_dir,
                    in_name,
                    spot_geno_df,
                    barcode_dict,
                    alpha,
                    thr_r2,
                    thr_prob,
                    thr_likelihood,
                    thr_vaf,
                    plot_supp,
                    fig_size,
                    method,
                    num_directions
                )
            ) as pool,
            open(out_spatial_features, "w") as f
        ):
            f.write("\t".join(colnames) + "\n")

            for values in tqdm(
                pool.imap(worker, mutation_identifier_list, chunksize=2),
                total=len(mutation_identifier_list),
                desc=f"spatial_feature {chunk}"
            ):
                if values:
                    f.write("\t".join(map(str, values)) + "\n")

        df = pd.read_csv(out_spatial_features, sep="\t", header=0)

        if not df.empty:
            df = df.set_index(["#chrom", "pos", "ref", "alt"])
            df.to_parquet(parquet_file, index=True, engine="pyarrow", compression="snappy")
        else:
            df.to_parquet(parquet_file, index=False, engine="pyarrow", compression="snappy")

        return {
            "chunk": chunk,
            "ind_geno_filter_mutation_list": mutation_list_file,
            "spot_geno_file": spot_genotype_file,
            "spatial_feature_txt": out_spatial_features,
            "spatial_feature_parquet": parquet_file,
        }

    def _run(self, context):
        inputs = self.get_inputs(context)

        manifest_file = inputs["genotype_results"]
        if not manifest_file or not os.path.exists(manifest_file):
            raise FileNotFoundError(f"spatial_feature input manifest not found: {manifest_file}")

        outputs = self.get_outputs(context)
        result_manifest = outputs["spatial_feature_results"]

        parameters = self.get_step_config()
        alpha = float(parameters["alpha"])
        thr_r2 = float(parameters["thr_r2"])
        thr_prob = float(parameters["thr_prob"])
        thr_likelihood = float(parameters["thr_likelihood"])
        thr_vaf = float(parameters["thr_vaf"])
        plot_supp = _str2bool(parameters["plot_supp"])
        fig_size = int(parameters["fig_size"])
        method = str(parameters["method"])
        num_directions = int(parameters["num_directions"])

        tissue_positions = self.config["tissue_position"]

        rows = load_manifest_tsv(manifest_file)
        if not rows:
            raise ValueError(f"No chunk records found in manifest: {manifest_file}")

        # 只保留有必要输入的 chunk
        valid_rows = []
        for row in rows:
            if not row:
                continue
            if not row.get("chunk"):
                continue
            if not row.get("spot_geno_file"):
                continue
            if not row.get("ind_geno_filter_mutation_list"):
                continue
            valid_rows.append(row)

        if not valid_rows:
            raise ValueError(f"No valid chunk rows found in manifest: {manifest_file}")

        max_workers = self.config.get("runtime", {}).get("max_parallel", self.threads)

        step_cfg = self.get_step_config()
        inner_processes = int(step_cfg.get("per_chunk_processes", 1))

        def chunk_worker(row: Dict[str, str]) -> Dict[str, str]:
            return self._run_one_chunk(
                row=row,
                tissue_positions=tissue_positions,
                alpha=alpha,
                thr_r2=thr_r2,
                thr_prob=thr_prob,
                thr_likelihood=thr_likelihood,
                thr_vaf=thr_vaf,
                plot_supp=plot_supp,
                fig_size=fig_size,
                method=method,
                num_directions=num_directions,
                inner_processes=inner_processes,
                context=context,
            )

        chunk_results = parallel_map(
            valid_rows,
            worker_fn=chunk_worker,
            max_workers=max_workers,
            desc=f"{self.name} parallel spatial_feature for sample={context.get('sample', 'unknown')}",
            raise_on_error=True,
        )

        save_manifest_tsv(chunk_results, result_manifest)

import os
from pathlib import Path

import pandas as pd

from SpaceTracer.utils.utils import as_int, as_float, as_str
from SpaceTracer.steps.base import BaseStep
from SpaceTracer.cores.cluster_domain_segementation import *
from SpaceTracer.utils.logger import get_logger

model_name = __name__
logger = get_logger(model_name)


class ClusterStep(BaseStep):

    def get_inputs(self, context):
        return {}

    def get_outputs(self, context):
        cluster = self.config["cluster"]

        # user provides a cluster file
        if isinstance(cluster, str):
            if Path(cluster).exists():
                return {'cluster_file': cluster}
            else:
                raise FileNotFoundError(
                    f'Provided cluster file does not exist: cluster={cluster}, seq_type={getattr(self, "seq_type", "unknown")}'
                )

        # auto clustering
        elif isinstance(cluster, int):
            if cluster == 0:
                return {'cluster_file': os.path.join(self.step_dir, "cluster.txt")}
            else:
                raise ValueError(
                    f'Invalid cluster input: cluster={cluster}. '
                    f'If cluster is int, only 0 is allowed for auto clustering.'
                )

        else:
            raise ValueError(
                f'Wrong cluster input: cluster={cluster}, seq_type={getattr(self, "seq_type", "unknown")}'
            )

    def optional_parameters(self):
        spaceranger_dir = self.config.get("spaceranger_dir", "")
        if spaceranger_dir is not None:
            if Path(spaceranger_dir).exists():
                self.spaceranger_dir = spaceranger_dir
                return True
            else:
                raise FileNotFoundError(f'Spaceranger_dir: {spaceranger_dir} not found!')
        else:
            return None

    def get_step_config(self):
        return self.config.get('steps', {}).get('cluster', {})

    def _resolve_spatial_file_name(self):
        """
        Resolve spatial file path for different spaceranger versions.
        Prefer tissue_positions_list.csv over tissue_positions.csv if both exist.
        """
        if self.seq_type=="visium":
            spatial_dir = Path(self.spaceranger_dir) / "spatial"
        elif self.seq_type=="visium-HD":
            bin_size=self.bin_size
            resolution_name = f"square_{bin_size:03d}um"
            spatial_dir = Path(self.spaceranger_dir)/ "binned_outputs" / resolution_name /"spatial"

        tissue_positions_list = spatial_dir / "tissue_positions_list.csv"
        tissue_positions = spatial_dir / "tissue_positions.csv"

        if tissue_positions_list.exists():
            return "spatial/tissue_positions_list.csv"
        elif tissue_positions.exists():
            return "spatial/tissue_positions.csv"
        else:
            raise FileNotFoundError(
                f"Neither 'tissue_positions_list.csv' nor 'tissue_positions.csv' was found in {spatial_dir}. "
                f"Please use the advanced config to provide the tissue position file path."
            )

    def _run(self, context):

        cluster = self.config["cluster"]
        barcode_mapping_parquet_file = self.config.get("barcode_mapping", "")

        outputs = self.get_outputs(context)
        out_cluster_file = outputs["cluster_file"]

        # Case 1: user provides cluster file directly
        if isinstance(cluster, str):
            logger.info(f"Using user-provided cluster file: {cluster}")
            return

        # Case 2: auto clustering
        if isinstance(cluster, int) and cluster == 0:
            self._cluster_spot(
                out_cluster_file=out_cluster_file,
                barcode_mapping_parquet_file=barcode_mapping_parquet_file
            )
            return

        raise ValueError(
            f'Unsupported cluster configuration: cluster={cluster}, seq_type={self.seq_type}'
        )

    def _cluster_spot(self, out_cluster_file, barcode_mapping_parquet_file):
        """
        Cluster the spots using default / SpaGCN / GraphST.

        Logic:
        1. method == "default" and seq_type in ["visium", "visium-HD"]
           -> directly extract default cluster from spaceranger outputs
        2. method in ["SpaGCN", "GraphST"]
           -> try to run; if failed, fallback to default for visium/visium-HD
        3. if seq_type not in ["visium", "visium-HD"] and no supported strategy
           -> create a pseudo empty cluster file
        4. invalid configuration -> raise
        """
        plot = True
        histology = True
        R_HOME = None
        refinement = True

        sample = self.sample
        sequence_type = self.seq_type
        outdir = self.step_dir
        parameters = self.get_step_config()

        method = as_str(parameters.get('method', 'default'))
        ncluster = as_int(parameters['ncluster'])
        init_method = as_str(parameters['init_method'])

        weight_histology = as_float(parameters['weight_histology'])
        spot_area = as_float(parameters['spot_area'])
        percentage = as_float(parameters['percentage'])
        tolerance = as_float(parameters['tol'])
        learning_rate = as_float(parameters['lr'])
        max_epochs_run = as_int(parameters['max_epochs'])

        distance_threshold = as_float(parameters['distance_threshold'])
        num_threshold = as_int(parameters['num_threshold'])
        min_samples = as_int(parameters['min_samples'])
        radius = as_float(parameters['radius'])

        graphst_tool = as_str(parameters['graphst_tool'])
        seed = as_int(parameters['seed'])

        supported_methods = {"default", "SpaGCN", "GraphST"}
        if method not in supported_methods:
            raise ValueError(
                f"Unsupported clustering method: method={method}, seq_type={self.seq_type}, "
                f"cluster={self.config.get('cluster')}"
            )

        # 1. default for visium / visium-HD
        if method == "default":
            if self.seq_type in ["visium", "visium-HD"]:
                logger.info(f'Using default clustering result for seq_type={self.seq_type}.')
                result = self._get_default_cluster_for_visium(
                    out_cluster_file,
                    barcode_mapping_parquet_file
                )
                if result is None:
                    logger.info("Default cluster not available, creating a pseudo empty cluster file.")
                    self._touch_empty_file(out_cluster_file)
                return out_cluster_file
            else:
                logger.info(
                    f'No default clustering strategy for seq_type={self.seq_type}. '
                    f'Creating a pseudo empty cluster file: {out_cluster_file}'
                )
                self._touch_empty_file(out_cluster_file)
                return out_cluster_file

        # 2. SpaGCN / GraphST
        try:
            self.h5_file_name = "raw_feature_bc_matrix.h5"
            self.image_file_name = "spatial/tissue_hires_image.png"
            has_spaceranger_dir=self.optional_parameters()
            if has_spaceranger_dir and self.seq_type=="visium":
                indir = self.config.get("spaceranger_dir", "")

                self.spatial_file_name = self._resolve_spatial_file_name()

                if method == "SpaGCN":
                    logger.info("Running clustering with SpaGCN.")
                    cluster_spot_SpaGCN(
                        indir, self.h5_file_name, self.spatial_file_name, self.image_file_name,
                        outdir, sample, ncluster,
                        plot=plot, init_cluster=init_method, s=weight_histology,
                        b=spot_area, histology=histology, p=percentage, seed=seed,
                        tol=tolerance, lr=learning_rate, max_epochs_run=max_epochs_run,
                        l_start=0.01, l_end=1000, l_tol=0.01, l_max_run=100,
                        res_start=0.7, res_step=0.1, res_tol=5e-3, res_lr=0.05, res_max_epochs=20,
                        out_cluster_file=out_cluster_file
                    )

                    cluster2domain(
                        outdir, sample, plot=plot, distance_threshold=distance_threshold,
                        min_samples=min_samples, num_threshold=num_threshold, shape="hexagon",
                        keep=False, out_cluster_file=out_cluster_file
                    )

                elif method == "GraphST":
                    logger.info("Running clustering with GraphST.")
                    cluster_spot_GraphST(
                        indir, outdir, sample, ncluster, type=sequence_type,
                        h5_file_name=self.h5_file_name,
                        R_HOME=R_HOME, radius=radius, tool=graphst_tool,
                        refinement=refinement, plot=plot,
                        out_cluster_file=out_cluster_file
                    )

                    if sequence_type == "visium":
                        cluster2domain(
                            outdir, sample, plot=plot, distance_threshold=distance_threshold,
                            min_samples=min_samples, num_threshold=num_threshold,
                            shape="hexagon", keep=False,
                            out_cluster_file=out_cluster_file
                        )

        except Exception as e:
            logger.warning(f"Something wrong during clustering: {e}")

            if self.seq_type in ["visium", "visium-HD"]:
                logger.info("Trying to fall back to default cluster result.")
                result = self._get_default_cluster_for_visium(
                    out_cluster_file,
                    barcode_mapping_parquet_file
                )
                if result is None:
                    logger.info("Default cluster not available after fallback, creating a pseudo empty cluster file.")
                    self._touch_empty_file(out_cluster_file)
            else:
                logger.info(
                    f'Clustering failed and no default fallback is supported for seq_type={self.seq_type}. '
                    f'Creating a pseudo empty cluster file: {out_cluster_file}'
                )
                self._touch_empty_file(out_cluster_file)

        if not os.path.exists(out_cluster_file):
            logger.info("Cluster output missing, creating empty cluster file.")
            self._touch_empty_file(out_cluster_file)

        return out_cluster_file

    def _get_default_cluster_for_visium(self, out_cluster_file, barcode_mapping_parquet_file):
        in_dir = self.config.get("spaceranger_dir", "")
        if not in_dir:
            logger.warning('"spaceranger_dir" is empty, could not find the raw cluster file!')
            return None

        if not os.path.exists(in_dir):
            logger.warning('"spaceranger_dir" does not exist, could not find the raw cluster file!')
            return None

        if self.seq_type == "visium":
            cluster_file = os.path.join(
                in_dir,
                "analysis/clustering/gene_expression_graphclust/clusters.csv"
            )

            if not os.path.exists(cluster_file):
                logger.warning(f'Cluster file not found: {cluster_file}')
                return None

            cluster_df = pd.read_csv(cluster_file)
            cluster_df.to_csv(out_cluster_file, sep="\t", index=False, header=False)
            return out_cluster_file

        elif self.seq_type == "visium-HD":
            bins = self.config.get("bin_size", 16)
            if bins==2:
                logger.warning("No cluster file found under square_002um, use square_016um by default.")
                bins=16
            resolution_name = f"square_{bins:03d}um"

            cluster_file = os.path.join(
                in_dir,
                "binned_outputs",
                resolution_name,
                "analysis",
                "clustering",
                "gene_expression_graphclust",
                "clusters.csv"
            )

            if not os.path.exists(cluster_file):
                logger.warning(f'Cluster file not found: {cluster_file}')
                return None

            cluster_df = pd.read_csv(cluster_file)
            cluster_df = cluster_df.rename(columns={"Barcode": resolution_name})

            if not os.path.exists(barcode_mapping_parquet_file):
                logger.info(
                    "Oops! No barcode_mapping file provided! The following steps will run with a pseudo cluster."
                )
                return None

            barcode_mapping_parquet_df = pd.read_parquet(barcode_mapping_parquet_file)

            merge_df = cluster_df.merge(
                barcode_mapping_parquet_df,
                on=resolution_name,
                how="inner"
            )

            required_cols = ["square_002um", "Cluster"]
            missing_cols = [col for col in required_cols if col not in merge_df.columns]
            if missing_cols:
                logger.warning(f"Missing columns in merged dataframe: {missing_cols}")
                return None

            merge_df[required_cols].to_csv(
                out_cluster_file,
                sep="\t",
                index=False,
                header=False
            )
            return out_cluster_file

        else:
            logger.warning(f"Unsupported seq_type: {self.seq_type}")
            return None

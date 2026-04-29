
import pysam
from pathlib import Path

from SpaceTracer.utils.utils import as_int,as_float,as_str
from SpaceTracer.steps.base import BaseStep
from SpaceTracer.cores.cluster_domain_segementation import *

from SpaceTracer.utils.logger import get_logger
model_name=__name__
logger = get_logger(model_name)

class ClusterStep(BaseStep):

    def get_inputs(self, context):
        inputs={}
        return inputs

    def get_outputs(self, context):
        cluster=self.config["cluster"] # not None value, this has been handled in ochestrator
        cell_num=self.config["cell_num"]
        # cell number
        if not isinstance(cell_num,int) and Path(cell_num).exists(): # provided a cell num list
            out_cell_num_file=cell_num
        elif isinstance(cell_num,int) and cell_num!=0:
            out_cell_num_file=cell_num # a constant int value will be used as cell number input.
        else:
            out_cell_num_file = os.path.join(self.work_dir, "cell_num.txt") # only for visium in this version
        
        # cluster
        if isinstance(cluster,int):
            if cluster==0: # 0 means we'll run cluster function.
                out_cluster_file = os.path.join(self.step_dir, "cluster.txt")
            else: 
                raise ValueError('Why you provide one int for cluster?')
        elif Path(cluster).exists(): 
            out_cluster_file=cluster
        else:
            raise ValueError(f'Wrong cluster input {cluster}')

        return {'cluster_file': out_cluster_file,
                'cell_num': out_cell_num_file}


    def optional_parameters(self, context):
        spaceranger_dir = self.config.get("spaceranger_dir","")
        if Path(spaceranger_dir).exists():
            self.spaceranger_dir=spaceranger_dir
        else:
            raise FileNotFoundError(f'Spaceranger_dir: {spaceranger_dir} not Found!')

    def get_step_config(self):
        return self.config.get('steps', {}).get('cluster', {})
    

    def _run(self,context):
        self.h5_file_name = "raw_feature_bc_matrix.h5"
        self.spatial_file_name = "spatial/tissue_positions.csv"
        self.image_file_name = "spatial/tissue_hires_image.png"
        
        # get optional parameters
        self.optional_parameters(context)

        # run
        cluster=self.config["cluster"] 
        cell_num=self.config["cell_num"]
        cell_num_type=type(cell_num)

        # the running step is seperated from output step, due to the compleable situations
        out_cluster_file = os.path.join(self.step_dir, "cluster.txt")
        out_cell_num_file = os.path.join(self.work_dir, "cell_num.txt")

        if isinstance(cluster,int):
            if cluster==0: # 0 means we'll run cluster function.
                self._cluster_spot(out_cluster_file,out_cell_num_file)
            else: # other int means the data was cell/sub-cell level.
                raise ValueError('Why you provide one int for cluster?')

        elif Path(cluster).exists() : # also file is allowed
            if isinstance(cell_num,int):
                pass
            elif not Path(cell_num).exists() :
                cluster_df= pd.read_csv(cluster, sep="\t", header=None, names=['barcode', 'cluster'], na_values=[])
                cluster_df['cluster'] = cluster_df['cluster'].apply(lambda x: str(int(x)) if isinstance(x, float) and x.is_integer() 
                                                                else str(x) if pd.notnull(x) else "NA")
                self._get_umi_from_cluster(cluster_df,out_cell_num_file)
        else:
            raise ValueError(f'Wrong input for cluster: {cluster}')


    def _cluster_spot(self,out_cluster_file,out_cell_num_file):
        """
        Cluster the spots
        """
        indir=self.spaceranger_dir
        sample="Sample"
        plot=True
        histology=True
        sequence_type="visium"
        R_HOME=None
        refinement=True

        outdir=self.step_dir
        parameters=self.get_step_config()
        method = as_str(parameters['method'])
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


        # method=parameters['method']
        # ncluster=int(parameters['ncluster'])
        # init_method=parameters['init_method']
        # weight_histology=parameters['weight_histology']
        # spot_area=parameters['spot_area']
        # percentage=parameters['percentage']
        # tolerance=parameters['tol']
        # learning_rate=parameters['lr']
        # max_epochs_run=parameters['max_epochs']
        # distance_threshold=parameters['distance_threshold']
        # num_threshold=parameters['num_threshold']
        # min_samples=parameters['min_samples']
        # radius=parameters['radius']
        # graphst_tool=parameters['graphst_tool']

        # seed=parameters['seed']


        if method == "SpaGCN":
            # input files from indir
            # h5_file_name = "raw_feature_bc_matrix.h5"
            # spatial_file_name = "spatial/tissue_positions.csv"
            # image_file_name = "spatial/tissue_hires_image.png"
            # run SpaGCN
            cluster_spot_SpaGCN(indir, self.h5_file_name, self.spatial_file_name, self.image_file_name, outdir, sample, ncluster, \
                                plot=plot, init_cluster=init_method, s=weight_histology, b=spot_area, histology=histology, \
                                p=percentage, seed=seed, tol=tolerance, lr=learning_rate, max_epochs_run=max_epochs_run, \
                                l_start=0.01, l_end=1000, l_tol=0.01, l_max_run=100, \
                                res_start=0.7, res_step=0.1, res_tol=5e-3, res_lr=0.05, res_max_epochs=20, \
                                out_cluster_file=out_cluster_file,out_cell_num_file=out_cell_num_file)
            
            # divide seperate cluster into domains
            cluster2domain(outdir, sample, plot=plot, distance_threshold=distance_threshold, \
                        min_samples=min_samples, num_threshold=num_threshold, shape="hexagon", keep=False, \
                        out_cluster_file=out_cluster_file,out_cell_num_file=out_cell_num_file)
            
        if method == "GraphST":
            # run GraphST
            cluster_spot_GraphST(indir, outdir, sample, ncluster, type=sequence_type, \
                                h5_file_name=self.h5_file_name, \
                                R_HOME=R_HOME, radius=radius, tool=graphst_tool, refinement=refinement, plot=plot, \
                                out_cluster_file=out_cluster_file,out_cell_num_file=out_cell_num_file)

            if sequence_type == "visium":
                # divide seperate cluster into domains
                cluster2domain(outdir, sample, plot=plot, distance_threshold=distance_threshold, \
                            min_samples=min_samples, num_threshold=num_threshold, shape="hexagon", keep=False, \
                            out_cluster_file=out_cluster_file,out_cell_num_file=out_cell_num_file)
            
        return out_cluster_file,out_cell_num_file


    def _get_umi_from_cluster(self,cluster_df,out_cell_num_file):
        adata = read_10x_h5(os.path.join(self.spaceranger_dir, self.h5_file_name))
        # mark mitochondrial genes (human + mouse safe)
        adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
        # compute QC metrics (UMI counts)
        sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)
        # add spatial info
        spatial_file = os.path.join(self.spaceranger_dir, self.spatial_file_name)
        spatial = pd.read_csv(spatial_file,sep=",",header='infer',na_filter=False,index_col=0) 
        adata.obs["in_tissue"] = spatial['in_tissue']
        adata.obs["x_array"] = spatial['array_row']
        adata.obs["y_array"] = spatial['array_col']
        adata.obs["x_pixel"] = spatial['pxl_row_in_fullres']
        adata.obs["y_pixel"] = spatial['pxl_col_in_fullres']
        # Select captured samples
        adata = adata[adata.obs["in_tissue"]==1]
        adata.var_names = [i.upper() for i in list(adata.var_names)]
        adata.var["genename"] = adata.var.index.astype("str")
        adata.obs["UMI_counts"] = adata.obs["total_counts"]
        adata.obs["gene_counts"] = adata.obs["n_genes_by_counts"]
        sub_df=adata.obs[["UMI_counts"]]
        sub_df['barcode']=sub_df.index

        if cluster_df.empty:                        
            file_merged = pd.merge(cluster_df, sub_df, on='barcode')
        else:
            file_merged=sub_df.copy()
            file_merged['cluster']='bulk'

        cluster_sums = file_merged.groupby('cluster')['UMI_counts'].sum()
        max_cluster = cluster_sums.idxmax()
        max_cluster_data = file_merged[file_merged['cluster'] == max_cluster]
        median_nUMI = max_cluster_data['UMI_counts'].median()
        file_merged['calculated_cell_num'] = np.ceil(20 * file_merged['UMI_counts'] / median_nUMI)
        file_merged['refined_cell_num'] = np.where(file_merged['calculated_cell_num'] > 25, 25, file_merged['calculated_cell_num'])
        
        final_output_df = file_merged[['barcode', 'cluster', 'UMI_counts', 'refined_cell_num']]
        final_output_df.to_csv(out_cell_num_file, index=False,sep="\t")

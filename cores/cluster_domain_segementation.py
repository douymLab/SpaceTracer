import os
import warnings

os.environ["KMP_WARNINGS"] = "0"

warnings.filterwarnings(
    "ignore",
    message="pkg_resources is deprecated as an API.*"
)
warnings.filterwarnings(
    "ignore",
    message="This package has been superseded by the `leidenalg` package.*"
)

import torch
import pandas as pd
import numpy as np
import scanpy as sc
import math
import SpaGCN as spg
from scipy.sparse import issparse
import random, torch
import warnings
from sklearn.cluster import DBSCAN

from SpaceTracer.utils.utils import check_dir, check_file
warnings.filterwarnings("ignore")
import matplotlib.colors as clr
import matplotlib.pyplot as plt
import colorcet as cc
#In order to read in image data, we need to install some package. Here we recommend package "opencv"
#install opencv in python
#!pip3 install opencv-python
import cv2
from scanpy import read_10x_h5

from sklearn import metrics
import multiprocessing as mp
# pip install POT
from GraphST import GraphST
from GraphST.utils import clustering
import subprocess
import shutil

def cluster_spot_SpaGCN(data_dir, h5_file_name, spatial_file_name, image_file_name, output_dir, sample_name, n_clusters, \
                        plot=False, init_cluster="louvain", s=1, b=49, histology=True, p=0.5, seed=100, tol=5e-3, lr=0.05, max_epochs_run=200, \
                        l_start=0.01, l_end=1000, l_tol=0.01, l_max_run=100, \
                        res_start=0.7, res_step=0.1, res_tol=5e-3, res_lr=0.05, res_max_epochs=20,
                        out_cluster_file="",out_cell_num_file=""):
    """
    Use SpaGCN to cluster the spots
    Citation: Hu, J., Li, X., Coleman, K., Schroeder, A., Ma, N., Irwin, D. J., Lee, E. B., Shinohara, R. T., & Li, M. (2021). 
    SpaGCN: Integrating gene expression, spatial location and histology to identify spatial domains and spatially variable genes 
    by graph convolutional network. Nature Methods, 18(11), 1342-1351. https://doi.org/10.1038/s41592-021-01255-8.

    Inputs:
        data_dir - the direction saving the original h5 data file and histology image
        h5_file_name - the name of the h5 data file from 10X Visium
        spatial_file_name - the name of the spatial info
        image_file_name - the name of the histology image file
        output_dir - the direction saving the output figures
        sample_name - the sample name of the output files
        n_clusters - Number of spatial domains wanted
        plot - the Boolean variable of whether plotting the prediction and refined prediction (default = False)
        init_cluster - the method for initial clustering, or "kmeans" (default = "louvain")
        s - the weight given to histology when calculating Euclidean distance between every two spots. (default = 1)
            's = 1' means that the histology pixel intensity value has the same scale variance as the (x,y) coordinates, 
            whereas higher value of 's' indicates higher scale variance, hence, higher weight to histology, when calculating the Euclidean distance
        b - the area of each spot when extracting color intensity (default = 49)
        histology - whether use the histology image for calculating the adjacency matrix (default = True)
        p - Percentage of total expression contributed by neighborhoods (default = 0.5)
        seed - the random seed (default = 100)
        tol - the tolerance (default = 5e-3)
        lr - the learning rate (default = 0.05)
        max_epochs_run - the maximum number of epochs when running SpaGCN (default = 200)
        l_start, l_end, l_tol, l_max_run - the parameters for finding the value of l, which controls p (default = 0.01, 1000, 0.01, 100)
        res_start, res_step, res_tol, res_lr, res_max_epochs - the parameters for finding the resolution (deault = 0.7, 0.1, 5e-3, 0.05, 20)
        out_cluster_file - the output cluster file if provided
        out_cell_num_file - the output cell number file if provided
    """

    # Create plot folder
    res_dir = os.path.join(output_dir, "result")
    check_dir(res_dir)
    
    # Read in gene expression and spatial location
    h5ad_file = os.path.join(output_dir, sample_name+"_data.h5ad")
    check_h5ad = os.path.isfile(h5ad_file)
    if check_h5ad:
        adata=sc.read(h5ad_file)
        print("Read data from h5ad file")
    else:
        # Read original 10x_h5 data and save it to h5ad
        h5_file = os.path.join(data_dir, h5_file_name)
        adata = read_10x_h5(h5_file)
        # mark mitochondrial genes (human + mouse safe)
        adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
        # compute QC metrics (UMI counts)
        sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)
        # add spatial info
        spatial_file = os.path.join(data_dir, spatial_file_name)
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
        adata.write_h5ad(h5ad_file)
    
    # Set coordinates
    x_array = adata.obs["x_array"].tolist()
    y_array = adata.obs["y_array"].tolist()
    x_pixel = adata.obs["x_pixel"].tolist()
    y_pixel = adata.obs["y_pixel"].tolist()
    # Read in hitology image
    if histology:
        img_file = os.path.join(data_dir, image_file_name)
        img = cv2.imread(img_file)
        # Test coordinates on the image
        img_new = img.copy()
        for i in range(len(x_pixel)):
            x=x_pixel[i]
            y=y_pixel[i]
            img_new[int(x-20):int(x+20), int(y-20):int(y+20),:]=0
        map_file = os.path.join(output_dir, sample_name+"_map.jpg")
        cv2.imwrite(map_file, img_new)
    else:
        img = None

    # Calculate adjacent matrix
    adj_file = os.path.join(output_dir, sample_name+"_adj.csv")
    check_adj = os.path.isfile(adj_file)
    if check_adj:
        adj=np.loadtxt(adj_file, delimiter=',')
    else:
        adj = spg.calculate_adj_matrix(x=x_pixel,y=y_pixel, x_pixel=x_pixel, y_pixel=y_pixel, image=img, beta=b, alpha=s, histology=histology)
        # If histlogy image is not available, SpaGCN can calculate the adjacent matrix using the fnction below
        #adj = calculate_adj_matrix(x=x_pixel,y=y_pixel, histology=False)
        np.savetxt(adj_file, adj, delimiter=',')


    # Expression data preprocessing
    adata.var_names_make_unique()
    spg.prefilter_genes(adata,min_cells=3) # avoiding all genes are zeros
    spg.prefilter_specialgenes(adata)
    # Normalize and take log for UMI
    # sc.pp.normalize_per_cell(adata)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    # Find the l value given p, the parameter to control p
    l = spg.search_l(p, adj, start=l_start, end=l_end, tol=l_tol, max_run=l_max_run)


    # If the number of clusters known, we can use the spg.search_res() fnction to search for suitable resolution(optional)
    # Set seed
    r_seed = t_seed = n_seed = seed
    # Search for suitable resolution in the initial Louvain's Clustering methods
    res = spg.search_res(adata, adj, l, n_clusters, start=res_start, step=res_step, tol=res_tol, lr=res_lr, max_epochs=res_max_epochs, r_seed=r_seed, t_seed=t_seed, n_seed=n_seed)


    # Run SpaGCN
    clf = spg.SpaGCN()
    clf.set_l(l)
    # Set seed
    random.seed(r_seed)
    torch.manual_seed(t_seed)
    np.random.seed(n_seed)
    # Run
    clf.train(adata, adj, init_spa=True, init=init_cluster, res=res, tol=tol, lr=lr, max_epochs=max_epochs_run)
    y_pred, prob=clf.predict()
    adata.obs["pred"]= y_pred
    adata.obs["pred"]=adata.obs["pred"].astype('category')
    # Do cluster refinement(optional)
    # shape="hexagon" for Visium data, "square" for ST data.
    adj_2d = spg.calculate_adj_matrix(x=x_array,y=y_array, histology=False)
    refined_pred = spg.refine(sample_id=adata.obs.index.tolist(), pred=adata.obs["pred"].tolist(), dis=adj_2d, shape="hexagon")
    adata.obs["refined_pred"] = refined_pred
    adata.obs["refined_pred"] = adata.obs["refined_pred"].astype('category')
    # Save results
    results_file = os.path.join(res_dir, sample_name+"_results.h5ad")
    adata.write_h5ad(results_file)
    # Save cluster results
    # out_cluster_file = os.path.join(res_dir, sample_name+"_cluster.txt")

    # add UMI counts
    adata.obs["cluster"] = adata.obs["refined_pred"]
    adata.obs["UMI_counts"] = adata.obs["total_counts"]
    adata.obs["gene_counts"] = adata.obs["n_genes_by_counts"]
    adata.obs['mt_counts'] = adata.obs["total_counts_mt"]
    cluster_df = adata.obs[["cluster", "UMI_counts", "gene_counts", "mt_counts"]]
    cluster_df.to_csv(out_cluster_file, sep="\t", index=True, na_rep='NA', header=True)

    cluster_sums = cluster_df.groupby('cluster')['UMI_counts'].sum()
    max_cluster = cluster_sums.idxmax()
    max_cluster_data = cluster_df[cluster_df['cluster'] == max_cluster]
    median_nUMI = max_cluster_data['UMI_counts'].median()
    cluster_df['calculated_cell_num'] = np.ceil(20 * cluster_df['UMI_counts'] / median_nUMI)
    cluster_df['refined_cell_num'] = np.where(cluster_df['calculated_cell_num'] > 25, 25, cluster_df['calculated_cell_num'])
        
    # out_cell_num_file = os.path.join(res_dir, sample_name+"_cell_num.txt")
    final_output_df = cluster_df[['cluster', 'UMI_counts', 'refined_cell_num']]
    final_output_df.to_csv(out_cell_num_file, index=True,sep="\t")

    # Plot the results
    if plot:
        # Set colors used
        plot_color=["#F56867","#FEB915","#C798EE","#59BE86","#7495D3","#D1D1D1","#6D1A9C","#15821E","#3A84E6","#997273","#787878","#DB4C6C","#9E7A7A","#554236","#AF5F3C","#93796C","#F9BD3F","#DAB370","#877F6C","#268785"]
        # Plot spatial domains
        domains="pred"
        num_celltype=len(adata.obs[domains].unique())
        adata.uns[domains+"_colors"]=list(plot_color[:num_celltype])
        ax=sc.pl.scatter(adata,alpha=1,x="y_pixel",y="x_pixel",color=domains,title=sample_name+" "+domains,color_map=plot_color,show=False,size=100000/adata.shape[0])
        ax.set_aspect('equal', 'box')
        ax.axes.invert_yaxis()
        pred_fig = os.path.join(res_dir, sample_name+"_pred.png")
        plt.savefig(pred_fig, dpi=600)
        plt.close()

        # Plot refined spatial domains
        domains="refined_pred"
        num_celltype=len(adata.obs[domains].unique())
        adata.uns[domains+"_colors"]=list(plot_color[:num_celltype])
        ax=sc.pl.scatter(adata,alpha=1,x="y_pixel",y="x_pixel",color=domains,title=sample_name+" "+domains,color_map=plot_color,show=False,size=100000/adata.shape[0])
        ax.set_aspect('equal', 'box')
        ax.axes.invert_yaxis()
        refined_pred_fig = os.path.join(res_dir, sample_name+"_refined_pred.png")
        plt.savefig(refined_pred_fig, dpi=600)
        plt.close()



def cluster_spot_GraphST(data_dir, output_dir, sample_name, n_clusters, type="Visium", \
                         h5_file_name='filtered_feature_bc_matrix.h5', \
                         R_HOME=None, radius=6, tool='louvain', refinement=True, plot=False, \
                         out_cluster_file="",out_cell_num_file=""):
    """
    Use GraphST to cluster the spots
    Citation: Long, Y., Ang, K. S. & Li, M., et al. (2023). 
    Spatially informed clustering, integration, and deconvolution of spatial transcriptomics with GraphST. 
    Nature Communications, 14(1), 1155. https://doi.org/10.1038/s41467-023-36796-3

    Inputs:
        data_dir - the direction saving the original h5 data file and histology image
        output_dir - the direction saving the output figures
        sample_name - the sample name of the output files
        n_clusters - Number of spatial domains wanted
        type - the data type ('Visium' or 'Stereo', default = 'Visium')
        h5_file_name - the name of the h5 data file from 10X Visium
                       (default = 'filtered_feature_bc_matrix.h5')
        R_HOME - the R installation path (default =  None)
        radius - the number of neighbors considered during refinement (default = 6)
        tool - clustering method ('mclust', 'leiden', and 'louvain', default = 'louvain')
        refinement - Boolean variable of whether refine the clustering result (refinement = True)
        plot - the Boolean variable of whether plotting the prediction and refined prediction (default = False)
        out_cluster_file - the output cluster file if provided
        out_cell_num_file - the output cell number file if provided
    """

    # Create plot folder
    res_dir = os.path.join(output_dir, "result")
    check_dir(res_dir)

    # ===============================================
    # Preparation
    # ===============================================
    # Run device, by default, the package is implemented on 'cpu'. We recommend using GPU.
    device = torch.device('cuda:1' if torch.cuda.is_available() else 'cpu')
    # the location of R, which is necessary for mclust algorithm.
    def get_r_home():
        """Get the R installation path"""
        try:
            # command to execute R script that returns the R installation path
            command = 'Rscript -e "cat(R.home())"'
            # execute the command
            r_home = subprocess.check_output(command, shell=True, universal_newlines=True)
            return r_home.strip()
        except subprocess.CalledProcessError as e:
            print(f"Error: {e}")
            return None
        
    # Get R installation path with command
    if (R_HOME is None) and (tool=='mclust'):
        R_HOME = get_r_home()
    if (R_HOME is None) and (tool=='mclust'):
        tool = 'louvain'
        print("Note: mclust need R installation path, use louvain for clustering instead")
    elif (R_HOME is not None) and (tool=='mclust'):
        os.environ['R_HOME'] = R_HOME

    # check the name for tissue position (different version output from Visium) 
    position_file = os.path.join(data_dir, 'spatial/tissue_positions.csv')
    position_list_file = os.path.join(data_dir, 'spatial/tissue_positions_list.csv')
    if type == "Visium":
        # copy correct file
        if check_file(position_list_file):
            print("Read tissue_position_list.csv file")
        elif check_file(position_file):
            shutil.copy(position_file, position_list_file)
            print("Copy to tissue_position_list.csv file")
        else:
            print(f"No tissue position file exist.")


    # ===============================================
    # Train with GraphST
    # ===============================================
    # Read in gene expression and spatial location
    h5ad_file = os.path.join(output_dir, sample_name+"_embed.h5ad")
    check_h5ad = os.path.isfile(h5ad_file)
    if check_h5ad:
        adata=sc.read(h5ad_file)
        print("Read embeded data from h5ad file directly")
    else:
        print("Learn the representations with GraphST")
        # read data
        if type == "Visium":
            adata = sc.read_visium(data_dir, count_file=h5_file_name, load_images=True)
            # mark mitochondrial genes (human + mouse safe)
            adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
            # compute QC metrics (UMI counts)
            sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)
        else:
            h5_file = os.path.join(data_dir, h5_file_name)
            adata = sc.read_h5ad(h5_file)
            # compute QC metrics (UMI counts)
            adata.obs["n_genes_by_counts"] = np.asarray(adata.X.sum(axis=1)).ravel()
            adata.obs["total_counts"] = np.asarray((adata.X > 0).sum(axis=1)).ravel()
            adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
            adata.obs["total_counts_mt"] = np.asarray(adata.X[:, adata.var["mt"]].sum(axis=1)).ravel()

        adata.var_names_make_unique()
        # add spatial location
        if type == "Visium":
            spatial = pd.read_csv(position_file, sep=",", header='infer', na_filter=False, index_col=0) 
            adata.obs["x_array"] = spatial['array_row']
            adata.obs["y_array"] = spatial['array_col']
            adata.obs["x_pixel"] = spatial['pxl_row_in_fullres']
            adata.obs["y_pixel"] = spatial['pxl_col_in_fullres']
        else:
            spatial_fold = os.path.join(data_dir, 'spatial')
            check_dir(spatial_fold)
            # write spatial location file for further analysis
            spatial = pd.DataFrame()
            spatial['barcode'] = adata.obs.index
            spatial['in_tissue'] = 1
            spatial['array_row'] = adata.obsm['spatial'][:, 1]
            spatial['array_col'] = adata.obsm['spatial'][:, 0]
            spatial['pxl_row_in_fullres'] = adata.obsm['spatial'][:, 1]
            spatial['pxl_col_in_fullres'] = adata.obsm['spatial'][:, 0]
            spatial.to_csv(position_file, index=False, sep=',', header=True, na_rep='NA')
            # add info into obs for plot
            adata.obs['x_pixel'] = adata.obsm['spatial'][:, 1]
            adata.obs['y_pixel'] = adata.obsm['spatial'][:, 0]
            adata.obs['x_array'] = adata.obsm['spatial'][:, 1]
            adata.obs['y_array'] = adata.obsm['spatial'][:, 0]
        # check the value type
        adata.obsm['spatial'] = adata.obsm['spatial'].astype(float)
        # define model
        if type == "Visium":
            model = GraphST.GraphST(adata, device=device)
        else:
            model = GraphST.GraphST(adata, datatype='Stereo', device=device)
            refinement = False  # Stereo-seq do not need refinement
        # train model
        adata = model.train()
        # save embeded data
        adata.write_h5ad(h5ad_file)

    # ===============================================
    # Cluster and Refinement
    # ===============================================
    # clustering
    if tool == 'mclust':
        clustering(adata, n_clusters, radius=radius, method=tool, refinement=refinement)
    elif tool in ['leiden', 'louvain']:
        clustering(adata, n_clusters, radius=radius, method=tool, start=0.1, end=2.0, increment=0.01, refinement=refinement)
    # rename the cluster result columns
    if type == "Visium":
        adata.obs.rename(columns={
                        tool: 'pred',
                        'domain': 'refined_pred'
                        }, inplace=True)
    # Save results
    results_file = os.path.join(res_dir, sample_name+"_results.h5ad")
    adata.write_h5ad(results_file)
    # Save cluster results
    # out_cluster_file = os.path.join(res_dir, sample_name+"_cluster.txt")
    # add UMI counts
    adata.obs["UMI_counts"] = adata.obs["total_counts"]
    adata.obs["gene_counts"] = adata.obs["n_genes_by_counts"]
    adata.obs['mt_counts'] = adata.obs["total_counts_mt"]
    if type == "Visium":
        adata.obs["cluster"] = adata.obs["refined_pred"]
    else:
        adata.obs["cluster"] = adata.obs["domain"]
    cluster_df = adata.obs[["cluster", "UMI_counts", "gene_counts", "mt_counts"]]
    cluster_df.to_csv(out_cluster_file, sep="\t", index=True, na_rep='NA', header=True)

    cluster_sums = cluster_df.groupby('cluster')['UMI_counts'].sum()
    max_cluster = cluster_sums.idxmax()
    max_cluster_data = cluster_df[cluster_df['cluster'] == max_cluster]
    median_nUMI = max_cluster_data['UMI_counts'].median()
    cluster_df['calculated_cell_num'] = np.ceil(20 * cluster_df['UMI_counts'] / median_nUMI)
    cluster_df['refined_cell_num'] = np.where(cluster_df['calculated_cell_num'] > 25, 25, cluster_df['calculated_cell_num'])
        
    # out_cell_num_file = os.path.join(res_dir, sample_name+"_cell_num.txt")
    final_output_df = cluster_df[['cluster', 'UMI_counts', 'refined_cell_num']]
    final_output_df.to_csv(out_cell_num_file, index=True,sep="\t")
    
    # Plot the results
    if plot and type == "Visium":
        # Set colors used
        plot_color=["#F56867","#FEB915","#C798EE","#59BE86","#7495D3","#D1D1D1","#6D1A9C",
                    "#15821E","#3A84E6","#997273","#787878","#DB4C6C","#9E7A7A","#554236",
                    "#AF5F3C","#93796C","#F9BD3F","#DAB370","#877F6C","#268785"]
        # Plot spatial domains
        domains="pred"
        num_celltype=len(adata.obs[domains].unique())
        adata.uns[domains+"_colors"]=list(plot_color[:num_celltype])
        ax=sc.pl.scatter(adata,alpha=1,x="y_pixel",y="x_pixel",color=domains,title=sample_name+" "+domains,color_map=plot_color,show=False,size=100000/adata.shape[0])
        ax.set_aspect('equal', 'box')
        ax.axes.invert_yaxis()
        pred_fig = os.path.join(res_dir, sample_name+"_pred.png")
        plt.savefig(pred_fig, dpi=600)
        plt.close()

        # Plot refined spatial domains
        domains="refined_pred"
        num_celltype=len(adata.obs[domains].unique())
        adata.uns[domains+"_colors"]=list(plot_color[:num_celltype])
        ax=sc.pl.scatter(adata,alpha=1,x="y_pixel",y="x_pixel",color=domains,title=sample_name+" "+domains,color_map=plot_color,show=False,size=100000/adata.shape[0])
        ax.set_aspect('equal', 'box')
        ax.axes.invert_yaxis()
        refined_pred_fig = os.path.join(res_dir, sample_name+"_refined_pred.png")
        plt.savefig(refined_pred_fig, dpi=600)
        plt.close()

    if plot and type == "Stereo":
        # Set colors used
        plot_color=["#F56867","#556B2F","#C798EE","#59BE86","#006400","#8470FF",
                    "#CD69C9","#EE7621","#B22222","#FFD700","#CD5555","#DB4C6C",
                    "#8B658B","#1E90FF","#AF5F3C","#CAFF70", "#F9BD3F","#DAB370",
                    "#877F6C","#268785", '#82EF2D', '#B4EEB4']
        # Plot refined spatial domains only as they are the same
        domains="domain"
        num_celltype=len(adata.obs[domains].unique())
        if num_celltype > 22:
            plot_color = cc.glasbey_light
        adata.uns[domains+"_colors"]=list(plot_color[:num_celltype])
        ax=sc.pl.scatter(adata,alpha=1,x="y_pixel",y="x_pixel",color=domains,title=sample_name+" "+domains,color_map=plot_color,show=False,size=100000/adata.shape[0])
        ax.set_aspect('equal', 'box')
        ax.axes.invert_yaxis()
        refined_pred_fig = os.path.join(res_dir, sample_name+"_domain.png")
        plt.savefig(refined_pred_fig, dpi=600)
        plt.close()



def cluster2domain(output_dir, sample_name, plot=False, distance_threshold=2, min_samples=1, num_threshold=10, shape="hexagon", keep=False, \
    out_cluster_file="",out_cell_num_file=""):
    """
    Segment the seperated clusters into domains

    Inputs:
        output_dir - the direction saving the output figures and results
        sample_name - the sample name of the output files
        plot - the Boolean variable of whether plotting the domain scatter plot (default = False)
        distance_threshold - the maximum threshold for the distance between two spots which could be included in one domain (default = 2)
        min_samples - the minimum number of spots in one domain when seperation (default = 1)
        num_threshold - the minimum number of spots in one domain that would be left (default = 10)
        shape - the shape of the cluster location  (default = "hexagon")
        keep - whether keep the small domain if their result is still itself
    """

    # read and format the result data
    res_dir = os.path.join(output_dir, "result")
    results_file = os.path.join(res_dir, sample_name+"_results.h5ad")
    adata=sc.read(results_file)
    adata.obs.rename(columns={'pred': 'cluster', 'refined_pred': 'refined_cluster'}, inplace=True)

    # get the spots distance
    dis = spg.calculate_adj_matrix(x=adata.obs['x_array'],y=adata.obs['y_array'], histology=False)
    sample_id = adata.obs.index.tolist()
    dis_df = pd.DataFrame(dis, index=sample_id, columns=sample_id)

    # determine the neighbor number
    if shape=="hexagon":
        num_nbs=6 
    elif shape=="square":
        num_nbs=4
    else:
        print("Shape not recongized, shape='hexagon' for Visium data, 'square' for ST data.")

    # ===============================================
    # Seperate into domains
    # ===============================================
    # subdivide seperated clusters
    def subdivide_cluster(df_cluster):
        """segment seperated cluster into domains"""
        # eps is the maximum distance between two samples for them to be considered as in the same neighborhood
        dbscan = DBSCAN(eps=distance_threshold, min_samples=min_samples).fit(df_cluster[['x_array', 'y_array']])
        df_cluster['sub_cluster'] = df_cluster['refined_cluster'].astype(str) + '_' + dbscan.labels_.astype(str)

        # count the number of points in each cluster
        cluster_counts = df_cluster['sub_cluster'].value_counts()
        small_clusters = cluster_counts[cluster_counts < num_threshold].index
        # get the distance between the spots in the cluster
        dis_cluster_df = dis_df.loc[df_cluster['index'], df_cluster['index']]

        # process each small cluster
        for sub_cluster in small_clusters:
            # the points in the sub cluster
            small_cluster_points = df_cluster[df_cluster['sub_cluster'] == sub_cluster]['index']
            # aggregate neighbor information for the entire small cluster
            neighbor_labels = pd.Series(dtype='int')

            for point in small_cluster_points:
                # find the neighbors of this point
                dis_tmp = dis_cluster_df.loc[point, :].sort_values()
                nbs = dis_tmp.iloc[1:num_nbs+1]  # exclude the point itself
                nbs_pred = df_cluster[df_cluster['index'].isin(nbs.index)]['sub_cluster']
                neighbor_labels = pd.concat([neighbor_labels, nbs_pred])

            # determine the new label for the small clusters
            if keep:
                # determine the most common label among neighbors
                most_common_label = neighbor_labels.value_counts().idxmax()
            else:
                # determine the most common label except the original one
                neighbor_labels_without_original = neighbor_labels[neighbor_labels != sub_cluster]
                if not neighbor_labels_without_original.empty:
                    most_common_label = neighbor_labels_without_original.value_counts().idxmax()
                else:
                    most_common_label = sub_cluster
            # update the lables
            for point in small_cluster_points:
                df_cluster.loc[df_cluster['index'] == point, 'sub_cluster'] = most_common_label
        return df_cluster

    # divide the cluster into subclusters
    df = adata.obs
    df = df.reset_index()
    df_domains = df.groupby('refined_cluster').apply(subdivide_cluster)

    # reset the spot barcode as index
    df_domains.reset_index(drop=True, inplace=True)
    df_domains = df_domains.set_index('index')

    # save results
    adata.obs = df_domains
    adata.obs["sub_cluster"] = adata.obs["sub_cluster"].astype('category')
    adata.write_h5ad(results_file)
    # save cluster results
    # cluster_file = os.path.join(res_dir, sample_name+"_cluster.txt")
    # add UMI counts
    adata.obs["cluster"] = adata.obs["sub_cluster"]
    adata.obs["UMI_counts"] = adata.obs["total_counts"]
    adata.obs["gene_counts"] = adata.obs["n_genes_by_counts"]
    adata.obs['mt_counts'] = adata.obs["total_counts_mt"]
    cluster_df = adata.obs[["cluster", "UMI_counts", "gene_counts", "mt_counts"]]
    cluster_df.to_csv(out_cluster_file, sep="\t", index=True, na_rep='NA', header=True)
    
    cluster_sums = cluster_df.groupby('cluster')['UMI_counts'].sum()
    max_cluster = cluster_sums.idxmax()
    max_cluster_data = cluster_df[cluster_df['cluster'] == max_cluster]
    median_nUMI = max_cluster_data['UMI_counts'].median()
    cluster_df['calculated_cell_num'] = np.ceil(20 * cluster_df['UMI_counts'] / median_nUMI)
    cluster_df['refined_cell_num'] = np.where(cluster_df['calculated_cell_num'] > 25, 25, cluster_df['calculated_cell_num'])
        
    # out_cell_num_file = os.path.join(res_dir, sample_name+"_cell_num.txt")
    final_output_df = cluster_df[['cluster', 'UMI_counts', 'refined_cell_num']]
    final_output_df.to_csv(out_cell_num_file, index=True,sep="\t")
    # Plot the results
    if plot:
        # Set colors used
        # plot_color=["#F56867","#FEB915","#C798EE","#59BE86","#7495D3","#D1D1D1","#6D1A9C",
        #             "#15821E","#3A84E6","#997273","#787878","#DB4C6C","#9E7A7A","#554236",
        #             "#AF5F3C","#93796C","#F9BD3F","#DAB370","#877F6C","#268785"]
        # plot_color=["#F56867","#556B2F","#C798EE","#59BE86","#006400","#8470FF",
        #             "#CD69C9","#EE7621","#B22222","#FFD700","#CD5555","#DB4C6C",
        #             "#8B658B","#1E90FF","#AF5F3C","#CAFF70", "#F9BD3F","#DAB370",
        #             "#877F6C","#268785", '#82EF2D', '#B4EEB4']
        plot_color = cc.glasbey_light
        # Plot spatial domains
        domains="sub_cluster"
        num_celltype=len(adata.obs[domains].unique())
        adata.uns[domains+"_colors"]=list(plot_color[:num_celltype])
        ax=sc.pl.scatter(adata,alpha=1,x="y_pixel",y="x_pixel",color=domains,title=sample_name+" domain",color_map=plot_color,show=False,size=100000/adata.shape[0])
        ax.set_aspect('equal', 'box')
        ax.axes.invert_yaxis()
        pred_fig = os.path.join(res_dir, sample_name+"_domain.png")
        plt.savefig(pred_fig, dpi=600)
        plt.close()



def tissue_separate(data_dir, h5_file_name, spatial_file_name, output_dir, sample_name, n_tissues, \
                    distance_threshold=2, min_samples=1, shape="hexagon", plot=False):
    """
    Separate the tissues if there are more than one tissue in the slice

    Inputs:
        data_dir - the direction saving the original h5 data file and histology image
        h5_file_name - the name of the h5 data file from 10X Visium
        spatial_file_name - the name of the spatial info
        output_dir - the direction saving the output figures
        sample_name - the sample name of the output files
        n_tissues - Number of tissues the slice contained
        distance_threshold - the maximum threshold for the distance between two spots which could be included in one domain (default = 2)
        min_samples - the minimum number of spots in one domain when seperation (default = 1)
        shape - the shape of the cluster location  (default = "hexagon")
        plot - the Boolean variable of whether plotting the prediction and refined prediction (default = False)
    """

    # check plot folder
    res_dir = os.path.join(output_dir, "result")
    check_dir(res_dir)
    
    # check whether we have already format the location data into h5ad file
    # Read in gene expression and spatial location
    h5ad_file = os.path.join(output_dir, sample_name+"_data.h5ad")
    results_file = os.path.join(res_dir, sample_name+"_results.h5ad")
    check_h5ad = os.path.isfile(h5ad_file)
    check_result = os.path.isfile(results_file)
    if check_result:
        adata = sc.read(results_file)
    elif check_h5ad:
        adata=sc.read(h5ad_file)
    else:
        # Read original 10x_h5 data and save it to h5ad
        h5_file = os.path.join(data_dir, h5_file_name)
        adata = read_10x_h5(h5_file)
        # mark mitochondrial genes (human + mouse safe)
        adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
        # compute QC metrics (UMI counts)
        sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)
        # add spatial info
        spatial_file = os.path.join(data_dir, spatial_file_name)
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
        adata.write_h5ad(h5ad_file)
    
    # get the spots distance
    dis = spg.calculate_adj_matrix(x=adata.obs['x_array'],y=adata.obs['y_array'], histology=False)
    sample_id = adata.obs.index.tolist()
    dis_df = pd.DataFrame(dis, index=sample_id, columns=sample_id)

    # determine the neighbor number
    if shape=="hexagon":
        num_nbs=6 
    elif shape=="square":
        num_nbs=4
    else:
        print("Shape not recongized, shape='hexagon' for Visium data, 'square' for ST data.")
    
    # save the data frame
    df = adata.obs
    df = df.reset_index()

    # eps is the maximum distance between two samples for them to be considered as in the same neighborhood
    dbscan = DBSCAN(eps=distance_threshold, min_samples=min_samples).fit(df[['x_array', 'y_array']])
    df['tissue'] = dbscan.labels_.astype(str)

    # count the number of points in each cluster
    tissue_counts = df['tissue'].value_counts()
    tissue_counts_sorted = tissue_counts.sort_values(ascending=True)
    # get clusters excluding the largest two
    small_tissue_list = tissue_counts_sorted[:-n_tissues].index

    # process each small cluster
    for small_tissue in small_tissue_list:
        # the points in the small tissue
        small_tissue_points = df[df['tissue'] == small_tissue]['index']
        # aggregate neighbor information for the entire small cluster
        neighbor_labels = pd.Series(dtype='int')

        for point in small_tissue_points:
            # find the neighbors of this point
            dis_tmp = dis_df.loc[point, :].sort_values()
            nbs = dis_tmp.iloc[1:num_nbs+1]  # exclude the point itself
            nbs_pred = df[df['index'].isin(nbs.index)]['tissue']
            neighbor_labels = pd.concat([neighbor_labels, nbs_pred])

        # determine the new label for the small clusters
        neighbor_labels_without_original = neighbor_labels[neighbor_labels != small_tissue]
        if not neighbor_labels_without_original.empty:
            most_common_label = neighbor_labels_without_original.value_counts().idxmax()
        else:
            most_common_label = "-1"
        # update the lables
        for point in small_tissue_points:
            df.loc[df['index'] == point, 'tissue'] = most_common_label

    # reset the spot barcode as index
    df.reset_index(drop=True, inplace=True)
    df = df.set_index('index')

    # save results
    adata.obs = df
    adata.obs["tissue"] = adata.obs["tissue"].astype('category')
    tissue_h5ad_file = os.path.join(res_dir, sample_name+"_tissue.h5ad")
    adata.write_h5ad(tissue_h5ad_file)
    # divide the adata into several tissues and save
    tissue_labels = adata.obs['tissue'].unique()
    sample_tissue_list = []
    for tissue in tissue_labels:
        tissue_split_adata = adata[adata.obs['tissue'] == tissue].copy()
        sample_tissue_name = sample_name+"_"+tissue
        tissue_split_h5ad_file = os.path.join(output_dir, sample_tissue_name+"_data.h5ad")
        tissue_split_adata.write_h5ad(tissue_split_h5ad_file)
        sample_tissue_list.append(sample_tissue_name)
    # save tissue results
    tissue_file = os.path.join(res_dir, sample_name+"_tissue.txt")
    adata.obs['tissue'].to_csv(tissue_file, sep="\t", index=True, na_rep='NA', header=False)


    # Plot the results
    if plot:
        # Set colors used
        plot_color=["#F56867","#FEB915","#C798EE","#59BE86","#7495D3","#D1D1D1","#6D1A9C","#15821E","#3A84E6","#997273","#787878","#DB4C6C","#9E7A7A","#554236","#AF5F3C","#93796C","#F9BD3F","#DAB370","#877F6C","#268785"]
        # Plot spatial domains
        domains="tissue"
        num_celltype=len(adata.obs[domains].unique())
        adata.uns[domains+"_colors"]=list(plot_color[:num_celltype])
        ax=sc.pl.scatter(adata,alpha=1,x="y_pixel",y="x_pixel",color=domains,title=sample_name+" "+domains,color_map=plot_color,show=False,size=100000/adata.shape[0])
        ax.set_aspect('equal', 'box')
        ax.axes.invert_yaxis()
        tissue_fig = os.path.join(res_dir, sample_name+"_tissue.png")
        plt.savefig(tissue_fig, dpi=600)
        plt.close()

    return sample_tissue_list

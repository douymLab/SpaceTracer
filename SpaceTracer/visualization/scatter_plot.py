
import os

from matplotlib import pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import numpy as np
import pandas as pd
import colorcet as cc
from sklearn.cluster import DBSCAN
from scipy.spatial import ConvexHull
from matplotlib.patches import Polygon
from SpaceTracer.utils.utils import check_dir


    
class plot_scatter_mutant:
    def __init__(self,plot_dir,tissue_position_file,barcode_dir,point_size,input,fig_size=None,in_name=""):
        self.in_name=in_name
        self.plot_dir=plot_dir # output dir
        self.barcode_dir=barcode_dir # the path save barcode
        self.tissue_position_file=tissue_position_file
        check_dir(self.plot_dir)
        self.input=input
        self.fig_size=fig_size

        # get the spatial scatter figure range
        barcode_df = pd.read_csv(self.tissue_position_file)
        min_array_row, max_array_row = barcode_df['array_row'].agg(['min', 'max'])
        min_array_col, max_array_col = barcode_df['array_col'].agg(['min', 'max'])
        
        # calculate proper point size
        array_size_row = max_array_row - min_array_row  # range of array_row
        array_size_col = max_array_col - min_array_col  # range of array_col
        array_area = array_size_row * array_size_col

        if self.point_size is not None:
            # use given point size
            point_size = int(self.point_size)
        else:
            # total number of points
            num_points = len(barcode_df)
            # define a base size for points
            base_point_size = 1.5
            # calculate point size based on point density (inverse of the density)
            point_density = num_points / array_area if array_area != 0 else 1
            point_size = round((base_point_size / point_density)**2, 1)
            # clamping the point size to reasonable limits
            min_point_size = 5
            max_point_size = 50
            point_size = np.clip(point_size, min_point_size, max_point_size)

        self.point_size=point_size
        self.min_array_row=min_array_row
        self.max_array_row=max_array_row
        self.min_array_col=min_array_col
        self.max_array_col=max_array_col
        
    def plot(self,input):
        # a list of mutations
        if is_file(input):
            with open(input, 'r') as file:
                for line in file:
                    mutation_name = line.strip()
                    sample_name = mutation_name
                    barcode_mutinfo_file_name=os.path.join(self.barcode_dir, mutation_name+'.barcode.mutinfo.txt')
                    if os.path.exists(barcode_mutinfo_file_name):
                        barcode_mutinfo_df_colnames=["barcode_name", "cluster","in_tissue","pos_x","pos_y","depth","vaf","mut_prob","l_mosaic", \
                                                    "mutant_allele_num","has_mutant_allele","high_prob","high_likelihood","mutated"]
                        barcode_mutinfo_df=pd.read_csv(barcode_mutinfo_file_name, sep='\t', header=None, names=barcode_mutinfo_df_colnames, comment = "#")
                        # print(mutation_name)
                        plot_scatter_spatial(df=barcode_mutinfo_df, output_dir=self.plot_dir, sample_name=sample_name, column="mut_prob", note="prob", \
                                    fig_size=self.fig_size, xmin=self.min_array_row, xmax=self.max_array_row, ymin=self.min_array_col, ymax=self.max_array_col, point_size=self.point_size)
                        plot_scatter_spatial(df=barcode_mutinfo_df, output_dir=self.plot_dir, sample_name=sample_name, column="l_mosaic", note="likelihood", \
                                    fig_size=self.fig_size, xmin=self.min_array_row, xmax=self.max_array_row, ymin=self.min_array_col, ymax=self.max_array_col, point_size=self.point_size)
                    else:
                        print(mutation_name, "does not have spatial info")
        else:
            mutation_name = input
            # sample name
            sample_name = mutation_name
            # plot
            barcode_mutinfo_file_name=os.path.join(self.barcode_dir, mutation_name+'.barcode.mutinfo.txt')
            barcode_mutinfo_df_colnames=["barcode_name", "cluster","in_tissue","pos_x","pos_y","depth","vaf","mut_prob","l_mosaic", \
                                        "mutant_allele_num","has_mutant_allele","high_prob","high_likelihood","mutated"]
            barcode_mutinfo_df=pd.read_csv(barcode_mutinfo_file_name, sep='\t', header=None, names=barcode_mutinfo_df_colnames, comment = "#")
            # print(mutation_name)
            plot_scatter_spatial(df=barcode_mutinfo_df, output_dir=self.plot_dir, sample_name=sample_name, column="mut_prob", note="prob", \
                                fig_size=self.fig_size, xmin=self.min_array_row, xmax=self.max_array_row, ymin=self.min_array_col, ymax=self.max_array_col, point_size=self.point_size)
            plot_scatter_spatial(df=barcode_mutinfo_df, output_dir=self.plot_dir, sample_name=sample_name, column="l_mosaic", note="likelihood", \
                                fig_size=self.fig_size, xmin=self.min_array_row, xmax=self.max_array_row, ymin=self.min_array_col, ymax=self.max_array_col, point_size=self.point_size)


def is_file(string):
    if os.path.exists(string):
        return True
    else:
        return False
    

def barcode_to_location(args):
    """
    Get the spatial location from the spots barcode
    """
    # set output file
    if args.output is not None:
        output_file = args.output
    else:
        output_file = os.path.join(args.outdir, "tissue_positions.csv")
    # run
    spatial_location_from_barcode(args.cluster, output_file)


def spatial_location_from_barcode(cluster_file, output_file):
    """
    Get the spatial location of the spot from spot barcodes (especially for Stereo-seq data)
    for example, the spot barcode is 'x_y', then extract the spatial location x and y

    Inputs:
        cluster_file - the file of the clusters
        output_file - the file saving the results
    """
    # read the clsuter file
    df = pd.read_csv(cluster_file, sep="\t", header=None, names=['barcode', 'cluster'], na_values=[])
    df['cluster'] = df['cluster'].apply(lambda x: str(int(x)) if isinstance(x, float) and x.is_integer() 
                                                        else str(x) if pd.notnull(x) else "NA")

    # Process the columns to create new required columns
    df['in_tissue'] = 1
    df['array_row'] = df['barcode'].apply(lambda x: x.split('_')[0])
    df['array_col'] = df['barcode'].apply(lambda x: x.split('_')[1])
    df['pxl_row_in_fullres'] = df['array_row']
    df['pxl_col_in_fullres'] = df['array_col']

    # Select and reorder columns as specified
    df = df[['barcode', 'in_tissue', 'array_row', 'array_col', 'pxl_row_in_fullres', 'pxl_col_in_fullres']]

    # Write to a new file
    df.to_csv(output_file, index=False, sep=",")


def plot_scatter_spatial(df, output_dir, sample_name, column="mut_prob", note="prob", fig_size=5, \
                            xmin=0, xmax=128, ymin=0, ymax=128, point_size=6, margin_padding=0.05):
    """
    scatter plot for the mutation probability (or mosaic likelihood) and mutant allele appearance
    on the original scale with cluster
    """
    plt.figure(figsize=(fig_size+2.5, fig_size+0.5))
    # sns.scatterplot(data=df, x='pos_x', y='pos_y', hue='mutated', s=6, edgecolor='none')
    # plt.legend(title="Mutated")

    # set cluster color
    # plot_color=["#F56867","#FEB915","#C798EE","#59BE86","#7495D3","#D1D1D1","#6D1A9C","#15821E","#3A84E6","#997273","#787878","#DB4C6C","#9E7A7A","#554236","#AF5F3C","#93796C","#F9BD3F","#DAB370","#877F6C","#268785"]
    plot_color = cc.glasbey_light
    # cluster_num = int(max(df['cluster']))
    cluster_num = len(df['cluster'].unique())
    # cluster_num = df['cluster'].nunique()
    cluster_colors = dict(zip(df['cluster'].unique(), plot_color[0:cluster_num]))
    # edge_colors = df['cluster'].map(cluster_colors)
    patches = [Patch(color=color, label=cluster) for cluster, color in cluster_colors.items()]

    # plot polygons around sub-clusters within each main cluster
    def are_points_collinear(points):
        if len(points) < 3:
            return True
        A = points[1:] - points[:-1]
        return np.linalg.matrix_rank(A) < 2

    for cluster, color in cluster_colors.items():
        points = df[df['cluster'] == cluster][['pos_y', 'pos_x']].to_numpy()
        if len(points) >= 3:  # need at least 3 points for convex hull
            # apply DBSCAN to identify sub-clusters
            dbscan = DBSCAN(eps=1.5, min_samples=2)  # adjust parameters as needed
            sub_clusters = dbscan.fit_predict(points)
            # print(f"Cluster {cluster}, DBSCAN labels: {set(sub_clusters)}")  # Print the unique labels
            for sub_cluster in set(sub_clusters):
                if sub_cluster != -1:  # -1 indicates noise points in DBSCAN
                    sub_points = points[sub_clusters == sub_cluster]
                    if len(sub_points) >= 3 and not are_points_collinear(sub_points):
                        hull = ConvexHull(sub_points)
                        # polygon = Polygon(sub_points[hull.vertices], edgecolor=color, fill=None)
                        polygon = Polygon(sub_points[hull.vertices], edgecolor=None, facecolor=color, alpha=0.1)
                        plt.gca().add_patch(polygon)

    # define a color map for mutant probability
    cmap_name = 'mut_prob_color'
    colors = [(0, '#1E466E'), (0.1, '#376795'), (0.2, '#528FAD'), (0.3, '#72BCD5'), (0.4, '#AADCE0'), \
              (0.5, '#FFE6B7'), (0.8, '#FFD06F'), (0.9, '#F7AA58'), (0.95, '#EF8A47'), (1, '#E76254')]
    cmap = LinearSegmentedColormap.from_list(cmap_name, colors, N=100)
    # set the linewidth according to the point size
    linewidth = 0.08 * point_size
    # scatter plot (or s=6, linewidth=0.5)
    sc = plt.scatter(df['pos_y'], df['pos_x'], c=df[column], cmap=cmap, s=point_size, marker='o', vmin=0, vmax=1, \
                        edgecolors=['#E76254' if vaf > 0 else '#1E466E' for vaf in df['vaf']], linewidth=linewidth)
    
    # control the limits and margins for the figure
    margins_y = margin_padding * (ymax-ymin)
    margins_x = margin_padding * (xmax-xmin)
    plt.xlim(ymin-margins_y, ymax+margins_y)
    plt.ylim(xmin-margins_x, xmax+margins_x)
    # plt.margins(x=margins_x, y=margins_y)
    # add a colorbar with specific ticks
    cb = plt.colorbar(sc)
    legend_label = 'Mosaic Mutation Probability' if column=="mut_prob" else 'Mosaic Likelihood'
    cb.set_label(legend_label)
    cb.set_ticks([0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1])
    cb.set_ticklabels(['0', '0.1', '0.2', '0.3', '0.4', '0.5', '0.6', '0.7', '0.8', '0.9', '1'])

    plt.title(sample_name+' '+note+' scatter plot')
    plt.xlabel('y-axis')
    plt.ylabel('x-axis')
    # plt.axis('equal')
    plt.gca().invert_yaxis()
    # add legend
    # plt.legend(handles=patches, loc='upper right', bbox_to_anchor=(1.4, 1))
    legend1 = plt.legend(handles=patches, loc='upper right', bbox_to_anchor=(1.4, 1))
    plt.gca().add_artist(legend1)
    # legend_elements = [
        # Line2D([0], [0], marker='o', color='w', label='vaf>0', markerfacecolor='#FFE6B7', markersize=6.5, markeredgewidth=1, markeredgecolor='#E76254'),
        # Line2D([0], [0], marker='o', color='w', label='vaf=0', markerfacecolor='#FFE6B7', markersize=6.5, markeredgewidth=1, markeredgecolor='#1E466E')
    # ]
    legend_elements = [
        Line2D([0], [0], marker='o', color='w', label='vaf>0', markerfacecolor='none', markersize=6.5, markeredgewidth=1, markeredgecolor='#E76254'),
        Line2D([0], [0], marker='o', color='w', label='vaf=0', markerfacecolor='none', markersize=6.5, markeredgewidth=1, markeredgecolor='#1E466E')
    ]
    plt.legend(handles=legend_elements, loc='lower right', bbox_to_anchor=(1.4, 1))
    plt.tight_layout()
    output_scatter = os.path.join(output_dir, sample_name+"_"+note+"_scatter.png")
    plt.savefig(output_scatter)
    plt.close()


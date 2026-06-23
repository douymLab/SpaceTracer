import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from scipy import stats
from scipy.stats import ks_2samp, wilcoxon
import seaborn as sns # cdf
# for moran's I index
import geopandas as gpd
from pysal.lib import weights  # Spatial weights
#import libpysal as lps
from esda.moran import Moran
# from splot.esda import plot_moran
# pca and standardize
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from scipy.spatial import ConvexHull
from sklearn.cluster import DBSCAN
from matplotlib.patches import Polygon, Patch
from matplotlib.lines import Line2D
import colorcet as cc


def get_spatial_features(input_file, output_dir, sample_name, alpha=0.05, thr_r2=0.5, \
                         thr_prob=0.5, thr_likelihood=0.5, thr_vaf=0, plot_supp=False, fig_size=5, \
                         fig_xmin=0, fig_xmax=128, fig_ymin=0, fig_ymax=128, point_size=6, \
                         method="LDA", num_directions=8, n_quantile=100):
    """
    Inputs:
        input_file - the direction and the name of the input file, which is the spot mutation info
        output_dir - the direction saving the output figures
        sample_name - the sample name of the output files
        alpha - the significance level for the tests (default=0.05)
        thr_prob - the threshold of mosaic mutation probability (default=0.9)
            the spots with higher probability than the threshold are recognized as having mosaic mutations
        thr_likelihood - the threshold of mosaic likelihood (default=0.9)
            the spots with higher mosaic likelihood than the threshold are recognized as having mosaic mutations
        thr_vaf - the threshold of mutant allele frequency (default=0)
            the spots with higher mutant AF than the threshold are recognized as having mosaic mutations
        thr_r2 - the threshold of the R-square used to distinguish early-embryonic somatic mutations (default=0.5)
        plot_supp - Boolean variable of whether plotting the supplementary plots (default=False), including
            the scatter plot between mutant probability versus vaf
            the bar plot for the cluster mutated spots count
            the standardized density plot of the depth and mutant allele number
            the QQ-plot between the depth and the mutant allele number
            the linear regression plot between the depth and the mutant allele number
            when the test result is significant
            the CDF plots for each KS-test and the Moran's I plot for each Moran's I test (if uncomment)
        fig_size - the height and the width of the square plot, the height of the rectangle plot (default=5)
        fig_xmin, fig_xmax, fig_ymin, fig_ymax - the range of the spatial scatter plot figure (default=0,128,0,128)
        point_size - the point size used in the spatial scatter plot (default=6)
        method - the method used to do the dimension-reduction on which we perform the KS-test (default="LDA"):
            LDA: use Linear Discriminant Analysis to find the most-different direction
            PCA: apply the PCA method to find the most-variant directions
            OPT: test number of directions to find the smallest KS-test p-value by optimization
        num_directions - the number of directions to test in the optimization method (OPT) (default=8)
        n_quantile - the number of quantiles plotted in the QQ-plot (default=100)
        
    Outputs:[int(test_sig), int(early_embryonic), int(late_mutation), \
            min_ks_stat, min_ks_pval, min_moran_stat, min_moran_pval, \
            mut_rate, mut_rate_prob, mut_rate_likelihood, mut_rate_vaf, \
            all_spots_vaf_mean, all_spots_vaf_max, mut_spots_vaf_mean, mut_spots_vaf_median, \
            num_spots, num_mut_spots, r_squared, wilcoxon_stat, wilcoxon_pval, wilcoxon_rbc]
        test_sig - True if pass the spatial test (including the early somatic mutation and the late, very late somatic mutation).
            i.e. Either for mutation probability or mosaic likelihood or mutant allele frequency, one of the KS test or global Moran's I test is passed
        early_embryonic, late_mutation - True if the site is identified as early embryonic or late somatic mutation, else False
        min_ks_stat, min_ks_pval - the statistic and the minimum p-value for the KS-test whether on the PC1/PC2 direction or 
                                      the most variant direction or the candidate direction lists
        min_moran_stat, min_moran_pval - the statistic and the p-value for the moran I's index versus the normal distribution
        mut_rate, mut_rate_prob, mut_rate_likelihood, mut_rate_vaf - the proportion of the spots that contain candidate mutations with
                                                                     high mutant probability or high mosaic likelihood or high vaf
        all_spots_vaf_mean, all_spots_vaf_max - the average / maximum spot mutant allele frequency
        mut_spots_vaf_mean, mut_spots_vaf_median - the average / median spot mutant allele frequency for spots with alternative alleles
        num_spots, num_mut_spots - the total number of spots and the number of spots with alternative alleles
        r_squared - the R-squared for the linear regression between the standardized depth and mutant allele number
        wilcoxon_stat, wilcoxon_pval, wilcoxon_rbc - the statistic, p-value and the rbc value of the Wilcoxon test 
                                                     between the standardized depth and mutant allele number
    """
    # read data and remove the rows with NA
    input = open(input_file,"r")
    df = pd.read_csv(input, sep="\t", header=None, names=['barcode_name', 'cluster', 'in_tissue', 'pos_x', 'pos_y', 'depth', 'vaf', 'mut_prob', 'l_mosaic'])
    df = df.dropna()
    
    # make sure the cluster column are factors
    df['cluster'] = df['cluster'].astype('category') if df['cluster'].dtype=='object' else df['cluster'].astype('int').astype('category')
    # identify if the spot is potential mutated
    df['mutant_allele_num'] = (np.ceil(df['depth'] * df['vaf'])).astype(int)
    df['has_mutant_allele'] = (df['vaf']>thr_vaf).astype(int)
    df['high_prob'] = (df['mut_prob']>thr_prob).astype(int)
    df['high_likelihood'] = (df['l_mosaic']>thr_likelihood).astype(int)
    df['mutated'] = df['has_mutant_allele'] & df['high_prob']
    # output the dataframe
    input.close()
    df.to_csv(input_file, sep='\t', index=False, header=False)

    # calculate mutation rate
    mut_rate = len(df[df['mutated']==1]) / len(df) if len(df)>0 else "NA"
    mut_rate_prob = len(df[df['high_prob']==1]) / len(df) if len(df)>0 else "NA"
    mut_rate_likelihood = len(df[df['high_likelihood']==1]) / len(df) if len(df)>0 else "NA"
    mut_rate_vaf = len(df[df['has_mutant_allele']==1]) / len(df) if len(df)>0 else "NA"
    # mut_rate_ratio = len(df[df['high_prob']==1]) / len(df[df['has_mutant_allele']==1]) if len(df)>0 else "NA"
    # mut_rate_ratio = len(df[df['high_likelihood']==1]) / len(df[df['has_mutant_allele']==1]) if len(df)>0 else "NA"

    # other features
    all_spots_vaf_max = df['vaf'].max()
    all_spots_vaf_mean = df['vaf'].mean()
    mut_spots_vaf_mean = df.loc[df["vaf"] > 0, "vaf"].mean()
    mut_spots_vaf_median = df.loc[df["vaf"] > 0, "vaf"].median()
    num_spots = len(df)
    num_mut_spots = (df["vaf"] > 0).sum()
    
    # identify output for NA sites  
    r_squared = "NA"
    wilcoxon_stat = "NA"
    wilcoxon_pval = "NA"
    wilcoxon_rbc = "NA"
    min_moran_stat = "NA"
    min_moran_pval = "NA"
    output_null = [0]*3 + ["NA"]*4 + [mut_rate, mut_rate_prob, mut_rate_likelihood, mut_rate_vaf, \
                                      all_spots_vaf_mean, all_spots_vaf_max, \
                                      mut_spots_vaf_mean, mut_spots_vaf_median, num_spots, num_mut_spots, \
                                      r_squared, wilcoxon_stat, wilcoxon_pval, wilcoxon_rbc]

    # standardize the two variables
    depth_data = pd.DataFrame({
        'depth': df['depth'],
        'mutant_allele_number': df['mutant_allele_num']
    })
    # consider the all NA depth data situation
    if depth_data.empty:
        return output_null
    standard_scaler = StandardScaler()
    depth_std = standard_scaler.fit_transform(depth_data)
    depth_std = pd.DataFrame(depth_std, columns=['depth', 'mutant_allele_num'])
    diff = depth_std['depth'] - depth_std['mutant_allele_num']
    # calculate R-quare and paired wilcoxon test for the linear regression between depth and mutant allele number
    if len(depth_std['depth'].unique()) > 1:
        _, _, r_value, _, _ = stats.linregress(depth_std['depth'], depth_std['mutant_allele_num'])
        r_squared = r_value**2
    else:
        r_squared = 1
    if len(diff.unique()) > 1:
        wilcoxon_stat, wilcoxon_pval = wilcoxon(depth_std['depth'], depth_std['mutant_allele_num'])
        wilcoxon_rbc = calculate_rbc_for_paired_wilcoxon(depth_std['depth'], depth_std['mutant_allele_num'])
    else:
        wilcoxon_stat = 0
        wilcoxon_pval = 1
        wilcoxon_rbc = 0

    # check if the site caused by allele dropout
    no_alleledrop = (r_squared < thr_r2) or (wilcoxon_pval < alpha)

    # if less than 2 spots has no NA mutated value
    # then the mutated spots are uniformally distributed so not pass the tests
    if len(df) < 2:
        return output_null
    # if all spots are mutated or not mutated
    # then the mutated spots are uniformally distributed so not pass the tests
    if (len(df[df['mutated']==1])==0 or len(df[df['mutated']==0])==0):
        return output_null
    # if only one mutated spot, then the site is a very late mutation
    elif len(df[df['mutated']==1])==1:
        output_verylate = [1,0,1] + ["NA"]*4 + [mut_rate, mut_rate_prob, mut_rate_likelihood, mut_rate_vaf, \
                                                all_spots_vaf_mean, all_spots_vaf_max, \
                                                mut_spots_vaf_mean, mut_spots_vaf_median, num_spots, num_mut_spots, \
                                                r_squared, wilcoxon_stat, wilcoxon_pval, wilcoxon_rbc]
        return output_verylate
    

    # =============================================================================
    # Global spatial tests for early mutations
    # =============================================================================
    # spatial tests for mutation probability from posterior (and have mutant alleles)
    early_mut = False
    if df['mutated'].isnull().any():
        print("There are missing values at the site: ", sample_name)
    elif np.isinf(df['mutated']).any():
        print("There are infinite values in the data at the site: ", sample_name)
    else:
        # try:
        _, ks_pass, min_ks_stat, min_ks_pval = spatial_kstest(df=df, col='mutated', alpha=alpha, method=method, num_directions=num_directions, \
                                                plot_supp=plot_supp, output_dir=output_dir, sample_name=sample_name, fig_size=fig_size, \
                                                note="")
        moran_pass, min_moran_stat, min_moran_pval = spatial_moran(df=df, col='mutated', alpha=alpha, \
                                                plot_supp=plot_supp, output_dir=output_dir, sample_name=sample_name, fig_size=fig_size, \
                                                note="")
        early_mut = True if (ks_pass or moran_pass) else False
        # except ZeroDivisionError as e:
        #     print("Division by zero error:", e)
        #     print("At the site: ", sample_name)
        # except Exception as e:
        #     print("An error occurred:", e)
        #     print("At the site: ", sample_name)


    # =============================================================================
    # Cluster-level spatial features
    # =============================================================================
    # group cluster info
    df_cluster = df.groupby(['cluster'], observed=True).agg({
        'barcode_name': 'count',
        'depth': 'sum',
        'mutant_allele_num': 'sum',
        'has_mutant_allele': 'sum',
        'high_prob': 'sum',
        'high_likelihood': 'sum',
        'mutated': 'sum',
        'vaf': 'mean'
    }).reset_index()
    df_cluster.rename(columns={'barcode_name': 'spot_count'}, inplace=True)
    # calculate cluster mutated spot proportion
    df_cluster['vaf_unnom'] = (df_cluster['mutant_allele_num'] / df_cluster['depth']).replace([np.inf, -np.inf], np.nan)
    df_cluster['mutant_prop'] = (df_cluster['mutated'] / df_cluster['spot_count']).replace([np.inf, -np.inf], np.nan)
    df_cluster['mutant_prop_prob'] = (df_cluster['high_prob'] / df_cluster['spot_count']).replace([np.inf, -np.inf], np.nan)
    df_cluster['mutant_prop_likelihood'] = (df_cluster['high_likelihood'] / df_cluster['spot_count']).replace([np.inf, -np.inf], np.nan)
    df_cluster['mutant_prop_vaf'] = (df_cluster['has_mutant_allele'] / df_cluster['spot_count']).replace([np.inf, -np.inf], np.nan)

    
    # =============================================================================
    # Conclusion and Plot
    # =============================================================================
    # identify all mutation types
    # test_sig = (early_mut or late_mut or verylate_mut) and no_alleledrop
    test_sig = early_mut and no_alleledrop

    # distinguish early embryonic mutation
    early_embryonic = False
    late_mutation = False
    if test_sig and r_squared > thr_r2:
        early_embryonic = True
    elif test_sig and r_squared <= thr_r2:
        late_mutation = True

    # scatter plot on the original scale with cluster
    # if plot:
    #     plot_scatter_spatial(df=df, output_dir=output_dir, sample_name=sample_name, column="mut_prob", note="prob", \
    #                          fig_size=fig_size, xmin=fig_xmin, xmax=fig_xmax, ymin=fig_ymin, ymax=fig_ymax, point_size=point_size)
    #     plot_scatter_spatial(df=df, output_dir=output_dir, sample_name=sample_name, column="l_mosaic", note="likelihood", \
    #                          fig_size=fig_size, xmin=fig_xmin, xmax=fig_xmax, ymin=fig_ymin, ymax=fig_ymax, point_size=point_size)
    if plot_supp:
        plot_scatter_spatial(df=df, output_dir=output_dir, sample_name=sample_name, column="mut_prob", note="prob", \
                             fig_size=fig_size, xmin=fig_xmin, xmax=fig_xmax, ymin=fig_ymin, ymax=fig_ymax, point_size=point_size)
        plot_scatter_spatial(df=df, output_dir=output_dir, sample_name=sample_name, column="l_mosaic", note="likelihood", \
                            fig_size=fig_size, xmin=fig_xmin, xmax=fig_xmax, ymin=fig_ymin, ymax=fig_ymax, point_size=point_size)

        plot_supp_spatial(df=df, df_cluster=df_cluster, output_dir=output_dir, sample_name=sample_name)
        plot_supp_depth(depth_std=depth_std, output_dir=output_dir, sample_name=sample_name, fig_size=fig_size)
        plot_qq_depth(depth_std=depth_std, output_dir=output_dir, sample_name=sample_name, num=n_quantile)
        plot_reg_depth(depth_std=depth_std, r_squared=r_squared, output_dir=output_dir, sample_name=sample_name, fig_size=fig_size)
        plot_box_depth(depth_std=depth_std, wilcoxon_pval=wilcoxon_pval, output_dir=output_dir, sample_name=sample_name, fig_size=fig_size)

    # outputs
    output = [int(test_sig), int(early_embryonic), int(late_mutation), \
            min_ks_stat, min_ks_pval, min_moran_stat, min_moran_pval, \
            mut_rate, mut_rate_prob, mut_rate_likelihood, mut_rate_vaf, \
            all_spots_vaf_mean, all_spots_vaf_max, mut_spots_vaf_mean, mut_spots_vaf_median, \
            num_spots, num_mut_spots, r_squared, wilcoxon_stat, wilcoxon_pval, wilcoxon_rbc]
    return output



def spatial_kstest(df, col="mutated", alpha=0.05, method="LDA", num_directions=8, \
                   plot_supp=False, output_dir="./", sample_name=None, fig_size=5, note=None):
    """
    Perform the Kolmogorov-Smirnov (KS) test to find whether the points are distributed uniformly 
    on the projected directions for the spots in one cluster.

    df - the dataframe we analyze with
    col - the column with Boolean variable
    """
    # get position matrix
    X = np.array([df['pos_x'], df['pos_y']]).transpose()
    # =============================================================================
    # KS-test with PCA
    # =============================================================================
    if method=="PCA":
        # normalise
        sc = StandardScaler()
        X_norm = sc.fit_transform(X)
        # pca
        pca = PCA(n_components = 2)
        X_pca = pca.fit_transform(X_norm)
        
        # find the positions of the spots with/without mutantations
        pc1_1 = X_pca[:,0][df[col]==1]
        pc1_0 = X_pca[:,0][df[col]==0]
        pc2_1 = X_pca[:,1][df[col]==1]
        pc2_0 = X_pca[:,1][df[col]==0]

        # perform Kolmogorov-Smirnov test
        ks_pc1 = ks_2samp(pc1_0, pc1_1)
        ks_pc2 = ks_2samp(pc2_0, pc2_1)
        # find the test with the minimal p-value
        min_test = min([ks_pc1, ks_pc2], key=lambda x: x.pvalue)
        # extract the minimum p-value and the corresponding statistic
        min_ks_pval = min_test.pvalue
        min_ks_stat = min_test.statistic
        ks_pass = True if (min_ks_pval<alpha) else False
        # X_proj = X_pca[:,0] if ks_pc1.pvalue<ks_pc2.pvalue else X_pca[:,1]
        # df['proj'] = X_pca[:,0] if ks_pc1.pvalue<ks_pc2.pvalue else X_pca[:,1]

    # =============================================================================
    # KS-test among several directions
    # =============================================================================
    if method=="OPT":
        # function to perform KS test for a given direction
        def perform_ks(direction, X, mutated):
            # project points onto a given direction
            projections = np.dot(X, direction)
            projection_0 = projections[mutated == 0]
            projection_1 = projections[mutated == 1]
            ks = ks_2samp(projection_0, projection_1)
            return ks

        # generate directions to sample
        thetas = np.linspace(0, np.pi, num_directions, endpoint=False)
        directions = np.array([np.cos(thetas), np.sin(thetas)]).T

        # compute the KS p-values for each direction
        ks_lists = [perform_ks(d, X, df[col]) for d in directions]
        # find the test with the minimal p-value
        min_test = min(ks_lists, key=lambda x: x.pvalue)
        # extract the minimum p-value and the corresponding statistic
        min_ks_pval = min_test.pvalue
        min_ks_stat = min_test.statistic

        ks_pass = True if min_ks_pval<alpha else False

        # # find the most-different direction
        # index = list(ks_lists).index(min_test)
        # diff_direction = directions[index]
        # X_proj = np.dot(X, diff_direction)
        # proj_1 = X_proj[df[col]==1]
        # proj_0 = X_proj[df[col]==0]
        # # df['proj'] = X_proj

    # =============================================================================
    # KS-test with LDA
    # =============================================================================
    if method=="LDA":
        # initialize the LDA model
        lda = LinearDiscriminantAnalysis(n_components=1)
        # apply LDA for dimensionality reduction
        X_proj = lda.fit_transform(X, df[col])
        # df['proj'] = X_proj

        # if the mutant and non-mutant spots has the same mean, return False
        if X_proj.size == 0 or X_proj.shape[1] == 0:
            return df, False, "NA", "NA"

        # find the positions of the spots with/without mutantations
        proj_1 = X_proj[df[col]==1]
        proj_0 = X_proj[df[col]==0]

        # perform Kolmogorov-Smirnov test
        proj_0 = proj_0.ravel()
        proj_1 = proj_1.ravel()
        ks_proj = ks_2samp(proj_0, proj_1)
        min_ks_pval = ks_proj.pvalue
        min_ks_stat = ks_proj.statistic
        ks_pass = True if (min_ks_pval<alpha) else False

    # =============================================================================
    # Plot CDFs
    # =============================================================================
    if plot_supp:
        if method=="PCA":
            plt.figure(figsize=(2*fig_size, fig_size))
            sns.ecdfplot(pc1_1, color='#ff7f0e', label='AF=1')
            sns.ecdfplot(pc1_0, color='#1f77b4', label='AF=0')
            plt.title(sample_name+note+' CDF on PC1')
            plt.legend()
            plt.xlabel("PC1")
            output_cdf1 = os.path.join(output_dir, sample_name+note+"_cdf1.png")
            plt.savefig(output_cdf1)
            plt.close()

            plt.figure(figsize=(2*fig_size, fig_size))
            sns.ecdfplot(pc1_1, color='#ff7f0e', label='AF=1')
            sns.ecdfplot(pc2_0, color='#1f77b4', label='AF=0')
            plt.title(sample_name+note+' CDF on PC2')
            plt.legend()
            plt.xlabel("PC2")
            output_cdf2 = os.path.join(output_dir, sample_name+note+"_cdf2.png")
            plt.savefig(output_cdf2)
            plt.close()
        else:
            plt.figure(figsize=(2*fig_size, fig_size))
            sns.ecdfplot(proj_1, color='#ff7f0e', label='AF=1')
            sns.ecdfplot(proj_0, color='#1f77b4', label='AF=0')
            plt.title(sample_name+note+' CDF on most-diff direction')
            plt.legend()
            plt.xlabel("direction")
            output_cdf = os.path.join(output_dir, sample_name+note+"_cdf.png")
            plt.savefig(output_cdf)
            plt.close()

    return df, ks_pass, min_ks_stat, min_ks_pval


def spatial_moran(df, col="mutated", alpha=0.05, \
                  plot_supp=False, output_dir="./", sample_name=None, fig_size=5, note=None):
    """
    Perform Moran's I index to find whether the spots with mutation are locally or globally distributed
    for one cluster
    """
    # keep only valid rows for spatial test and rebuild a clean index
    # so libpysal does not fail on sparse/misaligned indices in multiprocessing.
    moran_df = df[['pos_x', 'pos_y', col]].dropna().copy()
    if len(moran_df) < 2 or moran_df[col].nunique() <= 1:
        return False, np.nan, 1.0

    moran_df = moran_df.reset_index(drop=True)
    gdf = gpd.GeoDataFrame(moran_df, geometry=gpd.points_from_xy(moran_df['pos_x'], moran_df['pos_y']))

    try:
        # use a rebuilt 0..n-1 index to avoid KeyError from neighbor/cardinality lookup
        w = weights.Queen.from_dataframe(gdf, use_index=False)
        mi = Moran(moran_df[col], w)
    except Exception as queen_err:
        try:
            # fallback for degenerate geometries / disconnected queen graph
            k = min(4, len(moran_df) - 1)
            w = weights.KNN.from_dataframe(gdf, k=k)
            mi = Moran(moran_df[col], w)
            print(f"[WARN] Queen weights failed for {sample_name}; fallback to KNN (k={k}). Error: {queen_err}")
        except Exception as knn_err:
            print(f"[WARN] Moran test skipped for {sample_name}; Queen and KNN failed. "
                  f"Queen error: {queen_err}; KNN error: {knn_err}")
            return False, np.nan, 1.0
    # verify Moran's I results
    moran_pval = mi.p_rand
    moran_stat = mi.I
    
    moran_pass = True if moran_pval < alpha else False
    # print(mi.I) 
    # print(mi.p_rand)   # p-value for random distributed
    # print(mi.p_norm)   # p-value for normal distributed

    # # plot Moran's I test figure
    # if plot_supp:
    #     plot_moran(mi, zstandard=True, figsize=(2*fig_size,fig_size))
    #     output_moran = os.path.join(output_dir, sample_name+note+"_moran.png")
    #     plt.title(sample_name+note+" moran plot")
    #     plt.savefig(output_moran)
    #     plt.close()
    return moran_pass, moran_stat, moran_pval


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


def plot_supp_spatial(df, df_cluster, output_dir, sample_name):
    """
    Supplementary Plots for spatial test, including
    Scatter plot for vaf versus mutant probability
    Bar plot for mutated spot fraction for each cluster
    """
    # define a color map for mutant probability
    cmap_name = 'mut_prob_color'
    colors = [(0, '#1E466E'), (0.1, '#376795'), (0.2, '#528FAD'), (0.3, '#72BCD5'), (0.4, '#AADCE0'), \
              (0.5, '#FFE6B7'), (0.8, '#FFD06F'), (0.9, '#F7AA58'), (0.95, '#EF8A47'), (1, '#E76254')]
    cmap = LinearSegmentedColormap.from_list(cmap_name, colors, N=100)

    # vaf versus mutant probablity
    plt.scatter(df['vaf'], df['mut_prob'], c=df['mut_prob'], cmap=cmap)
    # plt.xscale('log')
    plt.title(sample_name + ' vaf vs mut_ptob')
    plt.xlabel('Mutant allele frequency')
    plt.ylabel('Mosaic mutation probability')
    plt.colorbar().set_label('mutant probability')
    output_vaf_mutprob = os.path.join(output_dir, sample_name+"_vaf_mutprob.png")
    plt.savefig(output_vaf_mutprob)
    plt.close()

    # mutated spot fraction for each cluster
    ax = df_cluster.plot(x='cluster', y=['mutant_prop_vaf', 'mutant_prop_prob'], kind='bar', color=['#1f77b4', '#ff7f0e'])
    # # adding the text labels above the bars
    # for container in ax.containers:
    #     ax.bar_label(container, fmt='%.2f')
    # custom legend labels
    legend_labels = ['mutant allele', 'mutation probability']
    handles = [plt.Rectangle((0,0),1,1, color=color) for color in ['#1f77b4', '#ff7f0e']]
    ax.legend(handles, legend_labels, title='Proportion Method')

    plt.xlabel('Cluster')
    plt.ylabel('Mutated spot fraction')
    plt.title(sample_name + ' mutated spot fraction')
    # plt.legend(title='Proportion')
    plt.xticks(rotation=0)  # Rotate x-axis labels to show them horizontally
    plt.tight_layout()  # Adjust layout to prevent clipping of labels
    output_mutprop = os.path.join(output_dir, sample_name+"_mutated_prop.png")
    plt.savefig(output_mutprop)
    plt.close()



def plot_supp_depth(depth_std, output_dir, sample_name, fig_size=5):
    """
    Plot the normalized density for depth and mutant allele number
    """
    plt.figure(figsize=(fig_size+3, fig_size))
    sns.kdeplot(depth_std['depth'], color="skyblue", label='depth (standardized)', fill=True)
    sns.kdeplot(depth_std['mutant_allele_num'], color="red", label='mutant allele number (standardized)', fill=True)
    plt.legend()
    plt.xlabel('Value')
    plt.ylabel('Density')
    plt.title(sample_name + ' standardized density plot')
    output_depth_density = os.path.join(output_dir, sample_name+"_depth_density.png")
    plt.savefig(output_depth_density)
    plt.close()


def plot_qq_depth(depth_std, output_dir, sample_name, num=100):
    """
    Plot the QQ-plot for the linear regression between normalised depth and mutant allele number
    """
    # get the values
    depth = depth_std['depth']
    mutant_num = depth_std['mutant_allele_num']
    # calculate quantiles
    quantiles = np.linspace(0, 1, num=num)  # adjust the quantile number
    q_depth = np.quantile(depth, quantiles)
    q_mutant = np.quantile(mutant_num, quantiles)
    # plot QQ-plot
    plt.figure(figsize=(6, 6))
    plt.plot(q_depth, q_mutant, marker='o', linestyle='', markersize=5)
    plt.plot([min(q_depth.min(), q_mutant.min()), max(q_depth.max(), q_mutant.max())], 
            [min(q_depth.min(), q_mutant.min()), max(q_depth.max(), q_mutant.max())], 
            color='r', linestyle='--')  # line y=x for reference
    plt.xlabel('Quantiles of depth')
    plt.ylabel('Quantiles of mutant allele number')
    plt.title(sample_name + ' QQ-plot')
    plt.grid(True)
    output_qqplot = os.path.join(output_dir, sample_name+"_qqplot.png")
    plt.savefig(output_qqplot)
    plt.close()


def plot_reg_depth(depth_std, r_squared, output_dir, sample_name, fig_size=5):
    """
    Plot the linear regression between the depth and the mutant allele number
    """
    plt.figure(figsize=(fig_size+3, fig_size))
    sns.regplot(x='depth', y='mutant_allele_num', data=depth_std, ci=95, line_kws={'label': f'$R^2$ = {r_squared:.2f}'})
    plt.legend()
    plt.xlabel('Depth (Standardized)')
    plt.ylabel('Mutant Allele Number (Standardized)')
    plt.title(sample_name + ' linear regression plot')
    output_depth_reg = os.path.join(output_dir, sample_name+"_depth_reg.png")
    plt.savefig(output_depth_reg)
    plt.close()


def plot_box_depth(depth_std, wilcoxon_pval, output_dir, sample_name, fig_size=5):
    """
    Plot the box plot for the depth and the mutant allele number
    """
    # prepare the data for boxplot
    box_df = pd.DataFrame({'Depth': depth_std['depth'], 'Mutant Allele Number':  depth_std['mutant_allele_num']})

    # plotting the boxplot with scatter points and lines connecting paired points
    plt.figure(figsize=(fig_size+3, fig_size))
    sns.boxplot(data=box_df)
    sns.stripplot(data=box_df, color='black', jitter=0, size=5)

    # adding lines between paired points
    for i in range(len(depth_std['depth'])):
        plt.plot([0, 1], [depth_std['depth'][i],  depth_std['mutant_allele_num'][i]], color='gray', linestyle='-', linewidth=0.5)

    # adding Wilcoxon p-value
    p_text = f'p-value: {wilcoxon_pval:.4f}'
    # plt.gcf().text(0.15, 0.95, p_text, fontsize=12, ha='left', va='top', 
    #                bbox=dict(facecolor='white', alpha=0.5))
    plt.text(0.5, max(max(depth_std['depth']), max(depth_std['mutant_allele_num'])) * 1, p_text, fontsize=12, ha='center')
    plt.title(sample_name + ' box plot')
    plt.ylabel('Normalized Values')
    output_depth_box = os.path.join(output_dir, sample_name+"_depth_box.png")
    plt.savefig(output_depth_box)
    plt.close()


def calculate_rbc_for_paired_wilcoxon(x, y):
    """ 
    Calculate the Rank-Biserial Correlation (RBC) value for the paired Wilcoxon signed-rank test.
    Parameters
    ----------
    x,y : array_like
        The data from the two samples.
    Returns
    -------
    rbc: float
        The Rank-Biserial Correlation (RBC) value, which is a measure of
        the strength and direction of the association between the two samples.
        It ranges from -1 to 1, where:
        - 1 indicates a perfect positive association,
        - -1 indicates a perfect negative association,
        - 0 indicates no association. 
    """
    # format as numpy arrays
    x, y = map(np.asarray, (x, y))
    # test whether the two arrays are of the same length
    if len(x) != len(y):
        raise ValueError("Input arrays must have the same length.")
    try:
        X = np.array([x, y]).T
        X_scaled = StandardScaler().fit_transform(X)
        diff = X_scaled[:, 0] - X_scaled[:, 1]
        non_zero = diff != 0
        if not np.any(non_zero):
            return 0
        abs_diff = np.abs(diff[non_zero])
        signs = np.sign(diff[non_zero])
        ranks = stats.rankdata(abs_diff)
        signed_ranks = ranks * signs
        W_plus = signed_ranks[signed_ranks > 0].sum()
        W_minus = -signed_ranks[signed_ranks < 0].sum()
        denom = W_plus + W_minus
        if denom == 0 or not np.isfinite(denom):
            return 0
        rbc = (W_plus - W_minus) / denom
        return 0 if not np.isfinite(rbc) else rbc
    except Exception:
        return 0


def handle_ind_geno_to_get_mutation_list(spot_geno_file):
    '''
    input:    
    #chrom  pos     strand      germline        mutant  cluster spot_number    consensus_read_count    l_germline      l_mosaic        max_spot_geno   G_spot_max      depth   vaf     p_mosaic
    chr12   52487210        +       A       T       ind       5      46,0,0,0        0.7874611700926878      0.21253882990731235     germline        0/0     46      0.0     0.3685326340225939
    
    output: get new allele info from the "ind_ref" and "ind_alt"
    list: [("chrX", "119811135", "T", "C,<*>"),...]
    '''
    mutation_identifier_list=[]
    f=open(spot_geno_file,"r")
    for line in f.readlines():
        if line[0]!="#" and line[0:5]!="chrom":
            s = line.strip().split()
            chrom=s[0];pos=s[1];ref=s[3];alt=s[4]
            mutation="_".join([str(chrom),str(pos),str(ref),str(alt)])
            if mutation not in mutation_identifier_list:
                mutation_identifier_list.append(mutation)

    return mutation_identifier_list


def read_spot_posterior_file(file):
    # df = pd.read_csv(file,sep='\t',header=None)
    colnames=["chr","pos","strand","germline","mutant","cluster","spot_barcode","consensus_read_count", \
              "l_germline", "l_mosaic", "max_spot_geno","G_spot_max","depth","vaf","p_mosaic"]
    df = pd.read_csv(file, sep='\t', header=None, names=colnames, comment = "#")
    # df.columns=colnames
    df['pos'] = df['pos'].astype(str)
    df['identifier'] = df['chr']+"_"+df['pos']+"_"+df['germline']+"_"+df['mutant'].astype(str)
    df = df.set_index('identifier')
    # spot_geno_file=BedTool(file)
    return df


def handle_per_line_for_sf(barcode_dir,output_dir,in_name,spot_geno_file,barcode_dict,alpha,thr_r2,thr_prob,thr_likelihood,thr_vaf, \
                    plot_supp,fig_size,method,num_directions,line):
    barcode_geno={}
    mutation_name=line
    mutation=line.split("_")
    barcode_mutinfo_file_name=os.path.join(barcode_dir,mutation_name+'.barcode.mutinfo.txt')
    barcode_mutinfo_file=open(barcode_mutinfo_file_name,"w")
    try:
        short_df = spot_geno_file.loc[[mutation_name]]
    except KeyError:
        return []

    required_cols = ["spot_barcode", "cluster", "depth", "vaf", "p_mosaic", "l_mosaic"]
    missing_cols = [c for c in required_cols if c not in short_df.columns]
    if missing_cols:
        raise KeyError(
            f"short_df missing required columns: {missing_cols}, columns={list(short_df.columns)}"
        )

    for row in short_df.itertuples(index=False):
        barcode_geno[row.spot_barcode] = (
            row.cluster,
            row.depth,
            row.vaf,
            row.p_mosaic,
            row.l_mosaic,
        )


    for barcode in barcode_dict.keys():
        if barcode not in barcode_geno:  # ← No need for .keys()
            continue
        
        cluster, depth, vaf, mut_prob, l_mosaic = barcode_geno[barcode]
        barcode_mutinfo_file.write(f'{barcode}\t{cluster}\t1\t{barcode_dict[barcode][0]}\t{barcode_dict[barcode][1]}\t{depth}\t{vaf}\t{mut_prob}\t{l_mosaic}\n')
    barcode_mutinfo_file.close()
    if in_name:
        barcode_mutinfo_file_name=os.path.join(barcode_dir,in_name+"."+mutation_name+'.barcode.mutinfo.txt')
    else:
        barcode_mutinfo_file_name=os.path.join(barcode_dir,mutation_name+'.barcode.mutinfo.txt')

    input_file=barcode_mutinfo_file_name
    if in_name:
        sample_name=in_name+"."+mutation_name
    else:
        sample_name=mutation_name
    
    try:
        res=get_spatial_features(input_file,output_dir, sample_name, alpha, thr_r2, \
                        thr_prob, thr_likelihood, thr_vaf, plot_supp, fig_size, method, num_directions)
    except:
        print("------------*******************------------input_file",input_file)
        raise
    if res[9] != "NA":
        out_line=list(mutation)+ [str(i) for i in res]
    else:
        out_line=[]

    # significant if pass at least one test
    # test_sig = res[0]

    # print("end:", mutation_name)
    return out_line


def merge_dataframes(df_list):
    if df_list:
        return pd.concat(df_list, ignore_index=False, axis=1)
    return pd.DataFrame()


def save_dataframe(df, filename):
    df.to_csv(filename, sep="\t", index=True, index_label='barcode_name', na_rep='NA')
    
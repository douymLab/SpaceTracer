import os
import random
import pandas as pd
import numpy as  np
from collections import Counter
from sklearn.preprocessing import OneHotEncoder, LabelEncoder
from sklearn.model_selection import train_test_split, cross_val_score, RandomizedSearchCV, GridSearchCV
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score, confusion_matrix, classification_report
from sklearn.preprocessing import StandardScaler
from imblearn.over_sampling import SMOTE
from hyperopt import hp, fmin, tpe, STATUS_OK, Trials
from scipy.stats import randint
from joblib import dump, load
import pickle
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
from wordcloud import WordCloud
import shap

from SpaceTracer.utils.utils import check_dir, check_file


def mutation_classification(input_file, output_dir, sample_name, model_dir="./", model_name="spatial_feature_preserved_model", random_seed=100, \
                            train=False, true_sites_file=None, validated_artifact_sites_file=None, phasable_artifact_sites_file=None, \
                            select_features=None, drop_features=None, no_spatial=False, phase_refine=False, \
                            true_to_artifact_ratio=None, encoder="label", save_models=True, \
                            plot=True, annotate_mosaic=True, annotate_outlier=False, n_features=20, \
                            save_pca=True, save_shap=False, use_lr=False, not_pred_het=True, transform_old_name=False, \
                            smote=True, tune="random_search", k_neighbors=4, sampling_strategy="auto", n_jobs=None, \
                            n_estimators=100, max_depth=None, min_samples_split=2):
    """
    Apply Random Forest and Logistic Regression method to classify somatic mutation from candidate sites

    Inputs:
        input_file - the direction and the name of the input file, which is the features file
        output_dir - the direction saving the outputs
        sample_name - the sample name of the output files
        model_dir - the directory of the saved trained models (default="./)
        model_name - the sample name of the saved trained model used (default="spatial_feature_preserved_model")
        random_seed - the random state used in the whole progress (default=100)
        train - Boolean variable shows whether training the models with input data or use exist models (default=False)
        true_sites_file - path to the file with true somatic sites in this data set which are manually checked (default=None)
        validated_artifact_sites_file - path to the file with validated artifact sites in this data set (default=None)
        phasable_artifact_sites_file - path to the file with phasable artifact sites in this data set, 
                                       if not given then use the sites with 'haplotype>3' (default=None)
        select_features - the features used to classificate the true somatic mutations, use all features if 'None' (default=None)
        drop_features - the features do not want to be used to classificate the true somatic mutations (default=None)
        no_spatial - Boolean variable of whether not using the spatial features in the model (default=False)
        phase_refine - whether doing the phase refinement (default=False)
                       True: use an additional random forest model to classificate the haplotype=3 phased sites into 
                             somatic mutation, heterozygous and artifact three groups and save the mutation sites
                       False: use the random forest model for the unphased sites with the features not related to the phasing
                              to classificate the haplotype=3 phased sites into somatic mutation and artifact groups
        true_to_artifact_ratio - the ratio of the true somatic mutations to artifacts in the training set 
                                 (default=None, not control the ratio)
        encoder - the encode method for object columns (default="label")
                  "label": use LabelEncoder
                  "onehot": use OneHotEncoder
        save_models - Boolean variable of whether saving trained models or not (default=True)
        plot - Boolean variable of whether plotting the feature importance figure and the PCA plots (default=True)
        annotate_mosaic - Boolean variable of whether annotating the mosaic sites in the training set PCA scatter plots (default=True)
        annotate_outlier - Boolean variable of whether annotating the outlier sites in all PCA scatter plots (default=False)
        n_features - the number of the most-important features from the random forest model used in the PCA projection (default=20)
        save_pca - Boolean variable of whether saving the PCA transformed value, together with the clean feature values (default=True)
        save_shap - Boolean variable of whether saving the SHAP values for the features (default=False)
        use_lr - Boolean variable of whether using logistic regression model to classify the somatic mutations (default=False)
        not_pred_het - Boolean variable of whether not predicting the heterozygous sites using logistic regression model (default=True)
        transform_old_name - Boolean variable of whether transforming the column names to new version before processing (default=False)
        
        # Following are the parameters used in the random forest and logistic regression model
        smote - whether use SMOTE (Synthetic Minority Over-sampling Technique) to over-sample the minority class 
                to treat the imbalance class values (default=True)
        tune - whether tuning the hyperparameters used in the random forest (n_estimators, max_depth, min_samples_split)
                "Bayesian_opt": Bayesian Optimization
                "random_search": Random Search
                "grid_search": Grid Search
                None: use the given parameters
                (default="random_search)
        k_neighbors - the nearest neighbors used to define the neighborhood of samples in SMOTE (default=4)
        sampling_strategy - sampling information to resample the data set (default="auto")
                'minority': resample only the minority class;
                'not minority': resample all classes but the minority class;
                'not majority': resample all classes but the majority class;
                'all': resample all classes;
                'auto': equivalent to 'not majority'.
        n_jobs - number of jobs to run in parallel (default=None) 
                None means 1 unless in a joblib.parallel_backend context. 
                -1 means using all processors.         
        n_estimators - the number of trees in the forest (default=100)
        max_depth - the maximum depth of the tree (default=None)
        min_samples_split - the minimum number of samples required to split an internal node (default=2)
    """
    # ====================================================================
    # Read Data
    # ====================================================================
    # read in feature info with warning bad rows (less or more columns)
    df = pd.read_parquet(input_file)
    df.replace({'no': None, 'NA': None, 'NaN': None}, inplace=True)
    # df = pd.read_csv(input_file, sep="\t", na_values=['no', 'NA', 'NaN'], index_col=False, on_bad_lines='warn')
    # convert tuple-like '("no",)' forms into NaN
    df.replace(r"^\('no',\)$", np.nan, regex=True, inplace=True)

    # # modify the column name to new version
    # if transform_old_name:
    #     print("Transform the column name to new version")
    #     df = rename_columns(df, old2new_mapping)

    # check index columns only appear once
    idx_cols = ["#chrom", "pos", "ref", "alt"]
    # drop duplicated genomic keys from regular columns if they are already index levels
    dup_index_cols = [c for c in idx_cols if (c in df.index.names) and (c in df.columns)]
    if dup_index_cols:
        df = df.drop(columns=dup_index_cols)

    # ensure sample is an index level
    if "sample" not in df.index.names:
        if all(c in df.index.names for c in idx_cols):
            df = (df.assign(sample=sample_name)
                    .set_index("sample", append=True)
                    .reorder_levels(["sample"] + [n for n in df.index.names if n != "sample"]))
        else:
            df = df.assign(sample=sample_name).set_index(["sample", *idx_cols])

    # keep the output features
    df_outputfeatures = df[['Filtration', 'consensus_UMI_count', 'consensus_alt_allele_count', 'AFind', 
                            'num_spots', 'num_mut_spots', 'DNAMutationType', 'RNAMutationType']].rename(columns={
                        'Filtration': 'FILTER',
                        'consensus_UMI_count': 'DP',
                        'consensus_alt_allele_count': 'ALT_DP',
                        'AFind': 'VAF',
                        'num_spots': 'NUM_SPOTS',
                        'num_mut_spots': 'NUM_MUT_SPOTS',
                        'DNAMutationType': 'DNA_MUT_TYPE',
                        'RNAMutationType': 'RNA_MUT_TYPE'})

    # ====================================================================
    # Data Format
    # ====================================================================
    # delete the rows if the likelihood columns and AF columns contain NA
    df = df.dropna(subset=['p_mosaic'])
    df = df.dropna(subset=['baseq_p_adj'])
    
    # all haplotype-related columns (keeping these as they are not in the mapping)
    haplotype_columns = ['phasing_most_phase_haplotype', 'phasing_nearest_phase_haplotype']
    # find the most frequent haplotype
    df['haplotype'] = df.apply(lambda row: merge_haplotype_columns(row, haplotype_columns), axis=1)
    # drop the original haplotype columns
    df = df.drop(haplotype_columns, axis=1)

    # set heterozygosity-related features
    het_features = ['AFind',
                    'mut_vs_nonmut_spots_KS_p',  # updated from 'KS_p'
                    'mut_spots_prop_by_vaf',      # updated from 'mut_rate_vaf'
                    'all_spots_vaf_mean',         # updated from 'mean_AFspot'
                    'mut_vs_nonmut_spots_MI_p',   # updated from 'MI_p'
                    'mut_spots_prop_by_probablity', # updated from 'mut_rate_prob'
                    'mismatches_p_adj']
    
    # ====================================================================
    # Remove the Not-used Features
    # ====================================================================
    if train:    
        # delete the gene function related columns and other non-reasonable features (updated with new column names)
        not_related_columns = ['major_read_strand', 'GCcontent', 'DNAMutationType', 'RNAMutationType', \
                               'num_mut_spots', 'num_spots', 'hFDR', 'falt', 'fref', \
                               'consensus_ref_allele_count', 'consensus_alt2_allele_count', 'consensus_alt_allele_count', \
                               'Filtration', 'editing_AtoG', 'editing_database', 'RNA_editing', \
                               'imprinted', 'ASE', 'hFDR', 'homopolymer', 'PON', 'cluster_event', \
                               '#chrom', 'pos', 'ref', 'alt', 'phasing_nearest_mut_origin', 'phasing_most_mut_origin']
        for col in not_related_columns:
            if col in df.columns:
                df = df.drop(col, axis=1)

        # drop some features which are not distinguishable for classification (updated with new column names)
        undistinguishable_columns = ['consensus_UMI_count', 'alt_read_number_perUMI_max', 'vaf', \
                                     'indel_proportion_for_site', 'mappabilityScore', 'cause_poly_alt']

        for col in undistinguishable_columns:
            if col in df.columns:
                df = df.drop(col, axis=1)

        # delete the columns which are the statistics for the features but relate to sample size (updated with new column names)
        removed_stats_columns = ['baseq_p', 'alt_baseq1b_p', 'querypos_p', 'leftpos_p', 'seqpos_p', 'mismatches_p', \
                                 'mapq_p', 'read_number_p', 'softclip_length_p', 'softclip_prop_p', 'strand_bias_p', \
                                 'reads_with_indel_p', 'multi_mapper_p', 'per_UMI_end_remove_clip_p']
        for col in removed_stats_columns:
            if col in df.columns:
                df = df.drop(col, axis=1)

    # Keep the features according to the model given
    else:
        # load the trained model feature names
        if phase_refine:
            rf_phased_feature_names_file = os.path.join(model_dir, model_name+"_rf_phased_feature_names.joblib")
            phased_feature_names = load(rf_phased_feature_names_file)
        if use_lr:
            lr_nohet_feature_names_file = os.path.join(model_dir, model_name+"_lr_nohet_feature_names.joblib")
            nohet_feature_names = load(lr_nohet_feature_names_file)
        else:
            rf_feature_names_file = os.path.join(model_dir, model_name+"_rf_feature_names.joblib")
            nohet_feature_names = load(rf_feature_names_file)
        # get the features kept
        features_kept = phased_feature_names if phase_refine else nohet_feature_names
        if not not_pred_het:
            features_kept = list(set(features_kept + het_features))
        # only keep the features we want
        df = df[features_kept + ['haplotype']]

    # ====================================================================
    # Treat Missing Values
    # ====================================================================
    # log transform all p-values
    pval_cols = df.filter(regex='_p$|p_adj').columns
    def log10_1p(x):
        return np.log10(1 + x)
    df[pval_cols] = log10_1p(df[pval_cols])

    # interpolate the NAs in the p-value and t-value columns as log(1)=0 and 0
    psval_cols = [col for col in df.columns if col.endswith('_p') or 
                  col.endswith('_s') or col.endswith('_odds') or col.endswith('_rbc')]
    for col in psval_cols:
        if col in df.columns:
            df[col] = df[col].fillna(0)

    # interpolate the NAs in the phasing-related columns but not p-values as 0
    phase_rel_cols = ['phasing_support_reads_prop_across_hSNPs', 'phasing_nearest_info_mutant_prop', \
                      'phasing_nearest_discordant_prop', 'phasing_nearest_phase_distance', 'phasing_most_info_mutant_prop', \
                      'phasing_most_discordant_prop', 'phasing_most_phase_distance']
    for col in phase_rel_cols:
        if col in df.columns:
            df[col] = df[col].fillna(0)
            
    # interpolate the other NAs with 0 (updated with new column names)
    other_num_cols = ['consensus_alt2_proportion', 'alt2_proportion_per_UMI', \
                      'alt_multi_map_prop', 'ref_multi_map_prop', 'multi_map_prop', \
                      'ref_softclip_prop', 'alt_softclip_prop', 'softclip_prop']
    
    for col in other_num_cols:
        if col in df.columns:
            df[col] = df[col].fillna(0)
        
    # delete the columns with all missing values
    df = df.dropna(axis=1, how='all')
    # drop the rows if still contains NAs
    df = df.dropna()

    # replace the infinity value to finite one
    very_large_number = 1e10
    df.replace(to_replace=np.inf, value=very_large_number, inplace=True)
    df.replace(to_replace=-np.inf, value=-very_large_number, inplace=True)

    # avoid duplicate rows in the data frame
    num_duplicates = df.index.duplicated().sum()
    if num_duplicates > 0:
        print(f"There exist {num_duplicates} duplicate rows, only keep the first appearance.")
        df = df[~df.index.duplicated(keep='first')]

    # ====================================================================
    # Sub Dataframe
    # ====================================================================
    # get the index for each haplotype
    het_index = df.index[df["haplotype"] == "haplo=2"]
    unphased_index = df.index[df["haplotype"] == "unphased"]
    phased_index = df.index[df["haplotype"] != "unphased"]
    haplo3_index = df.index[df["haplotype"] == "haplo=3"]
    haplo_gt3_index = df.index[df["haplotype"] == "haplo>3"]
    hyplotype_label = df['haplotype']
    # remove the hyplotype column
    df = df.drop('haplotype', axis=1)

    # for each non-numeric column, apply encoding
    obj_cols = df.select_dtypes(include=['object']).columns
    if encoder == "label":
        # use label encoding
        for column in obj_cols:
            le = LabelEncoder()
            df[column] = df[column].astype(str)
            df[column] = le.fit_transform(df[column])
            print(f"[LabelEncoder] {column}: {len(le.classes_)} classes")
    else:
        # use one-hot encoding
        print(f"[OneHotEncoder] Encoding columns: {list(obj_cols)}")
        encoder = OneHotEncoder(sparse=False, handle_unknown='ignore')
        encoded = encoder.fit_transform(df.select_dtypes(include=['object']))
        encoded = pd.DataFrame(encoded, columns=encoder.get_feature_names_out())
        # drop original non-numeric columns and concatenate the new DataFrame
        df = pd.concat([df.drop(df.select_dtypes(include=['object']), axis=1), encoded], axis=1)
    
    # sub-dataframe for het and mutation classification (updated with new column names)
    # features for classify heterozygous sites
    if not not_pred_het:
        df_het_features = df[het_features]
    
    # remove the spatial-related features if not using spatial features
    if no_spatial:
        spatial_features = ['pass_spatial_test', 'mut_vs_nonmut_spots_KS_p', 'mut_vs_nonmut_spots_KS_s', \
                            'mut_vs_nonmut_spots_MI_p', 'mut_vs_nonmut_spots_MI_s']
        for col in spatial_features:
            if col in df.columns:
                df = df.drop(col, axis=1)
    
    # get selected features for mutation classification
    if select_features is not None:
        phased_all_cols = list(set(select_features + phase_rel_cols))
        df = df[phased_all_cols]
    if drop_features is not None:
        for col in drop_features:
            if col in df.columns:
                df = df.drop(col, axis=1)

    # remove phase-related columns
    df_nophase = df.copy()
    for col in phase_rel_cols:
        if col in df_nophase.columns:
            df_nophase = df_nophase.drop(col, axis=1)
    # remove het-related columns
    het_used_cols = ['AFind', 'all_spots_vaf_mean']  # updated from 'mean_Afspot'
    for col in het_used_cols:
        if col in df_nophase.columns:
            df_nophase = df_nophase.drop(col, axis=1)

    # ====================================================================
    # Training Set Identify
    # ====================================================================
    # default empty index
    empty_index = pd.MultiIndex.from_tuples([], names=["sample", "#chrom", "pos", "ref", "alt"])
    if train:
        # get the index for artifacts
        validated_artifact_sites = load_or_default_sites(validated_artifact_sites_file, empty_index, 
                                                        "evaluated artifact", sample_name)
        phasable_artifact_sites = load_or_default_sites(phasable_artifact_sites_file, haplo_gt3_index,
                                                        "phasable artifact", sample_name)
        # check the sites given appeared in the data frame
        validated_artifact_sites = validated_artifact_sites.intersection(df.index)
        phasable_artifact_sites = phasable_artifact_sites.intersection(df.index)
        
        # check non-empty true mutation sites
        if true_sites_file is None:
            train = False
            true_sites = empty_index
            print("No true mutation sites given, use default trained model instead.")
        elif isinstance(true_sites_file, str):
            if check_file(true_sites_file):
                true_sites = load_site_index(true_sites_file, sample_name)
                print("Training models with the true sites given.")
            else:
                train = False
                true_sites = empty_index
                print(f"Can not find file {true_sites_file}, use default trained model instead.")
        else:
            raise TypeError("true_sites must be either None or a file path string.")
    else:
        true_sites = empty_index
        validated_artifact_sites = empty_index
        phasable_artifact_sites = empty_index

    # get the training set and the remaning candidate sets
    # convert to sets
    truesites_set = set(true_sites)
    validated_artifact_set = set(validated_artifact_sites)
    phasable_artifact_set = set(phasable_artifact_sites)
    het_set = set(het_index)
    haplo3_set = set(haplo3_index)
    phased_set = set(phased_index)
    unphased_set = set(unphased_index)

    # remove any overlap with true sites
    validated_artifact_set = validated_artifact_set - truesites_set
    phasable_artifact_set = phasable_artifact_set - truesites_set

    # get the phased true sites set
    phased_true_set = truesites_set & phased_set
    phased_true_index = pd.MultiIndex.from_tuples(
                        list(phased_true_set), names=df.index.names).intersection(df.index)
    true_sites = true_sites.intersection(df.index)

    # build final artifact set
    random.seed(random_seed)
    if true_to_artifact_ratio is None:
        selected_phasable_artifact_set = phasable_artifact_set
    else:
        target_artifact_count = int(true_to_artifact_ratio * len(true_sites))
        phasable_needed = target_artifact_count - len(validated_artifact_set)
        if phasable_needed <= 0:
            selected_phasable_artifact_set = set()
        elif phasable_needed >= len(phasable_artifact_set):
            selected_phasable_artifact_set = phasable_artifact_set
        else:
            selected_phasable_artifact_set = set(random.sample(list(phasable_artifact_set), phasable_needed))
    artifact_set = validated_artifact_set | selected_phasable_artifact_set
    artifact_index = pd.MultiIndex.from_tuples(list(artifact_set), names=df.index.names).intersection(df.index)

    # get the remaining candidate sets
    candidate_phased_set = haplo3_set - (truesites_set | artifact_set | het_set)
    candidate_unphased_set = unphased_set - (truesites_set | artifact_set | het_set)
    # get the candidate indices
    candidate_phased_index = pd.MultiIndex.from_tuples(list(candidate_phased_set), names=df.index.names).intersection(df.index)
    candidate_unphased_index = pd.MultiIndex.from_tuples(list(candidate_unphased_set), names=df.index.names).intersection(df.index)

    # print the number of each type site
    if train:
        print("=== Number of sites for each type ===")
        print("validated phased mosaic:", len(phased_true_index))
        print("validated mosaic:", len(true_sites))
        print("validated het:", len(het_index))
        print("validated artifact:", len(validated_artifact_set))
        print("phasable artifact (used):", len(selected_phasable_artifact_set))
        print("total artifact:", len(artifact_set))
        print("haplotype=3:", len(haplo3_index))
        print("candidate_phased:", len(candidate_phased_index))
        print("candidate_unphased:", len(candidate_unphased_index))
        print("=====================================")
    else:
        total_candidate_sites = len(candidate_phased_index) + len(candidate_unphased_index)
        print(f"Total candidate sites: {total_candidate_sites}")


    # ====================================================================
    # Training Models
    # ====================================================================
    if train:
        # saved path
        if save_models:
            model_dir = os.path.join(output_dir, "models")
            check_dir(model_dir)
            print(f"Save the trained models to {model_dir}")
        if plot:
            plot_dir = os.path.join(output_dir, "plots")
            check_dir(plot_dir)
            pca_data_dir = os.path.join(output_dir, "pca_data")
            check_dir(pca_data_dir)
            print(f"Save the plots to {plot_dir}")
            print(f"Save the PCA data to {pca_data_dir}")
        # ====================================================================
        # Phase Refinement Model
        # ====================================================================
        if phase_refine:
            # pair each index with its label
            phased_paired_indices = [(index, "mosaic") for index in phased_true_index] + \
                                    [(index, "artifact") for index in artifact_index] + \
                                    [(index, "het") for index in het_index]
            # shuffle the list to mix indices of different labels
            random.seed(random_seed)
            random.shuffle(phased_paired_indices)
            # generate training data
            X_phased_df = df.loc[[idx for idx, _ in phased_paired_indices]]
            y_phased_list = [label for _, label in phased_paired_indices]
            # convert to matrix
            X_phased = X_phased_df.to_numpy()
            y_phased = np.array(y_phased_list)
            
            # train the phase refinement random forest model
            print("Start training phasing refinement random forest model...")
            rf_phased, _, _, _ = random_forest(X_phased, y_phased, test_size=0, random_state=random_seed, n_labels=3, \
                                            smote=smote, tune=tune, sampling_strategy=sampling_strategy, n_jobs=n_jobs, \
                                            n_estimators=n_estimators, max_depth=max_depth, min_samples_split=min_samples_split)
            # save the model
            if save_models:
                rf_phased_file = os.path.join(model_dir, sample_name+"_rf_phased_model.joblib")
                dump(rf_phased, rf_phased_file)
                # save the feature names
                phased_feature_names = X_phased_df.columns.tolist()
                rf_phased_feature_names_file = os.path.join(model_dir, sample_name+"_rf_phased_feature_names.joblib")
                dump(phased_feature_names, rf_phased_feature_names_file)
                print("Save the phase refinement random forest model")
            # plots
            if plot:
                # plot the feature importance
                sorted_features_phased = feature_importamce_plot(rf_phased, X_phased_df, plot_dir, sample_name, \
                                                        title_note=" of Phasing Refinement", file_note="_phase")
                # perform PCA and plot
                pca_phased, top_features_phased, pca_phased_train_data = PCA_train(sorted_features_phased, X_phased_df, y_phased, plot_dir, sample_name, \
                                                    title_note=" of Phasing Refinement", file_note="_phase", n_features=n_features, fig_size=6, \
                                                    annotate_mosaic=annotate_mosaic, annotate_outlier=annotate_outlier)
                # save pca data
                pca_phased_train_file = os.path.join(pca_data_dir, sample_name+"_pca_phased_train_data.csv")
                pca_phased_train_data.to_csv(pca_phased_train_file, sep='\t', index=True, header=True)
                # save pca model
                if save_models:
                    pca_phased_file = os.path.join(model_dir, sample_name+"_pca_phased_model.pkl")
                    with open(pca_phased_file, 'wb') as file:
                        pickle.dump(pca_phased, file)
            #  add the labels to the data frame
            X_phased_df['evaluation'] = y_phased

        # ====================================================================
        # Logistic Regression Model For Het
        # ====================================================================
        if (len(het_index) > 0) and (not not_pred_het):
            # pair each index with its label
            all_paired_indices = [(index, "nothet") for index in true_sites] + \
                                [(index, "nothet") for index in artifact_index] + \
                                [(index, "het") for index in het_index]
            # shuffle the list to mix indices of different labels
            random.seed(random_seed)
            random.shuffle(all_paired_indices)
            # use part features to train the het
            X_all_df = df_het_features.loc[[idx for idx, _ in all_paired_indices]]
            y_all_list = [label for _, label in all_paired_indices]
            # convert to matrix
            X_all = X_all_df.to_numpy()
            y_all = np.array(y_all_list)

            # train the logistic regression model to seperate the het sites
            print("Start training heterozygous classification logistic regression model...")
            lr_model = logistic_regression(X_all, y_all, test_size=0, random_state=random_seed, k_neighbors=3, smote=False, \
                                           sampling_strategy=sampling_strategy, n_jobs=n_jobs)
            # save the model
            if save_models:
                lr_model_file = os.path.join(model_dir, sample_name+"_lr_model.joblib")
                dump(lr_model, lr_model_file)
                print("Save the heterozygous classification logistic regression model")
            #  add the labels to the data frame
            X_all_df['evaluation'] = y_all
            no_het_lr_model = False
        else:
            no_het_lr_model = True

        # # ====================================================================
        # # Mutaion Classification Random Forest Model
        # # ====================================================================
        # # pair each index with its label
        # nohet_paired_indices = [(index, "mosaic") for index in true_sites] + \
        #                         [(index, "artifact") for index in artifact_index]
        # # shuffle the list to mix indices of different labels
        # random.seed(random_seed)
        # random.shuffle(nohet_paired_indices)
        # # generate training data
        # X_nohet_df = df_nophase.loc[[idx for idx, _ in nohet_paired_indices]]
        # y_nohet_list = [label for _, label in nohet_paired_indices]
        # # convert to matrix
        # X_nohet = X_nohet_df.to_numpy()
        # y_nohet = np.array(y_nohet_list)

        # # train the mutation classification random forest model
        # print("Start training somatic mutation classification random forest model...")
        # rf, _, _, _ = random_forest(X_nohet, y_nohet, test_size=0, random_state=random_seed, n_labels=2, \
        #                             smote=smote, tune=tune, sampling_strategy=sampling_strategy, n_jobs=n_jobs, \
        #                             n_estimators=n_estimators, max_depth=max_depth, min_samples_split=min_samples_split)
        # # save the model
        # if save_models:
        #     rf_model_file = os.path.join(model_dir, sample_name+"_rf_model.joblib")
        #     dump(rf, rf_model_file)
        #     # save the feature names
        #     nohet_feature_names = X_nohet_df.columns.tolist()
        #     rf_feature_names_file = os.path.join(model_dir, sample_name+"_rf_feature_names.joblib")
        #     dump(nohet_feature_names, rf_feature_names_file)
        #     print("Save the somatic mutation classification random forest model")
        # # plots
        # if plot:
        #     # plot the feature importance
        #     sorted_features_nohet = feature_importamce_plot(rf, X_nohet_df, plot_dir, sample_name, 
        #                             title_note=" of Mutation Classification", file_note="_mutation")
        #     # perform PCA and plot
        #     pca_nohet, top_features_nohet, pca_nohet_train_data = PCA_train(sorted_features_nohet, X_nohet_df, y_nohet, plot_dir, sample_name, \
        #                                     title_note=" of Mutation Classification", file_note="_mutation", n_features=n_features, \
        #                                     fig_size=6, annotate_mosaic=annotate_mosaic, annotate_outlier=annotate_outlier)
        #     # save pca data
        #     pca_nohet_train_file = os.path.join(pca_data_dir, sample_name+"_pca_nohet_train_data.csv")
        #     pca_nohet_train_data.to_csv(pca_nohet_train_file, sep='\t', index=True, header=True)
        #     # save pca model
        #     if save_models:
        #         pca_nohet_file = os.path.join(model_dir, sample_name+"_pca_nohet_model.pkl")
        #         with open(pca_nohet_file, 'wb') as file:
        #             pickle.dump(pca_nohet, file)
        # #  add the labels to the data frame
        # X_nohet_df['evaluation'] = y_nohet          

        # ====================================================================
        # Mutaion Classification using Logistic Regression Model or Random Forest Model
        # ====================================================================
        # pair each index with its label
        nohet_paired_indices = [(index, "mosaic") for index in true_sites] + \
                                [(index, "artifact") for index in artifact_index]
        # shuffle the list to mix indices of different labels
        random.seed(random_seed)
        random.shuffle(nohet_paired_indices)
        # generate training data
        X_nohet_df = df_nophase.loc[[idx for idx, _ in nohet_paired_indices]]
        y_nohet_list = [label for _, label in nohet_paired_indices]
        # convert to matrix
        X_nohet = X_nohet_df.to_numpy()
        y_nohet = np.array(y_nohet_list)

        # choose models
        if use_lr:
            # train the logistic regression model to do mutation classification for the not-het sites
            print("Start training somatic mutation classification logistic regression model...")
            lr_nohet_model = logistic_regression(X_nohet, y_nohet, test_size=0, random_state=random_seed, 
                                                 penalty='l1', solver='liblinear', k_neighbors=3, smote=False, \
                                                 sampling_strategy=sampling_strategy, n_jobs=n_jobs)
            # get the non-zero coefficients
            lr_nohet_model_coef = lr_nohet_model.coef_[0]
            lr_nohet_model_features = X_nohet_df.columns[lr_nohet_model_coef != 0]
            # save the model
            if save_models:
                lr_nohet_model_file = os.path.join(model_dir, sample_name+"_lr_nohet_model.joblib")
                dump(lr_nohet_model, lr_nohet_model_file)
                # save the feature names
                nohet_feature_names = lr_nohet_model_features.tolist()
                lr_nohet_feature_names_file = os.path.join(model_dir, sample_name+"_lr_nohet_feature_names.joblib")
                dump(nohet_feature_names, lr_nohet_feature_names_file)
                print("Save the somatic mutation classification logistic regression model")
        else:
            # train the mutation classification random forest model
            print("Start training somatic mutation classification random forest model...")
            rf, _, _, _ = random_forest(X_nohet, y_nohet, test_size=0, random_state=random_seed, n_labels=2, \
                                        smote=smote, tune=tune, sampling_strategy=sampling_strategy, n_jobs=n_jobs, \
                                        n_estimators=n_estimators, max_depth=max_depth, min_samples_split=min_samples_split)
            # save the model
            if save_models:
                rf_model_file = os.path.join(model_dir, sample_name+"_rf_model.joblib")
                dump(rf, rf_model_file)
                # save the feature names
                nohet_feature_names = X_nohet_df.columns.tolist()
                rf_feature_names_file = os.path.join(model_dir, sample_name+"_rf_feature_names.joblib")
                dump(nohet_feature_names, rf_feature_names_file)
                print("Save the somatic mutation classification random forest model")
        # plots
        if plot:
            if use_lr:
                # plot the feature importance for logistic regression
                sorted_features_nohet = lr_model_feature_importance_plot(lr_nohet_model, X_nohet_df, plot_dir, sample_name, 
                                    title_note=" of Mutation Classification", file_note="_mutation")
            else:
                # plot the feature importance
                sorted_features_nohet = feature_importamce_plot(rf, X_nohet_df, plot_dir, sample_name, 
                                        title_note=" of Mutation Classification", file_note="_mutation")
            # perform PCA and plot
            pca_nohet, top_features_nohet, pca_nohet_train_data = PCA_train(sorted_features_nohet, X_nohet_df, y_nohet, plot_dir, sample_name, \
                                            title_note=" of Mutation Classification", file_note="_mutation", n_features=n_features, \
                                            fig_size=6, annotate_mosaic=annotate_mosaic, annotate_outlier=annotate_outlier)
            # save pca data
            pca_nohet_train_file = os.path.join(pca_data_dir, sample_name+"_pca_nohet_train_data.csv")
            pca_nohet_train_data.to_csv(pca_nohet_train_file, sep='\t', index=True, header=True)
            # save pca model
            if save_models:
                pca_nohet_file = os.path.join(model_dir, sample_name+"_pca_nohet_model.pkl")
                with open(pca_nohet_file, 'wb') as file:
                    pickle.dump(pca_nohet, file)
        #  add the labels to the data frame
        X_nohet_df['evaluation'] = y_nohet          

    else:
        # print(f"Load the saved trained models for {model_name}")
        if phase_refine:
            # get the model files
            rf_phased_file = os.path.join(model_dir, model_name+"_rf_phased_model.joblib")
            rf_phased_feature_names_file = os.path.join(model_dir, model_name+"_rf_phased_feature_names.joblib")
            # load the trained models
            rf_phased = load(rf_phased_file)
            phased_feature_names = load(rf_phased_feature_names_file)
            if plot:
                # get the pca model files
                pca_phased_file = os.path.join(model_dir, model_name+"_pca_phased_model.pkl")
                # load pca models
                with open(pca_phased_file, 'rb') as file:
                    pca_phased = pickle.load(file)
        # get the model files
        lr_model_file = os.path.join(model_dir, model_name+"_lr_model.joblib")
        if use_lr:
            lr_nohet_model_file = os.path.join(model_dir, model_name+"_lr_nohet_model.joblib")
            lr_nohet_feature_names_file = os.path.join(model_dir, model_name+"_lr_nohet_feature_names.joblib")
        else:
            rf_model_file = os.path.join(model_dir, model_name+"_rf_model.joblib")
            rf_feature_names_file = os.path.join(model_dir, model_name+"_rf_feature_names.joblib")
        # load the trained models
        if os.path.exists(lr_model_file):
            lr_model = load(lr_model_file)
            no_het_lr_model = False
        else:
            no_het_lr_model = True
        if use_lr:
            lr_nohet_model = load(lr_nohet_model_file)
            nohet_feature_names = load(lr_nohet_feature_names_file)
        else:
            rf = load(rf_model_file)
            nohet_feature_names = load(rf_feature_names_file)
        if plot:
            # get the pca model files
            pca_nohet_file = os.path.join(model_dir, model_name+"_pca_nohet_model.pkl")
            # load pca models
            with open(pca_nohet_file, 'rb') as file:
                pca_nohet = pickle.load(file)

    # create the folder to save SHAP values
    if save_shap and not use_lr:
        explainer = shap.TreeExplainer(rf)
        SHAP_data_dir = os.path.join(output_dir, "SHAP_values")
        check_dir(SHAP_data_dir)
        

    # ====================================================================
    # Predict
    # ====================================================================
    # predict the phsed haplotype=3 sites
    if len(candidate_phased_index) > 0:
        pred_phase = True
        if phase_refine:
            # using phase refinement
            candidate_phased_df = df.loc[candidate_phased_index]
            candidate_phased = candidate_phased_df.to_numpy()
            candidate_phased_pred = rf_phased.predict(candidate_phased)
            # get the number of each values
            phsed_pred_counts = Counter(candidate_phased_pred)
            # print("Phase refinement prediction result:")
            # print(phsed_pred_counts)
        else:
            # using random forest model for non-phased sites
            candidate_phased_df = df_nophase.loc[candidate_phased_index]
            candidate_phased = candidate_phased_df.to_numpy()
            if use_lr:
                candidate_phased_pred = lr_nohet_model.predict(candidate_phased)
            else:
                candidate_phased_pred = rf.predict(candidate_phased)
            # get the number of each values
            phsed_pred_counts = Counter(candidate_phased_pred)
            # print("Phased set somatic mutation prediction result:")
            # print(phsed_pred_counts)
            # calculate and save the SHAP values
            if save_shap and not use_lr:
                # calculate SHAP values
                shap_values_phased = explainer.shap_values(candidate_phased)
                # find the index of mosaic
                mosaic_class_index_phased = list(rf.classes_).index("mosaic")
                mosaic_shap_values_phased = shap_values_phased[:, :, mosaic_class_index_phased]
                # save SHAP values and expected values to a CSV file
                shap_phased_df = pd.DataFrame(mosaic_shap_values_phased, columns=candidate_phased_df.columns, index=candidate_phased_df.index)
                shap_values_phased_file = os.path.join(SHAP_data_dir, sample_name+"_shap_values_phased.csv")
                shap_phased_df.to_csv(shap_values_phased_file, sep='\t', index=True, header=True)
                # save the expected value for SHAP
                shap_expected_values_phased_file = os.path.join(SHAP_data_dir, sample_name+"_shap_expected_values_phased.csv")
                np.savetxt(shap_expected_values_phased_file, [explainer.expected_value[mosaic_class_index_phased]])
                # save training data for further analysis
                candidate_phased_file = os.path.join(SHAP_data_dir, sample_name+"_candidate_phased_data.csv")
                candidate_phased_df.to_csv(candidate_phased_file, sep='\t', index=True, header=True)

        # get predicted true sites
        candidate_phased_df['pred'] = candidate_phased_pred
        phased_pred_true = candidate_phased_df[candidate_phased_df['pred']=="mosaic"]
        phased_pred_true_barcode = phased_pred_true.index
        phased_pred_het = candidate_phased_df[candidate_phased_df['pred']=="het"]
        phased_pred_het_barcode = phased_pred_het.index
        # # save the true sites
        # results_dir = os.path.join(output_dir, "results")
        # check_dir(results_dir)
        # phased_pred_file = os.path.join(results_dir, sample_name+"_phased_pred_truesites.txt")
        # with open(phased_pred_file, 'w') as file:
        #     for item in phased_pred_true_barcode:
        #         file.write(f"{item}\n")
    else:
        pred_phase = False


    # predict the het sites
    # candidate_df = df.loc[candidate_index]
    if no_het_lr_model or not_pred_het:
        pred_nohet_barcode = candidate_unphased_index
    else:
        candidate_df = df_het_features.loc[candidate_unphased_index]
        candidate = candidate_df.to_numpy()
        # normalize
        scaler = StandardScaler()
        candidate = scaler.fit_transform(candidate)
        candidate_pred = lr_model.predict(candidate)
        # get the number of each values
        pred_counts = Counter(candidate_pred)
        # print("Heterozygous classification prediction result:")
        # print(pred_counts)
        # get predicted het and nohet sites
        candidate_df['pred'] = candidate_pred
        pred_het = candidate_df[candidate_df['pred']=="het"]
        pred_het_barcode = pred_het.index
        pred_nohet = candidate_df[candidate_df['pred']=="nothet"]
        pred_nohet_barcode = pred_nohet.index
        
    # predict the somatic mutations
    candidate_nohet_df = df_nophase.loc[pred_nohet_barcode]
    candidate_nohet = candidate_nohet_df.to_numpy()
    if use_lr:
        # normalize
        scaler = StandardScaler()
        candidate_nohet = scaler.fit_transform(candidate_nohet)
        candidate_nohet_pred = lr_nohet_model.predict(candidate_nohet)
    else:        
        candidate_nohet_pred = rf.predict(candidate_nohet)
        

    # get the number of each values
    nohet_pred_counts = Counter(candidate_nohet_pred)
    # print("Somatic mutation classification prediction result:")
    # print(nohet_pred_counts)
    # calculate and save the SHAP values
    if save_shap and not use_lr:
        # calculate SHAP values
        shap_values_nohet = explainer.shap_values(candidate_nohet)
        # find the index of mosaic
        mosaic_class_index_nohet = list(rf.classes_).index("mosaic")
        mosaic_shap_values_nohet = shap_values_nohet[:, :, mosaic_class_index_nohet]
        # save SHAP values and expected values to a CSV file
        shap_nohet_df = pd.DataFrame(mosaic_shap_values_nohet, columns=candidate_nohet_df.columns, index=candidate_nohet_df.index)
        shap_values_nohet_file = os.path.join(SHAP_data_dir, sample_name+"_shap_values_nohet.csv")
        shap_nohet_df.to_csv(shap_values_nohet_file, sep='\t', index=True, header=True)
        # save the expected value for SHAP
        shap_expected_values_nohet_file = os.path.join(SHAP_data_dir, sample_name+"_shap_expected_values_nohet.csv")
        np.savetxt(shap_expected_values_nohet_file, [explainer.expected_value[mosaic_class_index_nohet]])
        # save training data for further analysis
        candidate_nohet_file = os.path.join(SHAP_data_dir, sample_name+"_candidate_nohet_data.csv")
        candidate_nohet_df.to_csv(candidate_nohet_file, sep='\t', index=True, header=True)
    # get predicted true sites
    candidate_nohet_df['pred'] = candidate_nohet_pred
    nohet_pred_true = candidate_nohet_df[candidate_nohet_df['pred']=="mosaic"]
    nohet_pred_true_barcode = nohet_pred_true.index
    # # save predicted true sites in a file
    # results_dir = os.path.join(output_dir, "results")
    # check_dir(results_dir)
    # pred_file = os.path.join(results_dir, sample_name+"_pred_truesites.txt")
    # with open(pred_file, 'w') as file:
    #     for item in nohet_pred_true_barcode:
    #         file.write(f"{item}\n")

    # combine the predicted sites
    if pred_phase:
        total_pred_true_barcode = phased_pred_true_barcode.union(nohet_pred_true_barcode)
    else:
        total_pred_true_barcode = nohet_pred_true_barcode

    # save the total predicted true sites
    results_dir = os.path.join(output_dir, "results")
    check_dir(results_dir)

    # subset output
    df_output = df_outputfeatures.loc[total_pred_true_barcode].copy()
    df_output = df_output.reset_index()
    # drop sample column if present
    if "sample" in df_output.columns:
        df_output = df_output.drop(columns=["sample"])
    # rename index-derived columns to VCF-style names
    df_output = df_output.rename(
        columns={
            "#chrom": "#CHROM",
            "chrom": "#CHROM",
            "pos": "POS",
            "ref": "REF",
            "alt": "ALT",
        }
    )
    # chromosome-aware sorting
    if {"#CHROM", "POS"}.issubset(df_output.columns):
        df_output = df_output.copy()
        df_output["_chrom_sort_key"] = df_output["#CHROM"].map(chrom_sort_key)
        df_output = (
            df_output
            .sort_values(by=["_chrom_sort_key", "POS"])
            .drop(columns="_chrom_sort_key")
            .reset_index(drop=True)
        )

    # PASS-only sites
    df_output_pass = df_output[df_output["FILTER"] == "PASS"].copy()

    # output paths
    # parquet_output_file = os.path.join(results_dir, sample_name + "_total_pred_truesites.parquet")
    # tsv_output_file = os.path.join(results_dir, sample_name + "_total_pred_truesites.tsv")
    vcf_output_file = os.path.join(results_dir, sample_name + "_total_pred_truesites.vcf")
    # parquet_pass_output_file = os.path.join(results_dir, sample_name + "_total_pred_truesites_PASS.parquet")
    # tsv_pass_output_file = os.path.join(results_dir, sample_name + "_total_pred_truesites_PASS.tsv")
    vcf_pass_output_file = os.path.join(results_dir, sample_name + "_total_pred_truesites_PASS.vcf")
    pass_mutation_list_file = os.path.join(results_dir, sample_name + "_total_pred_truesites_PASS_mutation_list.txt")
    # save all sites
    # df_output.to_parquet(parquet_output_file, index=False)
    # df_output.to_csv(tsv_output_file, sep="\t", index=False, float_format="%.6g")
    write_simple_vcf(df_output, vcf_output_file, sample_name=sample_name)
    # save PASS-only sites
    # df_output_pass.to_parquet(parquet_pass_output_file, index=False)
    # df_output_pass.to_csv(tsv_pass_output_file, sep="\t", index=False, float_format="%.6g")
    write_simple_vcf(df_output_pass, vcf_pass_output_file, sample_name=sample_name)
    # save PASS-only mutation list in chrom_pos_ref_alt format
    with open(pass_mutation_list_file, "w") as file:
        for _, row in df_output_pass.iterrows():
            chrom = str(row["#CHROM"])
            pos = int(row["POS"]) if pd.notna(row["POS"]) else row["POS"]
            ref = str(row["REF"])
            alt = str(row["ALT"])
            file.write(f"{chrom}_{pos}_{ref}_{alt}\n")
    # print notes
    print(f"Total predicted true sites: {len(df_output)}")
    print(f"PASS predicted true sites: {len(df_output_pass)}")
    print("Saved output files:")
    print(f"  - All sites: {vcf_output_file}")
    print(f"  - PASS sites: {vcf_pass_output_file}")
    print(f"  - PASS mutation list: {pass_mutation_list_file}")


    # pca scatter plot
    if plot:
        plot_dir = os.path.join(output_dir, "plots")
        check_dir(plot_dir)
        pca_data_dir = os.path.join(output_dir, "pca_data")
        check_dir(pca_data_dir)
        if train:
            # plot for phasable sites
            if pred_phase and len(candidate_phased_index) > 1:
                if phase_refine:
                    pca_phased_pred_data = PCA_pred(X_phased_df, candidate_phased_df, plot_dir, sample_name, 
                                title_note=" of Phasing Refinement", file_note="_phase", 
                                fig_size=6, pca_train=pca_phased, top_features=top_features_phased,
                                annotate_outlier=annotate_outlier)
                else:
                    pca_phased_pred_data = PCA_pred(X_nohet_df, candidate_phased_df, plot_dir, sample_name, 
                                title_note=" of Phasing Refinement", file_note="_phase", 
                                fig_size=6, pca_train=pca_nohet, top_features=top_features_nohet)
            # plot for mutation prediction
            pca_nohet_pred_data = PCA_pred(X_nohet_df, candidate_nohet_df, plot_dir, sample_name, 
                        title_note=" of Mutation Classification", file_note="_mutation", 
                        fig_size=6, pca_train=pca_nohet, top_features=top_features_nohet)
        else:
            if phase_refine:
                # get feature importances
                rf_phased_importances = rf_phased.feature_importances_
                # sort the features and their importances
                rf_phsed_sorted_importances, rf_phased_sorted_features = zip(*sorted(zip(rf_phased_importances, phased_feature_names), reverse=True))
            # no het model
            if use_lr:
                lr_nohet_importances = np.abs(lr_nohet_model.coef_[0])
                lr_nohet_sorted_importances, lr_nohet_sorted_features = zip(*sorted(zip(lr_nohet_importances, nohet_feature_names), reverse=True))
            else:
                rf_importances = rf.feature_importances_
                rf_sorted_importances, rf_sorted_features = zip(*sorted(zip(rf_importances, nohet_feature_names), reverse=True))
            # get the top features
            if phase_refine:
                top_features_phased = list(rf_phased_sorted_features[:n_features])
            if use_lr:
                top_features_nohet = list(lr_nohet_sorted_features[:n_features])
            else:
                top_features_nohet = list(rf_sorted_features[:n_features])

            # plot PCA scatter plot
            if pred_phase and len(candidate_phased_index) > 1:
                if phase_refine:            
                    pca_phased_pred_data = PCA_pred(None, candidate_phased_df, plot_dir, sample_name, title_note=" of Phasing Refinement", file_note="_phase", fig_size=6,
                                top_features=top_features_phased)
                else:
                    pca_phased_pred_data = PCA_pred(None, candidate_phased_df, plot_dir, sample_name, title_note=" of Phasing Refinement", file_note="_phase", fig_size=6,
                                top_features=top_features_nohet)
            pca_nohet_pred_data = PCA_pred(None, candidate_nohet_df, plot_dir, sample_name, title_note=" of Mutation Classification", file_note="_mutation", fig_size=6,
                        top_features=top_features_nohet)
        
        # save pca data
        if save_pca:
            if pred_phase and len(candidate_phased_index) > 1:
                pca_phased_pred_file = os.path.join(pca_data_dir, sample_name+"_pca_phased_pred_data.csv")
                pca_phased_pred_data.to_csv(pca_phased_pred_file, sep='\t', index=True, header=True)
            pca_nohet_pred_file = os.path.join(pca_data_dir, sample_name+"_pca_nohet_pred_data.csv")
            pca_nohet_pred_data.to_csv(pca_nohet_pred_file, sep='\t', index=True, header=True)



    
# ====================================================================
# Supplementary Functions
# ====================================================================

def list2min(value):
    """find the minimum value in the comma-separated string"""
    # check if the input is a float number
    if isinstance(value, float):
        return value
    elif isinstance(value, int):
        return value
    else:
        return min(float(num) if num!="no" else np.nan for num in value.split(','))
    
def list2mean(value):
    """get the mean value in the comma-separated string"""
    if isinstance(value, float):
        return value
    elif isinstance(value, int):
        return value
    else:
        return np.mean([float(num) for num in value.split(',')])

def list2frequent(value):
    """find the most frequent element in the comma-separated string"""
    if isinstance(value, float):
        return value
    elif isinstance(value, int):
        return value
    else:
        return Counter(value.split(',')).most_common(1)[0][0]
    


def merge_haplotype_columns(row, haplotype_columns):
    """Merge the multiple haplotype-related columns to one haplotype columns"""
    # extract the haplotype-related values
    values = row[haplotype_columns].dropna().unique()  # remove NaN values and get unique values
    NA_name = "unphased"  #or np.nan
    
    # Rule 1: If all values are the same, then take the same value
    if len(values) == 1:
        return values[0]
    
    # Rule 2: If 'haplo=3' is one of the values, then keep 'haplo=3'
    if 'haplo=3' in values:
        return 'haplo=3'
    
    # Rule 3: Take the most frequent value among the specific columns
    if len(values) > 0:
        most_frequent = row[haplotype_columns].dropna().mode()
        return most_frequent.iloc[0] if not most_frequent.empty else NA_name
    else:
        return NA_name
    

def load_site_index(path, sample_name):
    """Read in sites file as multi-index data frame"""
    df_sites = pd.read_csv(path, sep="\t")
    if "sample" not in df_sites.columns:
        df_sites["sample"] = sample_name
    index_cols = ["sample", "#chrom", "pos", "ref", "alt"]
    df_sites = df_sites[index_cols]
    return pd.MultiIndex.from_frame(df_sites)

def load_or_default_sites(site_file, default_index, label, sample_name):
    """Load site file or return default index if not provided."""
    if site_file is None:
        print(f"No {label} sites provided. Using default {label} sites.")
        return default_index
    elif isinstance(site_file, str):
        if check_file(site_file):
            print(f"Training models with the {label} sites given.")
            return load_site_index(site_file, sample_name)
        else:
            print(f"Can not find file {site_file}, using default {label} sites.")
            return default_index
    else:
        raise TypeError(f"{label} must be None or a file path string.")
    


def random_forest(X, y, test_size=0.3, random_state=100, smote=True, tune="random_search", \
                  n_labels = 2, k_neighbors=4, sampling_strategy="auto", n_jobs=None, \
                  n_estimators=100, max_depth=None, min_samples_split=2):
    """
    Generate a random forest model for the classification of each type for the candidate mutation sets

    Inputs:
        X, y - the input data and labels (in the numpy format)
        test_size - the test set used when splitting the training and test set (default = 0.3)
        random_state - the random state for the train-test splitting and random forest model (default = 100)
        smote - whether use SMOTE (Synthetic Minority Over-sampling Technique) to over-sample the minority class 
                to treat the imbalance class values (default = True)
        tune - whether tuning the hyperparameters used in the random forest (n_estimators, max_depth, min_samples_split)
                "Bayesian_opt": Bayesian Optimization
                "random_search": Random Search
                "grid_search": Grid Search
                (default = "random_search)
        n_labels - the number of labels (choice = {2, 3}, default = 2)
        k_neighbors - the nearest neighbors used to define the neighborhood of samples in SMOTE (default = 4)
        sampling_strategy - sampling information to resample the data set (default = "auto")
                'minority': resample only the minority class;
                'not minority': resample all classes but the minority class;
                'not majority': resample all classes but the majority class;
                'all': resample all classes;
                'auto': equivalent to 'not majority'.
        n_job - number of jobs to run in parallel (default = None) 
                None means 1 unless in a joblib.parallel_backend context. 
                -1 means using all processors.         
        n_estimators - the number of trees in the forest (default = 100)
        max_depth - the maximum depth of the tree (default = None)
        min_samples_split - the minimum number of samples required to split an internal node (default = 2)
    """
    # split dataset into training set and test set
    if test_size != 0:
        X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=test_size, stratify=y, random_state=random_state)
    else:
        X_train = X
        y_train = y
    
    # ======================================================================
    # Over-sampling for the imbalance target values
    # ======================================================================
    # use SMOTE to over-sampling
    if smote:
        smote = SMOTE(random_state=random_state, sampling_strategy=sampling_strategy, k_neighbors=k_neighbors, n_jobs=n_jobs)
        X_train, y_train = smote.fit_resample(X_train, y_train)
        # get the number of each values
        print("Training set count:", dict(Counter(y_train)))


    # ======================================================================
    # Tune the hyperparameters
    # ======================================================================
    if tune == "Bayesian_opt":
        # define the search space of hyperparameters
        space = {
            'n_estimators': hp.choice('n_estimators', range(100, 500)),
            'max_depth': hp.choice('max_depth', [None, 10, 20, 30, 40]),
            'min_samples_split': hp.choice('min_samples_split', range(2, 11))
        }

        # objective function to minimize
        def objective(space):
            model = RandomForestClassifier(n_estimators=space['n_estimators'], max_depth=space['max_depth'], \
                                           min_samples_split=space['min_samples_split'])
            accuracy = cross_val_score(model, X_train, y_train, cv=5).mean()
            return {'loss': -accuracy, 'status': STATUS_OK}

        # run the algorithm
        trials = Trials()
        best = fmin(fn=objective, space=space, algo=tpe.suggest, max_evals=100, trials=trials)
        # print optimized hyperparameters
        print("Optimized hyperparameters:", best)

        # get the optimized parameter values respectively
        n_estimators_opt = range(100, 500)[best['n_estimators']]
        max_depth_opt = [None, 10, 20, 30, 40][best['max_depth']]
        min_samples_split_opt = range(2, 11)[best['min_samples_split']]
        # get the random forest model
        rf = RandomForestClassifier(n_estimators=n_estimators_opt, 
                                    max_depth=max_depth_opt, 
                                    min_samples_split=min_samples_split_opt,
                                    bootstrap=True, oob_score=True, 
                                    random_state=random_state)
    
    elif tune == "random_search":
        # define the parameter distribution
        param_dist = {
            'n_estimators': randint(100, 500),
            'max_depth': [None, 10, 20, 30, 40],
            'min_samples_split': randint(2, 11)
        }
        # initialize the classifier
        rf = RandomForestClassifier(bootstrap=True, oob_score=True, random_state=random_state)
        # initialize the Random Search model
        random_search = RandomizedSearchCV(estimator=rf, param_distributions=param_dist, n_iter=100, cv=5, verbose=0, \
                                           random_state=random_state, n_jobs=n_jobs)

        # fit the Random Search to the data
        random_search.fit(X_train, y_train)
        # get and print the optimized parameters
        best_params = random_search.best_params_
        print("Optimized hyperparameters:", best_params)
        # get the random forest model
        rf = RandomForestClassifier(**best_params, bootstrap=True, oob_score=True, random_state=random_state)
    
    elif tune == "grid_search":
        # define the parameter grid
        param_grid = {
            'n_estimators': [100, 200, 300],
            'max_depth': [None, 10, 20, 30],
            'min_samples_split': [2, 5, 10]
        }
        # initialize the classifier
        rf = RandomForestClassifier(bootstrap=True, oob_score=True, random_state=random_state)
        # initialize the Grid Search model
        grid_search = GridSearchCV(estimator=rf, param_grid=param_grid, cv=5, verbose=0, n_jobs=n_jobs)
        # fit the Grid Search to the data
        grid_search.fit(X_train, y_train)
        # get and print the optimized parameters
        best_params = grid_search.best_params_
        print("Optimized hyperparameters:", best_params)
        # get the random forest model
        rf = RandomForestClassifier(**best_params, bootstrap=True, oob_score=True, random_state=random_state)
    
    else:
        # create a Gaussian classifier with the default hyperparameters
        rf = RandomForestClassifier(n_estimators=n_estimators, max_depth=max_depth, min_samples_split=min_samples_split, \
                                    bootstrap=True, oob_score=True, random_state=random_state)
    

    # ======================================================================
    # Train the model
    # ======================================================================
    # train
    rf.fit(X_train, y_train)
    # print the average out-of-bag accuracy
    print("OOB Score:", rf.oob_score_)

    if test_size != 0:
        # predict for the test set
        y_pred = rf.predict(X_test)
        # predict for the test true mutation set
        true_indices = np.where(y_test == "mosaic")[0]
        X_test_true = X_test[true_indices, :]
        y_test_true = y_test[true_indices]
        y_pred_true = rf.predict(X_test_true)
        # generate confusion matrix
        if n_labels == 3:
            cm = confusion_matrix(y_test, y_pred, labels=['mosaic', 'het', 'artifact'])
        elif n_labels == 2:
            cm = confusion_matrix(y_test, y_pred, labels=['mosaic', 'artifact'])
        else:
            cm = confusion_matrix(y_test, y_pred)
        print(cm)
        # calculate the sensitivity
        test_sensitivity = cm[0, 0] / cm[0, :].sum()

        # print the test set accuracy
        test_accuracy = accuracy_score(y_test, y_pred)
        test_true_accuracy = accuracy_score(y_test_true, y_pred_true)
        print("Test Accuracy:", test_accuracy)
        print("Test True Mutation Accuracy:", test_true_accuracy)
        print(f"Test Sensitivity: {test_sensitivity}")
        
    else:
        test_accuracy, test_true_accuracy, test_sensitivity = None, None, None

    return rf, test_accuracy, test_true_accuracy, test_sensitivity



def logistic_regression(X, y, test_size=0.3, random_state=100, penalty='l2', solver='lbfgs', max_iter=1000, \
                        class_weight='balanced', smote=True, k_neighbors=4, sampling_strategy="auto", n_jobs=None):
    """
    Use Logistic regression method to seperate heterozygous sites from candidate sites

    Inputs:
        X, y - the input data and labels (in the numpy format)
        test_size - the test set used when splitting the training and test set (default = 0.3)
        random_state - the random state for the train-test splitting and random forest model (default = 100)
        penalty - the type of regularization to use (default = 'l2')
        solver - the algorithm to use in the optimization problem (default = 'lbfgs')
                 we would use penalty='l1', solver='liblinear' for non-het sites mutation classification
        max_iter - the maximum iteration number (default = 1000)
        class_weight - weights associated with classes (default = "balanced")
        smote - whether use SMOTE (Synthetic Minority Over-sampling Technique) to over-sample the minority class 
                to treat the imbalance class values (default = True)
        k_neighbors - the nearest neighbors used to define the neighborhood of samples in SMOTE (default = 4)
        sampling_strategy - sampling information to resample the data set (default = "auto")
                'minority': resample only the minority class;
                'not minority': resample all classes but the minority class;
                'not majority': resample all classes but the majority class;
                'all': resample all classes;
                'auto': equivalent to 'not majority'.
        n_job - number of jobs to run in parallel (default = None) 
                None means 1 unless in a joblib.parallel_backend context. 
                -1 means using all processors.
    """
    # split dataset into training set and test set
    if test_size != 0:
        X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=test_size, stratify=y, random_state=random_state)
    else:
        X_train = X
        y_train = y

    # ======================================================================
    # Over-sampling for the imbalance target values
    # ======================================================================
    # use SMOTE to over-sampling
    if smote:
        smote = SMOTE(random_state=random_state, sampling_strategy=sampling_strategy, k_neighbors=k_neighbors, n_jobs=n_jobs)
        X_train, y_train = smote.fit_resample(X_train, y_train)
        # get the number of each values
        print("Training set count:", dict(Counter(y_train)))

    # normalize the training set
    scaler = StandardScaler()
    X_train = scaler.fit_transform(X_train)
    if test_size != 0:
        X_test = scaler.transform(X_test)
    
    # ======================================================================
    # Train the model
    # ======================================================================
    # initialize the Logistic Regression model
    lr_model = LogisticRegression(random_state=random_state, penalty=penalty, solver=solver,
                                  max_iter=max_iter, class_weight=class_weight)
    # fit model with the training data
    lr_model.fit(X_train, y_train)

    if test_size != 0:
        # test accuracy
        y_pred = lr_model.predict(X_test)
        test_accuracy = accuracy_score(y_test, y_pred)

        # print result
        print(confusion_matrix(y_test, y_pred))
        print(classification_report(y_test, y_pred))
        print(f"Test Accuracy: {test_accuracy}")

    return lr_model



# from sklearn.model_selection import train_test_split
# from sklearn.linear_model import LogisticRegression
# from sklearn.metrics import classification_report, confusion_matrix, accuracy_score
# from imblearn.over_sampling import SMOTE
# from sklearn.preprocessing import StandardScaler

# def multinomial_logistic_regression(X, y, test_size=0.3, random_state=100, max_iter=1000, class_weight='balanced', \
#                                     smote=True, k_neighbors=4, sampling_strategy="auto", n_jobs=None):
#     """
#     Use Multinomial Logistic regression method to refine the phasing result

#     Inputs:
#         X, y - the input data and labels (in the numpy format)
#         test_size - the test set used when splitting the training and test set (default = 0.3)
#         random_state - the random state for the train-test splitting and random forest model (default = 100)
#         max_iter - the maximum iteration number (default = 1000)
#         class_weight - weights associated with classes (default = "balanced")
#         smote - whether use SMOTE (Synthetic Minority Over-sampling Technique) to over-sample the minority class 
#                 to treat the imbalance class values (default = True)
#         k_neighbors - the nearest neighbors used to define the neighborhood of samples in SMOTE (default = 4)
#         sampling_strategy - sampling information to resample the data set (default = "auto")
#                 'minority': resample only the minority class;
#                 'not minority': resample all classes but the minority class;
#                 'not majority': resample all classes but the majority class;
#                 'all': resample all classes;
#                 'auto': equivalent to 'not majority'.
#         n_job - number of jobs to run in parallel (default = None) 
#                 None means 1 unless in a joblib.parallel_backend context. 
#                 -1 means using all processors.
#     """
#     # split dataset into training set and test set
#     if test_size != 0:
#         X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=test_size, stratify=y, random_state=random_state)
#     else:
#         X_train = X
#         y_train = y


#     # ======================================================================
#     # Over-sampling for the imbalance target values
#     # ======================================================================
#     # use SMOTE to over-sampling
#     if smote:
#         smote = SMOTE(random_state=random_state, sampling_strategy=sampling_strategy, k_neighbors=k_neighbors, n_jobs=n_jobs)
#         X_train, y_train = smote.fit_resample(X_train, y_train)
#         # get the number of each values
#         print("Training set count:", dict(Counter(y_train)))

#     # normalize the training set
#     scaler = StandardScaler()
#     X_train = scaler.fit_transform(X_train)
#     if test_size != 0:
#         X_test = scaler.transform(X_test)
    
#     # ======================================================================
#     # Train the model
#     # ======================================================================
#     # initialize the Multinomial Logistic Regression model
#     mlr_model = LogisticRegression(random_state=random_state, max_iter=max_iter, class_weight=class_weight, \
#                                    multi_class='multinomial', solver='lbfgs')
#     # fit model with the training data
#     mlr_model.fit(X_train, y_train)

#     if test_size != 0:
#         # test accuracy
#         y_pred = mlr_model.predict(X_test)
#         test_accuracy = accuracy_score(y_test, y_pred)

#         # print result
#         print(confusion_matrix(y_test, y_pred))
#         print(classification_report(y_test, y_pred))
#         print(f"Test Accuracy: {test_accuracy}")

#     return mlr_model



def lr_model_feature_importance_plot(model, df, output_dir, sample_name, title_note="", file_note="", fig_size=10):
    """Plot the feature importance learnt from logistic regression model"""
    # get feature importances
    importances = np.abs(model.coef_[0])
    feature_names = df.columns
    # sort the features and their importances
    sorted_importances, sorted_features = zip(*sorted(zip(importances, feature_names), reverse=True))
    # plot feature importance bar plot
    plt.figure(figsize=(fig_size+2, fig_size))
    plt.title('Feature Importances' + title_note)
    plt.barh(range(len(sorted_importances)), sorted_importances, align='center')
    plt.yticks(range(len(sorted_importances)), sorted_features)
    plt.xlabel('Relative Importance')
    plt.gca().invert_yaxis()
    plt.tight_layout()
    fig_file = os.path.join(output_dir, sample_name+file_note+"logistic_reg_feature_importance.png")
    plt.savefig(fig_file)
    plt.close()

    return sorted_features


def feature_importamce_plot(model, df, output_dir, sample_name, title_note="", file_note="", fig_size=10):
    """Plot the feature importance learnt from random forest or other classification models"""
    # get feature importances
    importances = model.feature_importances_
    feature_names = df.columns
    # sort the features and their importances
    sorted_importances, sorted_features = zip(*sorted(zip(importances, feature_names), reverse=True))
    # plot feature importance bar plot
    plt.figure(figsize=(fig_size+2, fig_size))
    plt.title('Feature Importances' + title_note)
    plt.barh(range(len(sorted_importances)), sorted_importances, align='center')
    plt.yticks(range(len(sorted_importances)), sorted_features)
    plt.xlabel('Relative Importance')
    plt.gca().invert_yaxis()
    plt.tight_layout()
    fig_file = os.path.join(output_dir, sample_name+file_note+"_feature_importance.png")
    plt.savefig(fig_file)
    plt.close()

    # create a dict for feature importance
    feature_importance_dict = dict(zip(feature_names, importances))
    # check missing values and infinite values
    # feature_importance_dict = {k: (0 if np.isnan(v) else v) for k, v in feature_importance_dict.items()}
    feature_importance_dict = {k: v for k, v in feature_importance_dict.items() if not np.isnan(v)}
    feature_importance_dict = {k: v for k, v in feature_importance_dict.items() if np.isfinite(v)}
    # generate the word cloud
    wordcloud = WordCloud(width=800, height=400, background_color='white').generate_from_frequencies(feature_importance_dict)
    # plot the word cloud
    plt.figure(figsize=(fig_size*2, fig_size))
    plt.imshow(wordcloud, interpolation='bilinear')
    plt.axis('off')
    plt.title('Feature Importance Word Cloud' + title_note)
    feature_wordcloud_file = os.path.join(output_dir, sample_name+file_note+"_features_wordcloud.png")
    plt.savefig(feature_wordcloud_file)
    plt.close()
    return sorted_features


def pca_circle_plot(pca, feature_names, output_dir, sample_name, title_note="", file_note="", lim=1, figsize=7):
    """Plot features in the first two principal components as a circle plot"""
    # get components
    loadings = pca.components_.T
    fig, ax = plt.subplots(figsize=(figsize, figsize))
    # color pallete
    arrow_color = '#4DBAD6'
    text_color = '#E44A33'
    # plot the features
    for i, feature in enumerate(feature_names):
        plt.arrow(0, 0, loadings[i, 0], loadings[i, 1], 
                  color=arrow_color, alpha=0.5)
        plt.text(loadings[i, 0]*1.15, loadings[i, 1]*1.15, 
                 feature, color=text_color, ha='center', va='center')
    # set plot limits
    plt.xlim(-lim, lim)
    plt.ylim(-lim, lim)
    plt.axhline(0, color='gray', linestyle='--')
    plt.axvline(0, color='gray', linestyle='--')
    plt.xlabel('PC1')
    plt.ylabel('PC2')
    plt.title('PCA Componants Plot' + title_note)
    fig_file = os.path.join(output_dir, sample_name+file_note+"_features_circle.png")
    plt.savefig(fig_file)
    plt.close()


def PCA_train(sorted_features, X_df, y, output_dir, sample_name, title_note="", file_note="", 
                 n_features=10, fig_size=6, annotate_mosaic=True, annotate_outlier=False, thr_outlier=99.9):
    """Perform PCA with important features"""
    # get the important feature data frame
    features_list = list(sorted_features)
    if "context" in features_list:
        features_list.remove("context")
        sorted_features = tuple(features_list)
    top_features = list(sorted_features[:n_features])
    X_df_features = X_df[top_features]

    # modify the extreme p-values to a proper number
    p_value_cols = [col for col in X_df_features.columns if col.endswith('_p')]
    X_df_features = X_df_features.copy()  # make an explicit copy of the DataFrame
    # X_df_features.loc[:, p_value_cols] = X_df_features[p_value_cols].applymap(lambda x: -7 if x < -7 else x)
    X_df_features.loc[:, p_value_cols] = X_df_features[p_value_cols].clip(lower=-7)

    # standardize the features
    X_std = StandardScaler().fit_transform(X_df_features)

    # perform PCA
    pca = PCA(n_components=2)  # we choose 2 components for visualization purposes
    X_pca = pca.fit_transform(X_std)
    # plot features
    pca_circle_plot(pca, feature_names=X_df_features.columns, output_dir=output_dir, sample_name=sample_name, 
                    title_note=title_note, file_note=file_note, lim=0.85, figsize=9)

    # # create a data frame for the output
    # pca_df = pd.DataFrame(data = X_pca, columns = ['PC1', 'PC2'], index=X_df_features.index)
    # pca_df['evaluation'] = y
    # add PCA values to the original feature data frame
    df = X_df.copy()
    df['label'] = y
    df['type'] = 'evaluation'
    df[['PC1', 'PC2']] = X_pca
    # get the explained variance ratio
    explained_variance = pca.explained_variance_ratio_
    # calculate the distance from the origin (or centroid)
    distances = np.linalg.norm(X_pca, axis=1)
    # define a threshold (percentile) for outliers
    threshold = np.percentile(distances, thr_outlier)
    # identify outliers
    outliers = distances > threshold
    df['outlier'] = outliers

    # PCA plot
    color_palette = {
        'mosaic': '#fe4a49',
        'artifact': '#009fb7',
        'het': '#FF8C00'
    }
    # if n_labels==2:
    #     color_palette = ['#fe4a49', '#009fb7']
    # elif n_labels==3:
    #     color_palette = ['#fe4a49', '#009fb7', '#ffd166']
    plt.figure(figsize=(fig_size+2, fig_size))
    sns.scatterplot(data=df, x='PC1', y='PC2', hue='label', palette=color_palette, alpha=0.9, s=20)
    # add site info for true mutations
    for i, row in df.iterrows():
        if row['label'] == 'mosaic' and annotate_mosaic: 
            plt.annotate(row.name, 
                         (row['PC1'], row['PC2']),
                         textcoords="offset points",
                         xytext=(0,10),  # offset the label position
                         ha='center',
                         fontsize=9)
        if row['outlier'] and annotate_outlier:
            plt.annotate(row.name, 
                         (row['PC1'], row['PC2']),
                         textcoords="offset points",
                         xytext=(0,10),  # offset the label position
                         ha='center',
                         fontsize=9,
                         color="royalblue")
    plt.title('PCA Scatter Plot'+title_note)
    plt.xlabel(f'PC 1 ({explained_variance[0]*100:.2f}% variance)')
    plt.ylabel(f'PC 2 ({explained_variance[1]*100:.2f}% variance)')
    pca_scatter_plot = os.path.join(output_dir, sample_name+file_note+"_pca_scatter.png")
    plt.savefig(pca_scatter_plot)
    plt.close()

    return pca, top_features, df



def PCA_pred(train_df, pred_df, output_dir, sample_name, title_note="", file_note="", 
                fig_size=6, pca_train=None, top_features=None, annotate_outlier=False, thr_outlier=99.9):
    """Perform PCA with all features and plot for both training and testing data"""
    # combine the data frames
    pred_df.rename(columns={'pred':'label'}, inplace=True)
    pred_df['type'] = 'prediction'
    if train_df is not None:
        train_df.rename(columns={'evaluation':'label'}, inplace=True)
        train_df['type'] = 'evaluation'
        combined_df = pd.concat([train_df, pred_df], ignore_index=False)
    else:
        combined_df = pred_df

    # get the labels and the types
    df = combined_df.drop(columns=['label', 'type'])
    # select the features
    if top_features is not None:
        df = df[top_features]
    # modify the extreme p-values to a proper number
    p_value_cols = [col for col in df.columns if col.endswith('_p')]
    df = df.copy()
    # df.loc[:, p_value_cols] = df[p_value_cols].applymap(lambda x: -7 if x < -7 else x)
    df.loc[:, p_value_cols] = df[p_value_cols].clip(lower=-7)

    # standardize the features
    df_std = StandardScaler().fit_transform(df)
    # perform PCA
    if pca_train is not None:
        pca = pca_train
        df_pca = pca_train.transform(df_std)
    else:
        pca = PCA(n_components=2)  # we choose 2 components for visualization purposes
        df_pca = pca.fit_transform(df_std)
    # add the PCs to the combined dataframe
    combined_df[['PC1', 'PC2']] = df_pca
    # get the explained variance ratio
    explained_variance = pca.explained_variance_ratio_
    # calculate the distance from the origin (or centroid)
    distances = np.linalg.norm(df_pca, axis=1)
    # define a threshold (percentile) for outliers
    threshold = np.percentile(distances, thr_outlier)
    # identify outliers
    outliers = distances > threshold
    combined_df['outlier'] = outliers

    # plot features
    title_note = title_note + " with Prediction"
    file_note = file_note + "_pred"
    pca_circle_plot(pca, feature_names=df.columns, output_dir=output_dir, sample_name=sample_name, 
                    title_note=title_note, file_note=file_note, lim=0.55, figsize=9)
    
    # PCA plot
    color_palette = {
        'mosaic': '#fe4a49',
        'artifact': '#009fb7',
        'het': '#FF8C00'
    }
    # marker_styles = {'evaluation': 'o', 'prediction': '^'}
    # combined_df['marker'] = combined_df['type'].map(marker_styles)
    plt.figure(figsize=(fig_size+2, fig_size))
    if train_df is not None:
        sns.scatterplot(data=combined_df, x='PC1', y='PC2', hue='label', style='type', 
                        palette=color_palette, alpha=0.9, s=20, markers=['o', '^'])
    else:
        sns.scatterplot(data=combined_df, x='PC1', y='PC2', hue='label', style='type', 
                        palette=color_palette, alpha=0.9, s=20, markers='^')
    # add site info for outliers
    if annotate_outlier:
        for i, row in combined_df.iterrows():
            if row['outlier']: 
                plt.annotate(row.name, 
                            (row['PC1'], row['PC2']),
                            textcoords="offset points",
                            xytext=(0,10),  # offset the label position
                            ha='center',
                            fontsize=9,
                            color="royalblue")
    # for marker, subset in combined_df.groupby('marker'):
    #     sns.scatterplot(data=subset, x='PC1', y='PC2', hue='label', palette=color_palette, 
    #                     alpha=0.9, s=20, markers=marker, edgecolor='w')
    plt.title('PCA Scatter Plot'+title_note)
    plt.legend(loc='best')
    plt.xlabel(f'PC 1 ({explained_variance[0]*100:.2f}% variance)')
    plt.ylabel(f'PC 2 ({explained_variance[1]*100:.2f}% variance)')
    pca_scatter_plot = os.path.join(output_dir, sample_name+file_note+"_pca_scatter.png")
    plt.savefig(pca_scatter_plot)
    plt.close()

    return combined_df



def chrom_sort_key(chrom):
    """
    Sort chromosomes in biological order:
    chr1, chr2, ..., chr22, chrX, chrY, chrM/chrMT

    Unknown contigs are placed at the end.
    """
    chrom = str(chrom).strip()
    # remove common prefix
    if chrom.lower().startswith("chr"):
        chrom_core = chrom[3:]
    else:
        chrom_core = chrom

    chrom_upper = chrom_core.upper()
    special_map = {
        "X": 23,
        "Y": 24,
        "M": 25,
        "MT": 25,
    }

    if chrom_upper in special_map:
        return (0, special_map[chrom_upper])
    try:
        return (0, int(chrom_core))
    except ValueError:
        # place non-standard contigs after canonical chromosomes
        return (1, chrom_upper)


def format_info_value(val):
    """Format values safely for the VCF INFO field."""
    if pd.isna(val):
        return None
    if isinstance(val, bool):
        return "TRUE" if val else "FALSE"
    if isinstance(val, float):
        return f"{val:.6g}"
    val = str(val)
    val = val.replace(" ", "_").replace(";", ",")
    return val


def make_vcf_info(row, info_columns):
    """Build the INFO field from selected columns."""
    items = []
    for col in info_columns:
        if col in row.index:
            val = format_info_value(row[col])
            if val is not None:
                items.append(f"{col}={val}")
    return ";".join(items) if items else "."


def write_simple_vcf(df, output_file, sample_name=None, reference_name=None):
    """
    Write a simple site-level VCF from a DataFrame.

    Required columns:
        #CHROM, POS, REF, ALT, FILTER

    Optional INFO columns:
        DP, ALT_DP, VAF, NUM_SPOTS, NUM_MUT_SPOTS,
        DNA_MUT_TYPE, RNA_MUT_TYPE
    """
    # get basic info
    df = df.copy()
    required_cols = ["#CHROM", "POS", "REF", "ALT", "FILTER"]
    missing = [c for c in required_cols if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns for VCF export: {missing}")
    # combine info column
    info_columns = [
        "DP",
        "ALT_DP",
        "VAF",
        "NUM_SPOTS",
        "NUM_MUT_SPOTS",
        "DNA_MUT_TYPE",
        "RNA_MUT_TYPE",
    ]
    # write file
    with open(output_file, "w") as f:
        f.write("##fileformat=VCFv4.2\n")
        if sample_name is not None:
            f.write(f"##source={sample_name}\n")
        if reference_name is not None:
            f.write(f"##reference={reference_name}\n")
        # add chrom contig
        # contigs = sorted(df["#CHROM"].dropna().unique())
        contigs = sorted(df["#CHROM"].dropna().unique(), key=chrom_sort_key)
        for c in contigs:
            f.write(f"##contig=<ID={c}>\n")
        # add other title
        f.write('##INFO=<ID=DP,Number=1,Type=Integer,Description="Total depth">\n')
        f.write('##INFO=<ID=ALT_DP,Number=1,Type=Integer,Description="Mutant allele depth">\n')
        f.write('##INFO=<ID=VAF,Number=1,Type=Float,Description="Variant allele fraction">\n')
        f.write('##INFO=<ID=NUM_SPOTS,Number=1,Type=Integer,Description="Total number of spots">\n')
        f.write('##INFO=<ID=NUM_MUT_SPOTS,Number=1,Type=Integer,Description="Number of mutation-supporting spots">\n')
        f.write('##INFO=<ID=DNA_MUT_TYPE,Number=1,Type=String,Description="DNA mutation type">\n')
        f.write('##INFO=<ID=RNA_MUT_TYPE,Number=1,Type=String,Description="RNA mutation type">\n')
        f.write('##FILTER=<ID=PASS,Description="Passed all internal filters">\n')
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        # mutation sites
        for _, row in df.iterrows():
            chrom = row["#CHROM"]
            pos = int(row["POS"])
            ref = row["REF"]
            alt = row["ALT"]
            filt = row["FILTER"] if pd.notna(row["FILTER"]) else "."
            info = make_vcf_info(row, info_columns)

            f.write(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\t{filt}\t{info}\n")


# # column names corresponding to the features in the model
# old2new_mapping = {}

# # function to rename columns based on the mapping
# def rename_columns(df, name_mapping):
#     for col in df.columns:
#         if col in name_mapping:
#             # check if the new name is different, and if so, rename the column
#             if col != name_mapping[col]:
#                 df = df.rename(columns={col: name_mapping[col]})
#     return df

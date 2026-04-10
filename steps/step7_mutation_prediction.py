from SpaceTracer.cores.mutation_prediction import mutation_classification
import os

from SpaceTracer.steps.base import BaseStep
from SpaceTracer.utils.utils import check_dir, str2bool


def predict_mutation(args):
    """
    Predict somatic mutations with classification methods
    """
    check_dir(args.outdir)

class MutationPredictionStep(BaseStep):
    def get_inputs(self, context):
        inputs={
            'combine_feature': context.get('combine_feature')
        }
        return inputs

    def get_outputs(self, context):
        sample_name="Sample"
        vcf_output_file = os.path.join(self.step_dir, sample_name + "_total_pred_truesites.vcf")
        vcf_pass_output_file = os.path.join(self.step_dir, sample_name + "_total_pred_truesites_PASS.vcf")
        outputs={
            'raw_pred_vcf':vcf_output_file,
            'final_vcf':vcf_pass_output_file
        }
        return outputs

    def get_step_config(self):
        return self.config.get('steps', {}).get('mutation_prediction', {})

    def _run(self,context):
        inputs=self.get_inputs(context)

        parameter=self.get_step_config()
        model_dir=parameter["model_dir"]
        model_name=parameter["model_name"]
        random_seed=42
        train=False
        true_sites_file=None
        validated_artifact_sites_file=None
        phasable_artifact_sites_file=None
        select_features=None
        drop_features=None
        no_spatial=False
        phase_refine=False
        true_to_artifact_ratio=None
        encoder='label'
        save_models=True

        plot=True
        annotate_mosaic=True
        annotate_outlier=False
        n_features=20
        save_pca=True
        save_shap=True

        use_lr=False
        not_pred_het=True
        transform_old_name=False
        smote=True
        tune='random_search'
        k_neighbors=4
        sampling_strategy=False
        n_jobs=None
        n_estimators=100
        max_depth=None
        min_samples_split=2

        mutation_classification(inputs['combine_feature'], 
                self.step_dir, 
                "Sample",
                model_dir=model_dir, 
                model_name=model_name, 
                random_seed=random_seed,
                train=train, 
                true_sites_file=true_sites_file,
                validated_artifact_sites_file=validated_artifact_sites_file, 
                phasable_artifact_sites_file=phasable_artifact_sites_file, 
                select_features=select_features, 
                drop_features=drop_features, 
                no_spatial=no_spatial, 
                phase_refine=phase_refine, 
                true_to_artifact_ratio=true_to_artifact_ratio, 
                encoder=encoder, 
                save_models=save_models, 
                plot=plot, 
                annotate_mosaic=annotate_mosaic, 
                annotate_outlier=annotate_outlier, 
                n_features=n_features, 
                save_pca=save_pca, 
                save_shap=save_shap, 
                use_lr=use_lr, 
                not_pred_het=not_pred_het, 
                transform_old_name=transform_old_name, 
                smote=smote, 
                tune=tune, 
                k_neighbors=k_neighbors, 
                sampling_strategy=sampling_strategy, 
                n_jobs=n_jobs, 
                n_estimators=n_estimators,
                max_depth=max_depth, 
                min_samples_split=min_samples_split)
        


# if __name__ == '__main__':
#     ## parameters
#     parser = argparse.ArgumentParser()
#     subparsers = parser.add_subparsers(help='sub-command help')
#     # predict somatic mutations
#     parser_predict = subparsers.add_parser('predict', help='predict mutations') 
#     parser_predict.add_argument("--input","-i", required=True, help="the directory of features data input")
#     parser_predict.add_argument("--outdir",required=True, help="output dir")
#     parser_predict.add_argument("--sample","-s", required=True, help="sample name")
#     parser_predict.add_argument("--model_dir", required=False, default="./models_trained/spatial_feature_preserved_model", help="the directory of the trained models (default=./models_trained/spatial_feature_preserved_model)")
#     parser_predict.add_argument("--model_name", required=False, default="spatial_feature_preserved_model", help="the sample name of the trained models default=spatial_feature_preserved_model")
#     parser_predict.add_argument("--random_seed", required=False, default=100, type=int, help="the random seed (default=100)")
#     parser_predict.add_argument("--train", required=False, default=True, type=str2bool, choices=[True, False], help="Boolean variable for whether training the models or not (default=True)")
#     parser_predict.add_argument("--true_sites_file", required=False, default=None, help="path to the file with true somatic sites in this data set which are manually checked (default=None)")
#     parser_predict.add_argument("--validated_artifact_sites_file", required=False, default=None, help="path to the file with validated artifact sites in this data set (default=None)")
#     parser_predict.add_argument("--phasable_artifact_sites_file", required=False, default=None, help="path to the file with phasable artifact sites in this data set, if not given then use the sites with 'haplotype>3' (default=None)")
#     parser_predict.add_argument("--select_features", required=False, default=None, help="features used to classificate the true somatic mutations, use all features if 'None' (default=None)")
#     parser_predict.add_argument("--drop_features", required=False, default=None, help="the features do not want to be used to classificate the true somatic mutations (default=None)")
#     parser_predict.add_argument("--no_spatial", required=False, default=False, type=str2bool, choices=[True, False], help="Boolean variable of whether not using the spatial features in the model (default=False)")
#     parser_predict.add_argument("--phase_refine", required=False, default=False, type=str2bool, choices=[True, False], help="Boolean variable for whether use the phasing refinement model (default=False)")
#     parser_predict.add_argument("--true_to_artifact_ratio", required=False, default=None, type=float, help="the ratio of the true somatic mutations to artifacts in the training set (default=None, not control the ratio)")
#     parser_predict.add_argument("--encoder", required=False, default='label', choices=['label', 'onehot'], help="the encode method for object columns ('label' or 'onehot', default='label')")
#     parser_predict.add_argument("--save", required=False, default=True, type=str2bool, choices=[True, False], help="Boolean variable for whether saving the models (default=True)")
#     parser_predict.add_argument("--plot","-p", required=False, default=True, type=str2bool, choices=[True, False], help="Boolean variable for plotting the feature importances and the PCA figures (default=True)")
#     parser_predict.add_argument("--annotate_mosaic", required=False, default=True, type=str2bool, choices=[True, False], help="Boolean variable of whether annotating the mosaic sites in the training set PCA scatter plots (default=True)")
#     parser_predict.add_argument("--annotate_outlier", required=False, default=False, type=str2bool, choices=[True, False], help="Boolean variable of whether annotating the outlier sites in all PCA scatter plots (default=False)")
#     parser_predict.add_argument("--save_pca", required=False, default=True, type=str2bool, choices=[True, False], help="Boolean variable of whether saving the PCA transformed value, together with the clean feature values (default=True)")
#     parser_predict.add_argument("--save_shap", required=False, default=False, type=str2bool, choices=[True, False], help="Boolean variable of whether saving the SHAP values for the features (default=False)")
#     parser_predict.add_argument("--use_lr", required=False, default=False, type=str2bool, choices=[True, False], help="Boolean variable of whether using logistic regression model to classify the somatic mutations (default=False)")
#     parser_predict.add_argument("--not_pred_het", required=False, default=True, type=str2bool, choices=[True, False], help="Boolean variable of whether not predicting the heterozygous sites using logistic regression model (default=True)")
#     parser_predict.add_argument("--transform_old_name", required=False, default=False, type=str2bool, choices=[True, False], help="Boolean variable of whether transforming the column names to new version before processing (default=False)")
#     parser_predict.add_argument("--n_features", required=False, default=20, type=int, help="the number of the most-important features from the random forest model used in the PCA projection (default=20)")
#     parser_predict.add_argument("--smote", required=False, default=True, type=str2bool, choices=[True, False], help="whether use SMOTE to over-sample the minority class to treat the imbalance class values (default = True)")
#     parser_predict.add_argument("--tune", required=False, default='random_search', choices=['Bayesian_opt', 'random_search', 'grid_search', None], \
#                                 help="whether tuning the hyperparameters used in the random forest (choices=['Bayesian_opt', 'random_search', 'grid_search', None], default='random_search')")
#     parser_predict.add_argument("--k_neighbors", required=False, default=4, type=int, help="the nearest neighbors used to define the neighborhood of samples in SMOTE (default=4)")
#     parser_predict.add_argument("--sampling_strategy", required=False, default='auto', choices=['minority', 'not minority', 'not majority', 'all', 'auto'], \
#                                 help="sampling information to resample the data set (choices=['minority', 'not minority', 'not majority', 'all', 'auto'], default='auto')")
#     parser_predict.add_argument("--n_jobs", required=False, default=None, choices=[None, -1], help="number of jobs to run in parallel (default=None)")
#     parser_predict.add_argument("--n_estimators", required=False, default=100, type=int, help="the number of trees in the forest (default=100)")
#     parser_predict.add_argument("--max_depth", required=False, default=None, type=int, help="the maximum depth of the tree (default=None)")
#     parser_predict.add_argument("--min_samples_split", required=False, default=2, type=int, help="the minimum number of samples required to split an internal node (default=2)")    
#     parser_predict.set_defaults(func=predict_mutation)
#     args = parser.parse_args()
#     args.func(args)
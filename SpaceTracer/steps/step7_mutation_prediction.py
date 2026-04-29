from SpaceTracer.cores.mutation_prediction import mutation_classification
import os

from SpaceTracer.steps.base import BaseStep
from SpaceTracer.utils.utils import check_dir, str2bool

from importlib.resources import files
from pathlib import Path



def predict_mutation(args):
    """
    Predict somatic mutations with classification methods
    """
    check_dir(args.outdir)

class MutationPredictionStep(BaseStep):
    def get_inputs(self, context):
        inputs={
            'combine_feature_parquet': context.get('combine_feature_parquet')
        }
        return inputs

    def get_outputs(self, context):
        sample_name="Sample"
        vcf_output_file = os.path.join(self.step_dir, "results",sample_name + "_total_pred_truesites.vcf")
        vcf_pass_output_file = os.path.join(self.step_dir, "results", sample_name + "_total_pred_truesites_PASS.vcf")
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
        # model_dir=parameter["model_dir"]
        # model_name=parameter["model_name"]
        model_name = self.config["model_used"]
        model_dir = files("SpaceTracer").joinpath("models", model_name)

        random_seed=int(parameter['random_seed'])
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

        plot=bool(parameter['plot'])
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

        mutation_classification(inputs['combine_feature_parquet'], 
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
        

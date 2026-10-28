from SpaceTracer.cores.mutation_prediction import mutation_classification
import os

from SpaceTracer.steps.base import BaseStep
from SpaceTracer.utils.utils import check_dir

from SpaceTracer.utils.logger import get_logger

model_name=__name__
logger = get_logger(model_name)

def predict_mutation(args):
    """
    Predict somatic mutations with classification methods
    """
    check_dir(args.outdir)

class MutationPredictionStep(BaseStep):
    def _get_model_names(self):
        model_names = self.config.get("model_name")
        if not model_names:
            raise ValueError("model_name is required for mutation prediction")
        if isinstance(model_names, str):
            return [model_names]
        return list(model_names)

    def _paths_for_model(self, model_name):
        new_name = f"{self.sample}_{model_name}"
        results_dir = os.path.join(self.step_dir, "results")
        return {
            "new_name": new_name,
            "raw_pred_vcf": os.path.join(
                results_dir, new_name + "_total_pred_truesites.vcf"
            ),
            "final_vcf": os.path.join(
                results_dir, new_name + "_total_pred_truesites_PASS.vcf"
            ),
            "final_mutation_list": os.path.join(
                results_dir,
                new_name + "_total_pred_truesites_PASS_mutation_list.txt",
            ),
        }

    def get_inputs(self, context):
        inputs={
            'combine_feature_parquet': context.get('combine_feature_parquet')
        }
        return inputs

    def get_outputs(self, context):
        raw_pred_vcfs = []
        final_vcfs = []
        mutation_lists = []
        for model_name in self._get_model_names():
            paths = self._paths_for_model(model_name)
            raw_pred_vcfs.append(paths["raw_pred_vcf"])
            final_vcfs.append(paths["final_vcf"])
            mutation_lists.append(paths["final_mutation_list"])

        outputs={
            'raw_pred_vcf': raw_pred_vcfs,
            'final_vcf': final_vcfs,
            'final_mutation_list': mutation_lists,
        }
        return outputs

    def get_step_config(self):
        return self.config.get('steps', {}).get('mutation_prediction', {})

    def _run(self,context):
        inputs=self.get_inputs(context)

        parameter=self.get_step_config()
        model_dir = self.config["model_dir"]

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
        save_pca=False
        save_shap=False

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

        model_names = self._get_model_names()
        logger.info(f"Running mutation prediction for models: {model_names}")

        for model_name in model_names:
            paths = self._paths_for_model(model_name)
            logger.info(f"Predicting with model: {model_name}")
            mutation_classification(
                inputs['combine_feature_parquet'],
                self.step_dir,
                paths["new_name"],
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
                min_samples_split=min_samples_split,
            )

# Step 4: Mutation prediction

Script: `4_run_model_predict.py`

**What it does:** classify true somatic mutations versus artifacts using random forest models.

## How to run (template)

```bash
python 4_run_model_predict.py \
  --input <INPUT> \
  --outdir <OUTPUTDIRECTORY> \
  --outprefix <OUTPREFIX> \
  --model_dir <MODELDIRECTORY> \
  --model_name <MODELNAME> \
  --random_seed <RANDOMSEED> \
  --train <TRAIN> \
  --true_sites <TRUESITES> \
  --artifact_sites <ARTIFACTSITES> \
  --thr_altcount <THRALTCOUNT> \
  --thr_altSpotNum <THRALTSPOTNUMBER> \
  --subset <SUBSET> \
  --drop_subset <DROPSUBSET> \
  --hard_filter <HARDFILTER> \
  --phase_refine <PHASEREFINE> \
  --save <SAVE> \
  --plot <PLOT> \
  --n_features <FEATURENUMBER> \
  --tune <TUNE> \
  --k_neighbors <KNEIGHBORS>
```

## Parameters

| Parameter | Description |
| --- | --- |
| `--input` | Feature matrix (e.g. `*.features.add_hFDR.txt`) |
| `--outdir` | Output directory |
| `--outprefix` | Optional (default `sample`) |
| `--model_dir` | Optional (default `./models_trained/tumor_skin_model`) |
| `--model_name` | Optional (default `tumor_skin_model`) |
| `--random_seed` | Optional (default 100) |
| `--train` | Optional (`True`/`False`; default `FALSE`) |
| `--true_sites` | Optional. Confirmed somatic sites |
| `--artifact_sites` | Optional. Known artifact sites |
| `--thr_altcount` | Optional (default 5). Min alt alleles per site |
| `--thr_altSpotNum` | Optional. Min spots with alt (spot-specific mode) |
| `--subset` | Optional. Feature subset; `None` = all |
| `--drop_subset` | Optional. Features to exclude |
| `--hard_filter` | Optional (default `True`) |
| `--phase_refine` | Optional (default `True`) |
| `--save` | Optional (default `True`; training only) |
| `--plot` (`-p`) | Optional (default `True`) |
| `--n_features` | Optional (default 20). Top features for PCA |
| `--tune` | Optional: `Bayesian_opt`, `random_search`, `grid_search`, `None` |
| `--k_neighbors` | Optional (default 4). SMOTE neighbors |

## Example outputs

- `demo/output/predict/results/demo_phased_pred_truesites_filter.txt`
- `demo/output/predict/results/demo_phased_pred_truesites.txt`
- `demo/output/predict/results/demo_pred_truesites_filter.txt`
- `demo/output/predict/results/demo_pred_truesites.txt`
- `demo/output/predict/results/demo_total_pred_truesites.txt`

Main output before Step 5 cleanup:

- `demo/output/predict/results/demo_total_pred_truesites.txt`

## Demo command

```bash
python 4_run_model_predict.py \
  --input demo/output/features_dir/demo.features.add_hFDR.txt \
  --outdir demo/output/predict \
  --outprefix demo \
  --train FALSE \
  --thr_altcount 5 \
  --phase_refine FALSE \
  --save FALSE
```

This step may take about 12 seconds on the demo data.

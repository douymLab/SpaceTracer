# `cluster`

## Purpose

Provides cluster/domain assignments and cell-number information used by downstream genotype modeling.

## Upstream

None (DAG root step).

## Required config and inputs

- `steps.cluster.cluster_file` (optional existing file)
- `steps.cell_number` (integer or file path)
- `spaceranger_dir` (required when SpaceTracer needs to compute clusters internally)
- `sequence_type` (current implementation expects Visium when auto-clustering)

### Input interpretation

| Input/config key | Required | Interpretation |
| --- | --- | --- |
| `steps.cluster.cluster_file` | No | If provided and exists, clustering can be skipped and file is reused. |
| `steps.cell_number` | Conditional | Can be fixed integer/path; otherwise derived during preprocessing workflows. |
| `spaceranger_dir` | Conditional | Required when cluster must be computed from SpaceRanger outputs. |
| `sequence_type` | Yes | Defines data mode and auto-clustering expectations. |

## Parameters

From `steps.cluster`:

- `method`: clustering backend (for example `SpaGCN` or `GraphST`)
- `ncluster`, `init_method`
- `weight_histology`, `spot_area`, `percentage`
- `tol`, `lr`, `max_epochs`
- `distance_threshold`, `num_threshold`, `min_samples`, `radius`
- `graphst_tool`, `seed`

### Parameter interpretation highlights

| Parameter | Interpretation |
| --- | --- |
| `method` | Selects clustering backend (`SpaGCN`, `GraphST`, etc.). |
| `ncluster` | Target cluster count. |
| `weight_histology`, `spot_area`, `percentage` | Histology/spatial weighting controls in clustering objective. |
| `tol`, `lr`, `max_epochs` | Optimization convergence and learning-rate controls. |
| `distance_threshold`, `num_threshold`, `min_samples`, `radius` | Neighborhood density/smoothing behavior controls. |
| `seed` | Reproducibility control for stochastic components. |

## Outputs

Context keys:

- `cluster_file`: either provided file or generated `cluster.txt`
- `cell_num`: integer or generated/per-provided cell number file

Typical files:

- `output_dir/cluster/cluster.txt` (if computed)
- `output_dir/cell_num.txt` (if computed from Visium data)

## Tuning notes

- If `cluster_file` exists, this step can pass it through directly.
- If no cluster file is provided and `sequence_type` is Visium, clustering is computed from `spaceranger_dir`.

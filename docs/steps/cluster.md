# `cluster`

## Purpose

Provides cluster/domain assignments and cell-number information used by downstream genotype modeling.

## Upstream dependencies

None (DAG root step).

## Required config and inputs

- `steps.cluster.cluster_file` (optional existing file)
- `steps.cell_number` (integer or file path)
- `spaceranger_dir` (required when SpaceTracer needs to compute clusters internally)
- `sequence_type` (current implementation expects Visium when auto-clustering)

## Key step parameters

From `steps.cluster`:

- `method`: clustering backend (for example `SpaGCN` or `GraphST`)
- `ncluster`, `init_method`
- `weight_histology`, `spot_area`, `percentage`
- `tol`, `lr`, `max_epochs`
- `distance_threshold`, `num_threshold`, `min_samples`, `radius`
- `graphst_tool`, `seed`

## Outputs

Context keys:

- `cluster_file`: either provided file or generated `cluster.txt`
- `cell_num`: integer or generated/per-provided cell number file

Typical files:

- `output_dir/cluster/cluster.txt` (if computed)
- `output_dir/cell_num.txt` (if computed from Visium data)

## Notes

- If `cluster_file` exists, this step can pass it through directly.
- If no cluster file is provided and `sequence_type` is Visium, clustering is computed from `spaceranger_dir`.

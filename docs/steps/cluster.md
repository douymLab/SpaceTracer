# `cluster`

## Purpose

Provides cluster/domain assignments and cell-number information used by downstream genotype modeling.

## Upstream

None (DAG root step).

## Required config and inputs

- `cluster_file` (optional existing file)
- `cell_number` (integer or file path)
- `spaceranger_dir` (required when SpaceTracer needs to compute clusters internally)
- `sequence_type` (current implementation expects Visium when auto-clustering)

### Input interpretation

| Input/config key | Required | Interpretation |
| --- | --- | --- |
| `cluster_file` | No | If provided and exists, clustering is skipped and the file is reused. |
| `cell_number` | Conditional | Can be a fixed integer or file path; otherwise derived during preprocessing workflows. |
| `spaceranger_dir` | Conditional | Required when cluster must be computed from SpaceRanger outputs. |
| `sequence_type` | Yes | Defines data mode and auto-clustering expectations. |

## Parameters

From `steps.cluster`:

- `method`: clustering backend (`default`, `SpaGCN`, or `GraphST`)
- `ncluster`, `init_method`
- `weight_histology`, `spot_area`, `percentage`
- `tol`, `lr`, `max_epochs`
- `distance_threshold`, `num_threshold`, `min_samples`, `radius`
- `graphst_tool`, `seed`

### Parameter interpretation highlights

| Parameter | Interpretation |
| --- | --- |
| `method` | Selects clustering backend (`default`, `SpaGCN`, or `GraphST`). In the advanced config, set this with top-level `cluster_method`. |
| `ncluster` | Target cluster count. |
| `weight_histology`, `spot_area`, `percentage` | Histology/spatial weighting controls in clustering objective. |
| `tol`, `lr`, `max_epochs` | Optimization convergence and learning-rate controls. |
| `distance_threshold`, `num_threshold`, `min_samples`, `radius` | Neighborhood density/smoothing behavior controls. |
| `seed` | Reproducibility control for stochastic components. |

## Space Ranger 4.0.1 and `default` clustering

For Space Ranger 4.0.1 or cases where `SpaGCN`/`GraphST` cannot run, set:

```yaml
cluster_method: default
```

With `default`, SpaceTracer reads clustering results directly from the Space Ranger output instead of recomputing clusters. For Visium, it expects the graph-based clustering output under:

```text
analysis/clustering/gene_expression_graphclust/clusters.csv
```

SpaceTracer also handles the tissue-position filename change across Space Ranger versions by checking both:

- `spatial/tissue_positions_list.csv`
- `spatial/tissue_positions.csv`

For Visium-HD, `default` clustering uses binned outputs and may require `barcode_mappings.parquet` through `input_details.barcode_mapping`.

If `cluster_method` is `SpaGCN` or `GraphST` and the clustering backend fails for Visium/Visium-HD, SpaceTracer falls back to `default` clustering when possible.

## Outputs

Context keys:

- `cluster_file`: either provided file or generated `cluster.txt`
- `cell_num`: integer or generated/per-provided cell number file

Typical files:

- `output_dir/cluster/cluster.txt` (if computed)
- `output_dir/refined_cell_num.txt` (if computed from Visium data)

## Tuning notes

- If `cluster_file` exists, this step can pass it through directly.
- If no cluster file is provided and `sequence_type` is Visium, clustering is computed from `spaceranger_dir`.
- Use `cluster_method: default` when you want to reuse Space Ranger clustering, especially for Space Ranger 4.0.1 compatibility.

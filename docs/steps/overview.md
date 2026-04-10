# Step Reference Overview

This section documents each SpaceTracer pipeline step in practical terms:

- what the step does
- required upstream dependencies
- key input files and configuration parameters
- output files written to the run directory
- how to rerun/debug that step

## Active DAG Order

The default execution order is:

1. `cluster`
2. `bam_processing`
3. `mpileup`
4. `umi_combine`
5. `cell_num`
6. `prior`
7. `genotyping`
8. `spatial_feature`, `mappability_feature`, `read_feature`, `RNA_feature` (parallel branch)
9. `merge_feature`

`mutation_prediction` exists as an optional implementation step, but it is currently not enabled in the active orchestrator DAG by default.

## At-a-Glance Table

| Step | Depends on | Main output key(s) | Main parameter section |
| --- | --- | --- | --- |
| `cluster` | none | `cluster_file`, `cell_num` | `steps.cluster` |
| `bam_processing` | none | `in_bam`, `in_filter_bam` | `steps.bam_processing` |
| `mpileup` | `bam_processing` | `mpileup_file`, `filter_mpileup_file`, `db_path` | `steps.mpileup` |
| `umi_combine` | `mpileup` | `spot_count_file`, `error_count_file` | internal + `run.threads` |
| `cell_num` | `cluster`, `umi_combine` | `cell_num` | `steps.cell_number` |
| `prior` | `umi_combine` | `prior_file` | resource-driven |
| `genotyping` | `cluster`, `prior`, `cell_num` | `ind_geno_filter_file`, `spot_geno_file`, etc. | `steps.genotyping` |
| `spatial_feature` | `genotyping` | `spatial_feature` | `steps.spatial_feature` |
| `mappability_feature` | `genotyping` | `mappability_feature` | resource-driven |
| `read_feature` | `genotyping` | `read_feature` | `steps.read_feature` |
| `RNA_feature` | `genotyping` | `RNA_feature` | `steps.RNA_feature` |
| `merge_feature` | all feature branches | `combine_feature` | `steps.feature_filtration` |

## Important Conventions

- **Step outputs are written into context keys.** Downstream steps read these keys.
- **Step directories are created under `output_dir/<step_name>/`.**
- **Checkpoint files control resume behavior.** Use `--force` when you want a full recompute.
- **Validation defaults to path existence and non-empty outputs** (unless validation is skipped).

## Rerun Patterns

Run from a specific stage:

```bash
SpaceTracer run --config config.yaml --start-from genotyping
```

Recompute even if checkpoint says complete:

```bash
SpaceTracer run --config config.yaml --start-from genotyping --force
```

For practical recipes (for example, rerun only RNA or spatial branches), see [Single-Step Debug Cookbook](debug-cookbook.md).

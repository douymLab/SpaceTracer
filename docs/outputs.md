# Outputs and Troubleshooting

This page helps you interpret SpaceTracer outputs and debug common run failures.

## Output organization

SpaceTracer writes results under the configured `output_dir`.

Depending on your run scope and checkpoint state, you will see:

- step-level intermediate outputs
- checkpoint metadata for completed/failed steps
- merged feature outputs from `merge_feature`
- logs showing step order and runtime

## How to inspect a finished run

Start with three checks:

1. confirm all expected steps reached completion
2. inspect candidate counts after `genotyping`
3. confirm merged features were produced after `merge_feature`

If one of the feature branches is missing, inspect step dependencies and logs before rerunning.

## Typical failure modes

### Missing file paths during startup

**Symptom**
- run fails immediately with path-not-exist error

**Cause**
- one or more configured files under `input_details` or `resource_details` does not exist

**Fix**
- verify every absolute path in `config.yaml`
- rerun after correcting paths

### Cluster-related failure

**Symptom**
- error around `cluster` step

**Cause**
- `sequence_type` and cluster inputs are inconsistent

**Fix**
- for Visium, provide valid tissue positions and either:
  - provide `steps.cluster.cluster_file`, or
  - allow SpaceTracer to build cluster internally

### Validation failure in intermediate steps

**Symptom**
- run aborts with output validation failure

**Cause**
- step output missing/empty due to upstream data or parameter issue

**Fix**
- inspect the first failed step logs
- rerun from that step after fixing input/parameter issue:

```bash
SpaceTracer run --config config.yaml --start-from <failed_step> --force
```

### Unexpectedly skipped steps

**Symptom**
- some steps are skipped even though you changed parameters

**Cause**
- checkpoint state marks previous outputs as complete

**Fix**
- rerun with `--force` to recompute requested steps

## Debug workflow for stable reruns

Use this order for faster debugging:

1. fix paths and config schema first
2. rerun to the first failing step
3. validate intermediate output size/content
4. rerun downstream steps only

## Reproducibility tips

- keep one config file per analysis run
- archive logs and final merged outputs together
- record command line and commit version for each run

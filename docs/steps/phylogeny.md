# `phylogeny`

## Purpose

Builds a lineage tree from high-confidence mutation calls using the integrated PhyloSOLID workflow.

## Upstream

- `mutation_prediction`

## Code file

- `SpaceTracer/steps/step8_phylosolid.py`

## Runtime step name

- `phylogeny`

## What this step does

- reads the final PASS mutation list from `mutation_prediction`
- builds the mutation-by-spot matrix needed by PhyloSOLID
- runs the bundled PhyloSOLID preprocessing script with `Rscript`
- runs the bundled PhyloSOLID tree-building Python script
- writes logs and the final tree under `output_dir/phylogeny`

## Required inputs

| Input | Source | Required | Interpretation |
| --- | --- | --- | --- |
| `final_mutation_list` | `mutation_prediction` output | Yes | PASS mutation list used for tree building. |
| `in_filter_bam` | `bam_processing` output | Yes | Filtered BAM used to build the mutation matrix. |
| `tissue_position` | resolved from `spaceranger_dir` or `input_details.tissue_position` | Yes for automatic tree building | Barcode/tissue-position file used to map mutations to spatial spots. |
| `cell_info` | top-level `cell_info` or read-feature settings | No | Optional barcode-to-cell annotation file. |

## Skip behavior

The step is skipped when fewer than 3 PASS mutations are available, because there is not enough mutation signal to infer a useful tree.

If `tissue_position` cannot be resolved, the current implementation logs that automatic tree building is not supported for that run instead of building the tree.

## Dependencies

This step calls both Python and R code:

- `Rscript` must be available in the environment.
- The required R/Python dependencies should be installed from the SpaceTracer environment file.
- The bundled PhyloSOLID scripts are under `SpaceTracer/third_party/phylosolid`.

## Outputs

Typical outputs include:

- `output_dir/phylogeny/logs/<sample>_phylosolid_prepare.log`
- `output_dir/phylogeny/logs/<sample>_phylosolid_build_tree.log`
- `output_dir/phylogeny/tree/mutation_integrator/phylo/final_cleaned_M_full_basedPivots.filtered_sites_inferred.tree_scphylo.pdf`

## Running or skipping this step

A normal full run includes `phylogeny` after `mutation_prediction`:

```bash
spacetracer run --config config.yaml
```

To stop after VCF prediction and skip tree building:

```bash
spacetracer run --config config.yaml --stop-at mutation_prediction
```

To rerun only the tree-building step after mutation prediction outputs already exist:

```bash
spacetracer run --config config.yaml --start-from phylogeny --force
```

## Suggested practice

- Keep the PASS mutation list used for tree building with the final tree output.
- Check both PhyloSOLID logs if tree building fails.
- Interpret tree structure together with tissue spatial context and mutation evidence.

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

- reads every PASS mutation list produced by `mutation_prediction`
- writes a sorted unique merge of those lists to `output_dir/phylogeny/<sample>_merged_mutation_list.txt`
- builds a mutation-by-spot matrix for each list and for the merged list
- runs the bundled PhyloSOLID preprocessing script with `Rscript`
- runs the bundled PhyloSOLID tree-building Python script
- writes a separate tree directory for each list that contains at least 3 mutations

## Required inputs

| Input | Source | Required | Interpretation |
| --- | --- | --- | --- |
| `final_mutation_list` | `mutation_prediction` output | Yes | One PASS mutation list per `model_name`. Each list is used for its own tree, and all lists are also merged. |
| `in_filter_bam` | `bam_processing` output | Yes | Filtered BAM used to build the mutation matrix. |
| `tissue_position` | resolved from `spaceranger_dir` or `input_details.tissue_position` | Yes unless `target_barcodes` is provided | Barcode/tissue-position file used to map mutations to spatial spots. |
| `target_barcodes` | `steps.phylogeny.target_barcodes` | No | Optional file containing target barcodes for tree input preparation. |
| `cell_info` | top-level `cell_info` or read-feature settings | No | Optional barcode-to-cell annotation file. |

## Skip behavior

Each mutation list is handled on its own. A list with fewer than 3 mutations does not get a tree. If tree building fails for one list, that error is logged and the remaining lists continue.

If neither `tissue_position` nor `target_barcodes` can be resolved, the current implementation logs that automatic tree building is not supported for that run instead of building the tree.

`sequence_type: visium-HD` is accepted when building the mutation matrix. If `cell_info` is set, barcodes missing from that map are omitted from the matrix.

## Dependencies

This step calls both Python and R code:

- `Rscript` must be available in the environment.
- The required R/Python dependencies should be installed from the SpaceTracer environment file.
- The bundled PhyloSOLID scripts are under `SpaceTracer/third_party/phylosolid`.

## Outputs

For a per-model list, `<prefix>` is `<sample>_<model_name>`. For the merged list, `<prefix>` is `<sample>_merged_mutation_list.txt`.

Typical outputs include:

- `output_dir/phylogeny/<sample>_merged_mutation_list.txt`
- `output_dir/phylogeny/logs/<prefix>_phylosolid_prepare.log`
- `output_dir/phylogeny/logs/<prefix>_phylosolid_build_tree.log`
- `output_dir/phylogeny/<prefix>_tree/04_final_results/phylo/*_inferred.tree_scphylo.pdf`

Inside each `<prefix>_tree/` directory, PhyloSOLID also writes `01_germline_filter/`, `02_scaffold_builder/`, `03_mutation_integrator/`, `cv_search_summary.csv`, and `cv_selection_report.txt`.

While the mutation matrix is built, a mutation is left out when the number of spots with a non-missing genotype is less than or equal to `min_spot_number`. The pipeline default is `0`, which drops mutations with no covered spots. The matrix builder accepts `--min_spot_number` when run directly; a value above `0` also drops mutations present in only a few spots.

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

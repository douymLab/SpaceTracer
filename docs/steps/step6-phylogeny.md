# Step 6: Phylogeny

This optional stage reconstructs lineage relationships from high-confidence mutation calls.

## Input recommendation

Use curated loci after Step 5 filtering to reduce topology distortions from technical artifacts.

## Tooling

SpaceTracer does not enforce a single phylogeny backend. You can integrate external tools (for example project-standard SNV phylogeny packages) using SpaceTracer outputs as input.

## Suggested practice

- keep one reproducible export of loci/cell matrices used for tree building
- validate stability across filtering/stringency settings
- cross-check tree structure with tissue spatial context

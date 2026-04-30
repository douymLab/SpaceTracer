# Post-step: Phylogeny

This optional stage reconstructs lineage relationships from high-confidence mutation calls.

## Input recommendation

Use curated loci after post-step recurrent-artifact filtering to reduce topology distortions from technical artifacts.

### Input interpretation

| Input | Interpretation |
| --- | --- |
| curated mutation loci (post-filtering) | Candidate mutation set used to build lineage relationships. |
| cell/spot mutation matrix | Presence/absence or probabilistic mutation matrix required by external phylogeny tools. |

## Parameter interpretation

Phylogeny parameters are backend-specific (for example tree model, bootstrap count, distance metric, mutation filtering stringency). Document these in your project phylogeny script to keep lineage results reproducible.

## Tooling

SpaceTracer does not enforce a single phylogeny backend. You can integrate external tools (for example project-standard SNV phylogeny packages) using SpaceTracer outputs as input.

## Suggested practice

- keep one reproducible export of loci/cell matrices used for tree building
- validate stability across filtering/stringency settings
- cross-check tree structure with tissue spatial context

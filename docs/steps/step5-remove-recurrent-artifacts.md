# Step 5: Remove recurrent artifacts

This stage is a recommended post-processing layer to suppress recurrent technical artifacts before lineage interpretation.

## Typical filters

- population/common polymorphism resources
- panel-of-normals (PoN) style exclusion lists
- recurrent-position blacklists from internal cohorts

## Notes

- This step is often implemented as a project-specific filtering script after SpaceTracer core outputs are generated.
- Use the Step 3/4 outputs as input evidence and record all filter criteria for reproducibility.

## Related pages

- [Outputs](../outputs.md)
- [Single-Step Debug Cookbook](debug-cookbook.md)

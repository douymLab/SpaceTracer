# Step 5: Feature Extraction and Phasing

## Code files

- `SpaceTracer/steps/step5_spatial_feature.py`
- `SpaceTracer/steps/step5_mappability_feature.py`
- `SpaceTracer/steps/step5_read_feature.py`
- `SpaceTracer/steps/step5_RNA_level_feature.py`
- `SpaceTracer/steps/step5_phasing.py`

## Runtime step names

- `spatial_feature`
- `mappability_feature`
- `read_feature`
- `RNA_feature`
- `phasing`

## What this step does

Step 5 includes multiple parallel/related feature modules:

- `spatial_feature`: spatial-neighborhood evidence
- `mappability_feature`: mappability/confounder evidence
- `read_feature`: read-level quality/bias evidence
- `RNA_feature`: RNA-level context evidence
- `phasing`: phasing refinement before final merge/filter

## Detailed reference pages

- [spatial_feature](spatial-feature.md)
- [mappability_feature](mappability-feature.md)
- [read_feature](read-feature.md)
- [RNA_feature](rna-feature.md)
- [phasing](phasing.md)

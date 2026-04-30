# Step 4: Genotyping

## Code file

- `SpaceTracer/steps/step4_genotyping.py`

## Runtime step name

- `genotyping`

## What this step does

- combines Step 3 outputs (`umi_combine`, `cell_num`, `prior`)
- infers genotype evidence at individual/cluster/spot levels
- outputs genotype tables used by Step 5 feature branches

## Detailed reference page

- [genotyping](genotyping.md)

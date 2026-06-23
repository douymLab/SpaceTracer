import subprocess
import pandas as pd
import os
from pathlib import Path
import tempfile


def add_hFDR(df, error_profile_file, reference_error_profile, step_dir):
    sample_id = "Sample"
    rscript_path = str(Path(__file__).with_suffix('.R'))
    input_file = os.path.join(step_dir, "input.tsv")
    
    df_to_export = df.copy()
    df_to_export['Mutation_ID'] = df.index  
    df_to_export.to_csv(input_file, sep="\t", index=False)

    subprocess.run([
        "Rscript", rscript_path,
        sample_id, input_file, step_dir,
        error_profile_file, reference_error_profile
    ], check=True)

    result = pd.read_csv(
        os.path.join(step_dir, "features_with_hFDR.txt"),
        sep="\t"
    )
    
    required_cols = ["#chrom", "pos", "ref", "alt"]
    if not isinstance(result.index, pd.MultiIndex):
        if all(col in result.columns for col in required_cols):
            result = result.set_index(required_cols, drop=True)

    if len(df) != len(result):
        raise ValueError(f"(calculated hFDR df length:{len(result)}) is not equal to the length of input file:({len(df)}).")
    
    return result["refine_hFDR"].set_axis(df.index)
import subprocess
import pandas as pd
import os
from pathlib import Path
import tempfile


def add_hFDR(df, error_profile_file, reference_error_profile,step_dir):
    sample_id="Sample"
    rscript_path = str(Path(__file__).with_suffix('.R'))
    # def enrich(df: pd.DataFrame) -> pd.DataFrame:
    # with tempfile.TemporaryDirectory() as tmpdir:
    input_file = os.path.join(step_dir, "input.tsv")
    df.to_csv(input_file, sep="\t", index=False)

    subprocess.run([
        "Rscript", rscript_path,
        sample_id, input_file, step_dir,
        error_profile_file, reference_error_profile
    ], check=True)

    result = pd.read_csv(
        os.path.join(step_dir, "features_with_hFDR.txt"),
        sep="\t"
    )
    # keep_cols = [c for c in result.columns if c in df.columns] + ["refine_hFDR"]
    # short_result = result[keep_cols]
    # return short_result
    return result["refine_hFDR"].set_axis(df.index)


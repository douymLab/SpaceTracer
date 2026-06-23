import os
import pandas as pd

from SpaceTracer.cores.phasing_summary import build_phase_summary_df
from SpaceTracer.steps.base import BaseStep
from SpaceTracer.utils.read_files import load_parquet, load_text_file
from SpaceTracer.utils.utils import list2min, load_manifest_tsv
from SpaceTracer.utils.logger import get_logger

model_name=__name__
logger = get_logger(model_name)

def collect_files_from_manifest(manifest_file: str, column_name: str):
    rows = load_manifest_tsv(manifest_file)
    files = []
    for row in rows:
        if not row:
            continue
        f = row.get(column_name, "")
        if f:
            files.append(f)
    return files


def _ensure_variant_multiindex(df: pd.DataFrame) -> pd.DataFrame:
    """
    set index and remove (#chrom, pos, ref, alt)
    
    """
    if df is None:
        return pd.DataFrame()
        
    if df.empty:
        return df

    required_cols = ["#chrom", "pos", "ref", "alt"]

    df = df.copy()

    if not isinstance(df.index, pd.MultiIndex):
        if all(col in df.columns for col in required_cols):
            df = df.set_index(required_cols, drop=True)

    return df

def load_single_feature_file(path, sep="\t"):
    if not path:
        return pd.DataFrame()

    parquet_path = str(path).replace(".txt", ".parquet")

    df = None
    if os.path.exists(parquet_path):
        df = load_parquet(parquet_path)
    elif os.path.exists(path):
        df = load_text_file(path, sep=sep, header=0)
    else:
        return pd.DataFrame()

    if df is None or df.empty:
        return pd.DataFrame()

    df = _ensure_variant_multiindex(df)
    return df


def _merge_files(valid_files, output_file):
    if not valid_files:
        with open(output_file, "w"):
            pass
        return

    header_written = False

    with open(output_file, "w") as out_f:
        for file_path in valid_files:
            with open(file_path, "r") as in_f:
                for line_no, line in enumerate(in_f):
                    if not line.strip():
                        continue

                    if not header_written:
                        out_f.write(line.rstrip("\n") + "\n")
                        header_written = True
                        continue

                    if line_no == 0:
                        continue

                    out_f.write(line.rstrip("\n") + "\n")


def merge_feature_files_from_manifest(
    manifest_file: str,
    manifest_column: str,
    output_file: str,
    sep="\t",
    return_df=True,
):
   
    files = collect_files_from_manifest(manifest_file, manifest_column)

    valid_files = []
    for f in files:
        if not f or not os.path.exists(f):
            continue
        if os.path.getsize(f) == 0:
            print(f"Warning: Skipping empty file: {f}")
            continue
        valid_files.append(f)

    out_dir = os.path.dirname(output_file)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    parquet_file = output_file.replace(".txt", ".parquet")


    if not valid_files:
        empty_df = pd.DataFrame()
        empty_df.to_csv(output_file, sep=sep, index=False)
        try:
            empty_df.to_parquet(parquet_file, index=False)
        except Exception:
            pass
        return empty_df

    _merge_files(valid_files, output_file)

    merged_df = pd.read_csv(output_file, sep=sep)
    if merged_df is None or merged_df.empty:
        empty_df = pd.DataFrame()
        try:
            empty_df.to_parquet(parquet_file, index=False)
        except Exception:
            pass
        return empty_df

    merged_df = _ensure_variant_multiindex(merged_df)

    # save parquet 
    try:
        merged_df.to_parquet(parquet_file, index=False, engine="pyarrow", compression="snappy")
    except Exception as e:
        print(f"Warning: failed to write parquet cache: {e}")

    if return_df:
        return merged_df
    return pd.DataFrame()


class MergeFeatureStep(BaseStep):
    def get_inputs(self, context):
        return {
            "RNA_feature": context.get("RNA_feature", ""),
            "mappability_feature": context.get("mappability_feature", ""),
            "spatial_feature_results": (
                context.get("spatial_feature_results", "")
                or context.get("spatial_feature", "")
            ),
            "read_feature_results": (
                context.get("read_feature_results", "")
                or context.get("read_feature", "")
            ),
            "phasing_result": context.get("phasing_result"),
            "cluster_event_result": context.get("cluster_event_result")
        }

    def get_outputs(self, context):
        return {
            "combine_feature": os.path.join(self.work_dir, "all_feature.txt"),
            "combine_feature_parquet": os.path.join(self.work_dir, "all_feature.parquet"),
            "merged_spatial_feature": os.path.join(self.step_dir, "merged_spatial_feature.txt"),
            "merged_read_feature": os.path.join(self.step_dir, "merged_read_feature.txt"),
        }

    def get_step_config(self):
        return self.config.get("steps", {}).get("feature_filtration", {})

    def _load_single_feature(self, path):
        return load_single_feature_file(path, sep="\t")

    def _run(self, context):
        # auto_chrom=self.genome_details['chromosomes']['autosomes']
        # sex_chrom=self.genome_details['chromosomes']['sex_chromosomes']
        mitochondrial=self.genome_details['chromosomes']['mitochondrial']
        contigs=self.genome_details['chromosomes']['contigs']

        inputs = self.get_inputs(context)
        outputs = self.get_outputs(context)

        # 1. RNA_feature
        RNA_feature_df = self._load_single_feature(inputs["RNA_feature"]) \
            if inputs.get("RNA_feature") else pd.DataFrame()

        # 2. mappability_feature
        mappability_feature_df = self._load_single_feature(inputs["mappability_feature"]) \
            if inputs.get("mappability_feature") else pd.DataFrame()

        # 3. spatial_feature_results
        spatial_feature_df = pd.DataFrame()
        spatial_input = inputs.get("spatial_feature_results", "")
        print("******spatial_input",spatial_input)
        if spatial_input:
            if str(spatial_input).endswith(".tsv"):
                rows = load_manifest_tsv(spatial_input)
                if rows:
                    spatial_feature_df = merge_feature_files_from_manifest(
                        manifest_file=spatial_input,
                        manifest_column="spatial_feature_txt",
                        output_file=outputs["merged_spatial_feature"],
                        sep="\t",
                        return_df=True,
                    )
            else:
                logger.info("The spatial_feature_df is empty! Please check!")
                spatial_feature_df = self._load_single_feature(spatial_input)
        # 4. read_feature_results
        read_feature_df = pd.DataFrame()
        read_input = inputs.get("read_feature_results", "")
        if read_input:
            if str(read_input).endswith(".tsv"):
                rows = load_manifest_tsv(read_input)
                if rows:
                    read_feature_df = merge_feature_files_from_manifest(
                        manifest_file=read_input,
                        manifest_column="read_feature_txt",
                        output_file=outputs["merged_read_feature"],
                        sep="\t",
                        return_df=True,
                    )
            else:
                read_feature_df = self._load_single_feature(read_input)

        # 5. phasing_results
        phasing_result=inputs.get("phasing_result", "")
        if phasing_result:
            phase_df = pd.read_csv(phasing_result, sep="\t")
            phase_summary_df=build_phase_summary_df(phase_df)
        else:
            phase_summary_df=build_phase_summary_df(pd.DataFrame())
        
        phase_summary_df=phase_summary_df.rename(columns={
            "germline": "ref",
            "mutant": "alt",
        })

        #  6. cluster event
        cluster_event_result=inputs.get("cluster_event_result", "")
        if cluster_event_result:
            cluster_event_df = pd.read_csv(cluster_event_result, sep="\t")
            cluster_event_df = cluster_event_df.rename(columns={
                "#chr": "#chrom",
                "germline": "ref",
                "mutant": "alt",
            })

            key = ["#chrom", "pos", "ref", "alt"]
            if all(col in cluster_event_df.columns for col in key):
                dup_n = cluster_event_df.duplicated(subset=key).sum()
                if dup_n > 0:
                    # logger.warning(f"cluster_event_df has {dup_n} duplicated rows by variant key, dropping duplicates.")
                    cluster_event_df = cluster_event_df.drop_duplicates(subset=key, keep="first")
        else:
            cluster_event_df = pd.DataFrame()

        # reformat the index for each df
        RNA_feature_df = _ensure_variant_multiindex(RNA_feature_df)
        spatial_feature_df = _ensure_variant_multiindex(spatial_feature_df)
        read_feature_df = _ensure_variant_multiindex(read_feature_df)
        mappability_feature_df = _ensure_variant_multiindex(mappability_feature_df)
        phase_summary_df = _ensure_variant_multiindex(phase_summary_df)
        cluster_event_df = _ensure_variant_multiindex(cluster_event_df)


        output_feature = outputs["combine_feature"]

        short_spatial_feature_df = spatial_feature_df[[
            "pass_spatial_test",
            "mut_spots_vaf_mean", "mut_spots_vaf_median", "all_spots_vaf_mean",
            "num_spots", "num_mut_spots",
            "alt_vs_total_dp_r2", "alt_vs_total_dp_paired_wilcoxon_p", "alt_vs_total_dp_paired_wilcoxon_rbc",
            "mut_spots_prop_by_probablity", "mut_spots_prop_by_vaf",
            "mut_vs_nonmut_spots_KS_p", "mut_vs_nonmut_spots_KS_s",
            "mut_vs_nonmut_spots_MI_p", "mut_vs_nonmut_spots_MI_s"
        ]] if not spatial_feature_df.empty else pd.DataFrame()
        # combine by spatial feature
        if not short_spatial_feature_df.empty:
            merged_df = short_spatial_feature_df.copy()
        else:
            merged_df = pd.DataFrame()

        for name, df in [
            ("RNA_feature", RNA_feature_df),
            ("read_feature", read_feature_df),
            ("mappability_feature", mappability_feature_df),
            ("phasing", phase_summary_df)
        ]:
            if df is None:
                continue
            try:
                if not df.empty:
                    merged_df = merged_df.join(df, how="left")
                        
                else:
                    for col in df.columns:
                        if col not in merged_df.columns:
                            merged_df[col] = pd.NA
                
            except:
                print("******",name)
                print("******",df)
                print("******",merged_df)

        if not cluster_event_df.empty:
            merged_df = merged_df.join(cluster_event_df, how="left")
            merged_df["cluster_event"] = (
                merged_df["cluster_event"]
                .astype("boolean")
                .fillna(False)
                .astype(bool)
            )
        else:
            merged_df["cluster_event"] = False


        if "mappabilityScore" in merged_df.columns:
            merged_df["mappabilityScore"] = merged_df["mappabilityScore"].apply(list2min)

        if "falt" in merged_df.columns:
            merged_df["falt"] = merged_df["falt"].astype(float)

        thr_altcount = 3
        thr_popAF = 1e-4
        min_vaf = 0.01
        max_vaf = 0.5
        min_spot_num = 30
        if self.seq_type=="visium":
            min_alt_spot_num = 1
        else:
            min_alt_spot_num = 6
        
        filtrations_dict = self.get_step_config()

        CONDITION_GROUPS = {
            "ASE": ["ASE"],
            "hFDR": ["hFDR"],
            "imprinted": ["imprinted"],
            "homopolymer": ["homopolymer", "cause_poly_alt"],
            "PON": ["PON"],
            "RNA_editing": ["RNA_editing"],
            "ABNORMAL_MISMATCHES": ["HIGH_REF_MISMATCH", "HIGH_ALT_MISMATCH", "MISMATCH_BIAS"],
            "LOW_READ_DIVERSITY": ["LOW_READ_DIVERSITY"],
            "HIGH_MULTIPLE_MAPPING": ["HIGH_MULTIPLE_MAPPING"],
            "WIDE_DISTRIBUTION": ["HIGH_MUT_PROB"],
            "NEAR_READ_END": ["NEAR_READ_END1", "NEAR_READ_END2"],
            "LOW_MAPQ": ["LOW_MAPQ"],
            "LOW_BASEQ": ["LOW_BASEQ"],
            # "MAPPABILITY": ["MAPPABILITY"],
            "INDEL_PROPORTION": ["INDEL_PROPORTION"],
            "ALT_ALLELE_COUNT": ["ALT_ALLELE_COUNT"],
            "POPULATION_AF": ["POPULATION_AF"],
            "CLUSTER_EVENTS": ["CLUSTER_EVENTS","OTHER_CLUSTERED_NOISE"],
            "CONTIG": ["CONTIG"],
            "MITOCHONDRIA": ["MITOCHONDRIA"],
            "LOW_VAF": ["LOW_VAF"],
            "HIGH_VAF": ["HIGH_VAF"],
            "LOW_SPOT_NUM": ["LOW_SPOT_NUM"],
            "LOW_ALT_SPOT_NUM": ["LOW_ALT_SPOT_NUM"],
            "CLUSTERED_NOISE(RNA_editing)": ["CLUSTERED_NOISE(RNA_editing)"],
            "CLUSTERED_NOISE": ["CLUSTERED_NOISE"]
        }

        ALL_FILTER_CONDITIONS = {
            "LOW_UMI_CONSISTENCE": lambda df: pd.to_numeric(df["alt_UMI_avg_consistence"], errors="coerce") < 0.8,
            "HIGH_ALT_MISMATCH": lambda df: pd.to_numeric(df["alt_mismatches_mean"], errors="coerce") > 1.5,
            "HIGH_REF_MISMATCH": lambda df: pd.to_numeric(df["ref_mismatches_mean"], errors="coerce") > 1.5,
            "MISMATCH_BIAS": lambda df: (
                pd.to_numeric(df["alt_mismatches_mean"], errors="coerce") > 0.9
            ) & (
                pd.to_numeric(df["mismatches_p"], errors="coerce") < 0.01
            ),
            "LOW_READ_DIVERSITY": lambda df: pd.to_numeric(df["alt_querypos_num"], errors="coerce") <= 1,
            "HIGH_MULTIPLE_MAPPING": lambda df: pd.to_numeric(df["alt_multi_map_prop"], errors="coerce") > 0.2,
            "HIGH_MUT_PROB": lambda df: pd.to_numeric(df["mut_spots_prop_by_probablity"], errors="coerce") > 0.5,
            "NEAR_READ_END1": lambda df: (
                pd.to_numeric(df["per_alt_UMI_end_remove_clip_median"], errors="coerce") < 0.05
            ) | (
                pd.to_numeric(df["per_alt_UMI_end_remove_clip_median"], errors="coerce") > 0.95
            ),
            "NEAR_READ_END2": lambda df: (
                (pd.to_numeric(df["per_alt_UMI_end_remove_clip_median"], errors="coerce") < 0.1) |
                (pd.to_numeric(df["per_alt_UMI_end_remove_clip_median"], errors="coerce") > 0.9)
            ) & (
                pd.to_numeric(df["per_UMI_end_remove_clip_p"], errors="coerce") < 0.01
            ),
            "LOW_MAPQ": lambda df: pd.to_numeric(df["mapq255_prop"], errors="coerce") < 0.8,
            "LOW_BASEQ": lambda df: pd.to_numeric(df["alt_baseq_mean"], errors="coerce") < 20,
            "ASE": lambda df: df["ASE"] == True,
            "hFDR": lambda df: pd.to_numeric(df["hFDR"], errors="coerce") > 0.8,
            "imprinted": lambda df: df["imprinted"] == True,
            "cause_poly_alt": lambda df: df["cause_poly_alt"] == True,
            "homopolymer": lambda df: df["homopolymer"].notna(),
            "PON": lambda df: df["PON"] == True,
            "RNA_editing": lambda df: df["RNA_editing"] == True,
            # "RNA_editing": lambda df: df["editing_database"] == True,
            "CLUSTERED_NOISE(RNA_editing)": lambda df: df["AtoG_clustered_noise"] == True,
            "CLUSTERED_NOISE": lambda df: df["other_clustered_noise"] == True,
            # "OTHER_CLUSTERED_NOISE": lambda df: df["other_clustered_noise"] == True,
            # "MAPPABILITY": lambda df: df["mappabilityScore"] == 0,
            "INDEL_PROPORTION": lambda df: df["indel_proportion_for_site"] > 0.05,
            "ALT_ALLELE_COUNT": lambda df: df["consensus_alt_allele_count"] < thr_altcount,
            "POPULATION_AF": lambda df: df["falt"] > thr_popAF,
            "CLUSTER_EVENTS": lambda df: df["cluster_event"] == True,
            "CONTIG": lambda df: pd.Series(df.index.get_level_values("#chrom").isin(contigs),index=df.index),
            "MITOCHONDRIA": lambda df: pd.Series(df.index.get_level_values("#chrom").isin(mitochondrial),index=df.index),
            "LOW_VAF": lambda df: df["AFind"] < min_vaf,
            "HIGH_VAF": lambda df: df["AFind"] > max_vaf,
            "LOW_SPOT_NUM": lambda df: df["num_spots"] < min_spot_num,
            "LOW_ALT_SPOT_NUM": lambda df: df["num_mut_spots"] < min_alt_spot_num,
        }

        enabled_groups = [k for k, v in filtrations_dict.items() if v in [True,'true','TRUE','True',1]]

        if not enabled_groups:
            enabled_groups = list(CONDITION_GROUPS.keys())
        
        merged_df["Filtration"] = compute_group_filtration(
            merged_df,
            CONDITION_GROUPS,
            ALL_FILTER_CONDITIONS,
            enabled_groups
        )

        if not merged_df.empty:
            parquet_file = str(output_feature).replace(".txt", ".parquet")
            merged_df.to_parquet(parquet_file, index=True, engine="pyarrow", compression="snappy")

        # object_cols = merged_df.select_dtypes(include=["object"]).columns
        for col in merged_df.columns:
            merged_df[col] = merged_df[col].fillna("no")
            merged_df[col] = merged_df[col].replace("", "no")

        # Filtration fill
        if "Filtration" in merged_df.columns:
            merged_df["Filtration"] = merged_df["Filtration"].fillna("PASS").replace("", "PASS")

        merged_df.to_csv(output_feature, sep="\t", index=True)

        return {
            "combine_feature": output_feature
        }


def compute_group_filtration(df, condition_groups, all_filter_conditions, enabled_groups):
    group_masks = {}
    for group_key in enabled_groups:
        if group_key not in condition_groups:
            continue

        sub_keys = condition_groups[group_key]
        group_mask = None

        for sub_key in sub_keys:
            if sub_key not in all_filter_conditions:
                continue

            try:
                mask = all_filter_conditions[sub_key](df).fillna(False)
            except Exception:
                mask = pd.Series(False, index=df.index)

            group_mask = mask if group_mask is None else (group_mask | mask)

        if group_mask is not None:
            group_masks[group_key] = group_mask

    if not group_masks:
        return pd.Series(["PASS"] * len(df), index=df.index)

    flag_df = pd.DataFrame(group_masks, index=df.index)

    def _join_flags(row):
        hit = [name for name in flag_df.columns if row[name]]
        return ";".join(hit) if hit else "PASS"

    return flag_df.apply(_join_flags, axis=1)

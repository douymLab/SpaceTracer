import os
import pandas as pd

from SpaceTracer.steps.base import BaseStep
from SpaceTracer.utils.read_files import load_parquet, load_text_file
from SpaceTracer.utils.utils import list2min, load_manifest_tsv


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


def load_single_feature_file(path, sep="	"):
    """
    读取单个 feature 文件：
    - 优先读取同名 parquet
    - 否则读取 txt/csv
    并统一设置 MultiIndex
    """
    if not path:
        return pd.DataFrame()

    parquet_path = str(path).replace(".txt", ".parquet")
    if os.path.exists(parquet_path):
        df = load_parquet(parquet_path)
        return df

    if not os.path.exists(path):
        return pd.DataFrame()

    df = load_text_file(path, sep=sep, header=0)
    if df is None or df.empty:
        return pd.DataFrame()

    if not isinstance(df.index, pd.MultiIndex):
        required_cols = ["#chrom", "pos", "ref", "alt"]
        if all(col in df.columns for col in required_cols):
            df.index = pd.MultiIndex.from_arrays(
                [df["#chrom"], df["pos"], df["ref"], df["alt"]],
                names=["#chrom", "pos", "ref", "alt"]
            )
            # 保留原列也可以，不强制 drop
    return df


def merge_feature_files_from_manifest(
    manifest_file: str,
    manifest_column: str,
    output_file: str,
    sep="	"
):
    """
    从 manifest 中收集 feature 文件，逐个读取后按行合并。
    输出保存为 txt，同时如果需要也可写 parquet。
    """
    files = collect_files_from_manifest(manifest_file, manifest_column)
    files = [f for f in files if f]

    out_dir = os.path.dirname(output_file)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    if not files:
        empty_df = pd.DataFrame()
        empty_df.to_csv(output_file, sep=sep, index=False)
        return empty_df

    dfs = []
    for f in files:
        df = load_single_feature_file(f, sep=sep)
        if df is not None and not df.empty:
            dfs.append(df)

    if not dfs:
        empty_df = pd.DataFrame()
        empty_df.to_csv(output_file, sep=sep, index=False)
        return empty_df

    merged_df = pd.concat(dfs, axis=0)

    # 保存 txt
    merged_df.to_csv(output_file, sep=sep, index=True)

    # 如果是 MultiIndex，则顺手写 parquet 方便后续调试/复用
    parquet_file = str(output_file).replace(".txt", ".parquet")
    try:
        merged_df.to_parquet(parquet_file, index=True)
    except Exception:
        pass

    return merged_df


class MergeFeatureStep(BaseStep):
    def get_inputs(self, context):
        return {
            # 这两个是已经合并好的 sample-level 文件
            "RNA_feature": context.get("RNA_feature", ""),
            "mappability_feature": context.get("mappability_feature", ""),

            # 这两个现在是 chunk manifest
            "spatial_feature_results": (
                context.get("spatial_feature_results", "")
                or context.get("spatial_feature", "")
            ),
            "read_feature_results": (
                context.get("read_feature_results", "")
                or context.get("read_feature", "")
            ),
        }

    def get_outputs(self, context):
        return {
            "combine_feature": os.path.join(self.work_dir, "all_feature.txt"),
            "merged_spatial_feature": os.path.join(self.step_dir, "merged_spatial_feature.txt"),
            "merged_read_feature": os.path.join(self.step_dir, "merged_read_feature.txt"),
        }

    def get_step_config(self):
        return self.config.get("steps", {}).get("feature_filtration", {})

    def _load_single_feature(self, path):
        return load_single_feature_file(path, sep="	")

    def _run(self, context):
        inputs = self.get_inputs(context)
        outputs = self.get_outputs(context)

        # 1. RNA_feature: 已经是总文件，直接读
        RNA_feature_df = self._load_single_feature(inputs["RNA_feature"]) \
            if inputs.get("RNA_feature") else pd.DataFrame()

        # 2. mappability_feature: 已经是总文件，直接读
        mappability_feature_df = self._load_single_feature(inputs["mappability_feature"]) \
            if inputs.get("mappability_feature") else pd.DataFrame()

        # 3. spatial_feature_results: manifest，需要先合并
        spatial_feature_df = pd.DataFrame()
        spatial_input = inputs.get("spatial_feature_results", "")
        if spatial_input:
            if str(spatial_input).endswith(".tsv"):
                # 优先尝试 manifest 中的 spatial_feature 列
                # 如果你保存的是 spatial_feature_parquet 列，可以改成对应列名
                rows = load_manifest_tsv(spatial_input)
                spatial_column = None
                if rows:
                    sample_row = rows[0]
                    if "spatial_feature" in sample_row:
                        spatial_column = "spatial_feature"
                    elif "spatial_feature_parquet" in sample_row:
                        spatial_column = "spatial_feature_parquet"

                if spatial_column:
                    spatial_feature_df = merge_feature_files_from_manifest(
                        manifest_file=spatial_input,
                        manifest_column=spatial_column,
                        output_file=outputs["merged_spatial_feature"],
                        sep="	"
                    )
            else:
                spatial_feature_df = self._load_single_feature(spatial_input)

        # 4. read_feature_results: manifest，需要先合并
        read_feature_df = pd.DataFrame()
        read_input = inputs.get("read_feature_results", "")
        if read_input:
            if str(read_input).endswith(".tsv"):
                rows = load_manifest_tsv(read_input)
                read_column = None
                if rows:
                    sample_row = rows[0]
                    if "read_feature" in sample_row:
                        read_column = "read_feature"
                    elif "read_feature_parquet" in sample_row:
                        read_column = "read_feature_parquet"

                if read_column:
                    read_feature_df = merge_feature_files_from_manifest(
                        manifest_file=read_input,
                        manifest_column=read_column,
                        output_file=outputs["merged_read_feature"],
                        sep="	"
                    )
            else:
                read_feature_df = self._load_single_feature(read_input)

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

        dfs_to_merge = [
            df for df in [
                RNA_feature_df,
                short_spatial_feature_df,
                read_feature_df,
                mappability_feature_df
            ]
            if df is not None and not df.empty
        ]

        if dfs_to_merge:
            merged_df = pd.concat(dfs_to_merge, axis=1, join="outer")
        else:
            merged_df = pd.DataFrame()

        if "mappabilityScore" in merged_df.columns:
            merged_df["mappabilityScore"] = merged_df["mappabilityScore"].apply(list2min)

        thr_altcount = 3
        thr_popAF = 1e-4

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
            "MAPPABILITY": ["MAPPABILITY"],
            "INDEL_PROPORTION": ["INDEL_PROPORTION"],
            "ALT_ALLELE_COUNT": ["ALT_ALLELE_COUNT"],
            "POPULATION_AF": ["POPULATION_AF"],
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
            "ASE": lambda df: df["ASE"] == True,
            "hFDR": lambda df: pd.to_numeric(df["hFDR"], errors="coerce") > 0.8,
            "imprinted": lambda df: df["imprinted"] == True,
            "cause_poly_alt": lambda df: df["cause_poly_alt"] == True,
            "homopolymer": lambda df: df["homopolymer"].notna(),
            "PON": lambda df: df["PON"] == True,
            "RNA_editing": lambda df: df["RNA_editing"] == True,
            "MAPPABILITY": lambda df: df["mappabilityScore"] != 0,
            "INDEL_PROPORTION": lambda df: df["indel_proportion_for_site"] < 0.05,
            "ALT_ALLELE_COUNT": lambda df: df["consensus_alt_allele_count"] >= thr_altcount,
            "POPULATION_AF": lambda df: df["falt"] < thr_popAF,
        }

        filtrations_dict = self.get_step_config()

        if filtrations_dict.keys():
            enabled_groups = list(filtrations_dict.keys())
        else:
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

        for col in merged_df.columns:
            merged_df[col] = merged_df[col].fillna("no")
            merged_df[col] = merged_df[col].replace("", "no")

        merged_df.to_csv(output_feature, sep="	", index=True)

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

            mask = all_filter_conditions[sub_key](df).fillna(False)
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

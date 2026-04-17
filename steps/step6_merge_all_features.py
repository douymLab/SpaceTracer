import os
import pandas as pd
import subprocess
import shlex

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


def _ensure_variant_multiindex(df: pd.DataFrame) -> pd.DataFrame:
    """
    统一把 df 设置成 MultiIndex:
    (#chrom, pos, ref, alt)
    """
    if df is None or df.empty:
        return pd.DataFrame()

    required_cols = ["#chrom", "pos", "ref", "alt"]

    if not isinstance(df.index, pd.MultiIndex):
        if all(col in df.columns for col in required_cols):
            df = df.copy()
            df.index = pd.MultiIndex.from_arrays(
                [df["#chrom"], df["pos"], df["ref"], df["alt"]],
                names=["#chrom", "pos", "ref", "alt"]
            )
    return df


def load_single_feature_file(path, sep="\t"):
    """
    读取单个 feature 文件：
    - 优先读取同名 parquet
    - 否则读取 txt/csv
    - 无论哪种方式都统一设置 MultiIndex
    """
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


def _is_nonempty_file(path: str) -> bool:
    return bool(path) and os.path.exists(path) and os.path.getsize(path) > 0


def _newer_than_all(target_file: str, deps: list) -> bool:
    """
    target_file 是否新于所有依赖文件
    """
    if not os.path.exists(target_file):
        return False
    target_mtime = os.path.getmtime(target_file)
    for dep in deps:
        if dep and os.path.exists(dep) and os.path.getmtime(dep) > target_mtime:
            return False
    return True

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

def _awk_merge_files(valid_files, output_file):
    """
    用 awk 合并多个带表头 txt：
    - 保留第一个文件表头
    - 跳过后续文件表头
    - 去空行
    """
    if not valid_files:
        with open(output_file, "w"):
            pass
        return

    files_quoted = " ".join(shlex.quote(f) for f in valid_files)
    out_quoted = shlex.quote(output_file)

    awk_script = r'NR==1{print;next} FNR==1{next} $0 !~ /^[[:space:]]*$/ {print}'
    cmd = f"awk {shlex.quote(awk_script)} {files_quoted} > {out_quoted}"
    subprocess.run(["bash", "-lc", cmd], check=True,shell=True)


def merge_feature_files_from_manifest(
    manifest_file: str,
    manifest_column: str,
    output_file: str,
    sep="\t",
    return_df=True,
):
    """
    从 manifest 合并多个 chunk feature 文件：
    - 优先命中 parquet 缓存
    - 否则 awk 合并 txt（去空行、去重复表头）
    - 生成 parquet 缓存
    """
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

    # 先命中 parquet 缓存
    # deps = [manifest_file] + valid_files
    # if _newer_than_all(parquet_file, deps):
    #     if return_df:
    #         df = load_parquet(parquet_file)
    #         df = _ensure_variant_multiindex(df)
    #         return df
        # return pd.DataFrame()

    if not valid_files:
        empty_df = pd.DataFrame()
        empty_df.to_csv(output_file, sep=sep, index=False)
        try:
            empty_df.to_parquet(parquet_file, index=False)
        except Exception:
            pass
        return empty_df

    # 用 awk 快速合并 txt
    _merge_files(valid_files, output_file)

    # 读合并结果并补索引
    merged_df = pd.read_csv(output_file, sep=sep)
    if merged_df is None or merged_df.empty:
        empty_df = pd.DataFrame()
        try:
            empty_df.to_parquet(parquet_file, index=False)
        except Exception:
            pass
        return empty_df

    merged_df = _ensure_variant_multiindex(merged_df)

    # 写 parquet 缓存
    try:
        merged_df.to_parquet(parquet_file, index=False, engine="pyarrow", compression="snappy")
    except Exception as e:
        print(f"Warning: failed to write parquet cache: {e}")

    if return_df:
        return merged_df
    return pd.DataFrame()


def _log_df_index_info(name, df):
    if df is None or df.empty:
        print(f"[{name}] empty")
        return
    try:
        n_rows = len(df)
        n_unique = df.index.nunique()
        n_dup = n_rows - n_unique
        print(f"[{name}] rows={n_rows}, unique_index={n_unique}, duplicated_index={n_dup}")
    except Exception as e:
        print(f"[{name}] failed to inspect index: {e}")


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
        return load_single_feature_file(path, sep="\t")

    def _run(self, context):
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

        # 确保索引统一
        RNA_feature_df = _ensure_variant_multiindex(RNA_feature_df)
        spatial_feature_df = _ensure_variant_multiindex(spatial_feature_df)
        read_feature_df = _ensure_variant_multiindex(read_feature_df)
        mappability_feature_df = _ensure_variant_multiindex(mappability_feature_df)

        # 打印 index 唯一性，帮助排查 join 扩行
        _log_df_index_info("RNA_feature", RNA_feature_df)
        _log_df_index_info("spatial_feature", spatial_feature_df)
        _log_df_index_info("read_feature", read_feature_df)
        _log_df_index_info("mappability_feature", mappability_feature_df)

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

        # 以 spatial 为主表
        if not short_spatial_feature_df.empty:
            merged_df = short_spatial_feature_df.copy()
        else:
            merged_df = pd.DataFrame()

        # 按 spatial index 左连接
        for name, df in [
            ("RNA_feature", RNA_feature_df),
            ("read_feature", read_feature_df),
            ("mappability_feature", mappability_feature_df),
        ]:

            if df is not None and not df.empty:
                # 去掉和主表重复的列，避免 join 冲突
                duplicate_cols = [c for c in df.columns if c in merged_df.columns]
                if duplicate_cols:
                    print(f"[{name}] dropping duplicated columns before join: {duplicate_cols}")
                    df = df.drop(columns=duplicate_cols)

                merged_df = merged_df.join(df, how="left")

        print(f"[merged_df] rows={len(merged_df)}")

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

        # 只处理 object/string 列，避免把数值列全转成 object
        object_cols = merged_df.select_dtypes(include=["object"]).columns
        for col in object_cols:
            merged_df[col] = merged_df[col].fillna("no")
            merged_df[col] = merged_df[col].replace("", "no")

        # Filtration 单独补一下
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

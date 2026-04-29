import pandas as pd
import numpy as np


def _safe_int(x, default=0):
    try:
        return int(x)
    except Exception:
        return default


def _parse_detail_counts(details: str):
    """
    例如:
    A,G:29;A,T:17;G,G:2;G,T:4

    返回:
    values = [29,17,2,4]
    """
    if pd.isna(details) or not str(details).strip():
        return []

    values = []
    for item in str(details).split(";"):
        try:
            values.append(int(item.split(":")[1]))
        except Exception:
            continue
    return values


def _compute_discordant_prop(details: str):
    """
    复刻你原始代码逻辑：
        values = [int(item.split(":")[1]) for item in details.split(";")]
        pick_hSNP_value=values[0::2] if values[-2] < values[-1] else values[1::2]
        discordant=min(values[-2:])
        discordant_prop=int((pick_hSNP_value[-1]))/int(sum(pick_hSNP_value))

    注意：
    你原代码里 discordant 其实没被使用，只保留 discordant_prop
    """
    values = _parse_detail_counts(details)
    if len(values) < 2:
        return 0.0

    try:
        pick_hSNP_value = values[0::2] if values[-2] < values[-1] else values[1::2]
        if not pick_hSNP_value or sum(pick_hSNP_value) == 0:
            return 0.0
        return float(pick_hSNP_value[-1]) / float(sum(pick_hSNP_value))
    except Exception:
        return 0.0


def summarize_phase_group(phase_group: pd.DataFrame, phase_type: str = "combine") -> dict:
    """
    phase_group: 同一个 mutation 的多条 phasing 记录
    返回一行 summary dict
    """
    if phase_group.empty:
        return {}

    g = phase_group.copy()

    # 类型清理
    g["pos"] = pd.to_numeric(g["pos"], errors="coerce")
    g["total_count"] = pd.to_numeric(g["total_count"], errors="coerce")
    g["mut_allele"] = pd.to_numeric(g["mut_allele"], errors="coerce")

    # 从 new_mut_name 里取配对位点位置
    # 例如 chr1_9271218_G_T -> 9271218
    g["phase_target_pos"] = g["new_mut_name"].astype(str).str.split("_").str[1]
    g["phase_target_pos"] = pd.to_numeric(g["phase_target_pos"], errors="coerce")

    # 距离：完全复刻你原代码，不取绝对值
    g["phase_distance"] = g["phase_target_pos"] - g["pos"]

    # mutant_counts 就是 mut_allele
    mutant_counts = g["mut_allele"].fillna(0).astype(float).tolist()
    total_mutants = sum(mutant_counts)
    max_mutants = max(mutant_counts) if mutant_counts else 0.0

    support_reads_prop_across_hSNPs = (
        max_mutants / total_mutants if total_mutants > 0 else np.nan
    )

    # nearest: 与原逻辑一致，直接用 min(distance)，不是 abs(distance)
    nearest_idx = g["phase_distance"].idxmin()

    # most: mut_allele 最大
    most_idx = g["mut_allele"].fillna(-1).idxmax()

    row_nearest = g.loc[nearest_idx]
    row_most = g.loc[most_idx]

    def extract_row_features(row, prefix):
        total_count = row.get("total_count")
        mut_allele = row.get("mut_allele")

        info_mutant_prop = np.nan
        try:
            if pd.notna(total_count) and float(total_count) != 0:
                info_mutant_prop = float(mut_allele) / float(total_count)
        except Exception:
            pass

        return {
            f"{prefix}_phase_haplotype": row.get("haplo"),
            f"{prefix}_info_mutant_prop": info_mutant_prop,
            f"{prefix}_discordant_prop": _compute_discordant_prop(row.get("detail_count")),
            f"{prefix}_mut_origin": row.get("mut_origin"),
            f"{prefix}_phase_distance": row.get("phase_distance"),
        }

    summary = {
        "#chrom": g.iloc[0]["#chrom"],
        "pos": g.iloc[0]["pos"],
        "germline": g.iloc[0]["germline"],
        "mutant": g.iloc[0]["mutant"],
        f"{phase_type}_support_reads_prop_across_hSNPs": support_reads_prop_across_hSNPs,
    }

    summary.update(extract_row_features(row_nearest, f"{phase_type}_nearest"))
    summary.update(extract_row_features(row_most, f"{phase_type}_most"))

    return summary


def build_phase_summary_df(phase_df: pd.DataFrame, phase_type: str = "phasing") -> pd.DataFrame:
    if phase_df is None or phase_df.empty:
        columns = [
            "#chrom", "pos", "germline", "mutant",
            f"{phase_type}_support_reads_prop_across_hSNPs",
            f"{phase_type}_nearest_phase_haplotype",
            f"{phase_type}_nearest_info_mutant_prop",
            f"{phase_type}_nearest_discordant_prop",
            f"{phase_type}_nearest_mut_origin",
            f"{phase_type}_nearest_phase_distance",
            f"{phase_type}_most_phase_haplotype",
            f"{phase_type}_most_info_mutant_prop",
            f"{phase_type}_most_discordant_prop",
            f"{phase_type}_most_mut_origin",
            f"{phase_type}_most_phase_distance",
        ]
        return pd.DataFrame(columns=columns)

    # 一个 variant 一组
    group_cols = ["#chrom", "pos", "germline", "mutant"]

    summary_rows = []
    for _, g in phase_df.groupby(group_cols, dropna=False):
        summary_rows.append(summarize_phase_group(g, phase_type=phase_type))

    summary_df = pd.DataFrame(summary_rows)

    return summary_df

from collections import Counter
from locale import D_FMT
import os
import shutil
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
import gc

from SpaceTracer.steps.base import BaseStep
from SpaceTracer.cores.RNA_features_00_from_fasta import add_features_from_fasta
from SpaceTracer.cores.RNA_features_01_ASE import get_ase_germline_sites, intersect_somatic_with_ase
from SpaceTracer.cores.RNA_features_02_hFDR import add_hFDR
from SpaceTracer.cores.RNA_features_03_imprinted_editing_PON import add_col_from_bed, add_col_from_mutant, add_col_from_mutant_from_df
from SpaceTracer.utils.get_MutationType import get_mutation_type, reorder_mutation_df
from SpaceTracer.utils.logger import get_logger
from SpaceTracer.utils.utils import load_manifest_tsv

model_name = __name__
logger = get_logger(model_name)


def collect_files_from_manifest(manifest_file: str, column_name: str):
    rows = load_manifest_tsv(manifest_file)
    files = []
    for row in rows:
        if not row:
            continue
        file_path = row.get(column_name, "")
        if file_path:
            files.append(file_path)
    return files

def merge_parquet_files_from_list(parquet_files, output_file, compression="snappy"):
    parquet_files = [f for f in parquet_files if f and os.path.exists(f) and os.path.getsize(f) > 0]

    out_dir = os.path.dirname(output_file)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    if not parquet_files:
        raise ValueError(f"No valid parquet files to merge for {output_file}")

    writer = None
    try:
        for file in parquet_files:
            table = pq.read_table(file)

            if writer is None:
                writer = pq.ParquetWriter(
                    output_file,
                    table.schema,
                    compression=compression
                )

            writer.write_table(table)

            del table
            gc.collect()
    finally:
        if writer is not None:
            writer.close()

    return output_file

def detect_clustered_noise(df, window_size=100, min_mutations=3, homo_thresh=0.8):

    # 保存原始索引
    original_index = df.index
    
    # 如果 df 是 MultiIndex，先重命名索引的 levels 避免冲突
    if isinstance(df.index, pd.MultiIndex):
        # 重命名索引 levels
        new_index_names = []
        for i, name in enumerate(df.index.names):
            if name in df.columns:
                new_index_names.append(f"{name}_idx")
            else:
                new_index_names.append(name)
        df.index.names = new_index_names
    
    # 重置索引
    df_work = df.reset_index()
    
    # 添加行号
    df_work['_row_id'] = range(len(df_work))
    
    # 构建位置字典
    chrom_pos_dicts = {}
    for chrom, group in df_work.groupby('chrom'):
        pos_dict = group.groupby('pos').apply(
            lambda x: list(zip(x['ref'], x['alt'], x['_row_id']))
        ).to_dict()
        chrom_pos_dicts[chrom] = pos_dict
    
    bases = ['A', 'C', 'G', 'T']
    mutation_types = [(r, a) for r in bases for a in bases if r != a]
    
    atog_cluster_ids = set()
    other_cluster_ids = set()
    
    for ref_base, alt_base in mutation_types:
        is_atog_type = (ref_base == 'A' and alt_base == 'G') or (ref_base == 'T' and alt_base == 'C')
        
        for chrom, pos_dict in chrom_pos_dicts.items():
            unique_positions = sorted(pos_dict.keys())
            n_unique = len(unique_positions)
            
            i = 0
            while i < n_unique:
                current_pos = unique_positions[i]
                muts_at_i = pos_dict[current_pos]
                
                if not any(r == ref_base and a == alt_base for r, a, idx in muts_at_i):
                    i += 1
                    continue
                    
                end_pos = current_pos + window_size
                j = i
                unique_site_count = 0
                target_site_count = 0
                temp_row_ids = []
                
                while j < n_unique and unique_positions[j] <= end_pos:
                    muts_at_j = pos_dict[unique_positions[j]]
                    unique_site_count += 1
                    
                    has_target = False
                    for r, a, idx in muts_at_j:
                        if r == ref_base and a == alt_base:
                            has_target = True
                            temp_row_ids.append(idx)
                    
                    if has_target:
                        target_site_count += 1
                    j += 1
                    
                if unique_site_count >= min_mutations and (target_site_count / unique_site_count) >= homo_thresh:
                    if is_atog_type:
                        atog_cluster_ids.update(temp_row_ids)
                    else:
                        other_cluster_ids.update(temp_row_ids)
                    i = j
                else:
                    i += 1
    
    # 返回与输入df顺序一致的Series
    atog_result = pd.Series(df_work['_row_id'].isin(atog_cluster_ids).values, index=original_index)
    other_result = pd.Series(df_work['_row_id'].isin(other_cluster_ids).values, index=original_index)
    
    return atog_result, other_result


def merge_table_files_from_list(file_list, output_file, sep="\t"):
    file_list = [
        f for f in file_list
        if f and os.path.exists(f) and os.path.getsize(f) > 0
    ]

    out_dir = os.path.dirname(output_file)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    if not file_list:
        raise ValueError(f"No valid table files to merge for {output_file}")

    expected_header = None
    wrote_any = False

    with open(output_file, "w") as fout:
        for f in file_list:
            with open(f, "r") as fin:
                header = fin.readline()
                if not header:
                    continue

                if expected_header is None:
                    expected_header = header
                    fout.write(header)
                    wrote_any = True
                else:
                    if header != expected_header:
                        raise ValueError(f"Header mismatch in file: {f}")

                shutil.copyfileobj(fin, fout, length=1024 * 1024)

    if not wrote_any:
        with open(output_file, "w") as fout:
            fout.write("")

    return output_file


def _mask_placeholder_priors(priors: pd.DataFrame) -> pd.DataFrame:
    """Treat genotyping placeholders 1,1,1,1 and fillna 0,0,0,0 as missing."""
    fake_prior_mask = (
        (priors == 1.0).all(axis=1) | (priors == 0.0).all(axis=1)
    )
    priors = priors.copy()
    priors.loc[fake_prior_mask] = np.nan
    return priors


def _format_ref_prior(row) -> object:
    """Join germline allele priors; return np.nan if any allele prior is missing."""
    vals = [row[f"{b}_prior"] for b in row["germline"].split(",")]
    if any(pd.isna(v) for v in vals):
        return np.nan
    return ",".join(str(v) for v in vals)


def _prepare_ind_df(df: pd.DataFrame) -> pd.DataFrame:
    df["pos"] = df["pos"].astype(int)

    counts = df["consensus_read_count"].str.split(",", expand=True).astype(int)
    counts.columns = ["A_count", "T_count", "C_count", "G_count"]

    priors = df["prior_ATCG"].str.split(",", expand=True).astype(float)
    priors.columns = ["A_prior", "T_prior", "C_prior", "G_prior"]
    priors = _mask_placeholder_priors(priors)

    df = pd.concat([df, counts, priors], axis=1, copy=False)

    counts_np = counts.to_numpy()
    priors_np = priors.to_numpy()

    df["count"] = counts.sum(axis=1)

    def get_max_other_count(row):
        ref = row['germline'] 
        alt = row['mutant']
        
        other_bases = [base for base in ['A', 'T', 'C', 'G'] 
                    if base not in ref+alt]
        
        other_counts = [row[f'{base}_count'] for base in other_bases]
        
        return max(other_counts) if other_counts else 0

    df['alt2_count'] = df.apply(get_max_other_count, axis=1)
    
    base_to_idx = {'A': 0, 'T': 1, 'C': 2, 'G': 3}
    alt_idx = df["mutant"].map(base_to_idx).to_numpy()
    row_idx = range(len(df))
    df["alt_count"] = counts_np[row_idx, alt_idx]
    df["alt_prior"] = priors_np[row_idx, alt_idx]

    single_mask = ~df['germline'].str.contains(',')
    
    df['ref_count'] = df.apply(
        lambda r: sum(r[f"{b}_count"] for b in r['germline'].split(',')), axis=1
    )
    # if germline has 1 allele
    ref_idx = df.loc[single_mask, "germline"].map(base_to_idx).to_numpy()
    df.loc[single_mask, "ref_count"] = counts_np[df.index[single_mask], ref_idx]

    df['ref_prior'] = df.apply(_format_ref_prior, axis=1)
    
    del counts
    del priors
    del counts_np
    del priors_np
    return df


def _prepare_germline_df_for_feature(df: pd.DataFrame) -> pd.DataFrame:
    df["pos"] = df["pos"].astype(int)

    counts = df["germ_count"].str.split(",", expand=True).astype(int)
    counts.columns = ["A_count", "T_count", "C_count", "G_count"]

    priors = df["prior_ATCG"].str.split(",", expand=True).astype(float)
    priors.columns = ["A_prior", "T_prior", "C_prior", "G_prior"]
    priors = _mask_placeholder_priors(priors)

    df = pd.concat([df, counts, priors], axis=1, copy=False)

    counts_np = counts.to_numpy()
    priors_np = priors.to_numpy()

    df["count"] = counts.sum(axis=1)

    def get_max_other_count(row):
        ref = row['germline'] 
        alt = row['mutant']
        
        other_bases = [base for base in ['A', 'T', 'C', 'G'] 
                    if base not in ref+alt]
        
        other_counts = [row[f'{base}_count'] for base in other_bases]
        
        return max(other_counts) if other_counts else 0

    df['alt2_count'] = df.apply(get_max_other_count, axis=1)
    
    base_to_idx = {'A': 0, 'T': 1, 'C': 2, 'G': 3}
    alt_idx = df["mutant"].map(base_to_idx).to_numpy()
    row_idx = range(len(df))
    df["alt_count"] = counts_np[row_idx, alt_idx]
    df["alt_prior"] = priors_np[row_idx, alt_idx]

    single_mask = ~df['germline'].str.contains(',')
    
    df['ref_count'] = df.apply(
        lambda r: sum(r[f"{b}_count"] for b in r['germline'].split(',')), axis=1
    )
    # if germline has 1 allele
    ref_idx = df.loc[single_mask, "germline"].map(base_to_idx).to_numpy()
    df.loc[single_mask, "ref_count"] = counts_np[df.index[single_mask], ref_idx]

    df['ref_prior'] = df.apply(_format_ref_prior, axis=1)
    
    del counts
    del priors
    del counts_np
    del priors_np
    return df


class RNAFeatureStep(BaseStep):
    def get_inputs(self, context):
        return {
            # GenotypingStep chunk manifest
            "genotype_results": context.get("genotype_results", ""),

            # UMICombineStep error count list/manifest
            "error_count_results": (
                context.get("error_count_results", "")
                or context.get("error_count_file", "")
            ),
        }

    def get_outputs(self, context):
        # sample = context.get("sample", "sample")
        out_dir = os.path.join(self.step_dir)
        os.makedirs(out_dir, exist_ok=True)

        return {
            "RNA_feature": os.path.join(out_dir, "RNA_feature.txt"),
            "merged_ind_geno_filter_file": os.path.join(out_dir, "merged_ind_geno_filter_file.txt"),
            "merged_ind_geno_file": os.path.join(out_dir, "merged_ind_geno_file.txt"),
            "merged_germline_file": os.path.join(out_dir, "merged_germline.txt"),
            "merged_error_count_file": os.path.join(out_dir, "merged_error_count.parquet"),
        }

    def get_step_config(self):
        return self.config.get("steps", {}).get("RNA_feature", {})

    def _resolve_ind_geno_filter_file(self, inputs, outputs):
        genotype_results = inputs.get("genotype_results", "")
        if not genotype_results or not os.path.exists(genotype_results):
            raise FileNotFoundError(f"genotype_results not found: {genotype_results}")

        table_files = collect_files_from_manifest(genotype_results, "ind_geno_filter_file")
        if not table_files:
            raise ValueError(f"No ind_geno_filter_file found in manifest: {genotype_results}")

        table_files_raw = collect_files_from_manifest(genotype_results, "ind_geno_file")
        ind_geno_file=merge_table_files_from_list(
                            table_files_raw,
                            outputs["merged_ind_geno_file"],
                            sep="\t"
                        )
        ind_geno_filter_file=merge_table_files_from_list(
            table_files,
            outputs["merged_ind_geno_filter_file"],
            sep="\t"
        )
        return ind_geno_file,ind_geno_filter_file

    def _resolve_germline_file(self, inputs, outputs):
        genotype_results = inputs.get("genotype_results", "")
        if not genotype_results or not os.path.exists(genotype_results):
            raise FileNotFoundError(f"genotype_results not found: {genotype_results}")

        parquet_files = collect_files_from_manifest(genotype_results, "germline_file")
        if not parquet_files:
            raise ValueError(f"No germline_file found in manifest: {genotype_results}")

        return merge_table_files_from_list(
            parquet_files,
            outputs["merged_germline_file"],
            sep="\t"
        )

    def _resolve_error_count_file(self, inputs, outputs):
        error_count_results = inputs.get("error_count_results", "")
        if not error_count_results or not os.path.exists(error_count_results):
            raise FileNotFoundError(f"error_count_results not found: {error_count_results}")

        parquet_files = collect_files_from_manifest(error_count_results, "parquet_file")
        if not parquet_files:
            raise ValueError(f"No error_count_file found in manifest: {error_count_results}")

        return merge_parquet_files_from_list(
            parquet_files,
            outputs["merged_error_count_file"]
        )

    def _write_empty_outputs(self, outputs):
        empty_cols = [
            "major_read_strand",
            "consensus_UMI_count",
            "consensus_ref_allele_count",
            "consensus_alt_allele_count",
            "consensus_alt2_allele_count",
            "fref",
            "falt",
            "p_mosaic",
            "AFind",
            "DNAMutationType",
            "RNAMutationType",
            "GCcontent",
            "cause_poly_alt",
            "homopolymer",
            "editing_AtoG",
            "consensus_alt2_proportion",
            "ASE",
            "hFDR",
            "imprinted",
            "editing_database",
            "PON",
            "RNA_editing",
        ]
        empty_df = pd.DataFrame(columns=empty_cols)
        empty_df.index = pd.MultiIndex.from_arrays(
            [[], [], [], []],
            names=["#chrom", "pos", "ref", "alt"]
        )

        RNA_feature = outputs["RNA_feature"]
        empty_df.to_csv(RNA_feature, sep="\t", index=True, na_rep="NA")

        parquet_file = str(RNA_feature).replace(".txt", ".parquet")
        empty_df.to_parquet(parquet_file, index=True)

    def _run(self, context: dict):
        genome_details=self.config['genome_details']
        species=genome_details['species']
        # optional
        gene_bed = self.config.get("gene_bed","")
        dbsnp_vcf_file = self.config.get("dbsnp_vcf_file","")
        imprinted_bed = self.config.get("imprinted_bed","")
        editing_bed = self.config.get("editing_bed","")
        PON_file = self.config.get("PON_file","")

        fasta_file = self.config.get("genome_fasta")
        reference_error_profile = self.config.get("reference_error_profile")

        inputs = self.get_inputs(context)
        outputs = self.get_outputs(context)

        genotype_results = inputs["genotype_results"]
        error_count_results = inputs["error_count_results"]

        if not genotype_results or not os.path.exists(genotype_results):
            raise FileNotFoundError(f"genotype_results not found: {genotype_results}")
        if not error_count_results or not os.path.exists(error_count_results):
            raise FileNotFoundError(f"error_count_results not found: {error_count_results}")

        ind_geno_file,ind_geno_filter_file = self._resolve_ind_geno_filter_file(inputs, outputs)
        germline_file = self._resolve_germline_file(inputs, outputs)
        error_count_file = self._resolve_error_count_file(inputs, outputs)

        ase_germline_file = ""
        RNA_feature = outputs["RNA_feature"]

        parameters = self.get_step_config()
        count_threshold = int(parameters["min_count_for_germline"])
        prior_threshold = float(parameters["min_prior_for_germline"])
        default_range = int(parameters["default_range_of_gene"])
        p_threshold = float(parameters["p_threshold"])
        previous_base = int(parameters["previous_base"])

        # logger.info(f"Using merged ind_geno_filter_file: {ind_geno_filter_file}")
        # logger.info(f"Using merged germline_file: {germline_file}")
        # logger.info(f"Using merged error_count_file: {error_count_file}")
        import time
        t = time.time()
        df = pd.read_csv(ind_geno_filter_file, sep="\t")
        # logger.info(f"read merged ind file: {time.time()-t:.2f}s, n={len(df)}")
        if df.empty:
            logger.warning("Merged ind_geno_filter_file is empty, writing empty RNA feature outputs.")
            self._write_empty_outputs(outputs)
            del df
            gc.collect()
            return {
                "RNA_feature": RNA_feature
            }

        t = time.time()
        df = _prepare_ind_df(df)
        # logger.info(f"_prepare_ind_df: {time.time()-t:.2f}s, n={len(df)}")

        # rename
        if "germline" in df.columns and "ref" not in df.columns:
            df = df.rename(columns={"germline": "ref"})
        if "mutant" in df.columns and "alt" not in df.columns:
            df = df.rename(columns={"mutant": "alt"})
        if "#chrom" in df.columns:
            df = df.rename(columns={"#chrom": "chrom"})

        required_cols = ["chrom", "pos", "ref", "alt"]
        missing = [c for c in required_cols if c not in df.columns]
        if missing:
            raise ValueError(f"Merged ind_geno_filter_file missing required columns: {missing}")

        df = df.sort_values(by=['chrom', 'pos'])

        df.index = pd.MultiIndex.from_arrays(
            [df["chrom"], df["pos"], df["ref"], df["alt"]],
            names=["chrom", "pos", "ref", "alt"]
        )

        t = time.time()
        df[["DNAMutationType", "RNAMutationType", "GCcontent", "cause_poly_alt", "homopolymer"]] = \
            add_features_from_fasta(df, fasta_file, previous_base)
        df["editing_AtoG"] = df["RNAMutationType"].astype(str).str.contains("A>G", na=False)
        # logger.info(f"add_features_from_fasta: {time.time()-t:.2f}s")

        result_df = df[
            [
                "strand", "count", "ref_count", "alt_count", "alt2_count",
                "ref_prior", "alt_prior", "p_mosaic", 
                "DNAMutationType", "RNAMutationType", "GCcontent", "cause_poly_alt", "homopolymer",
                "editing_AtoG"
            ]
        ].copy()

        total_count = result_df["count"].astype(float)
        result_df["AFind"] = np.where(
            total_count > 0,
            result_df["alt_count"].astype(float) / total_count,
            np.nan,
        )
        result_df = result_df.rename(columns={
            "strand": "major_read_strand",
            "count": "consensus_UMI_count",
            "ref_count": "consensus_ref_allele_count",
            "alt_count": "consensus_alt_allele_count",
            "alt2_count": "consensus_alt2_allele_count",
            "ref_prior": "fref",
            "alt_prior": "falt"
        })

        total_umi = result_df["consensus_UMI_count"].astype(float)
        result_df["consensus_alt2_proportion"] = np.where(
            total_umi > 0,
            result_df["consensus_alt2_allele_count"].astype(float) / total_umi,
            np.nan,
        )

        # print(germline_file,
        #     dbsnp_vcf_file,
        #     ase_germline_file,
        #     gene_bed,
        #     count_threshold,
        #     prior_threshold,
        #     p_threshold,
        #     default_range)

        t = time.time()
        # ASE (not work for non-human species)
        if species == "human":
            ase_germline_df = get_ase_germline_sites(
                germline_file,
                dbsnp_vcf_file,
                ase_germline_file,
                gene_bed,
                count_threshold,
                prior_threshold,
                p_threshold,
                default_range
                )
            # logger.info(f"get_ase_germline_sites: {time.time()-t:.2f}s, n={len(ase_germline_df)}")

            t = time.time()
            result_df["ASE"] = intersect_somatic_with_ase(df, ase_germline_df, p_threshold)
            # logger.info(f"intersect_somatic_with_ase: {time.time()-t:.2f}s")
            del ase_germline_df
            gc.collect()
        else:
            result_df["ASE"] = "unknown"

        # hFDR (this version not support non-human species)
        if species == "human": 
            error_count_df = pd.read_parquet(error_count_file)
            if error_count_df.empty:
                logger.warning("Merged error_count_file is empty, hFDR may not be meaningful.")

            if "#chrom" in error_count_df.columns:
                error_count_df = error_count_df.rename(columns={"#chrom": "chrom"})

            mutation_type_df = get_mutation_type(error_count_df, fasta_file, "RNA")
            error_profile = reorder_mutation_df(Counter(mutation_type_df["RNAMutationType"]))
            error_profile_file = os.path.join(self.step_dir, "error_profile.txt")
            error_profile.to_csv(error_profile_file, sep="\t", index=False)

            del mutation_type_df
            del error_count_df
            gc.collect()

            result_df["hFDR"] = add_hFDR(df, error_profile_file, reference_error_profile, self.step_dir)
        else:
            result_df["hFDR"] = "unknown"

        if PON_file and os.path.exists(PON_file):
            PON_df = pd.read_csv(PON_file, sep="\t", header=None)

            if len(PON_df.columns) >= 5:
                PON_df = PON_df.iloc[:, :5]
                PON_df.columns = ['chrom', 'pos', 'ref', 'alt', 'sample']
                PON_df = PON_df[PON_df["sample"] != self.sample]
                PON_df = PON_df[['chrom', 'pos', 'ref', 'alt']]

            elif len(PON_df.columns) == 4:
                PON_df.columns = ['chrom', 'pos', 'ref', 'alt']

            else:
                raise ValueError(
                    f"PON file format error: expected at least 4 columns, got {len(PON_df.columns)}"
                )
        else:
            PON_df = pd.DataFrame()
        result_df["PON"] = add_col_from_mutant_from_df(df, PON_df)

        # others (if file is '', will return unknown)
        result_df["imprinted"] = add_col_from_bed(df, imprinted_bed)

        result_df["editing_database"] = add_col_from_mutant(df, editing_bed)

        # Clustered_Technical_Noise (containing RNA_editing_clustered and other_clustered_noise)
        window_size=100
        min_mutations=3
        homo_thresh=0.8
        atog_series, other_series = detect_clustered_noise(df,window_size,min_mutations,homo_thresh)

        result_df['AtoG_clustered_noise'] = atog_series
        result_df['other_clustered_noise'] = other_series

            
        del df
        gc.collect()
        
        result_df["RNA_editing"] = (
            (result_df["editing_database"] == True) |
            (result_df["editing_AtoG"] == True) 
        )

        # save: txt uses NA; parquet keeps np.nan / null
        result_df.index.names = ["#chrom", "pos", "ref", "alt"]
        result_df.to_csv(RNA_feature, sep="\t", index=True, na_rep="NA")

        parquet_file = str(RNA_feature).replace(".txt", ".parquet")
        result_df.to_parquet(parquet_file, index=True)

        return {
            "RNA_feature": RNA_feature
        }

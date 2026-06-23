import math
from functools import lru_cache
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
from functools import partial

import numpy as np
import pandas as pd
from scipy.special import comb

from SpaceTracer.utils.read_files import load_spot_count_data
from SpaceTracer.utils.utils import str2dict
from itertools import repeat
import gc

import pyarrow as pa
import pyarrow.parquet as pq


def rename_df(df):
    if '#chrom' in df.columns:
        df=df.rename(columns={'#chrom':'chrom'})
    return df

def _make_spot_geno_indexed_schema(parquet_cols):
    """
    Build a stable pyarrow schema with pandas metadata.

    The parquet file will store `identifier` as pandas index.
    All output columns are stored as string to keep schema stable across batches.
    """
    empty_df = pd.DataFrame(columns=parquet_cols)

    for c in parquet_cols:
        empty_df[c] = empty_df[c].astype("string")

    empty_df["identifier"] = pd.Series(dtype="string")
    empty_df = empty_df.set_index("identifier")

    table = pa.Table.from_pandas(empty_df, preserve_index=True)
    return table.schema


def _flush_spot_geno_rows_to_parquet(
    rows,
    columns,
    writer,
    out_file,
    schema,
    compression="snappy",
):
    """
    Flush buffered spot genotype rows to parquet.

    rows:
        List of spot genotype rows.
        Each row should match `columns`.

    columns:
        Parquet data columns.
        Usually:
        [
            "chr", "pos", "strand", "germline", "mutant", "cluster",
            "spot_barcode", "consensus_read_count", "l_germline", "l_mosaic",
            "max_spot_geno", "G_spot_max", "depth", "vaf", "p_mosaic"
        ]

    writer:
        Existing pq.ParquetWriter or None.

    Returns:
        writer
    """
    if not rows:
        return writer

    df = pd.DataFrame(rows, columns=columns)

    # Force stable schema across batches.
    # This avoids schema mismatch when one batch infers int/float and another infers string.
    for c in columns:
        df[c] = df[c].astype("string")

    df["identifier"] = (
        df["chr"].astype("string") + "_" +
        df["pos"].astype("string") + "_" +
        df["germline"].astype("string") + "_" +
        df["mutant"].astype("string")
    )

    df = df.set_index("identifier")

    table = pa.Table.from_pandas(
        df,
        schema=schema,
        preserve_index=True
    )

    if writer is None:
        writer = pq.ParquetWriter(
            str(out_file),
            schema,
            compression=compression
        )

    writer.write_table(table)

    return writer


def _write_empty_spot_geno_parquet(
    out_file,
    columns,
    schema,
    compression="snappy",
):
    """
    Write an empty spot genotype parquet with identifier as index.
    This prevents downstream FileNotFoundError when no valid spot genotype rows exist.
    """
    empty_df = pd.DataFrame(columns=columns)

    for c in columns:
        empty_df[c] = empty_df[c].astype("string")

    empty_df["identifier"] = pd.Series(dtype="string")
    empty_df = empty_df.set_index("identifier")

    table = pa.Table.from_pandas(
        empty_df,
        schema=schema,
        preserve_index=True
    )

    pq.write_table(
        table,
        str(out_file),
        compression=compression
    )


class SpotGenoCalculator:
    def __init__(self, bins, epsQ, thr_dp, pop_vaf, cell_num, max_workers=1):
        self.bins = bins
        self.epsQ = epsQ
        self.thr_dp = thr_dp
        self.pop_vaf = pop_vaf
        self.cell_num = cell_num
        self.max_workers = max_workers

    # def run_from_df(self, spot_count_df, ind_geno_df, cluster_df, cluster_vaf_df, output_file):
    #     raw_spot_count_df = spot_count_df.copy()
    #     raw_spot_count_df = raw_spot_count_df.fillna('NA')

    #     cluster_vaf_df = cluster_vaf_df.copy()
    #     cluster_vaf_df = cluster_vaf_df.drop(["spot_number", "consensus_read_count"], axis=1, errors='ignore')

    #     if self.bins is None:
    #         ind_geno_df = ind_geno_df.drop(["cluster", "spot_number", "consensus_read_count", "prior_ATCG"], axis=1, errors='ignore')
    #         spot_count_df = raw_spot_count_df
    #     else:
    #         bins = self.bins
    #         ind_geno_df = ind_geno_df.drop(["cluster", "spot_number", "consensus_read_count", "prior_ATCG"], axis=1, errors='ignore')
    #         spot_count_df = bin_df(raw_spot_count_df, bins)

    #     spot_count_df=rename_df(spot_count_df)
    #     ind_geno_df=rename_df(ind_geno_df)
    #     cluster_vaf_df=rename_df(cluster_vaf_df)
    #     if cluster_df.empty:
    #         count_cluster_df=spot_count_df.copy()
    #         count_cluster_df['cluster'] = "bulk"
    #     else:
    #         cluster_df=rename_df(cluster_df)
    #         count_cluster_df = pd.merge(spot_count_df, cluster_df, on='barcode', how='left')
            
    #     count_geno_join = pd.merge(count_cluster_df, ind_geno_df, on=['chrom', 'pos', 'strand'], how='inner')
    #     count_geno_join = count_geno_join.rename(columns={'vaf': 'ind_vaf'})
    #     count_geno_vaf_join = pd.merge(
    #         count_geno_join,
    #         cluster_vaf_df,
    #         on=['chrom', 'pos', 'strand', 'germline', 'mutant', 'cluster'],
    #         how='left'
    #     )
    #     count_geno_vaf_join = count_geno_vaf_join.rename(columns={'vaf': 'cluster_vaf'})

    #     use_cell_num_file = (not isinstance(self.cell_num, int)) and Path(self.cell_num).exists()
    #     if use_cell_num_file:
    #         cell_num_file = pd.read_csv(self.cell_num, sep="\t", header=0)
    #         if "index" in cell_num_file.columns:
    #             cell_num_file.rename(columns={"index": "barcode"}, inplace=True)
    #         cell_num_file = cell_num_file.drop(["cluster", "nUMI", "nREAD"], axis=1, errors='ignore')
    #         count_geno_vaf_join = pd.merge(count_geno_vaf_join, cell_num_file, on='barcode', how='left')

    #     out_colname = "\t".join([
    #         "#chrom", "pos", "strand", "germline", "mutant", "cluster", "spot_barcode", "consensus_read_count",
    #         "l_germline", "l_mosaic", "max_spot_geno", "G_spot_max", "depth", "vaf", "p_mosaic"
    #     ])

    #     mosaic_mask = count_geno_vaf_join["genotype"].eq("mosaic")
    #     mosaic_df = count_geno_vaf_join[mosaic_mask]
    #     if use_cell_num_file:
    #         cell_nums = mosaic_df.iloc[:, -1].to_numpy()
    #         mosaic_df_no_cell = mosaic_df.iloc[:, :-1]
    #     else:
    #         cell_nums = np.full(len(mosaic_df), self.cell_num)
    #         mosaic_df_no_cell = mosaic_df

    #     tuple_cols = list(mosaic_df_no_cell.columns)
    #     col_idx = {c: i for i, c in enumerate(tuple_cols)}
    #     row_records = list(mosaic_df_no_cell.itertuples(index=False, name=None))
    #     payloads = list(zip(row_records, cell_nums))

    #     worker_fn = partial(
    #         _spot_genotype_from_payload,
    #         col_idx=col_idx,
    #         epsQ=self.epsQ,
    #         thr_dp=self.thr_dp,
    #         pop_vaf=self.pop_vaf
    #     )
    #     if self.max_workers > 1 and len(payloads) > 1:
    #         with ThreadPoolExecutor(max_workers=self.max_workers) as executor:
    #             all_rows = list(executor.map(worker_fn, payloads))
    #     else:
    #         all_rows = [worker_fn(payload) for payload in payloads]

    #     out_rows = [spot_geno_info for spot_geno_info in all_rows if spot_geno_info[10] != "NA"]
    #     with open(output_file, 'w') as f:
    #         f.write(out_colname + "\n")
    #         for spot_geno_info in out_rows:
    #             f.write("\t".join(spot_geno_info) + "\n")

    #     # Build parquet directly to avoid re-reading large TSV.
    #     out_cols = [
    #         "chr", "pos", "strand", "germline", "mutant", "cluster",
    #         "spot_barcode", "consensus_read_count", "l_germline", "l_mosaic",
    #         "max_spot_geno", "G_spot_max", "depth", "vaf", "p_mosaic"
    #     ]
    #     df = pd.DataFrame(out_rows, columns=out_cols)
    #     if not df.empty:
    #         df['identifier'] = (
    #             df['chr'].astype(str) + '_' + df['pos'].astype(str) + '_' +
    #             df['germline'].astype(str) + '_' + df['mutant'].astype(str)
    #         )
    #         df = df.set_index('identifier')
    #     else:
    #         df = pd.DataFrame(columns=out_cols)
    #         df.index.name = "identifier"
    #     parquet_file = str(output_file).replace('.out', '.parquet')
    #     df.to_parquet(parquet_file, index=True, compression='snappy')

    # def run(self, spot_count_file, ind_geno_file, cluster_df, cluster_vaf_file, output_file):
    #     ind_geno_df = pd.read_csv(ind_geno_file, sep="\t")
    #     raw_spot_count_df = load_spot_count_data(spot_count_file)
    #     cluster_vaf_df = pd.read_csv(cluster_vaf_file, sep="\t")
    #     self.run_from_df(raw_spot_count_df, ind_geno_df, cluster_df, cluster_vaf_df, output_file)
    def run_from_df(
        self,
        spot_count_df,
        ind_geno_df,
        cluster_df,
        cluster_vaf_df,
        output_file
    ):
        """
        Buffered streaming spot genotype calculation.

        Main behavior:
        - Write spot genotype result directly to parquet.
        - Save `identifier` as pandas index in parquet.
        - Do not build full out_rows list.
        - Do not build full output DataFrame.
        - Do not start inner ThreadPoolExecutor by default.
        - Keep memory close to one parquet buffer plus merge intermediates.

        output_file can be:
            xxx/spot_geno.parquet  -> write this parquet
            xxx/spot_geno.out      -> write xxx/spot_geno.parquet by default
        """

        output_path = Path(output_file)

        if output_path.suffix == ".parquet":
            parquet_file = output_path
            tsv_file = output_path.with_suffix(".out")
        else:
            tsv_file = output_path
            parquet_file = output_path.with_suffix(".parquet")

        # Default production behavior.
        write_parquet = getattr(self, "write_parquet", True)
        write_tsv = getattr(self, "write_tsv", False)

        parquet_flush_rows = getattr(self, "parquet_flush_rows", 100_000)
        parquet_compression = getattr(self, "parquet_compression", "snappy")

        # Output columns.
        # TSV keeps old "#chrom".
        tsv_columns = [
            "#chrom", "pos", "strand", "germline", "mutant", "cluster",
            "spot_barcode", "consensus_read_count", "l_germline", "l_mosaic",
            "max_spot_geno", "G_spot_max", "depth", "vaf", "p_mosaic"
        ]

        # Parquet uses "chr", matching downstream loader.
        parquet_cols = [
            "chr", "pos", "strand", "germline", "mutant", "cluster",
            "spot_barcode", "consensus_read_count", "l_germline", "l_mosaic",
            "max_spot_geno", "G_spot_max", "depth", "vaf", "p_mosaic"
        ]

        schema = _make_spot_geno_indexed_schema(parquet_cols)

        parquet_writer = None
        parquet_buffer = []
        tsv_fh = None

        n_seen = 0
        n_written = 0

        try:
            # Remove stale output before writing.
            if write_parquet and parquet_file.exists():
                parquet_file.unlink()

            if write_tsv and tsv_file.exists():
                tsv_file.unlink()

            print(
                f"[PROFILE] run_from_df start | "
                f"parquet={parquet_file} write_parquet={write_parquet} | "
                f"tsv={tsv_file} write_tsv={write_tsv} | "
                f"parquet_flush_rows={parquet_flush_rows}",
                flush=True
            )

            # ------------------------------------------------------------------
            # 1. Prepare input DataFrames
            # ------------------------------------------------------------------
            # Note:
            # fillna returns a new DataFrame. This is unavoidable if bin_df or downstream
            # cannot handle NA. But we avoid extra .copy() after this.
            spot_count_df = spot_count_df.fillna("NA")

            # drop returns a new DataFrame view/copy depending on pandas internals,
            # but we do not explicitly copy.
            cluster_vaf_df = cluster_vaf_df.drop(
                ["spot_number", "consensus_read_count"],
                axis=1,
                errors="ignore"
            )

            ind_geno_df = ind_geno_df.drop(
                ["cluster", "spot_number", "consensus_read_count", "prior_ATCG"],
                axis=1,
                errors="ignore"
            )

            if self.bins is not None:
                spot_count_df = bin_df(spot_count_df, self.bins)

            spot_count_df = rename_df(spot_count_df)
            ind_geno_df = rename_df(ind_geno_df)
            cluster_vaf_df = rename_df(cluster_vaf_df)

            # ------------------------------------------------------------------
            # 2. Add cluster
            # ------------------------------------------------------------------
            if cluster_df.empty:
                # Avoid copy. This mutates local spot_count_df by adding cluster.
                # Since spot_count_df was already reassigned above, this should be safe.
                count_cluster_df = spot_count_df
                count_cluster_df["cluster"] = "bulk"
            else:
                cluster_df = rename_df(cluster_df)

                # If cluster_df barcode is expected unique, validate='m:1' is useful.
                # If your data can have duplicated barcode intentionally, remove validate.
                count_cluster_df = pd.merge(
                    spot_count_df,
                    cluster_df,
                    on="barcode",
                    how="left",
                    validate="m:1"
                )

            # spot_count_df is no longer needed separately.
            # If count_cluster_df is spot_count_df itself, del spot_count_df only removes the name.
            try:
                del spot_count_df
            except Exception:
                pass

            # ------------------------------------------------------------------
            # 3. Join individual genotype
            # ------------------------------------------------------------------
            count_geno_join = pd.merge(
                count_cluster_df,
                ind_geno_df,
                on=["chrom", "pos", "strand"],
                how="inner"
            )
            count_geno_join = count_geno_join.rename(columns={"vaf": "ind_vaf"})

            try:
                del count_cluster_df
                del ind_geno_df
            except Exception:
                pass

            # ------------------------------------------------------------------
            # 4. Join cluster VAF
            # ------------------------------------------------------------------
            count_geno_vaf_join = pd.merge(
                count_geno_join,
                cluster_vaf_df,
                on=["chrom", "pos", "strand", "germline", "mutant", "cluster"],
                how="left",
                validate="m:1"
            )
            count_geno_vaf_join = count_geno_vaf_join.rename(columns={"vaf": "cluster_vaf"})

            try:
                del count_geno_join
                del cluster_vaf_df
            except Exception:
                pass

            # ------------------------------------------------------------------
            # 5. Optional cell_num file
            # ------------------------------------------------------------------
            use_cell_num_file = False

            if not isinstance(self.cell_num, int):
                try:
                    use_cell_num_file = Path(self.cell_num).exists()
                except Exception:
                    use_cell_num_file = False

            if use_cell_num_file:
                cell_num_file = pd.read_csv(self.cell_num, sep="\t", header=0)

                if "index" in cell_num_file.columns:
                    cell_num_file.rename(columns={"index": "barcode"}, inplace=True)

                cell_num_file = cell_num_file.drop(
                    ["cluster", "nUMI", "nREAD"],
                    axis=1,
                    errors="ignore"
                )

                count_geno_vaf_join = pd.merge(
                    count_geno_vaf_join,
                    cell_num_file,
                    on="barcode",
                    how="left",
                    validate="m:1"
                )

                try:
                    del cell_num_file
                except Exception:
                    pass

            # ------------------------------------------------------------------
            # 6. Filter mosaic rows
            # ------------------------------------------------------------------
            # Boolean filtering generally creates a new DataFrame. But it lets us drop
            # the full joined DataFrame afterwards to reduce memory during streaming.
            mosaic_df = count_geno_vaf_join[
                count_geno_vaf_join["genotype"].eq("mosaic")
            ]

            try:
                del count_geno_vaf_join
            except Exception:
                pass

            # ------------------------------------------------------------------
            # 7. Prepare row iterator and worker
            # ------------------------------------------------------------------
            if use_cell_num_file:
                # Keep original logic: last column after merging cell_num_file is cell number.
                cell_nums = mosaic_df.iloc[:, -1].to_numpy()
                mosaic_df_no_cell = mosaic_df.iloc[:, :-1]
                cell_iter = iter(cell_nums)
            else:
                mosaic_df_no_cell = mosaic_df
                cell_iter = repeat(self.cell_num)

            tuple_cols = list(mosaic_df_no_cell.columns)
            col_idx = {c: i for i, c in enumerate(tuple_cols)}

            worker_fn = partial(
                _spot_genotype_from_payload,
                col_idx=col_idx,
                epsQ=self.epsQ,
                thr_dp=self.thr_dp,
                pop_vaf=self.pop_vaf
            )

            row_iter = zip(
                mosaic_df_no_cell.itertuples(index=False, name=None),
                cell_iter
            )

            # ------------------------------------------------------------------
            # 8. Open optional TSV
            # ------------------------------------------------------------------
            if write_tsv:
                tsv_fh = open(tsv_file, "w")
                tsv_fh.write("\t".join(tsv_columns) + "\n")

            # ------------------------------------------------------------------
            # 9. Buffered streaming calculation and parquet write
            # ------------------------------------------------------------------
            for row_tuple, cell_num in row_iter:
                n_seen += 1

                spot_geno_info = worker_fn((row_tuple, cell_num))

                # Original logic:
                # spot_geno_info[10] is max_spot_geno.
                # Skip invalid rows.
                if spot_geno_info[10] == "NA":
                    continue

                n_written += 1

                if write_tsv:
                    tsv_fh.write("\t".join(map(str, spot_geno_info)) + "\n")

                if write_parquet:
                    # spot_geno_info has same order as parquet_cols except old first col
                    # is chrom value. We name it "chr" in parquet columns.
                    parquet_buffer.append(spot_geno_info)

                    if len(parquet_buffer) >= parquet_flush_rows:
                        parquet_writer = _flush_spot_geno_rows_to_parquet(
                            rows=parquet_buffer,
                            columns=parquet_cols,
                            writer=parquet_writer,
                            out_file=parquet_file,
                            schema=schema,
                            compression=parquet_compression
                        )
                        parquet_buffer.clear()

            # ------------------------------------------------------------------
            # 10. Final flush
            # ------------------------------------------------------------------
            if write_parquet and parquet_buffer:
                parquet_writer = _flush_spot_geno_rows_to_parquet(
                    rows=parquet_buffer,
                    columns=parquet_cols,
                    writer=parquet_writer,
                    out_file=parquet_file,
                    schema=schema,
                    compression=parquet_compression
                )
                parquet_buffer.clear()

        finally:
            try:
                if tsv_fh is not None:
                    tsv_fh.close()
            except Exception:
                pass

            try:
                if parquet_writer is not None:
                    parquet_writer.close()
            except Exception:
                pass

            try:
                del parquet_buffer
            except Exception:
                pass

            gc.collect()

        # If no valid row was written, create an empty parquet.
        # This must happen after closing writer.
        if write_parquet and n_written == 0:
            _write_empty_spot_geno_parquet(
                out_file=parquet_file,
                columns=parquet_cols,
                schema=schema,
                compression=parquet_compression
            )


def bin_df(df, bins):
    # df = df.copy()      
    df['spot_barcode']=df['barcode']
    df[['x', 'y']] = df['barcode'].str.split('_', expand=True)
    df['x'] = df['x'].astype(np.int64)
    df['y'] = df['y'].astype(np.int64)
    df['new_x'] = np.ceil((df['x']) / bins) * bins
    df['new_y'] = np.ceil((df['y']) / bins) * bins
    df['new_x'] = df['new_x'].astype(str)
    df['new_y'] = df['new_y'].astype(str)
    df['barcode'] = df['new_x'] + "_" + df['new_y']
    df = df.drop(columns=['x', 'y', 'new_x', 'new_y'])
        
    return df


def spot_posterior(germline, mutant, cluster_vaf, qA, qT, qC, qG, cell_num=20, epsQ=20, thr_dp=1000, pop_vaf=1e-5):
    """
    Calculate the posterior probability of the four spot genotypes:
    refhom, althom, het amd mosaic

    Note: this function also treat following two situations:
    - if the spot does not contain any mutant allele but has high mosaic posterior(>0.5), 
      then mosaic_posterior=0.5, germline_likelihood=mosaic_likelihood=0.5
    - if the spot has mutant allele but the mosaic likelihood is low (<0.5),
      then germline_likelihood=mosaic_likelihood=0.5

    Inputs:
        germline, mutant - the germline alleles and mutant allele on the individual level
        cluster_vaf - the allele frequency (AF) of the mutant (alternative) allele on the cluster level
        qA, qT, qC, qG - the consensus read qualities for each nucleotide
        cell_num - the number of cell of the spot (default=20)
        epsQ - the threshold for the consensus read quality (default=20)
        thr_dp - the threshold for the total depth (default=1000)
                 when depth > thr_dp, downsample the allele numbers to let the depth = thr_dp
        pop_vaf - the probability of mutant allele in population (default=1e-5)

    Outputs: p, l_norm, [max_spot_geno, G_spot_max, depth, vaf, p_mosaic]
        p - the normalized array of the posterior probabilities for all spot genotypes
        l_norm - the normalized array of the likelihoods for all spot genotypes
        max_spot_geno - the spot genotype with the highest posterior probability
        G_spot_max - represent the spot genotype by 0 and 1
                     i.e. refhom(0/0), althom(1/1), het(0/1), mosaic(0/1)
        depth - the UMI count depth for the spot
        vaf - the mutant allele frequency
        p_mosaic - the probability of having mosaic mutation of the spot (mosaic posterior probability)
    """
    # delete the low quality consensus reads and the low AF alternative alleles
    # get the consensus read and quality dictionaries
    count_filter, q_filter = spot_filter(qA, qT, qC, qG, epsQ=epsQ, thr_dp=thr_dp)
    count_filter_nozero = {k:v for k,v  in count_filter.items() if v>0}
    allele_list_spot = list(count_filter_nozero.keys())
    depth = sum(count_filter.values())

    # get the spot reference allele list
    ind_ref_list = germline.split(",")
    ref_list = [ref for ref in ind_ref_list if ref in allele_list_spot]
    alt = mutant
    # get the spot allele list except the errors
    if mutant != ".":
        ind_allele_list = ind_ref_list + [mutant]
    else:
        ind_allele_list = ind_ref_list
    allele_list = [allele for allele in allele_list_spot if allele in ind_allele_list]
    n_allele = len(allele_list)

    # return NAs if the candidate alleles are none
    if n_allele == 0:
        return "NA", ["NA"]*2, ["NA", "NA", depth, "NA", "NA"]
    # set fake ref allele if no actual ref allele
    elif len(ref_list) == 0:
        ref_list = [ind_ref_list[0]]

    # get the alternative allele info
    genotype_list = ["germline", "mosaic"]
    alt_count = count_filter[alt]
    qalt_dict = q_filter[alt]
    # if no alternative allele, give alt allele count as 0
    if alt == ".":
        alt_count = 0

    ## find the highest genotype
    # calcualte the likelihoods for each genotype
    # and combine the likelihoods for each reference candidate
    l_germline = 0.0
    l_mosaic = 0.0
    for ref in ref_list:
        # likelihood
        ref_count = count_filter[ref]
        qref_dict = q_filter[ref]
        l_update = spot_likelihood(ref_count, alt_count, qref_dict, qalt_dict, cluster_vaf, pop_vaf)
        l_germline += l_update[0]
        l_mosaic += l_update[1]

    # normalized likelihood
    l = [l_germline, l_mosaic]
    s_likelihood = l_germline + l_mosaic
    s_likelihood = 1e-15 if s_likelihood==0 else s_likelihood # avoid all 0
    l_norm = [j/s_likelihood for j in l]
    
    ## prior probability
    if cluster_vaf < 0.5:
        prior_value = [(1-2*cluster_vaf)**cell_num, 1-(1-2*cluster_vaf)**cell_num]
    else:
        prior_value = [pop_vaf, 1-pop_vaf]
    # prior_value = [[1-2*mu, 0, 2*mu, 0], [0, 1-2*mu, 2*mu, 0], [0, 0, 1, 0], [(1-2*cluster_vaf)**cell_num, 0, 0, 1-(1-2*cluster_vaf)**cell_num]]
    # prior_dict = dict(zip(genotype_list, prior_value))
    
    ## calculate the posterior values
    # multiply prior
    posterior = [prior_value[0] * l_germline, prior_value[1] * l_mosaic]
    # normalise
    s_posterior = sum(posterior)
    s_posterior = 1e-15 if s_posterior==0 else s_posterior # avoid all 0
    p = [j/s_posterior for j in posterior]
    # transfer all nan probability to 0 to avoid 0 ind genotype
    p = np.nan_to_num(p)
    ## find the highest genotype
    max_spot_geno = "mosaic" if p[1] > p[0] else "germline"

    ## probability of mosaic mutation
    p_mosaic = p[1]
    ## calculate the mutant allele frequency at the spot level
    vaf = alt_count / depth
    ## represent it by 0 and 1
    G_list = ["0/0", "0/1"]
    G_dict = dict(zip(genotype_list, G_list))
    G_spot_max = G_dict[max_spot_geno]

    # balance the conflict situation
    if alt_count == 0:
        max_spot_geno = "germline"
        G_spot_max = "0/0"
        if p_mosaic > 0.5:
            p_mosaic = 0.5
        if l_norm[1] > l_norm[0]:
            l_norm = [0.5, 0.5]
    elif alt_count > 0 and l_norm[1] < 0.5:
        l_norm = [0.5, 0.5]
    
    ## output
    return p, l_norm, [max_spot_geno, G_spot_max, depth, vaf, p_mosaic]


def spot_filter(qA, qT, qC, qG, epsQ=20, thr_dp=1000):
    """
    Delete the low quality (quality < epsQ) consensus reads
    Also downsample the allele numbers if the depth is very high

    Outputs:
        count_dict: a dictionary for the number of consensus reads after filtering
        q_filter: a dictionary for the quality of consensus reads after filtering
    """
    nucleotide_list = ["A", "T", "C", "G"]
    # transfrom the string to dict
    q_list = [_cached_str2dict(_normalize_q_input(q)) for q in [qA, qT, qC, qG]]
    # delete the low quality consensus reads
    q_filter_list = [{k:v for k,v in q.items() if k>=epsQ} for q in q_list]

    # downsample if the depth is too high
    depth = sum([sum(q.values()) for q in q_filter_list])
    if depth > thr_dp:
        # scale for each quality dict
        scaling_factor = thr_dp / depth
        downsampled_q_list = [{key: int(round(value * scaling_factor)) for key, value in q.items()} for q in q_filter_list]
        q_filter_list = downsampled_q_list
    q_filter = dict(zip(nucleotide_list, q_filter_list))
    
    # calculate the numbers of the nucleotides
    count_list = [sum(q.values()) for q in q_filter_list]
    count_dict = dict(zip(nucleotide_list, count_list))

    return count_dict, q_filter


def spot_likelihood(ref_count, alt_count, qref_dict, qalt_dict, cluster_vaf, pop_vaf=1e-5):
    """
    Calculate the four genotype likelihoods with one ref allele for the spots
    genotype: refhom, althom, het amd mosaic
    here only outputs two likelihoods: refhom and mosaic, 
    cause we could see these two values as the likelihoods for germline genotype or mosaic mutation

    Input:
        ref_count, alt_count - the number of the consensus reads for the reference and alternative alleles
        qref_dict, qalt_dict - the dictionaries for the qualities of the consensus reads on the ref and alt alleles
        cluster_vaf - the allele frequency (AF) of the mutant (alternative) allele on the cluster level
        pop_vaf - the probability of mutant allele in population (default=1e-5)

    Output:
        the likelihood array for the four / two genotypes
    """
    ## depth
    depth = ref_count + alt_count
    # log combination
    log_comb = math.log10(comb(depth, alt_count, exact=True))

    ## refhom likelihood
    q_refhom = math.log10(1)
    if ref_count != 0:
        q_refhom = q_refhom + sum(math.log10(1-0.1**(i/10)) * qref_dict[i] for i in qref_dict.keys())
    if alt_count != 0:
        q_refhom = q_refhom + sum(math.log10(0.1**(i/10)) * qalt_dict[i] for i in qalt_dict.keys())
    l_refhom = log_comb + q_refhom
    l_refhom = 10 ** l_refhom
    
    ## mosaic likelihood
    # avoid cluster_vaf=0/1
    cluster_vaf = pop_vaf if cluster_vaf==0 else cluster_vaf
    cluster_vaf = 1-pop_vaf if cluster_vaf==1 else cluster_vaf
    # calculate
    r = 0
    if ref_count != 0:
        r = r + sum(0.1**(float(i)/10) * qref_dict[i] for i in qref_dict.keys())
    if alt_count != 0:
        r = r + sum((1-0.1**(float(i)/10)) * qalt_dict[i] for i in qalt_dict.keys())
    l_mosaic = log_comb + r*math.log10(cluster_vaf) + (depth-r)*math.log10(1-cluster_vaf)
    l_mosaic = 10 ** l_mosaic   

    # combine the likelihoods
    l = [l_refhom, l_mosaic]
    
    return l


# def spot_genotype(join_info, cell_num=20, epsQ=20, thr_dp=1000, pop_vaf=1e-5):
#     qA = join_info[8]
#     qT = join_info[9]
#     qC = join_info[10]
#     qG = join_info[11]
#     cluster = join_info[12]
#     germline = join_info[13]
#     mutant = join_info[14]
#     cluster_vaf = 0 if join_info[19] == "NA" else float(join_info[19])

#     _, l_norm, spot_geno = spot_posterior(
#         germline, mutant, cluster_vaf, qA, qT, qC, qG,
#         cell_num=cell_num, epsQ=epsQ, thr_dp=thr_dp, pop_vaf=pop_vaf
#     )

#     output = join_info[0:3] + [germline, mutant, str(cluster)] + [join_info[5] , join_info[7]] + l_norm + [str(i) for i in spot_geno]
#     return output

def spot_genotype(join_info, cell_num=20, epsQ=20, thr_dp=1000, pop_vaf=1e-5):
    """
    must contain: chrom,pos,strand,barcode,consensus_read_count,qA,qT,qC,qG,cluster,germline,mutant,cluster_vaf
    """
    # 支持传入 pd.Series / dict
    if hasattr(join_info, "to_dict"):
        row = join_info.to_dict()
    else:
        row = dict(join_info)

    qA = row.get("qA", "NA")
    qT = row.get("qT", "NA")
    qC = row.get("qC", "NA")
    qG = row.get("qG", "NA")

    cluster = row.get("cluster", "NA")
    germline = row.get("germline", "NA")
    mutant = row.get("mutant", "NA")

    cluster_vaf_raw = row.get("cluster_vaf", "NA")
    cluster_vaf = 0 if cluster_vaf_raw in ("NA", None, "") else float(cluster_vaf_raw)

    _, l_norm, spot_geno = spot_posterior(
        germline, mutant, cluster_vaf, qA, qT, qC, qG,
        cell_num=cell_num, epsQ=epsQ, thr_dp=thr_dp, pop_vaf=pop_vaf
    )

    # 输出字段顺序保持你原来的定义：
    # join_info[0:3] + [germline, mutant, cluster] + [barcode, consensus_read_count] + l_norm + spot_geno
    output = [
        str(row.get("chrom", "NA")),
        str(row.get("pos", "NA")),
        str(row.get("strand", "NA")),
        str(germline),
        str(mutant),
        str(cluster),
        str(row.get("barcode", "NA")),
        str(row.get("consensus_read_count", "NA")),
    ] + l_norm + [str(i) for i in spot_geno]

    return output


def _normalize_q_input(q):
    if q is None:
        return "NA"
    if isinstance(q, str):
        return q
    return str(q)


@lru_cache(maxsize=200000)
def _cached_str2dict(q):
    return str2dict(q)


def _tuple_get(row, col_idx, key, default="NA"):
    idx = col_idx.get(key)
    if idx is None:
        return default
    value = row[idx]
    if pd.isna(value):
        return default
    return value


def spot_genotype_fast(row, col_idx, cell_num=20, epsQ=20, thr_dp=1000, pop_vaf=1e-5):
    qA = _tuple_get(row, col_idx, "qA")
    qT = _tuple_get(row, col_idx, "qT")
    qC = _tuple_get(row, col_idx, "qC")
    qG = _tuple_get(row, col_idx, "qG")

    cluster = _tuple_get(row, col_idx, "cluster")
    germline = _tuple_get(row, col_idx, "germline")
    mutant = _tuple_get(row, col_idx, "mutant")
    cluster_vaf_raw = _tuple_get(row, col_idx, "cluster_vaf")

    try:
        cluster_vaf = 0 if cluster_vaf_raw in ("NA", None, "") else float(cluster_vaf_raw)
    except (TypeError, ValueError):
        cluster_vaf = 0

    _, l_norm, spot_geno = spot_posterior(
        str(germline), str(mutant), cluster_vaf, qA, qT, qC, qG,
        cell_num=cell_num, epsQ=epsQ, thr_dp=thr_dp, pop_vaf=pop_vaf
    )

    return [
        str(_tuple_get(row, col_idx, "chrom")),
        str(_tuple_get(row, col_idx, "pos")),
        str(_tuple_get(row, col_idx, "strand")),
        str(germline),
        str(mutant),
        str(cluster),
        str(_tuple_get(row, col_idx, "barcode")),
        str(_tuple_get(row, col_idx, "consensus_read_count")),
        str(l_norm[0]),
        str(l_norm[1]),
    ] + [str(i) for i in spot_geno]


def _spot_genotype_from_payload(payload, col_idx, epsQ=20, thr_dp=1000, pop_vaf=1e-5):
    row, cell_num = payload
    return spot_genotype_fast(
        row=row,
        col_idx=col_idx,
        cell_num=cell_num,
        epsQ=epsQ,
        thr_dp=thr_dp,
        pop_vaf=pop_vaf
    )
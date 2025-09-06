# calculate_prs_vds.py

import hail as hl
import gcsfs
import pandas as pd
import os
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

from .utils import calculate_effect_allele_count_na_hom_ref


# ---------------------------------------------------------------
# Timer helper
# ---------------------------------------------------------------
class StepTimer:
    def __init__(self, step_name, icon="⏩"):
        self.step_name = step_name
        self.icon = icon
        self.start = None

    def __enter__(self):
        print(f"{self.icon} {self.step_name}...")
        self.start = time.time()

    def __exit__(self, exc_type, exc_value, traceback):
        elapsed = round(time.time() - self.start, 2)
        print(f"✅ {self.step_name} done in {elapsed}s\n")


# ---------------------------------------------------------------
# Core single-run function
# ---------------------------------------------------------------
def _calculate_prs_single(vds, prs_df, prs_identifier, pgs_weight_path, output_path,
                          bucket=None, save_found_variants=False):
    """
    Run a single PRS calculation (no chunking).

    Parameters
    ----------
    vds : hail.vds.VariantDataset
        The Hail VariantDataset object to calculate PRS on.
    prs_df : pandas.DataFrame
        DataFrame containing the PRS weights and variant info. 
        Must include columns: 'variant_id', 'weight', 'contig', 'position', 
        'effect_allele', 'noneffect_allele'.
    prs_identifier : str
        Unique name/ID for this PRS run (used in output file naming).
    pgs_weight_path : str
        Path to the PRS weights file (CSV).
    output_path : str
        Local or bucket directory where results will be stored.
    bucket : str, optional
        Cloud bucket name (e.g. "gs://..."). If None, outputs are written locally.
    save_found_variants : bool, default=False
        If True, also save a CSV listing variants from the weight file that 
        were found in the dataset.

    Returns
    -------
    pandas.DataFrame
        DataFrame with per-sample PRS results:
        - sample_id
        - sum_weights (total PRS)
        - N_variants (number of variants used)
    """

    fs = gcsfs.GCSFileSystem()

    # Construct paths
    if bucket:
        PGS_path = f"{bucket}/{pgs_weight_path}"
        interval_fp = f"{bucket}/{output_path}/interval/{prs_identifier}_interval.tsv"
        hail_fp = f"{bucket}/{output_path}/hail/{prs_identifier}"
        gc_csv_fp = f"{bucket}/{output_path}/score/{prs_identifier}_scores.csv"
        gc_found_csv_fp = f"{bucket}/{output_path}/score/{prs_identifier}_found_in_aou.csv"
    else:
        PGS_path = pgs_weight_path
        interval_fp = f"{output_path}/interval/{prs_identifier}_interval.tsv"
        hail_fp = f"{output_path}/hail/{prs_identifier}"
        gc_csv_fp = f"{output_path}/score/{prs_identifier}_scores.csv"
        gc_found_csv_fp = f"{output_path}/score/{prs_identifier}_found_in_aou.csv"

    # Step 1. Save intervals
    with StepTimer("Saving intervals", "📦"):
        prs_df["end"] = prs_df["position"]
        interval_df = prs_df[["contig", "position", "end"]]
        if bucket:
            with fs.open(interval_fp, "w") as f:
                interval_df.to_csv(f, header=False, index=False, sep="\t")
        else:
            interval_df.to_csv(interval_fp, header=False, index=False, sep="\t")

    # Step 2. Filter VDS
    with StepTimer("Filtering VDS with intervals", "🔍"):
        prs_sites = hl.import_locus_intervals(interval_fp, reference_genome="GRCh38", skip_invalid_intervals=True)
        vds_prs = hl.vds.filter_intervals(vds, prs_sites, keep=True)

    # Step 3. Re-import PRS weights with schema
    with StepTimer("Re-importing PRS table with schema", "📖"):
        required = ["variant_id", "weight", "contig", "position", "effect_allele", "noneffect_allele"]
        missing = [c for c in required if c not in prs_df.columns]
        if missing:
            raise ValueError(f"PRS table missing required columns: {', '.join(missing)}")

        optional = [c for c in prs_df.columns if c not in required]
        col_types = {c: "str" for c in required}
        col_types.update({"weight": "float64", "position": "int32"})
        col_types.update({c: "str" for c in optional})

        prs_table = hl.import_table(PGS_path, types=col_types, delimiter=",")
        prs_table = prs_table.annotate(locus=hl.locus(prs_table.contig, prs_table.position))
        prs_table = prs_table.key_by("locus")

    # Step 4. Annotate MT
    with StepTimer("Annotating MT with PRS info", "📝"):
        mt = vds_prs.variant_data.annotate_rows(prs_info=prs_table[vds_prs.variant_data.locus])
        mt = mt.unfilter_entries()

    # Step 5. Effect allele counts
    with StepTimer("Calculating effect allele counts", "🧮"):
        eff = calculate_effect_allele_count_na_hom_ref(mt)
        mt = mt.annotate_entries(
            effect_allele_count=eff,
            weighted_count=eff * mt.prs_info["weight"]
        )

    # Step 6. Aggregate per sample
    with StepTimer("Summing weighted counts per sample", "➕"):
        mt = mt.annotate_cols(
            sum_weights=hl.agg.sum(mt.weighted_count),
            N_variants=hl.agg.count_where(hl.is_defined(mt.weighted_count))
        )

    # Step 7. Write scores
    with StepTimer("Writing PRS scores to Hail Table", "💾"):
        mt.key_cols_by().cols().write(hail_fp, overwrite=True)

    # Step 8. Export to CSV
    with StepTimer("Exporting PRS scores to CSV", "📊"):
        saved_mt = hl.read_table(hail_fp)
        saved_mt.export(gc_csv_fp, header=True, delimiter=",")

    # Step 9. Optional found variants
    if save_found_variants:
        with StepTimer("Extracting found variants", "🔎"):
            found = mt.filter_rows(hl.is_defined(mt.prs_info)).rows()
            found_df = found.select(found.prs_info).to_pandas()
            if bucket:
                with fs.open(gc_found_csv_fp, "w") as f:
                    found_df.to_csv(f, header=True, index=False, sep=",")
            else:
                found_df.to_csv(gc_found_csv_fp, header=True, index=False, sep=",")

    return pd.read_csv(gc_csv_fp)


# ---------------------------------------------------------------
# Wrapper with chunking
# ---------------------------------------------------------------
def calculate_prs_vds(vds, prs_identifier, pgs_weight_path, output_path,
                      bucket=None, save_found_variants=False,
                      chunk_size=None, max_workers=1, max_retries=2):
    """
    Calculate PRS from a VDS, with optional chunking and resume support.

    Parameters
    ----------
    vds : hail.vds.VariantDataset
        The Hail VariantDataset object to calculate PRS on.
    prs_identifier : str
        Unique name/ID for this PRS run (used in output file naming).
    pgs_weight_path : str
        Path to the PRS weights file (CSV).
    output_path : str
        Local or bucket directory where results will be stored.
    bucket : str, optional
        Cloud bucket name (e.g. "gs://..."). If None, outputs are written locally.
    save_found_variants : bool, default=False
        If True, save per-chunk CSVs with found variants and merge them.
    chunk_size : int, optional
        Number of variants per chunk. If None, run all variants in one job.
        Useful to avoid memory/timeouts with large scores.
    max_workers : int, default=1
        Number of parallel workers for chunked execution.
    max_retries : int, default=2
        How many times to retry a failed chunk before giving up.

    Returns
    -------
    pandas.DataFrame or None
        Merged PRS scores for all samples if successful.
        Returns None if all chunks fail.

    Notes
    -----
    - Resume logic:
        * If scores CSV already exists for a chunk, it is skipped.
        * If scores exist but found variants are missing, the chunk is recomputed.
        * If neither exists, the chunk is fully recomputed.
    - A final merged CSV is produced for both scores and found variants (if requested).
    
    """

    print("")
    print("##########################################")
    print("##                                      ##")
    print("##                AoUPRS                ##")
    print("##    A PRS Calculator for All of Us    ##")
    print("##         Author: Ahmed Khattab        ##")
    print("##           Scripps Research           ##")
    print("##                                      ##")
    print("##########################################")
    print("")

    import warnings
    warnings.simplefilter(action="ignore", category=pd.errors.SettingWithCopyWarning)

    fs = gcsfs.GCSFileSystem()

    # Ensure local dirs if no bucket
    if bucket is None:
        for sub in ["interval", "weights", "hail", "score"]:
            os.makedirs(f"{output_path}/{sub}", exist_ok=True)

    # Load weights into pandas
    print("[*] Reading PRS weights...")
    if bucket:
        with fs.open(f"{bucket}/{pgs_weight_path}", "rb") as f:
            prs_df = pd.read_csv(f)
    else:
        prs_df = pd.read_csv(pgs_weight_path)

    # If no chunking → run single-shot
    if chunk_size is None:
        return _calculate_prs_single(vds, prs_df, prs_identifier,
                                     pgs_weight_path, output_path,
                                     bucket, save_found_variants)

    # Otherwise: chunked run
    prs_df["chunk_id"] = prs_df.index // chunk_size
    chunk_ids = prs_df["chunk_id"].unique().tolist()
    print(f"[*] Splitting into {len(chunk_ids)} chunks (chunk_size={chunk_size})")

    all_scores, failed = [], []

    def run_with_retry(cid):
        sub_df = prs_df.loc[prs_df["chunk_id"] == cid].copy()

        # expected outputs
        score_fp = f"{output_path}/score/{prs_identifier}_chunk{cid}_scores.csv"
        found_fp = f"{output_path}/score/{prs_identifier}_chunk{cid}_found_in_aou.csv"
        hail_fp  = f"{output_path}/hail/{prs_identifier}_chunk{cid}"

        if bucket:
            score_exists = fs.exists(f"{bucket}/{score_fp}")
            found_exists = fs.exists(f"{bucket}/{found_fp}") if save_found_variants else True
        else:
            score_exists = os.path.exists(score_fp)
            found_exists = os.path.exists(found_fp) if save_found_variants else True

        # Case 1: both exist → skip
        if score_exists and found_exists:
            print(f"⏭️ [chunk {cid}] Skipping (already complete).")
            return pd.read_csv(f"{bucket}/{score_fp}" if bucket else score_fp)

        # Case 2: scores exist but found variants missing
        if score_exists and save_found_variants and not found_exists:
            print(f"🧩 [chunk {cid}] Scores exist but no found variants. Recomputing chunk to regenerate them...")
            # force full recompute so prs_info rows are available
            return _calculate_prs_single(
                vds, sub_df, f"{prs_identifier}_chunk{cid}",
                pgs_weight_path, output_path,
                bucket, save_found_variants
            )

        # Case 3: neither exists → run full
        for attempt in range(1, max_retries + 1):
            # Clean up stale hail dir if present but scores missing
            if bucket:
                hail_path = f"{bucket}/{hail_fp}"
                if fs.exists(hail_path):
                    print(f"🧹 [chunk {cid}] Removing stale Hail table before rerun")
                    fs.rm(hail_path, recursive=True)
            else:
                if os.path.exists(hail_fp):
                    print(f"🧹 [chunk {cid}] Removing stale Hail table before rerun")
                    import shutil
                    shutil.rmtree(hail_fp)

            with StepTimer(f"Chunk {cid} attempt {attempt}", "⏱️"):
                try:
                    return _calculate_prs_single(
                        vds, sub_df, f"{prs_identifier}_chunk{cid}",
                        pgs_weight_path, output_path,
                        bucket, save_found_variants
                    )
                except Exception as e:
                    print(f"💥 [chunk {cid}] Failed attempt {attempt}: {e}")
                    if attempt == max_retries:
                        print(f"💀 [chunk {cid}] Giving up after {max_retries} attempts.")
                        failed.append(cid)
                        return None

    # Run chunks in parallel
    with ThreadPoolExecutor(max_workers=max_workers) as ex:
        futures = {ex.submit(run_with_retry, cid): cid for cid in chunk_ids}
        for fut in as_completed(futures):
            result = fut.result()
            if result is not None:
                all_scores.append(result)

    # Merge results
    if all_scores:
        # --- merge scores ---
        final = pd.concat(all_scores)
        final = final.groupby("s").agg({
            "sum_weights": "sum",
            "N_variants": "sum"
        }).reset_index().rename(columns={"s": "sample_id"})

        final_fp = f"{output_path}/score/{prs_identifier}_final_scores.csv"
        if bucket:
            with fs.open(f"{bucket}/{final_fp}", "w") as f:
                final.to_csv(f, index=False)
        else:
            final.to_csv(final_fp, index=False)

        print(f"\n🎉 Final merged PRS scores → {final_fp}")

        # --- merge found variants if requested ---
        if save_found_variants:
            found_dfs = []
            for cid in chunk_ids:
                found_fp = f"{output_path}/score/{prs_identifier}_chunk{cid}_found_in_aou.csv"
                if bucket:
                    if fs.exists(f"{bucket}/{found_fp}"):
                        with fs.open(f"{bucket}/{found_fp}", "rb") as f:
                            found_dfs.append(pd.read_csv(f))
                else:
                    if os.path.exists(found_fp):
                        found_dfs.append(pd.read_csv(found_fp))

            if found_dfs:
                merged_found = pd.concat(found_dfs, ignore_index=True)
                merged_found_fp = f"{output_path}/score/{prs_identifier}_final_found_in_aou.csv"
                if bucket:
                    with fs.open(f"{bucket}/{merged_found_fp}", "w") as f:
                        merged_found.to_csv(f, index=False)
                else:
                    merged_found.to_csv(merged_found_fp, index=False)

                print(f"🔎 Final merged found variants → {merged_found_fp}")

        if failed:
            print(f"⚠️ Warning: These chunks failed: {failed}")
        return final

    print("💥 No scores produced (all chunks failed).")
    return None

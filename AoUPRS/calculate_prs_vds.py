# calculate_prs_vds.py

import hail as hl
import gcsfs
import pandas as pd
import os
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

from .utils import calculate_effect_allele_count_na_hom_ref


# ---------------------------------------------------------------
# Timer helper with icons
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
# Core single-run function (your original steps)
# ---------------------------------------------------------------
def _calculate_prs_single(vds, prs_df, prs_identifier, pgs_weight_path, output_path,
                          bucket=None, save_found_variants=False):
    """
    Run PRS calculation exactly as the old function (single set of variants).
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
# Public wrapper with chunking
# ---------------------------------------------------------------
def calculate_prs_vds(vds, prs_identifier, pgs_weight_path, output_path,
                      bucket=None, save_found_variants=False,
                      chunk_size=None, max_workers=1, max_retries=2):
    """
    Calculate PRS from a VDS.
    Uses the original single-run pipeline, with optional chunking of variants.
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
        sub_df = prs_df.loc[prs_df["chunk_id"] == cid].copy()  # ✅ no warning

        # chunk output file
        chunk_fp = f"{output_path}/score/{prs_identifier}_chunk{cid}_scores.csv"
        if bucket:
            chunk_exists = fs.exists(f"{bucket}/{chunk_fp}")
        else:
            chunk_exists = os.path.exists(chunk_fp)

        # if already exists, reload
        if chunk_exists:
            print(f"⏭️ [chunk {cid}] Skipping (already exists).")
            return pd.read_csv(f"{bucket}/{chunk_fp}" if bucket else chunk_fp)

        for attempt in range(1, max_retries + 1):
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

    with StepTimer("Total PRS run", "⏳"):  # ✅ grand total
        with ThreadPoolExecutor(max_workers=max_workers) as ex:
            futures = {ex.submit(run_with_retry, cid): cid for cid in chunk_ids}
            for fut in as_completed(futures):
                result = fut.result()
                if result is not None:
                    all_scores.append(result)

    # Merge results
    if all_scores:
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
        if failed:
            print(f"⚠️ Warning: These chunks failed: {failed}")
        return final

    print("💥 No scores produced (all chunks failed).")
    return None

# calculate_prs_mt.py

import hail as hl
import gcsfs
import pandas as pd

from .utils import calculate_effect_allele_count

def calculate_prs_mt(mt, prs_identifier, pgs_weight_path, output_path, bucket=None, save_found_variants=False):
    """
    Calculates the Polygenic Risk Score (PRS) for the given MatrixTable and saves the results.
    
    :param mt: Hail MatrixTable.
    :param prs_identifier: Identifier for the PRS.
    :param pgs_weight_path: Path to the PRS weights file.
    :param output_path: Path to save the output files.
    :param bucket: (Optional) GCS bucket name if reading/writing from/to GCS.
    :param save_found_variants: (Optional) If True, save the found variants information.
    :return: None
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
    print("******************************************")
    print("*                                        *")
    print("*       Ahoy, PRS treasures ahead!       *")
    print("*                                        *")
    print("******************************************")
    print("")    
    print("")
    print("<<<<<<<<<>>>>>>>>")
    print(f"    {prs_identifier}   ")
    print("<<<<<<<<<>>>>>>>>")
    print("")

    # Construct paths
    if bucket:
        PGS_path = f'{bucket}/{pgs_weight_path}'
        hail_fp = f'{bucket}/{output_path}/hail/'
        gc_csv_fp = f'{bucket}/{output_path}/score/{prs_identifier}_scores.csv'
        gc_found_csv_fp = f'{bucket}/{output_path}/score/{prs_identifier}_found_in_aou.csv'
    else:
        PGS_path = pgs_weight_path
        hail_fp = f'{output_path}/hail/'
        gc_csv_fp = f'{output_path}/score/{prs_identifier}_scores.csv'
        gc_found_csv_fp = f'{output_path}/score/{prs_identifier}_found_in_aou.csv'
    
    print("Reading PRS weight table...")
    try:
        if bucket:
            with gcsfs.GCSFileSystem().open(PGS_path, 'rb') as gcs_file:
                prs_df = pd.read_csv(gcs_file)
        else:
            prs_df = pd.read_csv(PGS_path)
    except FileNotFoundError:
        raise FileNotFoundError(f"The file {PGS_path} was not found. Please check the file path and try again.")
    except Exception as e:
        raise Exception(f"An error occurred: {e}")
    
    # Ensure the required columns are present
    required_columns = ["variant_id", "weight", "contig", "position", "effect_allele", "noneffect_allele"]
    missing_columns = [col for col in required_columns if col not in prs_df.columns]
    if missing_columns:
        raise ValueError(f"The PRS table is missing required columns: {', '.join(missing_columns)}")
    
    # Include only required columns and any additional columns present in the PRS table
    optional_columns = [col for col in prs_df.columns if col not in required_columns]
    column_types = {col: "tstr" for col in required_columns}
    column_types.update({"weight": "tfloat", "position": "tint32"})
    column_types.update({col: "tstr" for col in optional_columns})
    
    print("Re-importing PRS table with determined column types...")
    # Re-import the PRS table with the determined column types
    prs_table = hl.import_table(PGS_path, types=column_types, delimiter=',')
    prs_table = prs_table.annotate(locus=hl.locus(prs_table.contig, prs_table.position))
    prs_table = prs_table.key_by('locus')

    print("Filtering variants...")
    # Filter variants
    mt_prs = mt.semi_join_rows(prs_table)
    
    print("Annotating the MatrixTable with the PRS information...")
    # Annotate the MatrixTable with the PRS information and calculate effect allele count
    mt_prs = mt_prs.annotate_rows(prs_info=prs_table[mt_prs.locus])

    print("Calculating effect allele count and multiplying by variant weight...")
    # Calculate effect allele count and multiply by variant weight in a single step
    effect_allele_count_expr = calculate_effect_allele_count(mt_prs)
    mt_prs = mt_prs.annotate_entries(
        effect_allele_count=effect_allele_count_expr,
        weighted_count=effect_allele_count_expr * mt_prs.prs_info['weight'])

    print("Summing the weighted counts per sample and counting the number of variants with weights per sample...")
    # Sum the weighted counts per sample and count the number of variants with weights per sample
    mt_prs = mt_prs.annotate_cols(
        sum_weights=hl.agg.sum(mt_prs.weighted_count),
        N_variants=hl.agg.count_where(hl.is_defined(mt_prs.weighted_count)))
    
    # Optional: Extract and save found variants
    if save_found_variants:
        print("Extracting and saving found variants...")
        found_variants_table = mt_prs.filter_rows(hl.is_defined(mt_prs.prs_info)).rows()
        found_prs_info_df = found_variants_table.select(found_variants_table.prs_info).to_pandas()
        
        # Get the number of found variants
        number_of_found_variants = found_prs_info_df.shape[0]
        print(f"Number of found variants: {number_of_found_variants}")
        
        if bucket:
            with gcsfs.GCSFileSystem().open(gc_found_csv_fp, 'w') as gcs_file:
                found_prs_info_df.to_csv(gcs_file, header=True, index=False, sep=',')
        else:
            found_prs_info_df.to_csv(gc_found_csv_fp, header=True, index=False, sep=',')
        print(f"Found variants saved as: {gc_found_csv_fp}")

    print(f"Writing the PRS scores to a Hail Table at: {hail_fp}")
    # Write the PRS scores to a Hail Table
    mt_prs.key_cols_by().cols().write(hail_fp, overwrite=True)
    
    print("Exporting the Hail Table to a CSV file...")
    # Export the Hail Table to a CSV file
    saved_mt = hl.read_table(hail_fp)
    if bucket:
        saved_mt.export(gc_csv_fp, header=True, delimiter=',')
    else:
        saved_mt.export(gc_csv_fp, header=True, delimiter=',')
    print(f"PRS scores saved as: {gc_csv_fp}")
    
    return
import pandas as pd
import gcsfs

def prepare_prs_table(input_path, output_path, bucket=None):
    """
    Prepares the PRS table by reading it from the given path, modifying columns, and saving the modified table.
    
    The provided DataFrame must have 'chr' and 'bp' columns. If the values in the 'chr' column start with 'chr', 
    they are kept as is. If the values in the 'chr' column are numbers, 'chr' is added to the beginning of these values to match the AoU formatting.
    
    :param input_path: Path to the input PRS table file.
    :param output_path: Path to save the modified PRS table file.
    :param bucket: (Optional) GCS bucket name if reading/writing from/to GCS.
    :return: Path to the saved modified PRS table file.
    """
    print("")
    print("********************************************************")
    print("*                                                      *")
    print("*   Winter is coming... and so is your Weight table!   *")
    print("*                                                      *")
    print("********************************************************")
    print("")    
    
    # Read the PRS table
    if bucket:
        gcs_path = f'{bucket}/{input_path}'
        with gcsfs.GCSFileSystem().open(gcs_path, 'rb') as gcs_file:
            prs_df = pd.read_csv(gcs_file)
    else:
        prs_df = pd.read_csv(input_path)
    
    # Ensure the required columns are present
    if 'chr' not in prs_df.columns or 'bp' not in prs_df.columns:
        raise ValueError("The input DataFrame must have 'chr' and 'bp' columns.")
    
    # Modify columns to fit Hail
    prs_df['contig'] = prs_df['chr'].astype(str).apply(lambda x: x if x.startswith('chr') else f'chr{x}')
    prs_df['position'] = prs_df['bp']
    prs_df['variant_id'] = prs_df.apply(lambda row: f'{row["contig"]}:{row["position"]}', axis=1)
    
    # Save the modified table
    if bucket:
        hail_df_fp = f'{bucket}/{output_path}'
        with gcsfs.GCSFileSystem().open(hail_df_fp, 'w') as gcs_file:
            prs_df.to_csv(gcs_file, index=False)
    else:
        hail_df_fp = output_path
        prs_df.to_csv(output_path, index=False)
    
    # Get the number of variants in the modified table
    number_of_variants = prs_df.shape[0]
    print(f"Number of variants in the modified table: {number_of_variants}")
    
    print(f"Modified PRS table saved as: {hail_df_fp}")

    print("")
    print("********************************************************")
    print("*                                                      *")
    print("*       Your quest is nearly complete, brave one!      *")
    print("*                         BUT                          *")
    print("*        The PRS adventure is about to get epic!       *")
    print("*                                                      *")
    print("********************************************************")
    print("")    
    return 
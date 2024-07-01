# calculate_prs.py

import os
import pandas as pd
import hail as hl
import gcsfs
from .calculate_prs_mt import calculate_prs_mt
from .calculate_prs_vds import calculate_prs_vds

# Define your main logic here, including input parameters, paths, and function calls

if __name__ == "__main__":
    # Example usage
    mt = ...  # Your MatrixTable
    vds = ...  # Your VariantDataset
    prs_identifier = "example_prs"
    pgs_weight_path = "path/to/prs_weights.csv"
    output_path = "output_directory"
    bucket = "your_bucket_name"  # Optional if using cloud storage
    save_found_variants = True  # Optional
    
    # Example calls
    calculate_prs_mt(mt, prs_identifier, pgs_weight_path, output_path, bucket=bucket, save_found_variants=save_found_variants)
    calculate_prs_vds(vds, prs_identifier, pgs_weight_path, output_path, bucket=bucket, save_found_variants=save_found_variants)

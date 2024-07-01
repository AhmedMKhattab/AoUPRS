# AoUPRS

## Overview

AoUPRS is a Python module designed for calculating Polygenic Risk Scores (PRS) specific to the All of Us study. This tool leverages Hail, a scalable framework for exploring and analyzing genomic data, to provide efficient PRS calculations.

AoUPRS provides 2 different approaches for PRS calculation [Check hte publication for more details]:   

__Approach 1:__ Using Hail Dense MatrixTable (MT)  
__Approach 2:__ Using Hail Sparse Variant Dataset (VDS)   



## Installation

To install AoUPRS from GitHub, run the following command:

```bash
pip install AoUPRS
```
## Dependencies
AoUPRS requires the following Python packages:

- hail
- gcsfs
- pandas

These dependencies will be installed automatically when you install AoUPRS.


## Usage
1. __Setup your AoU cloud analysis environment by selecting the "Hail Genomic Analysis" environment and allocating the required resources.__
   
   __How to set up a Dataproc cluster:__
   - __Hail MT:__  
     Requires more resources and from our experience, you need to allocate 300 workers. It's expensive but you end up saving time and money because the kernel crashes with lower resources.
     
        __Cost when running:__ $72.91 per hour  
        __Main node:__ 4CPUs, 15GB RAM, 150 GB Disk   
        __Workers (300):__ 4CPUs, 15GB RAM, 150GB Disk   

    - __Hail VDS:__
      The defualt resources will mostly suffice but if you have a big score and want to run it faster, use preemptible workers which are a lot more cheaper.

        __Cost when running:__ $0.73 per hour  
        __Main node:__ 4CPUs, 15GB RAM, 150 GB Disk   
        __Workers (2):__ 4CPUs, 15GB RAM, 150GB Disk 

** AoUPRS gives you the option to save the output files locally or to the cloud. We recommend to always save to the cloud as the the local files will be deleted with the deletion of Hail environment.

2. __Importing the Packages__
   
    To use AoUPRS, first import the package:
```py
import AoUPRS
import os
import pandas as pd
import numpy as np
from datetime import datetime
import gcsfs
import glob
import hail as hl
```

3. __Initiate Hail__
```py
hl.init(tmp_dir='hail_temp/', default_reference='GRCh38')
```

4. __Define Bucket__
```py
bucket = os.getenv("WORKSPACE_BUCKET")
```

5. __Read Hail MT / VDS__
```py
# Hail MT
mt_wgs_path = os.getenv("WGS_ACAF_THRESHOLD_MULTI_HAIL_PATH")
mt = hl.read_matrix_table(mt_wgs_path)

# Hail VDS
vds_srwgs_path = os.getenv("WGS_VDS_PATH")
vds = hl.vds.read_vds(vds_srwgs_path)
```
6. __Drop Flagged srWGS samples__  
    AoU provides a table listing samples that are flagged as part of the sample outlier QC for the srWGS SNP and Indel joint callset.

    Read more: https://support.researchallofus.org/hc/en-us/articles/4614687617556-How-the-All-of-Us-Genomic-data-are-organized#h_01GY7QZR2QYFDKGK89TCHSJSA7

```py
# Read flagged samples
flagged_samples_path = "gs://fc-aou-datasets-controlled/v7/wgs/short_read/snpindel/aux/relatedness/relatedness_flagged_samples.tsv"

# Save flagged samples locally
!gsutil -u $$GOOGLE_PROJECT cat $flagged_samples_path > flagged_samples.cvs

# Import flagged samples into a hail table
flagged_samples = hl.import_table(flagged_samples_path, key='sample_id')

# Drop flagged sample from main Hail 
## If Hail MT
mt = mt.anti_join_cols(flagged_samples)

## If Hail VDS:
vds_no_flag = hl.vds.filter_samples(vds, flagged_samples, keep=False)
```
7. __Define the sample__
```py
# For MT:
## Convert the subset_sample_ids to a Python set
subset_sample_ids_set = set(map(str, sample_ids['person_id'].tolist()))
## Filter samples
mt = mt.filter_cols(hl.literal(subset_sample_ids_set).contains(mt.s))

# For VDS:
## Import the sample as a Hail table
sample_needed_ht = hl.import_table('sample_ids.csv', delimiter=',', key='person_id')
## Filter samples
vds_subset = hl.vds.filter_samples(vds_no_flag, sample_needed_ht, keep=True)
```
8. __Prepare PRS Weight Table__
```py
# Prepare PRS weight table using function 'prepare_prs_table'

AoUPRS.prepare_prs_table('PGS######_table.csv',
'PGS######_weight_table.csv', bucket=bucket)

# Read PRS weight table

with gcsfs.GCSFileSystem().open(f'{bucket}/prs_calculator_tutorial/prs_calculator_hail_vds/vat_check/PGS002774_weight_table.csv', 'rb') as gcs_file:
    PGS002774_weights_tabel = pd.read_csv(gcs_file)
```


9. __Calculate PRS__
```py
# Define paths

prs_identifier = 'PGS######'
pgs_weight_path = 'PGS######_weight_table.csv'
output_path = 'PGS######'

# Calculate PRS

## MT:
AoUPRS.calculate_prs_mt(mt, prs_identifier, pgs_weight_path, output_path, bucket=None, save_found_variants=False)

## VDS:
AoUPRS.calculate_prs_vds(vds_subset, prs_identifier, pgs_weight_path, output_path, bucket=bucket, save_found_variants=True)
```


## Example Notebooks
For a detailed example, refer to the provided Jupyter notebook in the notebooks directory [https://github.com/AhmedMKhattab/AoUPRS/tree/main/notebooks]. This notebook demonstrates how to use the AoUPRS package to calculate PRS step-by-step.

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Author
Ahmed Khattab
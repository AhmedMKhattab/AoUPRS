# AoUPRS

<p align="center">
  <img src="https://raw.githubusercontent.com/AhmedMKhattab/AoUPRS/main/AoUPRS/aou_logo_bigger.png" alt="AoUPRS Logo" width="300"/>
</p>



## 🌟 Overview
AoUPRS is a Python module designed for calculating **Polygenic Risk Scores (PRS)** on the **All of Us** dataset.  
It leverages [Hail](https://hail.is/), a scalable framework for genomic data, to provide efficient and cost-effective PRS calculations.  

AoUPRS provides two approaches:  
- **MatrixTable (MT)** – dense representation  
- **Variant Dataset (VDS)** – sparse representation (recommended for v8)  

📄 [Publication in BMC Genomics (2025)](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-025-11693-9)



## ⚠️ Dataset Compatibility (v7 & v8)
- AoUPRS supports both **v7** and **v8** Controlled Tier WGS datasets.  
- **Key change in v8**: the `GT` field was removed. AoUPRS now reconstructs genotypes using:
  - `LGT` (local genotype index)  
  - `LA` (local alleles array)  

This ensures PRS calculations remain accurate and efficient with the sparser v8 VDS format.   

  

## 🚀 Performance Update (v8)

### Why is PRS calculation slower in v8?
- **v7** ≈ 1.7B variants  
- **v8** ≈ **4.9B variants** (3× larger)  
- Bigger callset = heavier interval queries, more I/O, more shuffle overhead.  

### Chunking
- Large scores (>1M SNPs) **must be chunked**.  
- Example runtimes [8 CPU / 52 GB driver, 2/10 workers (4 CPU / 15 GB RAM / 150 GB disk)]:  
  - 50k SNPs → 7–8 min  
  - 100k SNPs → 15 min  
  - 150k SNPs → 52 min ⚠️ nonlinear slowdown  

👉 Best practice: `chunk_size=50000`  

### Workers
- Scaling up workers does not help.  
- Sweet spot: **10 preemptible workers**.  
- More workers = more shuffle overhead + stragglers.  

### Cost Example
- Master: 8 CPU / 52 GB RAM  
- Workers: 2/10 × (4 CPU / 15 GB RAM / 150 GB disk)  
- Cost: ~$1.95/hour running, ~$0.11/hour paused  
- ~1M SNP PRS (20 chunks) ≈ 3h wall time, ~$6  



## 🛠️ Resume & Checkpointing
AoUPRS now supports **chunked execution with resume**:
- Each chunk is saved immediately (`PGS######_chunkN.csv`).  
- If the environment crashes, rerun will skip completed chunks and continue.  
- At the end, all chunks are merged into the final PRS results file.  

👉 This makes long PRS runs on v8 robust and restartable.  




## 🔧 Installation

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


## 📘 Usage Guide

### 1. Setup Environment
Select **Hail Genomic Analysis** environment.  

#### Dataproc cluster options:
- **Hail MT**  
  - Requires many resources (not recommended for v8)  
  - ~300 workers  
  - **Cost**: ~$72.91/hour

        
- **Hail VDS (recommended)**  
  - Works reliably with modest resources  
  - Best setup for large PRS runs:  
    - **Master node**: 8 CPUs / 52 GB RAM  
    - **Workers**: 10 × (4 CPUs / 15 GB RAM / 150 GB disk, preemptible)  
  - **Cost**: ~$1.95/hour running, ~$0.11/hour paused  
  - Runtime example: ~1M SNP PRS (20 × 50k chunks) ≈ 3 hours wall time, ~$6 total  

👉 AoUPRS gives you the option to save the output files locally or to the cloud. We recommend always saving to the cloud as the local files will be deleted with the deletion of the Hail environment.

### 2. Query [VAT](https://support.researchallofus.org/hc/en-us/articles/4615256690836-Variant-Annotation-Table) (optional)
Before calculating PRS, you may want to restrict your weight table to variants that are actually present in the **All of Us Variant Annotation Table (VAT)**.  
This ensures you are only scoring variants found in the callset.  

📓 Example notebook:  
[Query VAT and filter PRS weights](https://github.com/AhmedMKhattab/AoUPRS/tree/main/notebooks/AoUPRS_hailvds_PGS000746_check_vat.ipynb)

👉 Skipping this step will still work, but may include variants not found in AoU.

### 3. Import Packages
   
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
### 4. Initiate Hail
```py
hl.init(tmp_dir='hail_temp/', default_reference='GRCh38')
```
### 5. Define Bucket
```py
bucket = os.getenv("WORKSPACE_BUCKET")
```
### 6. Read Hail MT / VDS

#### Hail MT
```py

mt_wgs_path = os.getenv("WGS_ACAF_THRESHOLD_MULTI_HAIL_PATH")
mt = hl.read_matrix_table(mt_wgs_path)
```

#### Hail VDS
```py


vds_srwgs_path = os.getenv("WGS_VDS_PATH")
vds = hl.vds.read_vds(vds_srwgs_path)
```
### 7. Drop Flagged srWGS samples (optional)
  AoU provides a table listing samples that are flagged as part of the sample outlier QC for the srWGS SNP and Indel joint callset.

  Read more: [How the All of Us Genomic data are organized](https://support.researchallofus.org/hc/en-us/articles/4614687617556-How-the-All-of-Us-Genomic-data-are-organized#h_01GY7QZR2QYFDKGK89TCHSJSA7)

```py
# Read flagged samples

flagged_samples_path = "gs://fc-aou-datasets-controlled/v7/wgs/short_read/snpindel/aux/relatedness/relatedness_flagged_samples.tsv"

# Save flagged samples locally

!gsutil -u $$GOOGLE_PROJECT cat $flagged_samples_path > flagged_samples.csv

# Import flagged samples into a hail table

flagged_samples = hl.import_table(flagged_samples_path, key='sample_id')

# Drop flagged sample from main Hail 

## If Hail MT
mt = mt.anti_join_cols(flagged_samples)

## If Hail VDS:
vds_no_flag = hl.vds.filter_samples(vds, flagged_samples, keep=False)
```

### 8. Filter to Your Samples

#### Hail MT
```py

## Convert the subset_sample_ids to a Python set
subset_sample_ids_set = set(map(str, sample_ids['person_id'].tolist()))
## Filter samples
mt = mt.filter_cols(hl.literal(subset_sample_ids_set).contains(mt.s))
```


#### Hail VDS:
```py
## Import the sample as a Hail table
sample_needed_ht = hl.import_table('sample_ids.csv', delimiter=',', key='person_id')
## Filter samples
vds_subset = hl.vds.filter_samples(vds_no_flag, sample_needed_ht, keep=True)
```
### 9. Prepare PRS Weight Table
     
   The weight table must have these columns:
   
   ["chr", "bp", "effect_allele", "noneffect_allele", "weight"]
    
   The table below shows an example of a PRS weight table
   
  | chr | bp        | effect_allele | noneffect_allele | weight   |
  |-----|-----------|---------------|------------------|----------|
  | 2   | 202881162 | C             | T                | 1.57E-01 |
  | 14  | 996676    | C             | T                | 6.77E-02 |
  | 2   | 202881162 | C             | T                | 1.57E-01 |
  | 14  | 99667605  | C             | T                | 6.77E-02 |
  | 6   | 12903725  | G             | A                | 1.13E-01 |
  | 13  | 110308365 | G             | A                | 6.77E-02 |

    
```py
# Prepare PRS weight table using function 'prepare_prs_table'

AoUPRS.prepare_prs_table('PGS######_table.csv',
'PGS######_weight_table.csv', bucket=bucket)

# Read PRS weight table

with gcsfs.GCSFileSystem().open('PGS######_weight_table.csv', 'rb') as gcs_file:
    PGS######_weights_table = pd.read_csv(gcs_file)
```


### 10. Calculate PRS


```py
prs_identifier = "PGS######"
pgs_weight_path = "PGS######_weight_table.csv"
output_path = "PGS######"
```
#### Hail MT
```py

AoUPRS.calculate_prs_mt(mt, prs_identifier,
                        pgs_weight_path, output_path,
                        bucket=None, save_found_variants=False)
```
#### Hail VDS
```py

AoUPRS.calculate_prs_vds(vds_subset, prs_identifier,
                         pgs_weight_path, output_path,
                         bucket=bucket, save_found_variants=True,
                         chunk_size=50000)  # ✅ recommended
```

## 📓 Example Notebooks

We provide ready-to-use Jupyter notebooks that demonstrate step-by-step how to run AoUPRS:

- [AoUPRS with Hail VDS](https://github.com/AhmedMKhattab/AoUPRS/tree/main/notebooks/AoUPRS_hailvds_PGS000746_check_vat.ipynb)  
  Example of filtering PRS weights using the Variant Annotation Table (VAT) and calculating scores.

- [Other AoUPRS notebooks](https://github.com/AhmedMKhattab/AoUPRS/tree/main/notebooks)  
  Full collection of usage examples, including MT and VDS approaches.

👉 You can also try AoUPRS directly on the All of Us Researcher Workbench:  
🔗 [Launch AoUPRS on the All of Us Researcher Workbench](https://workbench.researchallofus.org/workspaces/aou-rw-c346f546/aouprsacosteffectiveprscalculatorfortheallofusprogram/analysis)

You can explore, duplicate, and run the included notebooks.

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](https://github.com/AhmedMKhattab/AoUPRS/blob/main/LICENSE) file for details.

## 📚 Citation

If you use AoUPRS in your research, please cite:

> Khattab A, Chen S-F, Wineinger N, Torkamani A. **AoUPRS: A Cost-Effective and Versatile PRS Calculator for the All of Us Program**. *BMC Genomics*. 2025;26:412. https://doi.org/10.1186/s12864-025-11693-9

## Author
Ahmed Khattab

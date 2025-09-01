# AoUPRS

<p align="center">
  <img src="https://raw.githubusercontent.com/AhmedMKhattab/AoUPRS/main/AoUPRS/aou_logo_bigger.png" alt="AoUPRS Logo" width="300"/>
</p>

## Overview

AoUPRS is a Python module designed for calculating Polygenic Risk Scores (PRS) specific to the All of Us study. This tool leverages Hail, a scalable framework for exploring and analyzing genomic data, to provide efficient PRS calculations.

AoUPRS provides 2 different approaches for PRS calculation [[Check the publication](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-025-11693-9) for more details]:   

__Approach 1:__ Using Hail Dense MatrixTable (MT)  

__Approach 2:__ Using Hail Sparse Variant Dataset (VDS)   

## ⚠️ Dataset Compatibility Update (v8)

> 🧬 **AoUPRS now supports both v7 and v8 of the All of Us Controlled Tier WGS dataset.**  
>
> 🔄 **Key change in v8:** The `GT` field is **no longer present** in the VDS.  
> ✅ AoUPRS now infers genotype calls using the new fields:
> - `LGT` – local genotype index  
> - `LA` – local alleles array  
>
> This ensures seamless and cost-effective PRS calculation using the new, sparser v8 VDS format.  
>

# 🚀 AoUPRS v8 Performance Update

With the release of **All of Us WGS v8**, users will notice differences in runtime and stability compared to v7 when calculating PRS scores. This section summarizes the key updates, scaling behavior, and practical recommendations.

---

## 🔎 Why v8 is Slower
- **Variant count explosion**  
  - v7 VAT ≈ 1.7B variants  
  - v8 VAT ≈ **4.9B variants** (almost 3× larger)  
- **Impact**: Interval queries and filtering are heavier, requiring more I/O and shuffle.  
- Even with the same cluster resources, v8 chunks take **longer wall-clock time** compared to v7.

---

## ⚖️ Chunking Behavior
- Large scores (>1M variants) must be processed in chunks to avoid memory failures.  
- **Chunk size vs runtime (v8, 8 CPU / 52 GB driver, 10–50 workers):**
  - **50k variants** → ~7–8 minutes per chunk  
  - **100k variants** → ~15 minutes per chunk  
  - **150k variants** → ~52 minutes per chunk (nonlinear slowdown)  

👉 **Recommendation**: Use **50k chunks** for stability and predictability.  
Checkpointing ensures partial progress is saved if the environment resets.  

---

## ⚡ Worker Scaling
Experiments show that **adding more workers does not speed up runs** on v8:

- Increasing from **10 → 30 → 50 workers** actually introduced overhead (task scheduling, shuffle, stragglers).  
- Best performance and stability were achieved with **10 workers**.  
- Runtime per chunk stayed essentially the same, with fewer retries and less instability.

👉 **Recommendation**:  
- Use **10 preemptible workers** (4 CPUs, 15 GB RAM each) as the sweet spot.  
- More workers add cost without reducing runtime.

---

## 💰 Cost Notes
- Example setup:  
  - **Master node**: 8 CPUs / 52 GB RAM  
  - **Workers**: 10 × (4 CPUs, 15 GB RAM, 150 GB disk), preemptible  
  - **Cost**: ~$1.95/hour running, $0.11/hour paused  
- At ~7–8 min per 50k chunk, a 1M SNP PRS (~20 chunks) runs in ~3 hours wall time, costing ~ $6.

---

## ✅ Practical Tips
- Expect **slower runtimes on v8** than v7 for the same cluster size.  
- Use **50k chunks + checkpointing** for stability.  
- Stick to **10 workers** — scaling up won’t help.  
- Budget more compute time and cost for large PRS calculations.  

---

## 🛠️ Troubleshooting and Resume Support
AoUPRS code has been updated to support **chunked runs with checkpointing**:

- Each chunk is saved immediately as `{identifier}_chunkN.csv` in your output bucket.  
- If the AoU environment crashes mid-run, already completed chunks are **skipped automatically** on restart.  
- At the end, all chunks can be recombined into a single final file.  

👉 This makes long PRS runs on v8 more robust: you won’t lose hours of progress if the job is interrupted. 

---


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


## Usage
1. __Setup your AoU cloud analysis environment by selecting the "Hail Genomic Analysis" environment and allocating the required resources.__
   
   __How to set up a Dataproc cluster:__
   - __Hail MT:__  Requires more resources. From our experience, you need to allocate 300 workers. It's expensive but you end up saving time and money because the kernel crashes with lower resources.
     
        __Cost when running:__ $72.91 per hour  
        __Main node:__ 4CPUs, 15GB RAM, 150 GB Disk   
        __Workers (300):__ 4CPUs, 15GB RAM, 150GB Disk   

    - __Hail VDS:__
      The default resources will mostly suffice, but if you have a big score and want to run it faster, use preemptible workers which are much cheaper.

        __Cost when running:__ $0.73 per hour  
        __Main node:__ 4CPUs, 15GB RAM, 150 GB Disk   
        __Workers (2):__ 4CPUs, 15GB RAM, 150GB Disk 

** AoUPRS gives you the option to save the output files locally or to the cloud. We recommend always saving to the cloud as the local files will be deleted with the deletion of the Hail environment.

2. If you wish to query the [Variant Annotation Table](https://support.researchallofus.org/hc/en-us/articles/4615256690836-Variant-Annotation-Table) before calculating a PRS from Hail VDS to  include only variants present in the callset, follow this [notebook](https://github.com/AhmedMKhattab/AoUPRS/tree/main/notebooks/AoUPRS_hailvds_PGS000746_check_vat.ipynb).

3. __Importing the Packages__
   
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
4. __Initiate Hail__
```py
hl.init(tmp_dir='hail_temp/', default_reference='GRCh38')
```
5. __Define Bucket__
```py
bucket = os.getenv("WORKSPACE_BUCKET")
```
6. __Read Hail MT / VDS__

```py
# Hail MT

mt_wgs_path = os.getenv("WGS_ACAF_THRESHOLD_MULTI_HAIL_PATH")
mt = hl.read_matrix_table(mt_wgs_path)

# Hail VDS

vds_srwgs_path = os.getenv("WGS_VDS_PATH")
vds = hl.vds.read_vds(vds_srwgs_path)
```
7. __Drop Flagged srWGS samples__  
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

8. __Define the sample__
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
9. __Prepare PRS Weight Table__
     
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


10. __Calculate PRS__
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
For detailed examples, refer to the provided Jupyter notebooks in the [notebooks directory](https://workbench.researchallofus.org/workspaces/aou-rw-c346f546/aouprsacosteffectiveprscalculatorfortheallofusprogram/data)
. These notebooks demonstrate how to use the AoUPRS package to calculate PRS step-by-step.

## 🚀 Try it on the All of Us Researcher Workbench

This tool is live and fully executable in a public workspace:

🔗 [Launch AoUPRS on the All of Us Researcher Workbench](https://workbench.researchallofus.org/workspaces/aou-rw-c346f546/aouprsacosteffectiveprscalculatorfortheallofusprogram/analysis)

You can explore, duplicate, and run the included notebooks — no setup required.

## License

This project is licensed under the MIT License - see the [LICENSE](https://github.com/AhmedMKhattab/AoUPRS/blob/main/LICENSE) file for details.

## 📚 Citation

If you use AoUPRS in your research, please cite:

> Khattab A, Chen S-F, Wineinger N, Torkamani A. **AoUPRS: A Cost-Effective and Versatile PRS Calculator for the All of Us Program**. *BMC Genomics*. 2025;26:412. https://doi.org/10.1186/s12864-025-11693-9

## Author
Ahmed Khattab

# ClusterV: finding HIV subtypes from ONT sequencing data

[![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause) 


Contact: Ruibang Luo, Junhao Su  
Email: rbluo@cs.hku.hk, jhsu@cs.hku.hk  

----

## Introduction
ClusterV is a standalone pipeline for accurately identifying HIV subtypes from ONT sequencing data.  ClusterV takes the alignment BAM file, reference FastA file and target region BED file (defining the target regions) of the HIV data as input, and outputs all found subtypes with their abundance, alignment, variants (SNPs and INDELs), and the drug resistance reports. ClusterV consists of three major modules, 1) filtering and initialing, where we run quality control at alignment files; 2) subtypes finding, we iteratively run variant calling - hierarchical clustering to find all variants and subtypes; 3) HIV mutation finding, where we generate consensus from each subtype and get all HIV drug resistance mutations reports. 

<img src="./docs/github_wf.png" width = "800" alt="ClusterV wf">

---

## Contents

* [Introduction](#introduction)
* [Installation](#installation)
  + [Option 1. Docker pre-built image](#option-1-docker-pre-built-image)
  + [Option 2. Build an anaconda virtual environment](#option-2-build-an-anaconda-virtual-environment)
  + [Option 3. Bioconda](#option-3-bioconda)
  + [Option 4. Docker Dockerfile](#option-4-docker-dockerfile)
* [Usage](#usage)
* [Understand output file](/docs/output.md)
* [Using localized HIVdb](#using-localized-hivdb)

## Installation


### Option 2. Build an anaconda virtual environment

```
conda env create -f clusterV.yml
conda activate clusterV 

pypy3 -m ensurepip
pypy3 -m pip install --no-cache-dir intervaltree==3.0.2
```


## Usage

```bash

INPUT_BAM="{YOUR_INPUT_BAM}"        # Input BAM file
INPUT_REF="{YOUR_INPUT_FASTA}"      # Input reference Fasta file
INPUT_BED="{YOUR_INPUT_BED}"        # Input BED file for defining the target region
SAMPLE_ID="{YOUR_SAMPLE_ID}"        # Sample ID
OUTPUT_DIR="{YOUR_OUTPUT_DIR}"      # Output path

python cv.py ClusterV \
 --bam_fn ${INPUT_BAM} \
 --ref_fn ${INPUT_REF} \
 --bed_fn ${INPUT_BED} \
 --sample_id ${SAMPLE_ID} \
 --out_dir ${OUTPUT_DIR}
```


## Using localized HIVdb

ClusterV uses online API to query drug resistance mutations in default. If you wish to use the localized HIVdb for HIVdb, please launch the HIVdb's Sierra web service locally, please use the following setting to run the ClusterV program.

```
# Step 1, Start Sierra in local, as instructed in [this page](https://github.com/hivdb/sierra#start-sierra-with-docker).
docker pull hivdb/sierra:latest
docker run -it --publish=8111:8080 hivdb/sierra dev

# Step 2, Run ClusterV with addtional "--hivdb_url" setting.
python cv.py ClusterV \
--hivdb_url http://localhost:8111/sierra/rest/graphql \
--bam_fn ${INPUT_BAM} \
--ref_fn ${INPUT_REF} \
--bed_fn ${INPUT_BED} \
--sample_id ${SAMPLE_ID} \
--out_dir ${OUTPUT_DIR}
```

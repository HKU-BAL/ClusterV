# ClusterV: finding HIV subtypes from ONT sequencing data

[![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause) 


Contact: Ruibang Luo, Junhao Su  
Email: rbluo@cs.hku.hk, jhsu@cs.hku.hk  

----

## Introduction
ClusterV is a standalone pipeline for accurately identifying HIV subtypes from ONT sequencing data.  ClusterV takes the alignment BAM file, and target region BED files (defining the target regions) of the HIV data as input, and outputs all found subtypes with their abundance, alignment, variants (SNPs and INDELs), and the drug resistance reports. ClusterV consists of three major modules, 1) filtering and initialing, where we run quality control at alignment files; 2) subtypes finding, we iteratively run variant calling - hierarchical clustering to find all variants and subtypes; 3) HIV mutation finding, where we generate consensus from each subtype and get all HIV drug resistance mutations reports. 


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


## Usage


Option 2. 

conda env create -f clusterV.yml
pypy3 -m ensurepip
pypy3 -m pip install --no-cache-dir intervaltree==3.0.2

### General Usage

```bash

INPUT_BAM="{YOUR_INPUT_BAM}"        # Input BAM file
INPUT_REF="{YOUR_INPUT_FASTA}"      # Input reference Fasta file
INPUT_BED="{YOUR_INPUT_BED}"        # Input BED file for defining the target region
SAMPLE_ID="{YOUR_SAMPLE_ID}"        # Sample ID
OUTPUT_DIR="{YOUR_OUTPUT_DIR}"      # Output path
MIN_R_SUPPORT=50                    # minimal read support for forming a subtype
TOP_K_SUBTYPE=25                    # number of top K subtypes to predict

# STEP 1, run filtering and initializing
python cv.py run_initial_call --bam_fn ${INPUT_BAM} --ref_fn ${INPUT_REF} --bed_fn ${INPUT_BED} --sample_id ${SAMPLE_ID} --out_dir ${OUTPUT_DIR}


# STEP 2, run clustering
WORK_BAM=${OUTPUT_DIR}/${SAMPLE_ID}_f.bam 
WORK_VCF=${OUTPUT_DIR}/${SAMPLE_ID}.v/out.vcf
WORK_DIR=${OUTPUT_DIR}/${SAMPLE_ID}/clustering
python cv.py ClusterV --ref_fn ${INPUT_REF} --bed_fn ${INPUT_BED} --bam ${WORK_BAM} --vcf ${WORK_VCF} --out_dir ${WORK_DIR} --sample_id ${SAMPLE_ID} --top_k ${TOP_K_SUBTYPE} --n_min_supports ${MIN_R_SUPPORT}


# STEP 3, generate consesus and find drug resistance mutation
WORK_OUT=${OUTPUT_DIR}/consensus
WORK_TSV=${OUTPUT_DIR}/clustering/all_clusters_info.tsv
python ${CV_DIR}/cv.py get_consensus --ref_fn ${INPUT_REF} --bed_fn ${INPUT_BED} --tar_tsv ${WORK_TSV} --out_dir ${WORK_OUT}
```


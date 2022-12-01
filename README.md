# ClusterV: finding HIV subtypes from ONT sequencing data

[![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause) 


Contact: Ruibang Luo, Junhao Su  
Email: rbluo@cs.hku.hk, jhsu@cs.hku.hk  

----

## Introduction


---

## Contents

* [Introduction](#introduction)
* [Installation](#installation)
  + [Option 1. Docker pre-built image](#option-1-docker-pre-built-image)
  + [Option 2. Build an anaconda virtual environment](#option-2-build-an-anaconda-virtual-environment)
  + [Option 3. Bioconda](#option-3-bioconda)
  + [Option 4. Docker Dockerfile](#option-4-docker-dockerfile)
* [Usage](#usage)


## Usage


Option 2. 

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

# set working files
WORK_BAM=${OUTPUT_DIR}/${SAMPLE_ID}_f.bam 
WORK_VCF=${OUTPUT_DIR}/${SAMPLE_ID}.v/out.vcf
WORK_DIR=${OUTPUT_DIR}/${SAMPLE_ID}/clustering

# STEP 2, run clustering
python cv.py ClusterV --ref_fn ${INPUT_REF} --bed_fn ${INPUT_BED} --bam ${WORK_BAM} --vcf ${WORK_VCF} --out_dir ${WORK_DIR} --sample_id ${SAMPLE_ID} --top_k ${TOP_K_SUBTYPE} --n_min_supports ${MIN_R_SUPPORT}




```

## Understand output file

### The directory structure for the clusterV results is as below:

    .
    ├── ${OUTPUT_DIR}
    │   ├── consensus           		# output of clusterV
    │   │   ├── all_report.tsv          # drug resistance report for found subtypes
    │   │   ├── all_info.tsv            # raw statistics for all found subtypes, with abundance information
    │   │   ├── id_1_report.tsv			# HIVDB report for subtype 1
    │   │   ├── id_1_report.json		# [raw HIVDB](https://hivdb.stanford.edu/page/release-notes/#json.output) report output for subtype 1
    │   │   ├── id_1.fasta				# consensus for subtype 1
    │   │   ├── id_1                    # data directory for subtype 1
    │   │   │   ├── id_1.bam            # all read alignment in bam file for subtype 1
    │   │   │   ├── id_1_cs.bam         # single consensus alignment for subtype 1
    │   │   │   ├── id_1_cs.vcf         # all variants in consensus for subtype 1
    │   │   │   ├── id_1_low_af.vcf     # low AF variants in subtypes but not in consensus for subtype 1
    │   │   │   └── ...
    │   │   ├── id_2 					# data directory for subtype 2
    │   │   │   └── ...
    │   │   └── ...
    │   ├── clustering          		# working files for clusterV algorithm
    │   │   ├── all_clusters_info.tsv	# clustering results, subtypes, abundance, and files location information 
    │   │   ├── 1						# first iteration
    │   │   │   └── 1.tagged.bam  		# bam file with subtype labeling after 1st iteration                 
    │   │   │   └── 1.tagged.tag1.bam  	# subtype 1 after 1st iteration                 
    │   │   │   └── 1.tagged.tag2.bam  	# subtype 2 after 1st iteration                 
    │   │   │   └── ...                 
    │   │   ├── 1_1						# subtype after first iteration 
    │   │   │   └── ...                 
    │   │   └── ...             
    │   ├── ${SAMPLE_ID}_f.bam  # BAM file after initalize filtering
    │   └── ...                 
    └── ...

### drug resistance report for found subtypes

all drug resistance report is available at `all_report.tsv`

| id              | Sample Id                                                                                                                           |
|-----------------|-------------------------------------------------------------------------------------------------------------------------------------|
| gene            | [Gene name](https://hivdb.stanford.edu/page/release-notes/#drm.classification) from HIVDB                                           |
| drugClass       | Drug class name from HIVDB                                                                                                          |
| drugName        | Drug name from HIVDB                                                                                                                |
| drugScore       | [Drug score](https://hivdb.stanford.edu/page/release-notes/#drm.penalty.scores.and.resistance.interpretation) from HIVDB            |
| resistanceLevel | [Drug resistance level](https://hivdb.stanford.edu/page/release-notes/#drm.penalty.scores.and.resistance.interpretation) from HIVDB |
| subtype_ori     | Subtype raw ID when running clustering                                                                                              |
| subtype         | Subtype ID                                                                                                                          |
| abundance       | Subtype's abundance                                                                                                                 |
| is_in_consensus | Whether the drug resistance mutation found in consensus (1 if found in subtype, 0 if found in low AF set)                           |
| VAF             | variants AF in subtype                                                                                                              |
| mutationPos     | mutation genome position in HIV reference                                                                                           |
| mutation        | mutation name                                                                                                                       |
| mutationScore   | [mutation sscore](https://hivdb.stanford.edu/page/release-notes/#drm.penalty.scores.and.resistance.interpretation) in HIVDB         |
| mutationType    | [mutation type](https://hivdb.stanford.edu/page/release-notes/#drm.classification) in HIVDB                                         |
| comments        | [mutation comments](https://hivdb.stanford.edu/page/release-notes/#comments) from HIVDB                                             |



### statistics for all found subtypes, with abundance information

| column                                               | description                                                                                                                                |
|------------------------------------------------------|--------------------------------------------------------------------------------------------------------------------------------------------|
| _idx                                                 | subtype ID                                                                                                                                 |
| _vcf                                                 | raw VCF for subtypes during clustering                                                                                                     |
| _bam                                                 | raw BAM for the subtype during clustering                                                                                                  |
| _out_dir                                             | working directory for the subtype during clustering                                                                                        |
| _sample_id                                           | sample ID                                                                                                                                  |
| _sample_idx                                          | raw index for the subtype                                                                                                                  |
| is_check_clustering                                  | checking status for subtypes, whether try splitting in this subtype                                                                        |
| percentage                                           | subtype's abundance                                                                                                                        |
| _v_info:coverage;snp_c;indel_c;median_af;af1;af2;af3 | subtype's coverage, # of SNP, # of INDEL, the median of AF, # of raw AF in [0, 0.3], # of raw AF in [0.3, 0.7], # of raw AF in [0.7, 1.0]  |

all subtype information is available at `all_info.tsv`


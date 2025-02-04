# ClusterV: finding HIV quasispecies and drug resistance from ONT sequencing data

[![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause) 
[![docker](https://img.shields.io/badge/build-docker-brightgreen)](https://www.docker.com/)
[![ONT](https://img.shields.io/badge/Support-ONT-005c75)](https://nanoporetech.com/)




Contact: Ruibang Luo, Junhao Su  
Email: rbluo@cs.hku.hk, jhsu@cs.hku.hk  

----

## Introduction
ClusterV is a standalone pipeline for accurately identifying HIV quasispecies from ONT sequencing data.  ClusterV takes the alignment BAM file, reference FastA file and target region BED file (defining the target regions) of the HIV data as input, and outputs all found quasispecies with their abundance, alignment, variants (SNPs and INDELs), and the drug resistance reports. ClusterV consists of three major modules, 1) filtering and initialing, where we run quality control at alignment files; 2) quasispecies finding, we iteratively run variant calling - hierarchical clustering to find all variants and quasispecies; 3) HIV mutation finding, where we generate consensus from each quasispecies and get all HIV drug resistance / mutations reports. 

<img src="./docs/github_wf.png" width = "800" alt="ClusterV wf">

---

## Contents

* [Introduction](#introduction)
* [Latest Update](#latest-update)
* [Installation](#installation)
  + [Option 1. Docker pre-built image](#option-1-docker-pre-built-image)
  + [Option 2. Build an anaconda virtual environment](#option-2-build-an-anaconda-virtual-environment)
  + [Option 3. Docker Dockerfile](#option-3-docker-dockerfiles)
* [Usage](#usage)
  + [Parameters](/docs/options.md)
  + [Understand output files](/docs/output.md)
* [Quick demo](#quick-demo)
* [Using localized HIVdb](#using-localized-hivdb)
* [Changing pre-trained models of Clair3](#changing-pre-trained-models-of-clair3)
* [Acknowledgement](#acknowledgement)
* [Publication](#publication)


## Latest Update

*v1.3 (Dec 25, 2024)*: Support using [Clair3](https://github.com/HKU-BAL/Clair3) for variant calling.

*v1.2 (Nov 24, 2023)*: Update Flye [option](https://github.com/fenderglass/Flye/blob/flye/docs/USAGE.md) to support ONT Guppy5+ SUP or Q20 data.

*v1.1 (Dec 1, 2022)*: Initial release.

## Installation

### Option 1. Docker pre-built image
A pre-built docker image is available [here](https://hub.docker.com/repository/docker/hkubal/clusterv). With it you can run ClusterV using a single command.

Caution: Absolute path is needed for both `INPUT_DIR` and `OUTPUT_DIR`.

```
INPUT_DIR=`pwd`                                     # input path, e.g. /input
INPUT_REF=${INPUT_DIR}/HIV_1.fasta                  # change your reference file name here
INPUT_BED=${INPUT_DIR}/HIV_1_amplicon_region.bed    # change your bed file name here
INPUT_BAM=${INPUT_DIR}/mix_ds.bam                   # change your bam file name here
SAMPLE_ID="mix_ds"                                  # change your sample ID here
OUTPUT_DIR=`pwd`/out                                # output path, e.g. /output
THREADS=48                                          # threads, e.g. 48
# CLAIR3_MODEL_PATH="{YOUR_MODEL_PATH}"             # Optional: change to the absolute path of your clair3_model, defaulting to '../Clair3/models/ont'

docker run -it \
  -v ${INPUT_DIR}:${INPUT_DIR} \
  -v ${OUTPUT_DIR}:${OUTPUT_DIR} \
  # -v ${CLAIR3_MODEL_PATH}:${CLAIR3_MODEL_PATH} \  # Optional, delete if unchanged
  hkubal/clusterv:v1.3 \
  python /opt/bin/cv.py ClusterV \
  --ref_fn ${INPUT_REF} \
  --bed_fn ${INPUT_BED} \
  --bam_fn ${INPUT_BAM} \
  --sample_id ${SAMPLE_ID} \
  # --clair3_model_path ${CLAIR3_MODEL_PATH} \      # Optional, delete if unchanged
  --out_dir ${OUTPUT_DIR} \
  --threads ${THREADS}
```

### Option 2. Build an anaconda virtual environment

```
# get ClusterV code
git clone --recursive https://github.com/HKU-BAL/ClusterV.git
# if you see nothing in the Clair3 folder, please run:
# git submodule update --init --recursive
cd ClusterV

# create env
conda env create -f clusterV.yml
conda activate clusterV 

pypy3 -m ensurepip
pypy3 -m pip install mpmath==1.2.1
pypy3 -m pip install --no-cache-dir intervaltree==3.0.2

# prepare for Clair3
cd ./Clair3
make PREFIX=${CONDA_PREFIX}
# download pre-trained models to folder ./Clair3/models
mkdir models
wget http://www.bio8.cs.hku.hk/clair3/clair3_models/clair3_models.tar.gz 
tar -zxvf clair3_models.tar.gz -C ./models
cd ..

# run ClusterV
python cv.py ClusterV --help
```

### Option 3. Docker Dockerfiles
Building a docker image.

```
# clone ClusterV
git clone --recursive https://github.com/HKU-BAL/ClusterV.git
cd ClusterV

# build a docker image named hkubal/ClusterV:latest
# might require docker authentication to build docker image 
docker build -f ./Dockerfile -t hkubal/clusterv:v1.3 .

# run ClusterV docker image like 
docker run -it hkubal/clusterv:v1.3 python /opt/bin/cv.py --help
```



## Usage

For more parameters setting, please check the [parameters page](/docs/options.md).

For understanding the output files, please check the [output page](/docs/output.md).

Make sure you are in the clusterV environment by using `conda activate clusterV`, when using conda environment.

```bash
CV_DIR="{ClusterV repo path}"
INPUT_BAM="{YOUR_INPUT_BAM}"           # Input BAM file
INPUT_REF="{YOUR_INPUT_FASTA}"         # Input reference Fasta file
INPUT_BED="{YOUR_INPUT_BED}"           # Input BED file for defining the target region
SAMPLE_ID="{YOUR_SAMPLE_ID}"           # Sample ID
OUTPUT_DIR="{YOUR_OUTPUT_DIR}"         # Output path
CLAIR3_MODEL_PATH="{YOUR_MODEL_PATH}"  # Optional: change to the absolute path of your clair3_model, defaulting to '../Clair3/models/ont'

python ${CV_DIR}/cv.py ClusterV \
 --bam_fn ${INPUT_BAM} \
 --ref_fn ${INPUT_REF} \
 --bed_fn ${INPUT_BED} \
 --sample_id ${SAMPLE_ID} \
 --clair3_model_path ${CLAIR3_MODEL_PATH} \  # Optional, delete if unchanged
 --out_dir ${OUTPUT_DIR}
```


## Quick demo

Make sure you are in the clusterV environment by using `conda activate clusterV`, when using conda environment.

```
# step 1, download testing files
cd ClusterV
mkdir -p testing
(cd testing && wget http://www.bio8.cs.hku.hk/ClusterV/clusterv_testing.tar.gz && tar -xf clusterv_testing.tar.gz && rm clusterv_testing.tar.gz )

# step 2, run testing
TESTING_DIR=`pwd`/testing
INPUT_REF=${TESTING_DIR}/HIV_1.fasta
INPUT_BED=${TESTING_DIR}/HIV_1_amplicon_region.bed
INPUT_BAM=${TESTING_DIR}/mix_ds.bam
SAMPLE_ID="mix_ds"
OUTPUT_DIR=${TESTING_DIR}/out
CV_DIR=`pwd`

python ${CV_DIR}/cv.py ClusterV \
--bam_fn ${INPUT_BAM} \
--ref_fn ${INPUT_REF} \
--bed_fn ${INPUT_BED} \
--sample_id ${SAMPLE_ID} \
--out_dir ${OUTPUT_DIR}
```

## Using localized HIVdb

ClusterV uses online API to query drug resistance mutations in default. 

If you wish to use the [localized HIVdb](https://github.com/hivdb/sierra#start-sierra-with-docker) for HIV drug resistance mutation detection, please launch the HIVdb's [Sierra](https://github.com/hivdb/sierra#start-sierra-with-docker) web service locally, and use the following setting to run the ClusterV program.

```
# Step 1, Start Sierra in local, as instructed in https://github.com/hivdb/sierra#start-sierra-with-docker.
docker pull hivdb/sierra:latest
docker run -it --publish=8111:8080 hivdb/sierra dev

# Step 2, Run ClusterV with addtional "--hivdb_url" setting.
python ${CV_DIR}/cv.py ClusterV \
--hivdb_url http://localhost:8111/sierra/rest/graphql \
--bam_fn ${INPUT_BAM} \
--ref_fn ${INPUT_REF} \
--bed_fn ${INPUT_BED} \
--sample_id ${SAMPLE_ID} \
--out_dir ${OUTPUT_DIR}
```

## Changing pre-trained models of Clair3

You can modify the models used by Clair3 for variant calling by specifying `--clair3_model_path`. 

Caution: Absolute path is needed for `--clair3_model_path`.

In the clusterV environment, set `CLAIR3_MODEL_PATH` to the folder where your Clair3 models are located. When running `ClusterV`, add `--clair3_model_path ${CLAIR3_MODEL_PATH}`. By default, `--clair3_model_path` is set to `'../Clair3/models/ont'`, which refers to the `'ont'` models extracted in `Clair3/models/`. For usage, please refer to [Usage](#usage).

In a Docker image, HKU-provided models are located in `/opt/bin/Clair3/models/`. If you use the models in the Docker image, add `--clair3_model_path /opt/bin/Clair3/models/MODEL_NAME` to specify the model with `MODEL_NAME`. The available `MODEL_NAME`s are listed in [HKU-provided Models](#hku-provided-models). If you use models outside the image, you need to mount the folder containing the models using `-v` when running Docker, and then modify `--clair3_model_path`. For usage, please refer to [Option 1. Docker pre-built image](#option-1-docker-pre-built-image).

The following lists the available pre-trained models of Clair3. When specifying `--clair3_model_path`, be sure to modify `--platform` accordingly.

### HKU-provided Models

Download models from [here](http://www.bio8.cs.hku.hk/clair3/clair3_models/) or click on the links below.


|           Model name           |  Platform   | Option (`--platform`) |                       Training samples                       |   Date   |
| :----------------------------: | :---------: | :----------------------------------------------------------: | -------------------------------- | :------: |
|      r941_prom_sup_g5014 / ont_guppy5       |     ONT r9.4.1     |     `ont`     |                    HG002,4,5 (Guppy5_sup)                    | 20220112 |
|    r941_prom_hac_g360+g422 / ont     |     ONT r9.4.1    | `ont`    |                         HG001,2,4,5                          | 20210517 |
|       r941_prom_hac_g238 / ont_guppy2       |     ONT r9.4.1    | `ont`    |                         HG001,2,3,4                          | 20210627 |
|              hifi_revio              | PacBio HiFi Revio | `hifi` |                         HG002,4                         | 20230522 |
|             hifi_sequel2 / hifi             | PacBio HiFi Sequel II | `hifi` |                         HG001,2,4,5                          | 20210517 |
| ilmn | Illumina | `ilmn` | HG001,2,4,5 | 20210517 |

### ONT-provided Models

ONT provides models for some latest or specific chemistries and basecallers (including both Guppy and Dorado) through [Rerio](https://github.com/nanoporetech/rerio). These models are tested and supported by the ONT developers.

## Acknowledgement

ClusterV is designed and built by HKUCS BAL lab (HKUCS Bioinformatics Algorithm Lab).

We apply and thank the following tools during the development of ClusterV:

- [HIVdb](https://hivdb.stanford.edu/)
- [IGV](https://software.broadinstitute.org/software/igv/)
- [Flye](https://github.com/fenderglass/Flye)
- [minimap2](https://github.com/lh3/minimap2)
- [parallel](https://manpages.ubuntu.com/manpages/impish/man1/parallel.1.html)

## Publication

- Ng, Timothy Ting-Leung, Junhao Su, Hiu-Yin Lao, Wui-Wang Lui, Chloe Toi-Mei Chan, Amy Wing-Sze Leung, Stephanie Hoi-Ching Jim et al. "Long-Read Sequencing with Hierarchical Clustering for Antiretroviral Resistance Profiling of Mixed Human Immunodeficiency Virus Quasispecies." Clinical Chemistry 69, no. 10 (2023): 1174-1185.

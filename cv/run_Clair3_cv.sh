#!/usr/bin/env bash
# run_Clair3_cv.sh BAM SAMPLE_ID REF.fasta REF.bed

BAM_FILE_PATH=$1
SAMPLE_NAME=$2
REFERENCE_FASTA_FILE_PATH=$3
BED_F=$4
WORKING_DIRECTORY=$5
PARALLEL_THREADS=$6
CLAIR_PLATFORM=$7

echo "running Clair3"
echo "input bam: $BAM_FILE_PATH"
echo "sample ID: $SAMPLE_NAME"
echo "input ref: $REFERENCE_FASTA_FILE_PATH"
echo "input bed: $BED_F"
echo "parallel threads: $PARALLEL_THREADS"
echo "clair_platform: $CLAIR_PLATFORM" 
echo "working output dir: $WORKING_DIRECTORY"

DIR=$( dirname -- "$0"; )
#echo ${DIR}

CLAIR_PATH="${DIR}/../Clair3"
CLAIR_MODEL="${CLAIR_PATH}/models/${CLAIR_PLATFORM}"

mkdir -p ${WORKING_DIRECTORY}

${CLAIR_PATH}/run_clair3.sh \
  --bam_fn="${BAM_FILE_PATH}" \
  --ref_fn="${REFERENCE_FASTA_FILE_PATH}" \
  --bed_fn="${BED_F}" \
  --threads=16 \
  --platform="${CLAIR_PLATFORM}" \
  --model_path="${CLAIR_MODEL}" \
  --sample_name="${SAMPLE_NAME}" \
  --pileup_only \
  --chunk_size=200 \
  --output="${WORKING_DIRECTORY}"


echo "output vcf file at ${WORKING_DIRECTORY}/plieup.vcf.gz"



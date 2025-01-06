#!/usr/bin/env bash
# run_Clair3_cv.sh BAM SAMPLE_ID REF.fasta REF.bed

BAM_FILE_PATH=$1
SAMPLE_NAME=$2
REFERENCE_FASTA_FILE_PATH=$3
BED_F=$4
WORKING_DIRECTORY=$5
PARALLEL_THREADS=$6
CLAIR_PLATFORM=$7
CLAIR3_MODEL_PATH=$8
HAPLOID_PRECISE=$9
HAPLOID_SENSITIVE=${10}

DIR=$( dirname -- "$0"; )
# echo "DIR: ${DIR}"

echo "running Clair3"
echo "input bam: $BAM_FILE_PATH"
echo "sample ID: $SAMPLE_NAME"
echo "input ref: $REFERENCE_FASTA_FILE_PATH"
echo "input bed: $BED_F"
echo "parallel threads: $PARALLEL_THREADS"
echo "clair platform: $CLAIR_PLATFORM" 
# echo "HAPLOID_PRECISE: $HAPLOID_PRECISE"
# echo "HAPLOID_SENSITIVE: $HAPLOID_SENSITIVE"

if [[ $CLAIR3_MODEL_PATH == ../* ]]; then
    CLAIR3_MODEL_PATH=${CLAIR3_MODEL_PATH/../$(echo $DIR)\/..}
else
    CLAIR3_MODEL_PATH=$CLAIR3_MODEL_PATH
fi
echo "clair3 model path: $CLAIR3_MODEL_PATH"
echo "working output dir: $WORKING_DIRECTORY"
mkdir -p ${WORKING_DIRECTORY}


clair3_args=("--bam_fn=${BAM_FILE_PATH}" "--ref_fn=${REFERENCE_FASTA_FILE_PATH}" "--bed_fn=${BED_F}" \
"--threads=${PARALLEL_THREADS}" "--platform=${CLAIR_PLATFORM}" "--model_path=${CLAIR3_MODEL_PATH}" \
"--sample_name=${SAMPLE_NAME}" "--pileup_only" "--chunk_size=200" "--output=${WORKING_DIRECTORY}")

if [[ "$HAPLOID_PRECISE" == "True" ]]; then
    clair3_args+=("--haploid_precise")
fi

if [[ "$HAPLOID_SENSITIVE" == "True" ]]; then
    clair3_args+=("--haploid_sensitive")
fi

CLAIR_PATH="${DIR}/../Clair3"
${CLAIR_PATH}/run_clair3.sh "${clair3_args[@]}"


echo "output vcf file at ${WORKING_DIRECTORY}/plieup.vcf.gz"



#!/usr/bin/env bash
# run_Clair_ensemble_cv.sh BAM SAMPLE_ID REF.fasta REF.bed

BAM_FILE_PATH=$1
SAMPLE_NAME=$2
REFERENCE_FASTA_FILE_PATH=$3
BED_F=$4
WORKING_DIRECTORY=$5

#PARALLEL_THREADS="16"
PARALLEL_THREADS=$6
echo "running Clair-Ensemble"
echo "input bam: $BAM_FILE_PATH"
echo "sample ID: $SAMPLE_NAME"
echo "input ref: $REFERENCE_FASTA_FILE_PATH"
echo "input bed: $BED_F"
echo "parallel threads: $PARALLEL_THREADS"
echo "working output dir: $WORKING_DIRECTORY"

DIR=$( dirname -- "$0"; )
#echo ${DIR}

CLAIR_PATH=$DIR/../Clair-ensemble/
#"$CLAIR_PATH/model/model-000016"
CLAIR_MODELS=(
"$CLAIR_PATH/model/guppy5/model-000015"
)

#"$CLAIR_PATH/model/model"
CLAIR="$CLAIR_PATH/Clair.beta.ensemble.cpu/clair.py"
CLAIR_HOME_DIR=`dirname ${CLAIR}`
ENSEMBLE_CPP_EXECUTABLE="${CLAIR_HOME_DIR}/clair/ensemble"
RETRIES=3
MIN_AF=0.1
DCOV=(2000)

#set -x
SCRIPT_PATH=`readlink -f $0`
SCRIPT_DIR_PATH=`dirname $SCRIPT_PATH`
#set +x
NO_OF_CLAIR_MODELS=${#CLAIR_MODELS[@]}
DATE_TIME=`date "+%Y%m%d_%H%M%S"`
#WORKING_DIRECTORY="${ROOT_FOLDER_PATH}/`basename $1`_dir"
mkdir -p ${WORKING_DIRECTORY}
BAM_FILE_PATHS=("${BAM_FILE_PATH}")
INTERMEDIATE_OUTPUT_FOLDER="${WORKING_DIRECTORY}/tmp_output"
mkdir -p ${INTERMEDIATE_OUTPUT_FOLDER}
cd ${INTERMEDIATE_OUTPUT_FOLDER}
for i in "${!BAM_FILE_PATHS[@]}"
do
  INPUT_BAM_FILE_PATH="${BAM_FILE_PATHS[i]}"
  #SAMPLE_NAME="HIV"
  BAM_PREFIX=`printf "%02d" ${i}`

  for j in "${!CLAIR_MODELS[@]}"
  do
    CLAIR_MODEL="${CLAIR_MODELS[j]}"
    MODEL_PREFIX=`printf "%02d" ${j}`
    SCRIPT_OUTPUT_FOLDER="m${MODEL_PREFIX}_b${BAM_PREFIX}"

    mkdir -p ${SCRIPT_OUTPUT_FOLDER}
    OUTPUT_PREFIX="${SCRIPT_OUTPUT_FOLDER}/tmp"

    python ${CLAIR} callVarBamParallel \
    --chkpnt_fn "${CLAIR_MODEL}" \
    --ref_fn "${REFERENCE_FASTA_FILE_PATH}" \
    --bam_fn "${INPUT_BAM_FILE_PATH}" \
    --pysam_for_all_indel_bases \
    --output_for_ensemble \
    --refChunkSize "200" \
    --includingAllContigs \
    --threshold=$MIN_AF \
    --dcov="${DCOV[j]}" \
    --bed_fn=$BED_F \
    --sampleName "${SAMPLE_NAME}" \
    --output_prefix "${OUTPUT_PREFIX}" > ${SCRIPT_OUTPUT_FOLDER}/call.sh
    #--output_prefix "${OUTPUT_PREFIX}"
  done
done

( time cat */call.sh | parallel -j ${PARALLEL_THREADS} --retries ${RETRIES} -C ' ' --joblog ${INTERMEDIATE_OUTPUT_FOLDER}/1_out_call_pro.log ) |& tee ${INTERMEDIATE_OUTPUT_FOLDER}/log.call.txt

# Find incomplete VCF files and rerun them
for i in "${OUTPUT_PREFIX}".*.vcf; do if ! [ -z "$(tail -c 1 "$i")" ]; then echo "$i"; fi ; done | grep -f - m00_b00/call.sh | sh


FILES=(`ls m*_b00/*.vcf`)
#echo ${FILES[@]}
ENSEMBLE_OUTPUT_FOLDER="${INTERMEDIATE_OUTPUT_FOLDER}/ensemble"
mkdir -p ${ENSEMBLE_OUTPUT_FOLDER}
MININUM_NO_OF_VOTE_FOR_VARIANT="$(((${#BAM_FILE_PATHS[@]}*${#CLAIR_MODELS[@]}+2)/2))"
#rm ensemble_command.sh
for i in "${!FILES[@]}"
do
  TARGET_FILE_NAME="${FILES[i]:8}"
  CAT_COMMAND=""
  for j in "${!BAM_FILE_PATHS[@]}"
  do
    BAM_PREFIX=`printf "%02d" ${j}`
    for k in "${!CLAIR_MODELS[@]}"
    do
      MODEL_PREFIX=`printf "%02d" ${k}`
      FOLDER_NAME="m${MODEL_PREFIX}_b${BAM_PREFIX}"
      CAT_COMMAND="${CAT_COMMAND} ${FOLDER_NAME}/${TARGET_FILE_NAME}"
    done
  done


  echo "cat ${CAT_COMMAND:1} | ${ENSEMBLE_CPP_EXECUTABLE} ${MININUM_NO_OF_VOTE_FOR_VARIANT} > ${ENSEMBLE_OUTPUT_FOLDER}/${TARGET_FILE_NAME}" >> ensemble_command.sh
done
( time cat ensemble_command.sh | parallel -j${PARALLEL_THREADS} --retries ${RETRIES} -C ' ' --joblog ${INTERMEDIATE_OUTPUT_FOLDER}/2_out_ensemble.log ) |& tee ${INTERMEDIATE_OUTPUT_FOLDER}/log.ensemble.txt


VCF_OUTPUT_FOLDER="${WORKING_DIRECTORY}/output"
mkdir -p ${VCF_OUTPUT_FOLDER}
cd ${WORKING_DIRECTORY}
INPUT_FILES=(`ls tmp_output/ensemble/*.vcf`)

for i in "${!INPUT_FILES[@]}"
do
  FILE_NAME="${INPUT_FILES[i]:20}"
  echo "cat tmp_output/ensemble/${FILE_NAME} | \
  python ${CLAIR} call_var \
  --chkpnt_fn "${CLAIR_MODEL}" \
  --ref_fn "${REFERENCE_FASTA_FILE_PATH}" \
  --bam_fn "${BAM_FILE_PATH}" \
  --call_fn "${VCF_OUTPUT_FOLDER}/${FILE_NAME}" \
  --sampleName "${SAMPLE_NAME}" \
  --pysam_for_all_indel_bases \
  --input_probabilities" >> output.sh

 #--haploid \
done
( time cat output.sh | parallel -j${PARALLEL_THREADS} --retries ${RETRIES} -C ' ' --joblog ${INTERMEDIATE_OUTPUT_FOLDER}/3_out_callvar.log ) |& tee ${INTERMEDIATE_OUTPUT_FOLDER}/log.output.txt


cd ${WORKING_DIRECTORY}
#vcfcat ${VCF_OUTPUT_FOLDER}/*.vcf | sort -k1,1V -k2,2n > snp_and_indel.vcf
vcfcat ${VCF_OUTPUT_FOLDER}/*.vcf > all_cat.vcf
vcfsort all_cat.vcf > snp_and_indel.vcf
cat snp_and_indel.vcf | python ${CLAIR_PATH}/Clair.beta.ensemble.cpu/clair/post_processing/overlap_variant.py > snp_and_indel.filtered_ori.vcf
#cat snp_and_indel.vcf | python ${CLAIR_PATH}/Clair.beta.ensemble.cpu/clair/post_processing/overlap_variant.py | bcftools filter -e "FORMAT/GQ<50||FORMAT/AF<=${MIN_AF}" > snp_and_indel.filtered.vcf
#cp snp_and_indel.filtered.vcf ${SAMPLE_NAME}.vcf

_CV_PY=${DIR}/../cv.py 
#python ${_CV_PY} update_vcf_af --vcf_fn snp_and_indel_ori.vcf --bam "${INPUT_BAM_FILE_PATH}" --bed ${BED_F} --out_fn snp_and_indel.filtered_ori_u.vcf &> /dev/null
echo "python ${_CV_PY} update_vcf_af --vcf_fn ${WORKING_DIRECTORY}/snp_and_indel.filtered_ori.vcf --bam "${BAM_FILE_PATH}" --bed ${BED_F} --out_fn snp_and_indel.filtered_ori_u.vcf"
python ${_CV_PY} update_vcf_af --vcf_fn ${WORKING_DIRECTORY}/snp_and_indel.filtered_ori.vcf --bam "${BAM_FILE_PATH}" --bed ${BED_F} --out_fn snp_and_indel.filtered_ori_u.vcf &> _update_af_log
#bcftools filter -e "FORMAT/GQ<50||FORMAT/AF<=${MIN_AF}" snp_and_indel.filtered_ori_u.vcf > out.vcf
bcftools filter -e "FORMAT/GQ<50||FORMAT/AF<=${MIN_AF}" snp_and_indel.filtered_ori_u.vcf | bcftools norm --rm-dup all > out.vcf 
bcftools filter -e "FORMAT/GQ<748||FORMAT/AF<=${MIN_AF}" snp_and_indel.filtered_ori_u.vcf | bcftools norm --rm-dup all > out_1.vcf 
#rm all_cat.vcf
#rm snp_and_indel.vcf
#rm snp_and_indel.filtered_ori.vcf
#rm snp_and_indel.filtered_ori_u.vcf
echo "output vcf file at ${WORKING_DIRECTORY}/out.vcf"



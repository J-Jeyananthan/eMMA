#!/bin/bash
#$ -l tmem=32G
#$ -l h_vmem=32G
#$ -l h_rt=12:0:0
#$ -S /bin/bash
#$ -j y
#$ -N emma-dist-matrix
#$ -cwd
#$ -P cath
#$ -e /dev/null
#$ -o /dev/null
#$ -t 1-4

DATADIR=/SAN/orengolab/functional-families/janu/contrasted-ff/HUPs_data/3.20.20.70
PARENTDIR=/SAN/orengolab/functional-families/janu/contrasted-ff/HUPs_data
EMMA_EXP=${DATADIR}/eMMA_experiments
VENV=${PARENTDIR}/venv_emma/bin/activate

source ${VENV}

PROJECTS=(esm2 prostt5 esm2_contrastive prostt5_contrastive)
EMB_FILES=(
    "${DATADIR}/3.20.20.70_sc_reps_esm2.pt"
    "${DATADIR}/3.20.20.70_sc_reps_prostt5.pt"
    "${DATADIR}/esm2_contrastive_deduped/E_starting_contrastive.pt"
    "${DATADIR}/prostt5_contrastive/E_starting_contrastive.pt"
)

IDX=$((SGE_TASK_ID - 1))
PROJECT=${PROJECTS[$IDX]}
EMB_FILE=${EMB_FILES[$IDX]}

echo $(date) ${PROJECT} START

cath-emma-cli create-distance-matrix \
    --distance_source embeddings \
    --input_to_process ${EMB_FILE} \
    --embedding_distance euclidean \
    --matrix_output_file ${EMMA_EXP}/${PROJECT}/${PROJECT}_embedding_matrix.ssv \
    --labels_file ${DATADIR}/3.20.20.70_sc_reps_ids_list

echo $(date) ${PROJECT} END

#!/bin/bash
#$ -l tmem=16G
#$ -l h_vmem=16G
#$ -l h_rt=12:0:0
#$ -S /bin/bash
#$ -j y
#$ -N emma-run
#$ -cwd
#$ -P cath
#$ -e /dev/null
#$ -o /dev/null
#$ -t 1-4

DATADIR=/SAN/orengolab/functional-families/janu/contrasted-ff/HUPs_data/3.20.20.70
PARENTDIR=/SAN/orengolab/functional-families/janu/contrasted-ff/HUPs_data
EMMA_EXP=${DATADIR}/eMMA_experiments
PLENV=${PARENTDIR}/perl_local/bin
CODE_DIR=${PARENTDIR}/eMMA

source ~/.bash_profile

PROJECTS=(esm2 prostt5 esm2_contrastive prostt5_contrastive)
IDX=$((SGE_TASK_ID - 1))
PROJECT=${PROJECTS[$IDX]}

echo $(date) ${PROJECT} START

cd ${EMMA_EXP}/${PROJECT}

${PLENV}/perl ${CODE_DIR}/Cath-Gemma/script/prepare_research_data.pl \
    --local \
    --projects-list-file projects.txt \
    --output-root-dir . \
    --embs-file ${EMMA_EXP}/${PROJECT}/${PROJECT}_embedding_matrix.ssv \
    1> ${EMMA_EXP}/${PROJECT}/${PROJECT}.stdout \
    2> ${EMMA_EXP}/${PROJECT}/${PROJECT}.stderr

echo $(date) ${PROJECT} END

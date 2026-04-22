#!/bin/bash
# Run this interactively on the login node to prepare eMMA experiment directories for 3.40.50.620.

DATADIR=/SAN/orengolab/functional-families/janu/contrasted-ff/HUPs_data/HUPs_mmseqs_s90_3.40.50.620
PARENTDIR=/SAN/orengolab/functional-families/janu/contrasted-ff/HUPs_data
EMMA_EXP=${DATADIR}/eMMA_experiments
REPS_CSV=${DATADIR}/HUPs_mmseqs_s90_3.40.50.620_reps.csv
VENV=${PARENTDIR}/venv_emma/bin/activate

module load python/3.8.5
source ${VENV}

for PROJECT in esm2_contrastive prostt5_contrastive; do
    mkdir -p ${EMMA_EXP}/${PROJECT}/starting_clusters/${PROJECT}
    echo ${PROJECT} > ${EMMA_EXP}/${PROJECT}/projects.txt

    cath-emma-cli create-starting-clusters-from-centroids \
        --cluster_reps_file ${REPS_CSV} \
        --starting_clusters_dir ${EMMA_EXP}/${PROJECT}/starting_clusters/${PROJECT} \
        --cluster_mapping_file ${EMMA_EXP}/${PROJECT}/${PROJECT}_cluster_mapping.tsv
done

echo "Setup complete."

#!/bin/bash
# Run this interactively on the login node to prepare all 4 eMMA experiment directories.

DATADIR=/SAN/orengolab/functional-families/janu/contrasted-ff/HUPs_data/3.20.20.70
PARENTDIR=/SAN/orengolab/functional-families/janu/contrasted-ff/HUPs_data
EMMA_EXP=${DATADIR}/eMMA_experiments
VENV=${PARENTDIR}/venv_emma/bin/activate

module load python/3.8.5
source ${VENV}

# Create experiment directories
for PROJECT in esm2 prostt5 esm2_contrastive prostt5_contrastive; do
    mkdir -p ${EMMA_EXP}/${PROJECT}/starting_clusters/${PROJECT}
done

# Generate reps CSV from fasta (shared across all 4 runs)
cath-emma-cli convert-fasta-to-csv-for-embed \
    -i ${DATADIR}/3.20.20.70_sc_reps.fasta \
    -o ${DATADIR}/3.20.20.70_sc_reps.csv

# Extract IDs list (shared across all 4 runs)
awk -F ',' '{print $1}' ${DATADIR}/3.20.20.70_sc_reps.csv \
    > ${DATADIR}/3.20.20.70_sc_reps_ids_list

# Create starting clusters and projects.txt for each experiment
for PROJECT in esm2 prostt5 esm2_contrastive prostt5_contrastive; do
    echo ${PROJECT} > ${EMMA_EXP}/${PROJECT}/projects.txt

    cath-emma-cli create-starting-clusters-from-centroids \
        --cluster_reps_file ${DATADIR}/3.20.20.70_sc_reps.csv \
        --starting_clusters_dir ${EMMA_EXP}/${PROJECT}/starting_clusters/${PROJECT} \
        --cluster_mapping_file ${EMMA_EXP}/${PROJECT}/${PROJECT}_cluster_mapping.tsv
done

echo "Setup complete."

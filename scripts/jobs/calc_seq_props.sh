#!/bin/bash
#SBATCH --job-name=calc_seq_props_CL_SK_plant
#SBATCH --account=ac_stallerlab
#SBATCH --partition=savio3
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=1:00:00
#SBATCH --output=../logs/calc_seq_props_%j.log

module load anaconda3/2024.02-1-11.4
source activate /global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/conda/parrot

# # Run Python script with arguments
# python /global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/consensus_AD_predictor/scripts/calculate_properties.py \
#     --input /global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/consensus_AD_predictor/output/interpro_uniprot_evid_1_2.fasta \
#     --output /global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/consensus_AD_predictor/output/interpro_uniprot_evid_1_2_seq_props.csv

python /global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/consensus_AD_predictor/scripts/calculate_properties.py \
    --input /global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/consensus_AD_predictor/data/CL_SK_plant_input.fasta \
    --output /global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/consensus_AD_predictor/output/CL_SK_plant_seq_props.csv
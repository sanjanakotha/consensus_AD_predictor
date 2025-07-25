#!/bin/bash
#SBATCH --job-name=allpredictors_len5_test
#SBATCH --account=ac_stallerlab
#SBATCH --partition=savio3
#SBATCH --time=1:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sanjana.kotha@berkeley.edu
#SBATCH --output=logs/%x_%j.out

source activate sk_allpredictors
cd /global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/consensus_AD_predictor/scripts/
python run-all-predictors.py -f /global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/consensus_AD_predictor/data/all_predictors_test_50.fasta -o /global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/consensus_AD_predictor/output/all_predictors_test
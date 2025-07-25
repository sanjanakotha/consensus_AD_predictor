#!/bin/bash
# Job name:
#SBATCH --job-name=allpredictors_test
#
# Account:
#SBATCH --account=ac_stallerlab
#
# Partition:
#SBATCH --partition=savio3
#
# Wall clock limit:
#SBATCH --time=1:00:00
#
## Command(s) to run:
source activate sk_allpredictors

python /global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/consensus_AD_predictor/scripts/run-all-predictors.py -f global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/consensus_AD_predictor/data/test_long_TFs_TADA.fasta -o /global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/consensus_AD_predictor/output/all_predictors_test_long

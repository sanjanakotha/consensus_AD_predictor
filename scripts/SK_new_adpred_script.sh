#!/bin/bash
# Job name:
#SBATCH --job-name=adpred_yeast_TFs
#
# Account:
#SBATCH --account=ac_stallerlab
#
# Partition:
#SBATCH --partition=savio2
#
# Wall clock limit:
#SBATCH --time=72:00:00
#
## Command(s) to run:
source activate sk_allpredictors

python /global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/consensus_AD_predictor/scripts/run_adpred.py


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
#SBATCH --time=1:00:00
#
## Command(s) to run:
source activate /global/scratch/projects/fc_mvslab/conda/adpred

python /global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/misc_AD_predictors/adpred_files/use_ADPred.py -c ProteinWindowSeq -p yes /global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/consensus_AD_predictor/data/test_tiled_yeast_TFs_for_adpred.csv


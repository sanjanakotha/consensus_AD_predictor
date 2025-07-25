#!/bin/bash
# Job name:
#SBATCH --job-name=TADA_yeast_TFs
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

ml python/3.7; source activate /global/scratch/projects/fc_mvslab/conda/tada
cd /global/scratch/projects/fc_mvslab/predictors/TADA
python src/Predict.py /global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/consensus_AD_predictor/data/tiled_yeast_TFs_for_tada.csv

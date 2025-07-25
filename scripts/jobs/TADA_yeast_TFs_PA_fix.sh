#!/bin/bash
#SBATCH --job-name=TADA_PA_fix_yeast
#SBATCH --account=ac_stallerlab
#SBATCH --partition=savio3
#SBATCH --time=24:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sanjana.kotha@berkeley.edu

source activate sk_allpredictors
cd /global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/consensus_AD_predictor/scripts/
python run_TADA.py --input /global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/consensus_AD_predictor/data/yeast_TFs.fasta --output /global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/consensus_AD_predictor/output/tada_pa_fix_yeast
#!/bin/bash
#SBATCH --job-name=TADA_PA_fix_human
#SBATCH --account=ac_stallerlab
#SBATCH --partition=savio3
#SBATCH --time=72:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sanjana.kotha@berkeley.edu

source activate sk_allpredictors
cd /global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/consensus_AD_predictor/scripts/
python run_TADA.py --input /global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/consensus_AD_predictor/data/lambert_TFs_10-21-24_with_DBD_coords.fasta --output /global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/consensus_AD_predictor/output/tada_pa_fix_human
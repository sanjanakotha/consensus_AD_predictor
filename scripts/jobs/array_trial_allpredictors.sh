#!/bin/bash
#SBATCH --job-name=all_predictors_two_parts
#SBATCH --account=ac_stallerlab
#SBATCH --partition=savio3
#SBATCH --time=72:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sanjana.kotha@berkeley.edu
#SBATCH --array=0-1
#SBATCH --output=array_interpro_1_2_%A_task_%a.out

source activate sk_allpredictors
cd /global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/consensus_AD_predictor/scripts/
python /global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/consensus_AD_predictor/scripts/run-all-predictors.py -f /global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/consensus_AD_predictor/data/interpro_uniprot_evid_1_2_split/pt_$SLURM_ARRAY_TASK_ID.fasta -o /global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/consensus_AD_predictor/output/interpro_uniprot_evid_1_2_split_preds
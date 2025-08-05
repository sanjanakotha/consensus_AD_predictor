#!/bin/bash
# Job name:
#SBATCH --job-name=adpred_split_input
#
# Account:
#SBATCH --account=ac_stallerlab
#
# Partition:
#SBATCH --partition=savio3
#
# Wall clock limit:
#SBATCH --time=72:00:00
#
#SBATCH --output=../logs/%x_%j_task_%a.out
#SBATCH --array=77,87
## Command(s) to run:

echo "Starting job on $(hostname)"
echo "Loading modules..."
module load anaconda3/2024.02-1-11.4

echo "Activating environment..."
source activate /global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/conda/new_dask

echo "Python path: $(which python)"
echo "Python version: $(python --version)"

in_file="/global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/consensus_AD_predictor/data/interpro_uniprot_evid_1_2_split/pt_$SLURM_ARRAY_TASK_ID.fasta"
out_dir="/global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/consensus_AD_predictor/output/adpred_on_unfinished_parts/pt_$SLURM_ARRAY_TASK_ID"
script="/global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/consensus_AD_predictor/scripts/run_adpred.py"

echo "Script exists? $(ls -l $script)"
echo "Input exists? $(ls -l $in_file)"
echo "Output dir exists? $(ls -ld $out_dir)"

echo "Running: python $script -f $in_file -o $out_dir"
srun /global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/conda/new_dask/bin/python "$script" -f "$in_file" -o "$out_dir"

echo "Done."

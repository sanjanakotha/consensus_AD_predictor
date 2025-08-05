#!/bin/bash
#SBATCH --job-name=TADA_test
#SBATCH --account=ac_stallerlab
#SBATCH --partition=savio2_htc
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=1:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sanjana.kotha@berkeley.edu
#SBATCH --output=../logs/%x_%j.out

unset SLURM_MEM_PER_NODE

echo "Starting job on $(hostname)"
echo "Loading modules..."
module load anaconda3/2024.02-1-11.4

echo "Activating environment..."
source activate /global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/conda/new_dask

echo "Python path: $(which python)"
echo "Python version: $(python --version)"

in_file="/global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/consensus_AD_predictor/data/test_yeast_TFs_for_tada.fasta"
out_dir="/global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/consensus_AD_predictor/output/dask_scripts_test"
script="/global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/consensus_AD_predictor/scripts/run_TADA.py"

echo "Script exists? $(ls -l $script)"
echo "Input exists? $(ls -l $in_file)"
echo "Output dir exists? $(ls -ld $out_dir)"

echo "Running: python $script --input $in_file --output $out_dir"
srun /global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/conda/new_dask/bin/python "$script" --input "$in_file" --output "$out_dir"

echo "Done."

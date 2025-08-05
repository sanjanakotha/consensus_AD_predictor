#!/bin/bash
# Job name:
#SBATCH --job-name=dask_tada_test_TFs
#
# Account:
#SBATCH --account=fc_mvslab
#
# Partition:
#SBATCH --partition=savio3
#
# Number of nodes:
#SBATCH --nodes=1
#
# Number of tasks
#SBATCH --ntasks=7
#
# Processors per task:
#SBATCH --cpus-per-task=4
#
# Wall clock limit:
#SBATCH --time=1:00:00
#
#SBATCH --export=ALL
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sanjana.kotha@berkeley.edu
#SBATCH --output=../logs/%x_%j.out
## Command(s) to run:

echo "Starting job on $(hostname)"
echo "Loading modules..."
module load anaconda3/2024.02-1-11.4

echo "Activating environment..."
source activate /global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/conda/new_dask

echo "Python path: $(which python)"
echo "Python version: $(python --version)"

in_file="/global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/consensus_AD_predictor/data/test_yeast_TFs_for_tada.fasta"
out_dir="/global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/consensus_AD_predictor/output/dask_scripts_test"
script="/global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/consensus_AD_predictor/scripts/run_tada_dask.py"

echo "Script exists? $(ls -l $script)"
echo "Input exists? $(ls -l $in_file)"
echo "Output dir exists? $(ls -ld $out_dir)"

echo "Running: python $script -f $in_file -o $out_dir"
srun /global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/conda/new_dask/bin/python "$script" -f "$in_file" -o "$out_dir"

echo "Done."

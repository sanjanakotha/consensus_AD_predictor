#!/bin/bash
# Job name:
#SBATCH --job-name=tada_input_split
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
#SBATCH --ntasks=1
#
# Processors per task:
#SBATCH --cpus-per-task=32
#
# Wall clock limit:
#SBATCH --time=4:00:00
#
#SBATCH --output=../logs/%x_%j_%a.out
#SBATCH --array=136,137,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127

## Command(s) to run:

# module load openjdk/17.0.8.1_1-gcc-11.4.0

echo "Total CPUs: $SLURM_CPUS_ON_NODE"

echo "Starting job on $(hostname)"
echo "Loading modules..."
module load anaconda3/2024.02-1-11.4

echo "Activating environment..."
source activate /global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/conda/new_dask

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export BLIS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export TF_NUM_INTRAOP_THREADS=1
export TF_NUM_INTEROP_THREADS=1

in_file="/global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/consensus_AD_predictor/data/interpro_uniprot_evid_1_2_split/pt_$SLURM_ARRAY_TASK_ID.fasta"
out_dir="/global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/consensus_AD_predictor/output/interpro_uniprot_evid_1_2_dask"
script="/global/scratch/projects/fc_mvslab/predictors/TADA/run_tada_new_parallel.py"

echo "Running: python $script -f $in_file -o $out_dir"
python "$script" -f "$in_file" -o "$out_dir"
echo "Done."

#!/bin/bash
# Job name:
#SBATCH --job-name=parrot_kinetics
#
# Account:
#SBATCH --account=ac_stallerlab
#
# Partition:
#SBATCH --partition=savio2
#
# Wall clock limit:
#SBATCH --time=12:00:00
#
## Command(s) to run:
source activate /global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/conda/parrot

parrot-train '../output/kinetics_parrot_input.tsv' '../output/parrot_kinetics' -d 'sequence' -c 1 --include-figs

#!/bin/bash
# Job name:
#SBATCH --job-name=parrot_human_yeast
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

parrot-train '../output/human_yeast_aggreg_preds.tsv' '../output/parrot_human_yeast_aggreg' -d 'residues' -c 1 --include-figs

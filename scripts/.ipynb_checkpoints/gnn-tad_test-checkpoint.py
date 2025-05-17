import subprocess

command = """
source ~/.bashrc
conda activate gnntad
cd /global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/consensus_AD_predictor/models/GNN-TAD
python3 main.py --single_seq=IGIRTIVADVGISVPFVTIDVGVEEFYCMI --multi_mode=0 --mode=63 --type=1 --modelpath=best_models/model_53.pth.tar --threshold=0.85 --gpu=1
"""

subprocess.run(command, shell=True, executable="/bin/bash")

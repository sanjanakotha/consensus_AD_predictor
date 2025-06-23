# Load necessary modules
import numpy as np
import pandas as pd
import tensorflow as tf
from Bio import SeqIO
import glob, pickle, os, re, subprocess, csv, sys
os.chdir("/global/scratch/projects/fc_mvslab/predictors/TADA/")
from uuid import uuid4
import ipyparallel as ipp
sys.path.append(os.path.abspath('..'))
from TADA.src.Preprocessing_FIXED import scale_features_predict
from TADA.src.Preprocessing import create_features
from TADA.src.Predict_model_only import tada_predict
# Set Numpy to display floats with 3 decimal places
np.set_printoptions(formatter={'float': lambda x: "{0:0.3f}".format(x)})
import argparse

# Argument parsing
parser = argparse.ArgumentParser(description="Run TADA predictions on input FASTA.")
parser.add_argument("--input", "-i", required=True, help="Path to input FASTA file")
parser.add_argument("--output", "-o", required=True, help="Path to output directory")
args = parser.parse_args()

# Input/output handling
fasta_name = args.input
tmp_output = args.output

# # INPUT PATH TO FASTA FILE HERE
# fasta_name = '/global/scratch/projects/fc_mvslab/predictors/all_predictors/example/example.fasta'
# tmp_output = '/global/scratch/projects/fc_mvslab/predictors/TADA/tmp_output'

aa_lst = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']

# Set up SeqRecord objects and check that sequences are valid
#fasta_name_tmp = fasta_name.split('/')[-1].strip('.fasta').strip('.fa')
output_dir = tmp_output #+ '/' + fasta_name_tmp
subprocess.run(['mkdir', output_dir])
recs = list(SeqIO.parse(fasta_name, 'fasta'))
#fasta_name = fasta_name_tmp
for r in recs:
    sequence = str(r.seq)
    assert len(sequence) >= 40, 'Protein sequence must be at least 40 amino acids long' 
    assert len(set(sequence) - set(aa_lst)) == 0, 'No non-standard amino acids allowed'


# Defines the sequence window size and steps (stride length). Change values if needed.
SEQUENCE_WINDOW = 5
STEPS = 1
LENGTH = 40

seqs_40 = []
seqs_idx = 0
tada_out = dict()

#Iterate through all sequences
for r in recs:
    seq = str(r.seq)
    
    # Chop sequence into 40-mer tiles
    seqs = np.array([seq[i:i+40] for i in range(len(seq)-39)], dtype=object)
    seqs_40.append(seqs)
    
    # Record indices of tiles corresponding to the sequence
    tada_out[seq] = np.arange(seqs_idx, seqs_idx+len(seqs))
    seqs_idx += len(seqs)

# Flatten list of 40-mer tiles
seqs_40 = np.concatenate(seqs_40)

assert seqs_idx == len(seqs_40)

#Calculate features on 40-mer tiles
features = create_features(seqs_40, SEQUENCE_WINDOW, STEPS)
features_scaled = scale_features_predict(features)
pickle.dump(features_scaled, open(output_dir + '/features_scaled.pkl', 'wb'))
    
# Make TADA classification predictions
tada_preds = tada_predict(loaded_features=features_scaled)

# Use indices of tiles to retrieve predictions for each sequence
for seq, val in tada_out.items():
    preds = tada_preds[val]
    centers = np.arange(len(preds)) + 40/2
    tada_out[seq] = (centers, preds)

# Save predicted values to a csv file
with open(output_dir + '/TADA_preds.csv', 'w') as f:
    w = csv.writer(f)
    w.writerow(['sequence', 'tada_centers', 'tada_preds'])
    for r in recs:
        sequence = str(r.seq)
        tada_centers, tada_preds = tada_out[sequence]
        tada_preds = re.sub(r'\s+', ',', str(tada_preds))
        tada_centers = re.sub(r'\s+', ',', str(tada_centers))
        w.writerow([sequence, tada_centers, tada_preds])
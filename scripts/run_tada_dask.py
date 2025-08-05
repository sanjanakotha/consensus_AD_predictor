# Load necessary modules
import numpy as np
import pandas as pd
from Bio import SeqIO
import glob, pickle, os, re, subprocess, csv, sys, io
import torch
from sklearn import preprocessing
import argparse
from contextlib import redirect_stdout
import dask
from dask_jobqueue.slurm import SLURMRunner
from dask.distributed import Client, progress
from dask.diagnostics import ProgressBar
import sys

sys.path.append(os.path.abspath('/global/scratch/projects/fc_mvslab/predictors/TADA'))
from create_features_fast import calculate_properties

sys.path.append(os.path.abspath('/global/scratch/projects/fc_mvslab/predictors'))
from TADA.src.Preprocessing_FIXED import scale_features_predict
from tqdm import tqdm

sys.path.append(os.path.abspath('/global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/consensus_AD_predictor/scripts/'))
from create_features_fast_dask import create_features_faster


#### Setting up the parallel workers
cluster = SLURMRunner()
client = Client(cluster)
client.wait_for_workers(5)

# The goal of this is to use a individual paddle model for each worker
from dask.distributed import get_worker, WorkerPlugin


# This should prevent the model from being loaded every tme
class TadaPlugin(WorkerPlugin):
    def setup(self, worker):
        from TADA.src.Predict_model_only import load_model
        with redirect_stdout(io.StringIO()) as f:
            model = load_model()
        self.model = model

client.register_plugin(TadaPlugin(), name="tada")

def dask_create_features(tile):
    calculate_properties(tile)

# ---- Helper to tile, calculate, and return features ---- #
def wrapper_create_features(sequence):
    SEQUENCE_WINDOW = 5
    STEPS = 1
    LENGTH = 40
    PROPERTIES = 42

    all_features = []
    for i in range(0, len(sequence) - LENGTH + 1, STEPS):
        tile = sequence[i:i+LENGTH]
        features = calculate_properties(tile)
        all_features.append(features)
    return np.array(all_features)  # shape (36, 42)

def wrapper_pred(seq, features_scaled):
    """
    Runs TADA on the given sequence

    seq: Amino acid sequence (str)
    ss_dict: dictionary of secondary structures

.
    @returns tuple with aa sequence and predictions
    """
    from TADA.src.Predict_model_only import tada_predict

    worker = get_worker()
    model = worker.plugins['tada'].model

    # Trying to prevent TADA output from printing
    with redirect_stdout(io.StringIO()) as f:
        results = tada_predict(loaded_features=features_scaled, loaded_model=model)
    tada_centers = np.arange(len(results)) + 40/2
    
    return (seq, (tada_centers, results))

def main(fasta_name, output_dir):

    # These are used to run in parallel on savio --> Automatically detects available cpus
    # After running this, can use dask as normal

    print("Loading sequences")
    aa_lst = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
    
    fasta_name_tmp = fasta_name.split('/')[-1].removesuffix('.fasta').removesuffix('.fa')
    recs = list(SeqIO.parse(fasta_name, 'fasta'))
    
    seqs_to_keep = []
    for r in recs:
        seq_str = str(r.seq)

        # Filtering sequences with too few AAs and non-standard AAs
        if (len(seq_str) < 40):
            print("Removing", seq_str, "because it is too short")
            continue
        elif (len(set(seq_str) - set(aa_lst)) > 0):
            print("Removing",seq_str, "because it contains non-standard amino acids")
            continue # No need for 
        else:
            seqs_to_keep.append(r)

    recs = seqs_to_keep

    # Defines the sequence window size and steps (stride length). Change values if needed.
    SEQUENCE_WINDOW = 5
    STEPS = 1
    LENGTH = 40

    seqs_40 = []
    seqs_idx = 0
    tada_out = dict()

    #Iterate through all sequences
    print("Tiling sequences")
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
    print("Calculating sequence features for", len(seqs_40), "sequences (this may take a while)")
    start = time.time()
    i  = 0
    features = np.zeros((len(seqs_40), 36, 42))
    features = create_features_faster(seqs_40, SEQUENCE_WINDOW, STEPS)
    features_scaled = scale_features_predict(features)
    print("Done calculating sequence features")
    


    # # Step 1: Feature calculation
    # print("Creating features...")
    # futures = client.map(wrapper_create_features, sequences)
    # progress(futures)
    # features_result = client.gather(futures)  # list of (36, 42)
    # print("Created features.")

    # # Step 2: Feature scaling
    # print("Scaling features...")
    # features_batch = np.array(features_result)  # shape (N, 36, 42)
    # features_scaled_result = scale_features_predict(features_batch)  # same shape
    # print("Scaled features.")

    # scattered_ss_dict = client.scatter(ss_dict, broadcast=True)
    
    ### ------------- RUN TADA  -------------------####
    print("Running TADA") 
    
    futures = client.map(wrapper_pred, sequences, features_scaled_result) 
    progress(futures) # Adds a progress bar
    result = client.gather(futures)
    # Result is a list of tuples, turns into a dict
    tada_out = dict(result)

    print("Finished running TADA, writing results")
    # Write TADA results
    data = []
    
    for r in recs:
        sequence = str(r.seq)
        name = str(r.id)
        tada_centers, tada_preds = tada_out[sequence]
        tada_preds = "[" + ",".join(str(x) for x in tada_preds) + "]"
        tada_centers = "[" + ",".join(str(x) for x in tada_centers) + "]"
        data.append([name, sequence, tada_centers, tada_preds])
    out_df = pd.DataFrame(data, columns=['name', 'sequence', 'tada_centers','tada_preds'])
    out_df.to_csv(f"{output_dir}/{fasta_name_tmp}_tada_dask_preds.csv", encoding='utf-8')


    # These were printing annoying output
    with redirect_stdout(io.StringIO()) as f:
        client.shutdown()
        cluster.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", help="Input fasta file", type=str)
    parser.add_argument("-o", help="Output directory", type=str)
    args = parser.parse_args()
    main(args.f, args.o)

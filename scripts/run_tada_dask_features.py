# Load necessary modules
import numpy as np
import pandas as pd
import matplotlib
matplotlib.rcParams['figure.dpi'] = 100
matplotlib.rcParams['font.sans-serif'] = 'Arial'
from matplotlib import pyplot as plt
import matplotlib.backends.backend_pdf
import seaborn as sns
sns.set_style("white")
sns.set_context("notebook")
from Bio import SeqIO
import glob, pickle, os, re, subprocess, csv, sys
import argparse
import time

sys.path.append(os.path.abspath('/global/scratch/projects/fc_mvslab/predictors/'))
from TADA.src.Preprocessing_FIXED import scale_features_predict
from TADA.src.Predict_model_only import tada_predict

# Set Numpy to display floats with 3 decimal places
np.set_printoptions(formatter={'float': lambda x: "{0:0.3f}".format(x)})

from localcider.sequenceParameters import SequenceParameters
import alphaPredict as alpha
from tqdm import tqdm
from dask_jobqueue.slurm import SLURMRunner
from dask.distributed import Client, progress
#from dask.distributed import LocalCluster

from contextlib import redirect_stdout
import io

def create_features_faster(sequences, client, SEQUENCE_WINDOW = 5, STEPS = 1, LENGTH = 40, PROPERTIES = 42, pools=1):  
    #### Setting up the parallel workers
    futures = client.map(calculate_properties, sequences)
    progress(futures)
    features_result = client.gather(futures)
        # These were printing annoying output        
    return np.array(features_result)

# def create_features_faster(sequences, client, SEQUENCE_WINDOW=5, STEPS=1, batch_size=500):
#     # sequences: list or array of all sequences to process
#     # Break sequences into batches
#     batches = [sequences[i:i+batch_size] for i in range(0, len(sequences), batch_size)]

#     futures = client.map(calculate_properties_batch, batches)
#     progress(futures)
#     batch_results = client.gather(futures)

#     # Concatenate results from batches into one array
#     features_result = np.concatenate(batch_results, axis=0)
#     return features_result

# def calculate_properties_batch(seq_batch):
#     # seq_batch: list or array of sequences
#     # Returns a numpy array of shape (batch_size, ..., ...)
#     batch_features = [calculate_properties(seq) for seq in seq_batch]
#     return np.array(batch_features)
    
def calculate_properties(sequence, SEQUENCE_WINDOW = 5, STEPS = 1, LENGTH = 40, PROPERTIES = 42):
    #print("Calculating properties")
    SEQUENCE_LENGTH = len(sequence)
    SeqOb = SequenceParameters(sequence)
    
    # Break the sequence into smaller pieces of type "SequenceParameters"
    num_sub_seq = int((SEQUENCE_LENGTH-SEQUENCE_WINDOW)/STEPS+1)
    num_features = PROPERTIES
    data = np.zeros((num_sub_seq, num_features)) 
    amino_acids = 'RKDEQNHSTYCWMAILFVPG'

    data[:,0] = np.full(int((SEQUENCE_LENGTH-SEQUENCE_WINDOW)/STEPS+1), SeqOb.get_kappa()) # a parameter to describe the extent of charged amino acid mixing in a sequence
    data[:,1] = np.full(int((SEQUENCE_LENGTH-SEQUENCE_WINDOW)/STEPS+1), SeqOb.get_Omega())
    
    sub_seq = np.array([sequence[STEPS*j:(STEPS*j+SEQUENCE_WINDOW)] for j in range(int((SEQUENCE_LENGTH-SEQUENCE_WINDOW)/STEPS+1))])
    data[:,21] = np.array(list(map(lambda x: sum(alpha.predict(x))/len(x), sub_seq))) # secondary structure --> Why not do with the full sequence??


    sub_seq = np.array([SequenceParameters(sequence[STEPS*j:(STEPS*j+SEQUENCE_WINDOW)]) for j in range(num_sub_seq)])
    # return hydropathy
    data[:,2] = np.array(list(map(lambda x: x.get_mean_hydropathy(),sub_seq))) # hydropathy
    data[:,3] = np.array(list(map(lambda x: x.get_WW_hydropathy(),sub_seq))) # hydropathy_ww
    data[:,4] = np.array(list(map(lambda x: x.get_NCPR(),sub_seq))) # ncpr
    data[:,5] = np.array(list(map(lambda x: x.get_fraction_disorder_promoting(),sub_seq))) # promoting
    data[:,6] = np.array(list(map(lambda x: x.get_FCR(),sub_seq))) #fcr
    data[:,7] = np.array(list(map(lambda x: x.get_mean_net_charge(),sub_seq))) # charge
    data[:,8] = np.array(list(map(lambda x: x.get_fraction_negative(),sub_seq))) # negative
    data[:,9] = np.array(list(map(lambda x: x.get_fraction_positive(),sub_seq))) # positive

    for i in range(SEQUENCE_LENGTH):
        aa = sequence[i]
        if aa == 'I':
            for j in range(SEQUENCE_WINDOW):
                if (i-j) >= 0 and (i-j) < len(sub_seq):
                    aa_index = amino_acids.index(aa)
                    data[i-j,aa_index+22] += 1
                    data[i-j,10] += 1 # add 1 to aliphatics
                    data[i-j,12] += 1 # add 1 to branching
                    data[i-j, 17] += 1 # add 1 to hydrophobics
        elif aa == 'V':
            for j in range(SEQUENCE_WINDOW):
                if (i-j) >= 0 and (i-j) < len(sub_seq):
                    aa_index = amino_acids.index(aa)
                    data[i-j,aa_index+22] += 1
                    data[i-j,10] += 1 # add 1 to aliphatics
                    data[i-j,12] += 1 # add 1 to branching
                    data[i-j, 17] += 1 # add 1 to hydrophobics
        elif aa == 'L':
            for j in range(SEQUENCE_WINDOW):
                if (i-j) >= 0 and (i-j) < len(sub_seq):
                    aa_index = amino_acids.index(aa)
                    data[i-j,aa_index+22] += 1
                    data[i-j,10] += 1 # add 1 to aliphatics
                    data[i-j, 17] += 1 # add 1 to hydrophobics
        elif aa == 'A':
            for j in range(SEQUENCE_WINDOW):
                if (i-j) >= 0 and (i-j) < len(sub_seq):
                    aa_index = amino_acids.index(aa)
                    data[i-j,aa_index+22] += 1
                    data[i-j,10] += 1 # add 1 to aliphatics
                    data[i-j, 20] += 1 # add 1 to tinys
        elif aa == 'W':
            for j in range(SEQUENCE_WINDOW):
                if (i-j) >= 0 and (i-j) < len(sub_seq):
                    aa_index = amino_acids.index(aa)
                    data[i-j,aa_index+22] += 1
                    data[i-j, 11] += 1 # add 1 to aromatics 
                    data[i-j, 17] += 1 # add 1 to hydrophobics
        elif aa == 'F':
            for j in range(SEQUENCE_WINDOW):
                if (i-j) >= 0 and (i-j) < len(sub_seq):
                    aa_index = amino_acids.index(aa)
                    data[i-j,aa_index+22] += 1
                    data[i-j, 11] += 1 # add 1 to aromatics 
                    data[i-j, 17] += 1 # add 1 to hydrophobics
        elif aa == 'Y':
            for j in range(SEQUENCE_WINDOW):
                if (i-j) >= 0 and (i-j) < len(sub_seq):
                    aa_index = amino_acids.index(aa)
                    data[i-j,aa_index+22] += 1
                    data[i-j, 11] += 1 # add 1 to aromatics 
                    data[i-j, 15] += 1 # add 1 to phosphorylatables
                    data[i-j, 16] += 1 # add 1 to polars
        elif aa == 'T':
            for j in range(SEQUENCE_WINDOW):
                if (i-j) >= 0 and (i-j) < len(sub_seq):
                    aa_index = amino_acids.index(aa)
                    data[i-j,aa_index+22] += 1
                    data[i-j,12] += 1 # add 1 to branching
                    data[i-j, 15] += 1 # add 1 to phosphorylatables
        elif aa == 'K':
            for j in range(SEQUENCE_WINDOW):
                if (i-j) >= 0 and (i-j) < len(sub_seq):
                    aa_index = amino_acids.index(aa)
                    data[i-j,aa_index+22] += 1
                    data[i-j,13] += 1 # add 1 to charged
                    data[i-j, 16] += 1 # add 1 to polars
                    data[i-j, 18] += 1 # add 1 to positives
        elif aa == 'R':
            for j in range(SEQUENCE_WINDOW):
                if (i-j) >= 0 and (i-j) < len(sub_seq):
                    aa_index = amino_acids.index(aa)
                    data[i-j,aa_index+22] += 1
                    data[i-j,13] += 1 # add 1 to charged
                    data[i-j, 16] += 1 # add 1 to polars
                    data[i-j, 18] += 1 # add 1 to positives
        elif aa == 'H':
            for j in range(SEQUENCE_WINDOW):
                if (i-j) >= 0 and (i-j) < len(sub_seq):
                    aa_index = amino_acids.index(aa)
                    data[i-j,aa_index+22] += 1
                    data[i-j,13] += 1 # add 1 to charged
                    data[i-j, 18] += 1 # add 1 to positives
        elif aa == 'D':
            for j in range(SEQUENCE_WINDOW):
                if (i-j) >= 0 and (i-j) < len(sub_seq):
                    aa_index = amino_acids.index(aa)
                    data[i-j,aa_index+22] += 1
                    data[i-j,13] += 1 # add 1 to charged
                    data[i-j,14] += 1 # add 1 to negative
                    data[i-j, 16] += 1 # add 1 to polars
        elif aa == 'E':
            for j in range(SEQUENCE_WINDOW):
                if (i-j) >= 0 and (i-j) < len(sub_seq):
                    aa_index = amino_acids.index(aa)
                    data[i-j,aa_index+22] += 1
                    data[i-j,13] += 1 # add 1 to charged
                    data[i-j,14] += 1 # add 1 to negative
                    data[i-j, 16] += 1 # add 1 to polars
        elif aa == 'S':
            for j in range(SEQUENCE_WINDOW):
                if (i-j) >= 0 and (i-j) < len(sub_seq):
                    aa_index = amino_acids.index(aa)
                    data[i-j,aa_index+22] += 1
                    data[i-j, 15] += 1 # add 1 to phosphorylatables
                    data[i-j, 20] += 1 # add 1 to tinys
        elif aa == 'Q':
            for j in range(SEQUENCE_WINDOW):
                if (i-j) >= 0 and (i-j) < len(sub_seq):
                    aa_index = amino_acids.index(aa)
                    data[i-j,aa_index+22] += 1
                    data[i-j, 16] += 1 # add 1 to polars
        elif aa == 'N':
            for j in range(SEQUENCE_WINDOW):
                if (i-j) >= 0 and (i-j) < len(sub_seq):
                    aa_index = amino_acids.index(aa)
                    data[i-j,aa_index+22] += 1
                    data[i-j, 16] += 1 # add 1 to polars
        elif aa == 'C':
            for j in range(SEQUENCE_WINDOW):
                if (i-j) >= 0 and (i-j) < len(sub_seq):
                    aa_index = amino_acids.index(aa)
                    data[i-j,aa_index+22] += 1
                    data[i-j, 17] += 1 # add 1 to hydrophobics
                    data[i-j, 19] += 1 # add 1 to sulfur containing 
        elif aa == 'M':
            for j in range(SEQUENCE_WINDOW):
                if (i-j) >= 0 and (i-j) < len(sub_seq):
                    aa_index = amino_acids.index(aa)
                    data[i-j,aa_index+22] += 1
                    data[i-j, 17] += 1 # add 1 to hydrophobics
                    data[i-j, 19] += 1 # add 1 to sulfur containing 
        elif aa == 'G':
            for j in range(SEQUENCE_WINDOW):
                if (i-j) >= 0 and (i-j) < len(sub_seq):
                    aa_index = amino_acids.index(aa)
                    data[i-j,aa_index+22] += 1
                    data[i-j, 20] += 1 # add 1 to tinys
        elif aa == 'P':
            for j in range(SEQUENCE_WINDOW):
                if (i-j) >= 0 and (i-j) < len(sub_seq):
                    aa_index = amino_acids.index(aa)
                    data[i-j,aa_index+22] += 1
                    data[i-j, 20] += 1 # add 1 to tinys
        
    x = data.copy()
    shape = np.shape(x)
    
    return x


def run_tada(fasta_name, output_dir, client):
    # Provide the path to a fasta file containing your sequences

    aa_lst = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']

    # Set up SeqRecord objects and check that sequences are valid
    fasta_name_tmp = fasta_name.split('/')[-1].removesuffix('.fasta').removesuffix('.fa')
    os.makedirs(output_dir, exist_ok=True)

    recs = list(SeqIO.parse(fasta_name, 'fasta'))
    fasta_name = fasta_name_tmp

    seqs_to_keep = []
    for seq in recs:
        seq_str = (seq.seq)
        if (len(seq_str) < 39):
            print("Removing", seq_str, "because it is too short")
        elif (len(set(seq_str) - set(aa_lst)) > 0):
            print("Removing",seq_str, "because it contains non-standard amino acids")
        else:
            seqs_to_keep.append(seq)

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
        
    features = create_features_faster(seqs_40, client, SEQUENCE_WINDOW, STEPS)
            
    features_scaled = scale_features_predict(features)
    print("Done calculating sequence features")
    
    # Make TADA classification predictions
    print("Predicting activities")
    tada_preds = tada_predict(loaded_features=features_scaled)
    print("Done predicting activities")

    print("Saving predictions")
    # Use indices of tiles to retrieve predictions for each sequence
    for seq, val in tada_out.items():
        preds = tada_preds[val]
        centers = np.arange(len(preds)) + 40/2
        tada_out[seq] = (centers, preds)

    # Save predicted values to a csv file
    with open(output_dir + '/' + fasta_name_tmp + 'TADA_preds.csv', 'w+') as f:
        w = csv.writer(f)
        w.writerow(['sequence', 'tada_centers', 'tada_preds'])
        for r in recs:
            sequence = str(r.seq)
            tada_centers, tada_preds = tada_out[sequence]
            #tada_preds = re.sub(r'\s+', ',', str(tada_preds))
            #tada_centers = re.sub(r'\s+', ',', str(tada_centers))
            tada_preds = "[" + ",".join(str(x) for x in tada_preds) + "]"
            tada_centers = "[" + ",".join(str(x) for x in tada_centers) + "]"
            w.writerow([sequence, tada_centers, tada_preds])

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", help="Input fasta file", type=str)
    parser.add_argument("-o", help="Output directory", type=str)
    args = parser.parse_args()

    cluster = SLURMRunner()
    client = Client(cluster)
    client.wait_for_workers(5)   
    run_tada(args.f, args.o, client)

    with redirect_stdout(io.StringIO()) as f:
        client.shutdown()
        cluster.close()

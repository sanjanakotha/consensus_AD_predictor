# Load necessary modules
import warnings
warnings.filterwarnings('ignore')

import os
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'

import numpy as np
import pandas as pd
import matplotlib
from matplotlib import pyplot as plt
import matplotlib.backends.backend_pdf
from Bio import SeqIO
import glob, pickle, os, re, subprocess, csv, sys, io
from contextlib import redirect_stdout

from uuid import uuid4
import argparse
from tqdm import tqdm 
import time

sys.path.append(os.path.abspath('/global/scratch/projects/fc_mvslab/predictors'))
from adpred import ADpred as adp
from iupred3.iupred import get_iupred

import tensorflow as tf
from tensorflow.keras.models import Model
from tensorflow.keras.layers import Input, Conv2D, Flatten, Dense, Dropout, Activation
from tensorflow.keras import regularizers
from tensorflow.keras.activations import softplus
from tensorflow.keras import backend as K
from tensorflow.keras import metrics
from contextlib import redirect_stdout, redirect_stderr

# Limit tensorflow thread usage
import ipyparallel as ipp

# Amino acids list
aa_lst = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']

def return_adpred_ss(r, output_dir):
    sequence = r.seq
    sub_fasta_name = re.sub(r'\s+', '_', r.id)
    computed = f"{output_dir}/ss/{sub_fasta_name}/{sub_fasta_name}.ss2"
    df = pd.read_table(computed, skiprows=[0], sep='\s+', names=['position', 'residue', 'ss', 'c', 'h', 'e'])
    return ''.join(df.ss.values).replace('C','-')

def make_pred(r, output_dir):
    with redirect_stdout(io.StringIO()) as f:
        sequence = r.seq
        adp_prot = adp.protein('x', sequence)
        adp_prot.second_struct = return_adpred_ss(r, output_dir)
        adp_prot.predict()
        adpred_preds = adp_prot.predictions
    return (sequence, adpred_preds)

def main():
    # --- Argument parsing ---
    parser = argparse.ArgumentParser(description="Run ADpred pipeline on a FASTA file.")
    parser.add_argument('-f', '--fasta_path', type=str, required=True, help='Path to input FASTA file')
    parser.add_argument('-o', '--output_dir', type=str, required=True, help='Directory to write output files')
    args = parser.parse_args()

    fasta_path = args.fasta_path
    output_dir = args.output_dir

    print("Loading sequences...")
    recs = list(SeqIO.parse(fasta_path, 'fasta'))

    # Filter sequences
    seqs_to_keep = []
    for r in recs:
        seq_str = str(r.seq)
        if len(seq_str) < 30:
            print("Removing", seq_str, "because it is too short")
            continue
        elif len(set(seq_str) - set(aa_lst)) > 0:
            print("Removing", seq_str, "because it contains non-standard amino acids")
            continue
        else:
            seqs_to_keep.append(r)

    recs = seqs_to_keep
    print(f"Kept {len(recs)} sequences.")

    # Get output filename stem
    fasta_name_tmp = os.path.basename(fasta_path).removesuffix('.fasta').removesuffix('.fa')

    # Run predictions
    print("Running ADpred predictions...")
    results = []
    for r in tqdm(recs):
        try:
            results.append(make_pred(r, output_dir))
        except Exception as e:
            print(f"Failed to process sequence {r.id}: {e}")
            results.append((str(r.seq), []))

    results_dict = dict(results)

    # Write ADpred results
    print("Saving results...")
    data = []
    for r in recs:
        sequence = str(r.seq)
        name = str(r.id)
        adpred_preds = results_dict[sequence]
        adpred_preds = "[" + ",".join(str(x) for x in adpred_preds) + "]"
        data.append([name, sequence, adpred_preds])

    out_df = pd.DataFrame(data, columns=['name', 'sequence', 'adpred_preds'])
    out_csv_path = os.path.join(output_dir, f"{fasta_name_tmp}_ADpred_preds.csv")
    out_df.to_csv(out_csv_path, encoding='utf-8', index=False)
    print("Done.")

if __name__ == "__main__":
    main()

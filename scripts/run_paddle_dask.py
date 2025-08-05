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
from PADDLE import paddle
from iupred3.iupred import get_iupred

import dask
from dask_jobqueue.slurm import SLURMRunner
from dask.distributed import Client, progress
from dask.diagnostics import ProgressBar

# Limit tensorflow thread usage
import os
# os.environ["OMP_NUM_THREADS"] = "1"
# os.environ["TF_NUM_INTRAOP_THREADS"] = "1"
# os.environ["TF_NUM_INTEROP_THREADS"] = "1"

cluster = SLURMRunner()
client = Client(cluster)

client.wait_for_workers(5)

# The goal of this is to use a individual paddle model for each worker
from dask.distributed import get_worker, WorkerPlugin

class PaddlePlugin(WorkerPlugin):
    def setup(self, worker):
        from PADDLE import paddle
        # worker.extensions["paddle_model"] = paddle.PADDLE()
        # The other option
        self.paddle = paddle.PADDLE()

client.register_plugin(PaddlePlugin(), name="paddle")

# path to a local installation of psipred
local_psipred = '/global/scratch/projects/fc_mvslab/predictors/psipred/runpsipred' 

def get_psipred(seq, output_dir, fasta_name, computed=None):
    '''
    Runs PSIPRED on the provided sequence and returns secondary structure
    data for PADDLE and ADpred, viz
    (paddle_seq, paddle_helix, paddle_coil), 'adpred_ss_as_a_string'
    
    seq: provided sequence
    output_dir: path to output directory containing fasta file
    fasta_name: name of fasta file, minus the .fasta/.fa file ending
    computed: path to already computed .ss2 file to skip re-running PSIPRED
    ''' 
    if computed is not None:
        return paddle.read_ss2(computed)

    p = subprocess.run('cd "{}" && {} "{}".fa*'.format(output_dir, local_psipred, fasta_name), 
                       capture_output=True, shell=True, encoding='utf-8')
    rootname = re.search('Final output files: (.*).ss2 (.*).horiz', p.stdout).group(1).strip()
    
    paddle_ss = paddle.read_ss2(output_dir + '/' + rootname + '.ss2')

    return paddle_ss
    

def read_iupred_batch(recs, output_dir, fasta_name, computed_long=None, computed_short=None):
    '''
    Runs IUPred on the provided sequences and returns two dictionaries:
    {sequences: lists of long-style predicted disorder values} and
    {sequences: lists of short-style predicted disorder values}
    
    recs: list of Bio.SeqRecord objects containing the provided sequences
    output_dir: path to output directory
    fasta_name: path to fasta file
    computed_long: path to already computed long-style .dis file to skip re-running IUPred
    computed_short: path to already computed short-style .dis file to skip re-running IUPred
    ''' 
    if computed_long is not None and computed_short is not None:
        return read_iupred_dict(computed_long), read_iupred_dict(computed_short)
    
    with open(output_dir + '/ss/' + fasta_name + '_long.dis', 'w') as dis:
        for r in recs:
            dis.write('> ' + r.id + '\n')
            dis.write(get_iupred(str(r.seq), iupred_type='long'))
            dis.write('\n\n')
    with open(output_dir + '/ss/' + fasta_name + '_short.dis', 'w') as dis:
        for r in recs:
            dis.write('> ' + r.id + '\n')
            dis.write(get_iupred(str(r.seq), iupred_type='short'))
            dis.write('\n\n')
    
    long_out = read_iupred_dict(output_dir + '/ss/' + fasta_name + '_long.dis')
    short_out = read_iupred_dict(output_dir + '/ss/' + fasta_name + '_short.dis')
    
    return long_out, short_out


def read_iupred_dict(file):
    """
    Read disorder predictions produced by IUPRED2 from the .dis file
    Adapted from PADDLE's native read_iupred fn due to formatting
    Inputs:
        - file: Filename of the .dis file output by IUPRED2
    Returns:
        - iupred_dict: Dictionary of prot keys and dis values
        where
        - prot: Protein sequence (string)
        - dis:  List of predicted disorder (values between 0 and 1)
    """

    f = open(file, 'r')
    iupred_dict = dict()
    dis = []
    prot = ''

    for line in f:
        line = line.strip()
        if line.startswith('>'):
            iupred_dict[prot] = dis
            dis = []
            prot = ''
            continue
        if '#' in line or line == "":
            continue
        tokens = line.split()
        assert len(tokens) == 3, tokens

        prot += tokens[1]
        dis.append(float(tokens[2]))

    iupred_dict[prot] = dis
    iupred_dict.pop("")

    f.close()
    return iupred_dict

def wrapper_ss(r, output_dir):
    """
    Calculates the secondary structure for a given seq record

    r: SeqIO seq record

    @returns tuple with aa sequuence and list of ss predictions
    """
    sequence = r.seq
    sub_fasta_name = re.sub(r'\s+', '_', r.id)
    # print(output_dir)
    
    if os.path.exists("{}/ss/{}/{}.ss2".format(output_dir, sub_fasta_name, sub_fasta_name)):
        # print("SS already exists, calling get psipred")
        paddle_ss = get_psipred(sequence, output_dir=output_dir + '/' + sub_fasta_name, fasta_name=sub_fasta_name,
                                          computed="{}/ss/{}/{}.ss2".format(output_dir, sub_fasta_name, sub_fasta_name) )
    else:
        # print(r.id, "doesn't have ss prediction")
        paddle_ss = get_psipred(sequence, output_dir=output_dir + '/ss/' + sub_fasta_name, fasta_name=sub_fasta_name)
        
    return paddle_ss
    
def wrapper_pred(seq, ss, dis_short, dis_long):
    # Load pre-computed secondary structure predicted by PSIPRED V4.
    # pad = paddle.PADDLE()
    worker = get_worker()
    # pad = worker.extensions["paddle_model"]
    pad = worker.plugins['paddle'].paddle
    
    prot, helix, coil = ss[0], ss[1], ss[2]

    # Run predictions on all 53 amino acid long tiles across the protein.
    # This function requires matching protein sequence and secondary structure scores.
    # Returns a Numpy array of size (protein_length-52) which gives the
    # predicted activation Z-score for the 53aa tiles starting at positions
    # 1, 2, 3, ..., protein_length-52.
    # High-strength ADs can be called by finding >=5 consecutive positions with Z-score > 6.
    # Medium-strength ADs can be called by finding >=5 consecutive positions with Z-score > 4.
    paddle_preds = pad.predict_protein(prot, helix, coil, dis_short, dis_long)
    paddle_centers = np.arange(len(paddle_preds)) + (53+1)/2

    # paddle_preds = worker.paddle.predict_protein(prot, helix, coil, dis_short, dis_long)
    # paddle_centers = np.arange(len(paddle_preds)) + (53+1)/2
    
    return (seq, (paddle_centers, paddle_preds))
    
def main(fasta_name, output_dir):

    # These are used to run in parallel on savio --> Automatically detects available cpus
    # After running this, can use dask as normal

    print("Loading sequences")
    aa_lst = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
    
    fasta_name_tmp = fasta_name.split('/')[-1].removesuffix('.fasta').removesuffix('.fa')
    recs = list(SeqIO.parse(fasta_name, 'fasta'))
    fasta_name = fasta_name_tmp
    
    seqs_to_keep = []
    for r in recs:
        seq_str = str(r.seq)

        # Filtering sequences with too few AAs and non-standard AAs
        if (len(seq_str) < 53):
            print("Removing", seq_str, "because it is too short")
            continue
        elif (len(set(seq_str) - set(aa_lst)) > 0):
            print("Removing",seq_str, "because it contains non-standard amino acids")
            continue # No need for 
        else:
            seqs_to_keep.append(r)

        # Writes each sequence to its own fasta file
        sub_fasta_name = re.sub(r'\s+', '_', r.id)

        # make folder to store secondary structures
        if not os.path.exists(output_dir + '/ss/'):
            subprocess.run(['mkdir', output_dir + '/ss/'])
            
        if not os.path.exists(output_dir + '/ss/' + sub_fasta_name):
            subprocess.run(['mkdir', output_dir + '/ss/' + sub_fasta_name])
        with open(output_dir + '/ss/' + sub_fasta_name + '/' + sub_fasta_name + '.fasta', 'w') as fa:
            fa.write(seq_str)
    
    recs = seqs_to_keep
    
    print("Predicting SS") 
    futures = client.map(wrapper_ss, recs, output_dir=output_dir)    # Send many tasks
    progress(futures)
    ss_result = client.gather(futures)
    
    print("Completed SS predictions") 
    
    ### ------------- RUN ADPRED  -------------------####
    #Need a wrapper function because map() only operates on one argument
    print("Running PADDLE") 
    if os.path.exists("{}/ss/{}_long.dis".format(output_dir, fasta_name)):
        dis_long_dict, dis_short_dict = read_iupred_batch(recs,output_dir,fasta_name,
                                                          computed_long="{}/ss/{}_long.dis".format(output_dir, fasta_name),
                                                          computed_short="{}/ss/{}_short.dis".format(output_dir, fasta_name))
    else:
        dis_long_dict, dis_short_dict = read_iupred_batch(recs,output_dir,fasta_name)

    sequences = [str(r.seq) for r in recs]

    # Load pre-computed disorder predicted by IUPRED2, in both
    # the short and long modes.
    dis_shorts = [dis_short_dict[str(r.seq)] for r in recs] 
    dis_longs = [dis_long_dict[str(r.seq)] for r in recs] 
 
    
    # This is where the parallelization happens
    futures = client.map(wrapper_pred, sequences, ss_result, dis_shorts, dis_longs) 
    progress(futures) # Adds a progress bar

    # Waits until all the parallel tasks are finished
    result = client.gather(futures)

    # Result is a list of tuples, turns into a dict
    paddle_out = dict(result)

    print("Finished running PADDLE, writing results")
    # Write PADDLE results
    data = []
    for r in recs:
        sequence = str(r.seq)
        name = str(r.id)
        paddle_centers, paddle_preds = paddle_out[sequence]
        paddle_preds = "[" + ",".join(str(x) for x in paddle_preds) + "]"
        paddle_centers = "[" + ",".join(str(x) for x in paddle_centers) + "]"
        data.append([name, sequence, paddle_centers, paddle_preds])
    out_df = pd.DataFrame(data, columns=['name', 'sequence', 'paddle_centers','paddle_preds'])
    out_df.to_csv(f"{output_dir}/{fasta_name_tmp}_Paddle_dask_preds.csv", encoding='utf-8')

    # These were printing annoying output (still are) 
    with redirect_stdout(io.StringIO()) as f:
        client.shutdown()
        cluster.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", help="Input fasta file", type=str)
    parser.add_argument("-o", help="Output directory", type=str)
    args = parser.parse_args()
    main(args.f, args.o)


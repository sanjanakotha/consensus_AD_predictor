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
# from adpred import ADpred as adp
# from PADDLE import paddle
from iupred3.iupred import get_iupred

import dask
from dask_jobqueue.slurm import SLURMRunner
from dask.distributed import Client, progress, get_worker, WorkerPlugin
from dask.diagnostics import ProgressBar

import tensorflow as tf
from tensorflow.keras.models import Model
from tensorflow.keras.layers import Input, Conv2D, Flatten, Dense, Dropout, Activation
from tensorflow.keras import regularizers
from tensorflow.keras.activations import softplus
from tensorflow.keras import backend as K
from tensorflow.keras import metrics

# Limit tensorflow thread usage
import os
# os.environ["OMP_NUM_THREADS"] = "1"
# os.environ["TF_NUM_INTRAOP_THREADS"] = "1"
# os.environ["TF_NUM_INTEROP_THREADS"] = "1"

def get_model():
    """
    Function to reload the ADPred model. 
    Must be done once in each worker otherwise
    there are issues. 
    """
    K.clear_session()                                                           
    inputs = Input(shape=(30,23,1))                                             
    x = Conv2D(29, (6,23), activation=softplus)(inputs)                         
    x = Flatten()(x)                                                            
    x = Dense(100, activation=softplus, kernel_regularizer=regularizers.l2(0.001), name='dense_layer_1')(x)
    x = Dropout(0.5, name='dropout_1')(x)                                                         
    x = Dense(30, activation=softplus, kernel_regularizer=regularizers.l2(0.001), name='dense_layer_2')(x)
    x = Dropout(0.5, name='dropout_2')(x)                                                         
    x = Dense(1, name='dense_layer_3')(x)                                                             
    output = (Activation('sigmoid'))(x)           
    
    model = Model(inputs=inputs, outputs=output)                               
    model.compile(optimizer='adam', loss='binary_crossentropy', metrics=[auc]) 
    model.load_weights('/global/scratch/projects/fc_mvslab/predictors/adpred/run_w_dask/ADpred.h5')     
    return model

# Functions to reimplement ADpred model
def auc(y_true, y_pred):
    auc = metrics.auc(y_true, y_pred)[1]  # using defaults parameters --> num_thresholds=200
    K.get_session().run(tf.local_variables_initializer())
    return auc

def make_ohe(seq, struct):
    '''
        function returns the data in ohe shape. The columns correspond to the lexicon.
        INPUT: sequence. Sequence of amino acids or secondary structure (ss) elements.
               lexicon. Ordered list of all 20 amino acids or ss elements.
        OUTPUT: ohe_data (shape = (1, len(lexicon))
        e.g. of lexicon for ss: ["E","H","-"] --> beta, alpha, coil

        NOTE: This function can be vectorized since it will constitute a ufunc 
              and the result matrix should have a shape = (len(sequences), len(lexicon))
    '''
    # initialize tensors
    ohe_seq = np.zeros(shape=(len(seq), 20))
    ohe_ss = np.zeros(shape=(len(struct),3))
    
    aa = ['R','H','K','D','E','S','T','N','Q','A','V','L','I','M','F' ,'Y', 'W', 'C','G','P']
    ss = ['E','H','-'] # list of secondary structure elements

    # encode sequence and secondary structure
    for n in range(len(seq)):
        ohe_seq[n,aa.index(seq[n])] = 1
        ohe_ss[n, ss.index(struct[n])] = 1

    # join botho tensor 
    ohe = np.vstack([ohe_seq.T, ohe_ss.T]).T #.reshape(1,len(seq),23,1)

    return ohe

#### Setting up the parallel workers
cluster = SLURMRunner()
client = Client(cluster)

client.wait_for_workers(5)

# This should prevent the model from being loaded every tme
class ADpredPlugin(WorkerPlugin):
    def setup(self, worker):
        self.model = get_model()

client.register_plugin(ADpredPlugin(), name="adpred")


# path to a local installation of psipred
local_psipred = '/global/scratch/projects/fc_mvslab/predictors/psipred/runpsipred'

def predict(model, seq, struct=None):
    '''
    Copied from ADPred package
    ------------
    Assigns to each aminoacid of the sequence the probability of being in a 
    AD region.

    parameters
    ----------
      - sequence: sequence of amino acids (aa). 
        *** REMOVED PART ALLOW UNIPROT IDS ***
      - second_struct: sequence of secondary elements of each aa in sequence.

    returns
    -------
      - numpy array with probabilities of AD for each aa in the input sequence.
    '''
    
    # extend adapters for the extremes #####, unless is a 30mer (e.g. in saturated mutagenesis)
    seq = ''.join(['G']*15) + seq + ''.join(['G']*15)
    struct = ''.join(['-']*15) + struct + ''.join(['-']*15)
                                                                                
    # encode for keras and initialize results                                   
    ohe = make_ohe(seq,struct)                                                  
    # results = np.zeros(len(seq)-30)                                             
                                                                                
    # roll window of predictions --> Slow!!                                               
    # for n in range(results.shape[0]):                                           
    #     results[n] = model.predict(ohe[n:n+30].reshape(1,30,23,1))[0][0]     

    # BATCH METHOD
    X = np.stack([ohe[n:n+30].reshape(30, 23, 1) for n in range(len(seq) - 30)])
    results = model.predict(X, batch_size=128, verbose=0)

    return results[:, 0]
    
def adpred_read_ss2(ss2):
    '''
    Converts .ss2 data to horizontal secondary structure information (e.g., --HHHHH--HHHHH-HH)
    
    ss2: name of .ss2 file
    '''
    df = pd.read_table(ss2, skiprows=[0], sep='\s+', names=['position', 'residue', 'ss', 'c', 'h', 'e'])
    return ''.join(df.ss.values).replace('C','-')


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
        # print("Already computed ss")
        return adpred_read_ss2(computed)
    
    p = subprocess.run('cd "{}" && {} "{}".fa*'.format(output_dir, local_psipred, fasta_name), 
                       capture_output=True, shell=True, encoding='utf-8')
    rootname = re.search('Final output files: (.*).ss2 (.*).horiz', p.stdout).group(1).strip()
    
    # paddle_ss = paddle.read_ss2(output_dir + '/' + rootname + '.ss2')
    f = open(output_dir + '/' + rootname + '.horiz')
    adpred_ss = ''.join([i.group(1) for i in re.finditer('Pred: (.*)\n', ''.join(f))]).replace("C","-")

    return adpred_ss
    
def wrapper_ss(r, output_dir):
    sequence = r.seq
    sub_fasta_name = re.sub(r'\s+', '_', r.id)
    computed="{}/ss/{}/{}.ss2".format(output_dir, sub_fasta_name, sub_fasta_name)
    df = pd.read_table(computed, 
                       skiprows=[0], sep='\s+', names=['position', 'residue', 'ss', 'c', 'h', 'e'])
    return ''.join(df.ss.values).replace('C','-')
    
# def wrapper_ss(r, output_dir):
#     """
#     Calculates the secondary structure for a given seq record

#     r: SeqIO seq record

#     @returns tuple with aa sequuence and list of ss predictions
#     """
#     sequence = r.seq
#     sub_fasta_name = re.sub(r'\s+', '_', r.id)
#     # print(output_dir)
    
#     if os.path.exists("{}/ss/{}/{}.ss2".format(output_dir, sub_fasta_name, sub_fasta_name)):
#         adpred_ss = get_psipred(sequence, output_dir=output_dir + '/' + sub_fasta_name, fasta_name=sub_fasta_name,
#                                           computed="{}/ss/{}/{}.ss2".format(output_dir, sub_fasta_name, sub_fasta_name) )
#     else:
#         adpred_ss = get_psipred(sequence, output_dir=output_dir + '/ss/' + sub_fasta_name, fasta_name=sub_fasta_name)
#     # assert paddle_ss[0] == sequence
#     return adpred_ss

def wrapper_pred(seq, ss):
    """
    Runs ADPred on the given sequence

    seq: Amino acid sequence (str)
    ss_dict: dictionary of secondary structures

    @returns tuple with aa sequence and predictions
    """
    worker = get_worker()
    model = worker.plugins['adpred'].model

    print(seq)
    # Trying to prevent ADPred output from printing
    with redirect_stdout(io.StringIO()) as f:
        results = predict(model, seq, ss) 
    
    return (seq, results)
    
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
        if (len(seq_str) < 30):
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

    # ss_dict = dict(result)
    # print(ss_dict[list(ss_dict.keys())[0]])
    
    print("Completed SS predictions") 
    # scattered_ss_dict = client.scatter(ss_dict, broadcast=True)
    
    ### ------------- RUN ADPRED  -------------------####
    #Need a wrapper function because map() only operates on one argument
    print("Running ADPred") 

    sequences = [str(r.seq) for r in recs]
    
    # This is where the parallelization happens
    futures = client.map(wrapper_pred, sequences, ss_result) 
    progress(futures) # Adds a progress bar

    # Waits until all the parallel tasks are finished
    result = client.gather(futures)

    # Result is a list of tuples, turns into a dict
    adpred_out = dict(result)

    print("Finished running ADPred, writing results")
    # Write ADPred results
    data = []
    for r in recs:
        sequence = str(r.seq)
        name = str(r.id)
        adpred_preds = adpred_out[sequence] 
        adpred_preds = "[" + ",".join(str(x) for x in adpred_preds) + "]"
        data.append([name, sequence, adpred_preds])
    out_df = pd.DataFrame(data, columns=['name', 'sequence', 'adpred_preds'])
    out_df.to_csv(f"{output_dir}/{fasta_name_tmp}_ADpred_preds_dask.csv", encoding='utf-8')

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


import tensorflow as tf 
tf.config.list_physical_devices('GPU')
# Load necessary modules
import numpy as np
import pandas as pd
from Bio import SeqIO
import glob, pickle, os, re, subprocess, csv, sys
from uuid import uuid4
import ipyparallel as ipp

os.chdir("/global/scratch/projects/fc_mvslab/predictors/all_predictors/")

sys.path.append(os.path.abspath('..'))
from PADDLE import paddle
from iupred3.iupred import get_iupred

import torch
from torch.utils.data import TensorDataset, DataLoader
from sklearn import preprocessing
from pickle import dump
import pytorch_lightning as pl
from actpred.models import ActCNNSystem
from actpred.utils import get_threshold, get_stratified_split

from scipy.stats import spearmanr
from pytorch_lightning.callbacks import ModelCheckpoint, EarlyStopping
from pytorch_lightning.loggers import CSVLogger, WandbLogger
import argparse
import time

overall_start = time.time()

def main(fasta_name, output_dir):
    # Set Numpy to display floats with 3 decimal places
    np.set_printoptions(formatter={'float': lambda x: "{0:0.3f}".format(x)})
    
    ipp.__version__   # check that at least version 7
    mycluster = ipp.Cluster(n = int(os.getenv('SLURM_CPUS_ON_NODE')))
    mycluster.start_cluster_sync()
    
    c = mycluster.connect_client_sync()
    c.wait_for_engines(n = int(os.getenv('SLURM_CPUS_ON_NODE')))
    c.ids
    
    dview = c[:]
    dview.block = True

    # Load necessary modules on workers
    dview.execute('import sys, os')
    dview.execute('os.chdir("/global/scratch/projects/fc_mvslab/predictors/all_predictors/")')
    dview.execute('sys.path.append(os.path.abspath(".."))')
    dview.execute('import pandas as pd')
    dview.execute('import numpy as np')
    dview.execute('import re, subprocess')
    dview.execute('import Bio')
    dview.execute('from PADDLE import paddle')
    dview.execute('from adpred import ADpred as adp')
    dview.execute('import torch')
    dview.execute('from TADA.src.Preprocessing_FIXED import scale_features_predict')
    dview.execute('from TADA.src.Preprocessing import create_features')
    dview.execute('from TADA.src.Predict_model_only import tada_predict')
    dview.execute('from TADA.src.Predict_model_only import load_model')
    
    # Provide the path to a fasta file containing your sequences
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

    recs = seqs_to_keep
    
    # path to a local installation of psipred
    local_psipred = '/global/scratch/projects/fc_mvslab/predictors/psipred/runpsipred'
    
    def get_psipred(seq, output_dir=output_dir, fasta_name=fasta_name, computed=None):
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
        print(p.stdout)
        rootname = re.search('Final output files: (.*).ss2 (.*).horiz', p.stdout).group(1).strip()
        
        paddle_ss = paddle.read_ss2(output_dir + '/' + rootname + '.ss2')
        f = open(output_dir + '/' + rootname + '.horiz')
        adpred_ss = ''.join([i.group(1) for i in re.finditer('Pred: (.*)\n', ''.join(f))]).replace("C","-")
        
        
        return paddle_ss
    
    dview = c[:]
    dview.block = True
    
    # Push variables and functions to workers
    var = dict( output_dir=output_dir, local_psipred=local_psipred, fasta_name=fasta_name,get_psipred=get_psipred )
    dview.push(var)

    
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
            paddle_ss = get_psipred(sequence, output_dir=output_dir + '/' + sub_fasta_name, fasta_name=sub_fasta_name, computed="{}/ss/{}/{}.ss2".format(output_dir, sub_fasta_name, sub_fasta_name) )
        else:
            # print(r.id, "doesn't have ss prediction")
            paddle_ss = get_psipred(sequence, output_dir=output_dir + '/ss/' + sub_fasta_name, fasta_name=sub_fasta_name)
        
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

    print("GETTING IUPRED")
    dis_long_dict, dis_short_dict = read_iupred_batch(recs,
                                                      output_dir,
                                                      fasta_name,
                                                      computed_long="{}/{}_long.dis".format(output_dir, fasta_name),
                                                      computed_short="{}/{}_short.dis".format(output_dir,fasta_name))
    ss_dict = {str(r.seq): wrapper_ss(r, output_dir) for r in recs}

    print(f"DONE GETTING IUPRED, {round(time.time()-overall_start, 2)}s elapsed\n")
    
    ########
    # Run predictions using PADDLE. This network requires predicted secondary
    # structure (from PSIPRED) and predicted disorder (from IUPRED2, in both the
    # short and long modes) as input in addition to the protein sequence. This
    # model is the most accurate and should be used for predicted ADs in wild-type
    # proteins.
    ########
    
    # Load PADDLE models. There are 10 in total, and their individual predictions
    # are averaged to obtain the final result.
    dview.execute('pad=paddle.PADDLE()')
    
    # Push PSIPRED and IUPred output to workers
    dview.push(dict( ss_dict=ss_dict, dis_long_dict=dis_long_dict, dis_short_dict=dis_short_dict ))
    
    # Need a wrapper function because map() only operates on one argument
    
    def wrapper_pred(r):
        sub_fasta_name = re.sub(r'\s+', '_', r.id)
        seq = str(r.seq)
        
        paddle_ss = ss_dict[seq]
        prot, helix, coil = paddle_ss[0], paddle_ss[1], paddle_ss[2]
    
        
        # Run predictions on all 53 amino acid long tiles across the protein.
        # This function requires matching protein sequence and secondary structure scores.
        # Returns a Numpy array of size (protein_length-52) which gives the
        # predicted activation Z-score for the 53aa tiles starting at positions
        # 1, 2, 3, ..., protein_length-52.
        # High-strength ADs can be called by finding >=5 consecutive positions with Z-score > 6.
        # Medium-strength ADs can be called by finding >=5 consecutive positions with Z-score > 4.
        
        paddle_preds = pad.predict_protein(prot, helix, coil, dis_short, dis_long)
        paddle_centers = np.arange(len(paddle_preds)) + (53+1)/2
        return (seq, (paddle_centers, paddle_preds))
    
    print("RUNNING PADDLE")
    # Run a parallel map, executing the wrapper function on indices 0,...,n-1
    lview = c.load_balanced_view()
    # Cause execution on main process to wait while tasks sent to workers finish
    lview.block = True 
    out = lview.map(wrapper_pred, recs)   # Run calculation in parallel
    paddle_out = dict(out)
    print(f"DONE RUNNING PADDLE, {round(time.time()-overall_start, 2)}s elapsed\n")
    
    
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
    out_df.to_csv(f"{output_dir}/{fasta_name_tmp}_Paddle_preds.csv", encoding='utf-8')
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", help="Input fasta file", type=str)
    parser.add_argument("-o", help="Output directory", type=str)
    args = parser.parse_args()
    main(args.f, args.o)
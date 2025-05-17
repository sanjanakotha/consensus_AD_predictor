# Provide the path to a fasta file containing your sequences
# INPUT PATH TO FASTA FILE HERE
fasta_name = '/global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/consensus_AD_predictor/data/yeast_TFs.fasta'
tmp_output = '/global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/consensus_AD_predictor/output/adpred_yeast'

# Python
import os
os.chdir("/global/scratch/projects/fc_mvslab/predictors/adpred/")

# Load necessary modules
import numpy as np
import pandas as pd
import tensorflow as tf
import matplotlib
from matplotlib import pyplot as plt
import matplotlib.backends.backend_pdf
from Bio import SeqIO
import glob, pickle, os, re, subprocess, csv, sys
from uuid import uuid4
import ipyparallel as ipp

#sys.path.append(os.path.abspath('..'))
from adpred import ADpred as adp


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
dview.execute('import pandas as pd')
dview.execute('import numpy as np')
dview.execute('import re, subprocess')
dview.execute('import Bio')
dview.execute('from adpred import ADpred as adp')

aa_lst = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']

# Set up directories for storing output, including PSIPRED and IUPred predictions
fasta_name_tmp = fasta_name.split('/')[-1].strip('.fasta').strip('.fa')
output_dir = tmp_output + '/' + fasta_name_tmp
subprocess.run(['mkdir', output_dir])
recs = list(SeqIO.parse(fasta_name, 'fasta'))
fasta_name = fasta_name_tmp
for r in recs:
    sequence = str(r.seq)
    # ADpred uses tiles of 30aa, but will pad shorter tiles, including those at the N-terminal end of the sequence
    #assert len(sequence) >= 30, 'Protein sequence must be at least 30 amino acids long'   
    assert len(set(sequence) - set(aa_lst)) == 0, 'No non-standard amino acids allowed'
    sub_fasta_name = re.sub(r'\s+', '_', r.id)
    subprocess.run(['mkdir', output_dir + '/' + sub_fasta_name])
    with open(output_dir + '/' + sub_fasta_name + '/' + sub_fasta_name + '.fasta', 'w') as fa:
        fa.write(sequence)

# path to a local installation of psipred
local_psipred = '/global/scratch/projects/fc_mvslab/predictors/psipred/runpsipred'


def adpred_read_ss2(ss2):
    '''
    Converts .ss2 data to horizontal secondary structure information (e.g., --HHHHH--HHHHH-HH)
    
    ss2: name of .ss2 file
    '''
    df = pd.read_table(ss2, skiprows=[0], sep='\s+', names=['position', 'residue', 'ss', 'c', 'h', 'e'])
    return ''.join(df.ss.values).replace('C','-')

def get_psipred(seq, output_dir=output_dir, fasta_name=fasta_name, computed=None):
    '''
    Runs PSIPRED on the provided sequence and returns secondary structure
    data for ADpred, viz 'adpred_ss_as_a_string'
    
    seq: provided sequence
    output_dir: path to output directory containing fasta file
    fasta_name: name of fasta file, minus the .fasta/.fa file ending
    computed: path to already computed .ss2 file to skip re-running PSIPRED
    ''' 
    if computed is not None:
        return adpred_read_ss2(computed)
    
    p = subprocess.run('cd "{}" && {} "{}".fa*'.format(output_dir, local_psipred, fasta_name), 
                       capture_output=True, shell=True, encoding='utf-8')
    rootname = re.search('Final output files: (.*).ss2 (.*).horiz', p.stdout).group(1).strip()
    
    f = open(output_dir + '/' + rootname + '.horiz')
    adpred_ss = ''.join([i.group(1) for i in re.finditer('Pred: (.*)\n', ''.join(f))]).replace("C","-")

    return adpred_ss

dview = c[:]
dview.block = True

# Push variables and functions to workers
var = dict( output_dir=output_dir, local_psipred=local_psipred, fasta_name=fasta_name, adpred_read_ss2=adpred_read_ss2, get_psipred=get_psipred )
dview.push(var)

# Need a wrapper function because map() only operates on one argument
def wrapper(r):
    sequence = str(r.seq)
    sub_fasta_name = re.sub(r'\s+', '_', r.id)
    adpred_ss = get_psipred(sequence, output_dir=output_dir + '/' + sub_fasta_name, fasta_name=sub_fasta_name,)
                                    #computed="{}/{}/{}.ss2".format(output_dir, sub_fasta_name, sub_fasta_name) )
    assert len(adpred_ss) == len(sequence)
    return (sequence, adpred_ss)

# Run a parallel map, executing the wrapper function on indices 0,...,n-1
lview = c.load_balanced_view()
# Cause execution on main process to wait while tasks sent to workers finish
lview.block = True 
out = lview.map(wrapper, recs)   # Run calculation in parallel
ss_dict = dict(out)


# Push PSIPRED output to workers
dview.push(dict( ss_dict=ss_dict ))

# Need a wrapper function because map() only operates on one argument
def wrapper(r):
    seq = str(r.seq)
    prot = adp.protein('x', seq)
    
    # Load horizontal secondary structure data
    prot.second_struct = ss_dict[seq]
    
    # Run ADpred
    prot.predict()
    return (prot.sequence, prot.predictions)


# Run a parallel map, executing the wrapper function on indices 0,...,n-1
lview = c.load_balanced_view()
# Cause execution on main process to wait while tasks sent to workers finish
lview.block = True 
out = lview.map(wrapper, recs)   # Run calculation in parallel
adpred_out = dict(out)

pdf = matplotlib.backends.backend_pdf.PdfPages(output_dir + "/ADpred_traces.pdf")

for r in recs:
    sequence = str(r.seq)
    adpred_preds = adpred_out[sequence]

    fig, ax = plt.subplots(figsize=(10,5))

    # Plot activity trace and activity threshold
    ax.plot(np.arange(len(sequence)), adpred_preds, c='orange', label='ADPred')
    ax.plot([0, len(sequence)-1], [.8, .8], '-', c='orange')

    ax.set_title(r.id)
    ax.legend()
    ax.set_ylim(0, 1)
    ax.margins(x=0)

    ax.set_xlabel('Center position of tile')
    ax.set_ylabel('Activity scores, normalized')

    # Plot protein sequence on a secondary x-axis
    # Comment out this block of code if it looks too messy for you
    ax2 = ax.twiny()
    ax2.set_xticks(np.arange(len(sequence)))
    ax2.set_xticklabels([aa for aa in sequence], fontsize=3)

    pdf.savefig()
    plt.show()
    
pdf.close()

# Save predicted values to a csv file
data = []
for r in recs:
    sequence = str(r.seq)
    adpred_preds = adpred_out[sequence]
    
    data.append([sequence, adpred_preds.tolist()])

out_df = pd.DataFrame(data, columns=['sequence', 'adpred_preds'])
out_df.to_csv(output_dir + "/ADpred_preds.csv", encoding='utf-8')



import multiprocessing
import numpy as np
import time
from localcider.sequenceParameters import SequenceParameters
import alphaPredict as alpha
from tqdm import tqdm
from dask_jobqueue.slurm import SLURMRunner
from dask.distributed import Client, progress
from contextlib import redirect_stdout
import io

def create_features_faster(sequences, SEQUENCE_WINDOW = 5, STEPS = 1, LENGTH = 40, PROPERTIES = 42, pools=1):  
    #### Setting up the parallel workers
    cluster = SLURMRunner()
    client = Client(cluster)
    client.wait_for_workers(5)    
    futures = client.map(calculate_properties, sequences)
    progress(futures)
    features_result = client.gather(futures)
        # These were printing annoying output        
    return np.array(features_result)

def calculate_properties(sequence, SEQUENCE_WINDOW = 5, STEPS = 1, LENGTH = 40, PROPERTIES = 42):
    #print("Calculating properties")
    SEQUENCE_LENGTH = len(sequence)
    SeqOb = SequenceParameters(sequence)
    
    # Break the sequence into smaller pieces of type "SequenceParameters"
    num_sub_seq = int(SEQUENCE_LENGTH-SEQUENCE_WINDOW/STEPS+1)
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



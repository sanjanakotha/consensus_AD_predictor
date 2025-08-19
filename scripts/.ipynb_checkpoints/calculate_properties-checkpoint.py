import pandas as pd
from joblib import Parallel, delayed
from tqdm import tqdm
import metapredict
import protfasta
import os
import argparse

# Parse command-line arguments
parser = argparse.ArgumentParser(description="Calculate protein properties and disorder")
parser.add_argument("--input", "-i", required=True, help="Path to input FASTA file")
parser.add_argument("--output", "-o", required=True, help="Path to output CSV file")
args = parser.parse_args()

# Read sequences
seq_dict = protfasta.read_fasta(args.input, invalid_sequence_action='remove')
seq_df = pd.DataFrame({"id": seq_dict.keys(), "seq": seq_dict.values()})

# Calculate net charge
seq_df["net charge"] = seq_df["seq"].str.count('R') + seq_df["seq"].str.count('K') - \
                       seq_df["seq"].str.count('D') - seq_df["seq"].str.count('E')
print("Calculated charge.")

# Calculate hydrophobic proportion
seq_df["hydrophobics_prop"] = (
    seq_df["seq"].str.count('W') +
    seq_df["seq"].str.count('F') +
    seq_df["seq"].str.count('Y') +
    seq_df["seq"].str.count('L')
) / seq_df["seq"].str.len()
print("Calcuated hydrophobics.")

# Define function for joblib
def calc_percent_disorder(seq):
    return metapredict.percent_disorder(seq)

# Use SLURM-aware number of cores
n_cores = int(os.environ.get("SLURM_CPUS_PER_TASK", 1))
print(f"Using {n_cores} cores for parallel processing")

# Calculate disorder in parallel
seq_df["percent_disorder"] = list(
    Parallel(n_jobs=n_cores)(
        delayed(calc_percent_disorder)(seq)
        for seq in tqdm(seq_df["seq"], desc="Calculating disorder")
    )
)
print("Calculated disorder.")

# Save results
seq_df.to_csv(args.output, index=False)
print(f"Results saved to {args.output}")

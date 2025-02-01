import numpy as np
import pandas as pd

def generate_random_sequences(num_seqs=4, seq_length=7):
    """
    Generates a list of random sequences composed of the base pairs A, C, G, T.
    
    Parameters:
        num_seqs (int): The number of sequences to generate.
        seq_length (int): The length of each sequence.
        
    Returns:
        list: A list of randomly generated sequences (as strings).
    """
    nucleotides = ['A', 'C', 'G', 'T']
    # Leverage numpy's random.choice to generate each sequence
    sequences = [''.join(np.random.choice(nucleotides, size=seq_length)) for _ in range(num_seqs)]
    return sequences

# Generate random sequences for testing
random_sequences = generate_random_sequences(num_seqs=10, seq_length=7)
print("Randomly Generated Sequences:")
print(random_sequences)
print("\n---------------------\n")

# Compute Counts and Frequencies from the sequences

# Step 1: initialise nucleotide positions
nucleotides = ['A', 'C', 'G', 'T']
motif_length = len(random_sequences[0])  # assume all sequences have the same length

# Step 2: compute Counts Matrix
counts_matrix = {nuc: [0] * motif_length for nuc in nucleotides}

for seq in random_sequences:
    for pos, base in enumerate(seq):
        counts_matrix[base][pos] += 1

# Convert counts_matrix to a Pandas DataFrame
counts_df = pd.DataFrame(counts_matrix, index=[f'Pos {i}' for i in range(motif_length)])
print("Counts Matrix:")
print(counts_df)
print("\n---------------------\n")

# Step 3: compute Frequency Matrix (applying Pseudocount of 0.5)
pseudocount = 0.5

# Add pseudocount to each count and update the total count per position accordingly.
# Since there are 4 nucleotides, we add pseudocount * 4 to the sum.
total_counts = counts_df.sum(axis=1) + (len(nucleotides) * pseudocount)

frequency_matrix = {
    nuc: [(counts_df[nuc].iloc[i] + pseudocount) / total_counts.iloc[i] for i in range(motif_length)]
    for nuc in nucleotides
}

# Convert frequency_matrix to a Pandas DataFrame
frequency_df = pd.DataFrame(frequency_matrix, index=[f'Pos {i}' for i in range(motif_length)])
print("Frequency Matrix with Pseudocounts:")
print(frequency_df)
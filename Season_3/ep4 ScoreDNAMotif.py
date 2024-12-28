def Count(motif):
    """
    Counts the frequency of each nucleotide (A, C, G, T) at each position in the motif.
    """
    count = {nt: [0] * len(motif[0]) for nt in "ACGT"}
    for sequence in motif:
        for i, nucleotide in enumerate(sequence):
            count[nucleotide][i] += 1
    return count

def Consensus(motif):
    """
    Finds the consensus sequence for a given motif.
    """
    k = len(motif[0])  # Length of sequences
    count = Count(motif)  # Get nucleotide counts
    consensus = ""  # Initialize consensus sequence
    
    for j in range(k):
        max_count = 0
        frequent_symbol = ""
        for symbol in "ACGT":
            if count[symbol][j] > max_count:
                max_count = count[symbol][j]
                frequent_symbol = symbol
        consensus += frequent_symbol  # Append the most frequent nucleotide
    
    return consensus

def Score(motif):
    """
    Calculates the score of a motif based on its consensus sequence.
    
    Args:
        motif (list): A list of DNA sequences of the same length.
        
    Returns:
        int: The total score (number of mismatches with the consensus sequence).
    """
    consensus = Consensus(motif)  # Find consensus sequence
    print("Consensus Sequence:", consensus)
    k = len(motif[0])  # Length of each sequence
    score = 0  # Initialize score
    
    # Calculate mismatches at each position
    for j in range(k):
        for sequence in motif:
            if sequence[j] != consensus[j]:
                score += 1  # Increment score for mismatches
    
    return score

# Example motif input
motif = [
    'AACGTA',
    'CCCGTT',
    'CACCTT',
    'GGATTA',
    'TTCCGG'
]

# Calculate and display the score of the motif
motif_score = Score(motif)
print("Motif Score:", motif_score)

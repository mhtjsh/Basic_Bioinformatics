def Count(motif):
    """
    Counts the frequency of each nucleotide (A, C, G, T) at each position in the motif.
    
    Args:
        motif (list): A list of DNA sequences of the same length.
        
    Returns:
        dict: A dictionary where keys are nucleotides and values are lists of frequencies.
    """
    # Initialize the count dictionary with zeros for each nucleotide
    count = {nt: [0] * len(motif[0]) for nt in "ACGT"}
    
    # Iterate through each sequence in the motif
    for sequence in motif:
        for i, nucleotide in enumerate(sequence):
            # Increment the count for the corresponding nucleotide and position
            count[nucleotide][i] += 1
    
    return count

# Example motif input
motif = [
    'AACGTA',
    'CCCGTT',
    'CACCTT',
    'GGATTA',
    'TTCCGG'
]

# Calculate and display nucleotide counts
counts = Count(motif)
print("Nucleotide Frequency Count:\n", counts)

def Count(motif):
    """
    Counts the frequency of each nucleotide (A, C, G, T) at each position in the motif.
    
    Args:
        motif (list): A list of DNA sequences of the same length.
        
    Returns:
        dict: A dictionary where keys are nucleotides and values are lists of frequencies.
    """
    count = {nt: [0] * len(motif[0]) for nt in "ACGT"}
    for sequence in motif:
        for i, nucleotide in enumerate(sequence):
            count[nucleotide][i] += 1
    return count

def Profile(motif):
    """
    Calculates the profile matrix for a given motif by converting nucleotide counts to probabilities.
    
    Args:
        motif (list): A list of DNA sequences of the same length.
        
    Returns:
        dict: A dictionary where keys are nucleotides and values are lists of probabilities.
    """
    t = len(motif)  # Total number of sequences
    k = len(motif[0])  # Length of each sequence
    profile = {}  # Initialize profile matrix
    
    # Get nucleotide counts using the Count function
    count_matrix = Count(motif)
    
    # Convert counts to probabilities
    for nt in "ACGT":
        profile[nt] = [count / t for count in count_matrix[nt]]
    
    return profile

# Example motif input
motif = [
    'AACGTA',
    'CCCGTT',
    'CACCTT',
    'GGATTA',
    'TTCCGG'
]

# Calculate and display the profile matrix
profile_matrix = Profile(motif)
print("Profile Matrix:\n", profile_matrix)

def HammingDistance(seq1, seq2):
    """
    Calculate the Hamming Distance between two DNA sequences of equal length.

    Args:
        seq1 (str): The first DNA sequence.
        seq2 (str): The second DNA sequence.

    Returns:
        int: The number of differing positions between the two sequences.
    """
    # Initialize a counter for differing positions
    counter = 0

    # Compare characters at each position in both sequences
    for i in range(len(seq1)):
        if seq1[i] != seq2[i]:
            counter += 1  # Increment counter for each mismatch

    return counter

# Example Execution
seq1 = "CTTGAAGTGGACCTCTAGTTCCTCTACAAAGAACAGGTTGACCTGTCGCGAAG"
seq2 = "ATGCCTTACCTAGATGCAATGACGGACGTATTCCTTTTGCCTCAACGGCTCCT"
print(HammingDistance(seq1, seq2))

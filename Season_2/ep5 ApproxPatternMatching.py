def HammingDistance(seq1, seq2):
    """
    Calculate the Hamming Distance between two DNA sequences of equal length.

    Args:
        seq1 (str): The first DNA sequence.
        seq2 (str): The second DNA sequence.

    Returns:
        int: The number of differing positions between the two sequences.
    """
    counter = 0
    for i in range(len(seq1)):
        if seq1[i] != seq2[i]:
            counter += 1
    return counter

def ApproximatePatternMatching(Text, Pattern, d):
    """
    Find all positions where a pattern matches the text with up to 'd' mismatches.

    Args:
        Text (str): The DNA sequence.
        Pattern (str): The pattern to search for.
        d (int): Maximum allowed mismatches.

    Returns:
        list: Positions where the pattern approximately matches the text.
    """
    positions = []  # Store matching positions
    n = len(Text)   # Length of the DNA sequence
    m = len(Pattern)  # Length of the pattern

    # Check all substrings of Text with length equal to the pattern
    for i in range(n - m + 1):
        # Extract the substring and calculate Hamming distance
        if HammingDistance(Text[i:i + m], Pattern) <= d:
            positions.append(i)  # Add position if within mismatch threshold

    return positions

# Example Execution
Text = "ATGCATATGACTACTAGATACTGATACTGATACATA"
Pattern = "ATAG"
d = 1
print(ApproximatePatternMatching(Text, Pattern, d))

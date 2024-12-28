import random

def RandomizedMotifSearch(DNA, k, t):
    """
    Perform Randomized Motif Search to find motifs in DNA sequences.

    Args:
        DNA (list): A list of DNA sequences.
        k (int): Length of the motif.
        t (int): Number of sequences in DNA.

    Returns:
        list: The best motifs found.
    """
    # Initialize motifs randomly
    motifs = RandomMotifs(DNA, k)
    best_motifs = motifs

    while True:
        profile = ProfileWithPseudocounts(motifs)
        motifs = MotifsFromProfile(DNA, k, profile)
        if Score(motifs) < Score(best_motifs):
            best_motifs = motifs
        else:
            return best_motifs


def RandomMotifs(DNA, k):
    """
    Generate a random set of motifs from DNA sequences.
    """
    t = len(DNA)
    motifs = []
    for sequence in DNA:
        start = random.randint(0, len(sequence) - k)
        motifs.append(sequence[start:start + k])
    return motifs


def ProfileWithPseudocounts(motifs):
    """
    Generate a profile matrix with pseudocounts for a set of motifs.
    """
    t = len(motifs) + 4  # Add pseudocounts for A, C, G, T
    k = len(motifs[0])
    profile = {nt: [1] * k for nt in "ACGT"}  # Initialize with pseudocounts

    for motif in motifs:
        for i, nt in enumerate(motif):
            profile[nt][i] += 1

    # Normalize by the total count
    for nt in profile:
        profile[nt] = [count / t for count in profile[nt]]

    return profile


def MotifsFromProfile(DNA, k, profile):
    """
    Generate motifs from a profile matrix for a set of DNA sequences.
    """
    return [ProfileMostProbablePattern(seq, k, profile) for seq in DNA]


def ProfileMostProbablePattern(Text, k, profile):
    """
    Find the most probable k-mer in a DNA sequence given a profile.
    """
    max_prob = -1
    most_probable = ""
    for i in range(len(Text) - k + 1):
        pattern = Text[i:i + k]
        prob = 1
        for j in range(k):
            prob *= profile[pattern[j]][j]
        if prob > max_prob:
            max_prob = prob
            most_probable = pattern
    return most_probable


def Score(motifs):
    """
    Calculate the score of a set of motifs based on the consensus sequence.
    """
    consensus = Consensus(motifs)
    score = 0
    for j in range(len(motifs[0])):
        score += sum(1 for motif in motifs if motif[j] != consensus[j])
    return score


def Consensus(motifs):
    """
    Find the consensus sequence for a set of motifs.
    """
    k = len(motifs[0])
    count = {nt: [0] * k for nt in "ACGT"}
    for motif in motifs:
        for i, nt in enumerate(motif):
            count[nt][i] += 1
    consensus = ""
    for j in range(k):
        consensus += max("ACGT", key=lambda nt: count[nt][j])
    return consensus


# Example Input
DNA = [
    "CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA",
    "GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG",
    "TAGTACCGAGACCGAAAGAAGTATACAGGCGT",
    "TAGATCAAGTTTCAGGTGCACGTCGGTGAACC",
    "AATCCACCAGCTCCACGTGCAATGTTGGCCTA"
]
k = 8
t = len(DNA)

# Run the Randomized Motif Search
result = RandomizedMotifSearch(DNA, k, t)
print("Best Motifs:", result)

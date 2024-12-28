def GreedyMotifSearch(DNA, k, t):
    """
    Perform Greedy Motif Search to find the best motifs in DNA sequences.

    Args:
        DNA (list): A list of DNA sequences.
        k (int): Length of the motif.
        t (int): Number of sequences in DNA.

    Returns:
        list: The best motifs found in the DNA sequences.
    """
    n = len(DNA[0])  # Length of each DNA sequence
    best_motifs = [DNA[i][0:k] for i in range(t)]  # Initialize best motifs with first k-mers

    # Iterate through all possible k-mers in the first sequence
    for i in range(n - k + 1):
        motifs = [DNA[0][i:i + k]]  # Start with a single k-mer
        for j in range(1, t):  # Build motifs for the remaining sequences
            profile = Profile(motifs[:j])  # Calculate the profile matrix for current motifs
            motifs.append(ProfileMostProbableKmer(DNA[j], k, profile))  # Find most probable k-mer
        if Score(motifs) < Score(best_motifs):  # Update best motifs if a better score is found
            best_motifs = motifs

    return best_motifs


# Example Functions Required
def Profile(motifs):
    """
    Calculate the profile matrix for a set of motifs.
    """
    t = len(motifs)
    k = len(motifs[0])
    profile = {nt: [0] * k for nt in "ACGT"}
    for motif in motifs:
        for i, nt in enumerate(motif):
            profile[nt][i] += 1
    for nt in "ACGT":
        profile[nt] = [x / t for x in profile[nt]]
    return profile


def ProfileMostProbableKmer(Text, k, Profile):
    """
    Find the most probable k-mer in a sequence given a profile matrix.
    """
    max_probability = -1
    most_probable_kmer = ''
    for i in range(len(Text) - k + 1):
        kmer = Text[i:i + k]
        probability = Pr(kmer, Profile)
        if probability > max_probability:
            max_probability = probability
            most_probable_kmer = kmer
    return most_probable_kmer


def Pr(Text, Profile):
    """
    Calculate the probability of a k-mer given a profile matrix.
    """
    p = 1
    for i in range(len(Text)):
        p *= Profile[Text[i]][i]
    return p


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
    "GGCGTTCAGGCA",
    "AAGAATCAGTCA",
    "CAAGGAGTTCGC",
    "CACGTCAATCAC",
    "CAATAATATTCG"
]
k = 3
t = 5

# Find the best motifs
result = GreedyMotifSearch(DNA, k, t)
print("Best Motifs:", result)

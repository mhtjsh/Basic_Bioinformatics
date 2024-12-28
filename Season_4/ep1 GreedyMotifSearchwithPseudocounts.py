def GreedyMotifSearchWithPseudocounts(DNA, k, t):
    """
    Perform Greedy Motif Search with Pseudocounts.

    Args:
        DNA (list): A list of DNA sequences.
        k (int): Length of the motif.
        t (int): Number of sequences in DNA.

    Returns:
        list: The best motifs found in the DNA sequences.
    """
    n = len(DNA[0])
    best_motifs = [DNA[i][:k] for i in range(t)]

    for i in range(n - k + 1):
        motifs = [DNA[0][i:i + k]]
        for j in range(1, t):
            profile = ProfileWithPseudocounts(motifs)
            motifs.append(ProfileMostProbableKmer(DNA[j], k, profile))
        if Score(motifs) < Score(best_motifs):
            best_motifs = motifs

    return best_motifs


# Supporting Functions
def ProfileWithPseudocounts(motifs):
    """
    Calculate the profile matrix with pseudocounts.
    """
    t = len(motifs)
    k = len(motifs[0])
    profile = {nt: [1] * k for nt in "ACGT"}  # Initialize pseudocounts as 1
    for motif in motifs:
        for i, nt in enumerate(motif):
            profile[nt][i] += 1
    for nt in "ACGT":
        profile[nt] = [x / (t + 4) for x in profile[nt]]  # Normalize with pseudocounts
    return profile


def ProfileMostProbableKmer(Text, k, profile):
    """
    Find the most probable k-mer in a sequence given a profile matrix.
    """
    max_probability = -1
    most_probable_kmer = ''
    for i in range(len(Text) - k + 1):
        kmer = Text[i:i + k]
        probability = Pr(kmer, profile)
        if probability > max_probability:
            max_probability = probability
            most_probable_kmer = kmer
    return most_probable_kmer


def Pr(Text, profile):
    """
    Calculate the probability of a k-mer given a profile matrix.
    """
    p = 1
    for i in range(len(Text)):
        p *= profile[Text[i]][i]
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

#Example input
DNA = [
    "GGCGTTCAGGCA",
    "AAGAATCAGTCA",
    "CAAGGAGTTCGC",
    "CACGTCAATCAC",
    "CAATAATATTCG"
]
k = 3
t = 5

result = GreedyMotifSearchWithPseudocounts(DNA, k, t)
print("Best Motifs:", result)

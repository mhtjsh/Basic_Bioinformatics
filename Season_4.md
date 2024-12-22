## 🧬 Season 4 - Colouring the Building an Final touch-ups
Now the Series finale is this particular season which leads to end of this particualar git repo and in turn ends your initial journey of understanding very major Bioinformatics algorithms as well.

### 🐍 Ep 1 GreedyMotifSearch with Pseudocounts
🧠 What is GreedyMotifSearchWithPseudocounts?

Motif finding is a critical task in bioinformatics to identify conserved patterns in DNA sequences. This variant of the GreedyMotifSearch algorithm incorporates pseudocounts, ensuring better performance on datasets with low motif representation by avoiding zero probabilities in profile matrices.


🧠 Why PseudoCounts?

They stabilize the profile matrix, preventing a probability of zero for any nucleotide, which might lead to the exclusion of possible motifs.
- Objective: To find motifs across sequences that minimize the score, accounting for both observed and potential variability.

🛠️ Program Code
```Python
def GreedyMotifSearchWithPseudocounts(Dna, k, t):
    """
    Perform Greedy Motif Search with Pseudocounts.

    Args:
        Dna (list): A list of DNA sequences.
        k (int): Length of the motif.
        t (int): Number of sequences in Dna.

    Returns:
        list: The best motifs found in the DNA sequences.
    """
    n = len(Dna[0])
    best_motifs = [Dna[i][:k] for i in range(t)]

    for i in range(n - k + 1):
        motifs = [Dna[0][i:i + k]]
        for j in range(1, t):
            profile = ProfileWithPseudocounts(motifs)
            motifs.append(ProfileMostProbableKmer(Dna[j], k, profile))
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
Dna = [
    "GGCGTTCAGGCA",
    "AAGAATCAGTCA",
    "CAAGGAGTTCGC",
    "CACGTCAATCAC",
    "CAATAATATTCG"
]
k = 3
t = 5

result = GreedyMotifSearchWithPseudocounts(Dna, k, t)
print("Best Motifs:", result)
```
Output of which looks something like this:
```
Best Motifs: ['GGC', 'AAG', 'CAA', 'CAC', 'CAA']
```

💡 Key concept of ``GreedyMotifSearch with Pseudocounts``
1. Input
 - ``DNA`` : list of DNA sequences
 - ``k`` : length of the motif to search for
 - ``t`` : number of DNA sequences in Dna
2. Randomized Initialization: Selects random k-mers (substrings of length k) from each DNA sequence as initial motifs
3. Iterative Refinement:
 -Score: Measures dissimilarity from the consensus sequence (lower is better)
4. Updates the best motifs if the new score is lower
5. Stops when no improvement is observed
6. Output is observed 

### 🐍 Ep 2 Randomized Motif Search
🧠 What is Randomized Motif Search?
- **Randomized Selection**: Unlike deterministic algorithms, this starts by randomly selecting motifs
- **Iterative Refinement**: The algorithm iteratively refines motifs by building profiles and scoring
- **Objective**: Find motifs that minimize the Score function (how dissimilar the motifs are to the consensus sequence
- **Stochastic Nature**: This randomness helps explore more possible motif combinations, potentially avoiding local optima

🛠️ Program Code
```Python
import random

def RandomizedMotifSearch(Dna, k, t):
    """
    Perform Randomized Motif Search to find motifs in DNA sequences.

    Args:
        Dna (list): A list of DNA sequences.
        k (int): Length of the motif.
        t (int): Number of sequences in Dna.

    Returns:
        list: The best motifs found.
    """
    # Initialize motifs randomly
    motifs = RandomMotifs(Dna, k)
    best_motifs = motifs

    while True:
        profile = ProfileWithPseudocounts(motifs)
        motifs = MotifsFromProfile(Dna, k, profile)
        if Score(motifs) < Score(best_motifs):
            best_motifs = motifs
        else:
            return best_motifs


def RandomMotifs(Dna, k):
    """
    Generate a random set of motifs from DNA sequences.
    """
    t = len(Dna)
    motifs = []
    for sequence in Dna:
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


def MotifsFromProfile(Dna, k, profile):
    """
    Generate motifs from a profile matrix for a set of DNA sequences.
    """
    return [ProfileMostProbablePattern(seq, k, profile) for seq in Dna]


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
Dna = [
    "CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA",
    "GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG",
    "TAGTACCGAGACCGAAAGAAGTATACAGGCGT",
    "TAGATCAAGTTTCAGGTGCACGTCGGTGAACC",
    "AATCCACCAGCTCCACGTGCAATGTTGGCCTA"
]
k = 8
t = len(Dna)

# Run the Randomized Motif Search
result = RandomizedMotifSearch(Dna, k, t)
print("Best Motifs:", result)
```
Output which looks something like this:
```
Best Motifs: ['TCTCGGGG', 'CCAAGGTG', 'TACAGGCG', 'TTCAGGTG', 'TCCACGTG']
```
💡 Key concept of ``RandomizedMotifSearch``
1. Input:
 - ``DNA`` - list of DNA sequences, each of the same length
 - ``k`` - length of the motifs (k-mers) to search for
 - ``t`` - number of DNA sequences provided
2. Random Motif Initialization: Selects one random k-mer from each sequence in ``Dna`` as the initial motif set using the ``RandomMotifs`` function
3. Profile Matrix Construction (with Pseudocounts):
 - Generates a profile matrix based on the frequency of nucleotides in the motifs, adding pseudocounts (+1 for stability)
 - Normalizes frequencies to probabilities for each position in the motif
4. Refining Motifs:
 - Finds the most probable k-mer in each sequence according to the current profile matrix using the ``ProfileMostProbablePattern`` function
 - Replaces the existing motifs with these updated k-mers.
5. Scoring:
 - Calculates the score of the motif set using the ``Score`` function, which measures mismatches from the consensus sequence
 - If the new score is better, update the best motifs
6. Stops iterating when no improvement in the score is observed and the output is printed


## 🧬 Season 3 - Adding Concrete slurry in the Skeletal Structure
Now season 3 and the upcoming last season involves a bit of Umbrella type of Code where I will walk you through each and every type of functions separately and then finally after using that small function codes we will make something very meaningful out of it. 

Therefore, I encourage you to just stick with me and just enjoy the miracle or magic that you are about to experience at the end during this journey.

Now to understand the further episodes (Ep) I would urge you to take a look at _NumPy_ an interesting and very efficient Python library to sort and visualize the Array is been used. I would urge you to go through the [_NumPy_ Documentation](https://numpy.org/doc/) a bit and this will help you to understand all the further episodes with a breeze.

### 🐍 Ep 0 - Basic Numpy Array understanding from DNA Motif Matrix
A motif matrix represents aligned DNA sequences. Each row corresponds to a sequence, and each column represents a position in the alignment.

🛠️ Program Code
```Python
import numpy as np

# Define a motif matrix as a list of lists
motif = [
    ['A', 'T', 'G', 'C'],
    ['A', 'G', 'T', 'G'],
    ['T', 'C', 'C', 'A'],
    ['G', 'T', 'A', 'G'],
    ['C', 'G', 'A', 'T']
]

# Convert the list of lists into a NumPy array
motif_array = np.array(motif)

# Extract the third column (index 2) using NumPy slicing
print("Extracted Column (using NumPy):", motif_array[:, 2])

# Alternative method: Extract the third column using a loop
column = []  # Initialize an empty list
for row in motif_array:
    column.append(row[2])  # Append the element at index 2 in each row

# Print the manually extracted column
print("Extracted Column (using Loop):", column)

# Display the original motif array
print("Motif Array:\n", motif_array)
```

Output of which looks something like this:
```
Extracted Column (using NumPy): ['G' 'T' 'C' 'A' 'A']
Extracted Column (using Loop): ['G', 'T', 'C', 'A', 'A']
Motif Array:
 [['A' 'T' 'G' 'C']
  ['A' 'G' 'T' 'G']
  ['T' 'C' 'C' 'A']
  ['G' 'T' 'A' 'G']
  ['C' 'G' 'A' 'T']]
```

### 🐍 Ep 1 - Counting Nucleotide Frequencies in DNA Motifs
This code demonstrates how to calculate nucleotide frequencies(``A``, ``T``, ``G`` and ``C``) (counts) for each position in a set of DNA motifs. The Count function is a foundational step in generating position weight matrices (PWMs) and understanding sequence conservation.

🛠️ Program Code
```Python
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
```

Output of which looks something like this:
```
Nucleotide Frequency Count:
 {'A': [1, 2, 1, 0, 0, 2],
  'C': [2, 1, 4, 2, 0, 0],
  'G': [1, 1, 0, 2, 1, 1],
  'T': [1, 1, 0, 1, 4, 2]}
```

💡 Key Concept of ``CountMotif``
1. A list of aligned DNA sequences (``motif``), where all sequences have the same length.
2. Initialize Count Dictionary where keys are ``A``, ``C``, ``G``, and ``T``, and values are lists of zeros. The length of each list matches the length of the sequences.
```Python
count = {nt: [0] * len(motif[0]) for nt in "ACGT"}
```
3. Update counts by iterating through each sequence and each position in the sequence.
4. As the output, a dictionary of nucleotide frequencies is built at every position.


### 🐍 Ep 2 - Profile Matrix from DNA Motifs
The Profile function calculates the probabilities of each nucleotide at every position in a set of aligned DNA sequences (motifs). This extends the Count function by converting counts into probabilities, enabling applications like motif scoring and PWM generation.

It is derived by **Dividing the counts of each nucleotide by the total number of sequences**.

🛠️ Program Code
```Python
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
```

Output of which looks something like this:
```
Profile Matrix:
 {'A': [0.2, 0.4, 0.2, 0.0, 0.0, 0.4],
  'C': [0.4, 0.2, 0.8, 0.4, 0.0, 0.0],
  'G': [0.2, 0.2, 0.0, 0.4, 0.2, 0.2],
  'T': [0.2, 0.2, 0.0, 0.2, 0.8, 0.4]}
```

💡 Key Concept of ``ProfileMotif``
1. A list of aligned DNA sequences (``motif``) with equal lengths as an input.
2. **Count Frequencies**: Use the Count function to calculate nucleotide counts at each position.
3. **Convert to Probabilities**: For each nucleotide, divide the count at each position by the total number of sequences.
4. A dictionary showing the probabilities of each nucleotide at every position is shown as output.

### 🐍 Ep 3 - Consensus Sequence from DNA Motifs
The Consensus function determines the most frequent nucleotide at each position in a motif. It uses the nucleotide Count function and constructs a single sequence, called the consensus sequence, which represents the most likely sequence in the motif.

 At each position, the **nucleotide with the highest count** is chosen, making it the most likely nucleotide at that position.

 🛠️ Program Code
 ```Python
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

def Consensus(motif):
    """
    Finds the consensus sequence for a given motif.
    
    Args:
        motif (list): A list of DNA sequences of the same length.
        
    Returns:
        str: The consensus sequence.
    """
    k = len(motif[0])  # Length of sequences
    count = Count(motif)  # Get nucleotide counts
    consensus = ""  # Initialize consensus sequence
    
    # Find the most frequent nucleotide at each position
    for j in range(k):
        max_count = 0
        frequent_symbol = ""
        for symbol in "ACGT":
            if count[symbol][j] > max_count:
                max_count = count[symbol][j]
                frequent_symbol = symbol
        consensus += frequent_symbol  # Append the most frequent nucleotide
    
    return consensus

# Example motif input
motif = [
    'AACGTA',
    'CCCGTT',
    'CACCTT',
    'GGATTA',
    'TTCCGG'
]

# Calculate and display the consensus sequence
consensus_sequence = Consensus(motif)
print("Consensus Sequence:", consensus_sequence)
```
Output of which looks like ``Consensus Sequence: CACCTA``

💡 Key Concept of ``ConsensusMotif``
1. A list of aligned DNA sequences (``motif``) of equal length as input.
2. Use the ``Count`` function to calculate the frequencies of each nucleotide at every position.
3. Identify the nucleotide with the highest frequency at each position and append it to the consensus sequence.
4. A single DNA sequence representing the consensus sequence as the output.

### 🐍 Ep 4 - Scoring DNA Motifs

Scoring evaluates how much each sequence in the motif deviates from the consensus sequence. It is calculated as the number of mismatches between each sequence and the consensus sequence at every position. This helps assess the level of conservation within the motif:
- **Low Score**: High conservation (the motif sequences closely match the consensus).
- **High Score**: Low conservation (the motif sequences vary significantly from the consensus).

🛠️ Program Code
```Python
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
```
Output of which looks something like this:
```
Consensus Sequence: CACCTA
Motif Score: 14
```
💡 Key concept of ``ScoringDNA`` is:
1. A list of DNA sequences (``motif``) of the same length as input
2. Compute the consensus sequence of the motif using the ``Consensus`` function
3. Count mismatches at every position and sum them up to calculate the score
4. The total number of mismatches (score) between the motif and its consensus sequence as output.

### 🐍 Ep 5 - Profile Most Probable K-mer (Similar to ``ConsensusMotif``)

The ProfileMostProbableKmer function is a utility for identifying the most likely k-mer (substring of length k) in a given DNA sequence (Text) based on a given profile matrix. 

This is especially useful in motif discovery, where we find conserved sequences in DNA.

🛠️ Program Code
```Python
def Pr(Text, Profile):
    """
    Calculate the probability of a given k-mer Text based on a profile matrix.

    Args:
        Text (str): A k-mer (substring of DNA sequence).
        Profile (dict): A profile matrix where keys are nucleotides and values 
                        are lists of probabilities for each position.

    Returns:
        float: Probability of Text being generated by Profile.
    """
    p = 1  # Initialize probability
    for i in range(len(Text)):
        p *= Profile[Text[i]][i]  # Multiply probabilities for each nucleotide
    return p


def ProfileMostProbableKmer(Text, k, Profile):
    """
    Find the most probable k-mer in Text based on a given profile matrix.

    Args:
        Text (str): DNA sequence.
        k (int): Length of the k-mer.
        Profile (dict): A profile matrix where keys are nucleotides and values 
                        are lists of probabilities for each position.

    Returns:
        str: The most probable k-mer in Text based on Profile.
    """
    max_probability = -1  # Initialize max probability as a very low value
    most_probable_kmer = ''  # Initialize the most probable k-mer

    # Slide through the Text to consider all possible k-mers
    for i in range(len(Text) - k + 1):
        kmer = Text[i:i + k]  # Extract k-mer
        probability = Pr(kmer, Profile)  # Calculate its probability

        # Update the most probable k-mer if a higher probability is found
        if probability > max_probability:
            max_probability = probability
            most_probable_kmer = kmer

    return most_probable_kmer


# Example Input
Text = "ACGTACGTGACG"
k = 3
Profile = {
    'A': [0.2, 0.4, 0.3],
    'C': [0.4, 0.3, 0.1],
    'G': [0.3, 0.1, 0.5],
    'T': [0.1, 0.2, 0.1]
}

# Find and display the most probable k-mer
result = ProfileMostProbableKmer(Text, k, Profile)
print("Most Probable k-mer:", result)
```
Output of which looks somehthing like this ``Most Probable k-mer: ACG``

💡 Key concept of ``ProfileMostProbableKmer``
1. Input:
 - ``Text``: A DNA sequence in which we search for the most probable k-mer.
 - ``k``: Length of the k-mer.
 - ``Profile``: A profile matrix that contains probabilities for each nucleotide at each position
2. Use the ``Pr`` function to calculate the probability of each k-mer in ``Text`` being generated by ``Profile``.
3. Keep track of the k-mer with the highest probability.
4. The k-mer with the highest probability according to the profile matrix.

### 🐍 Ep 6 - Greedy Motif Search
🧠 What is ``GreedyMotifSearch``?

- **Motifs** : Short, conserved subsequences in DNA found across multiple sequences (e.g., regulatory sites).
- **Greedy Search** : The algorithm builds motifs iteratively by considering one sequence at a time, refining the motif set based on a scoring function.
- **Objective** : Minimize the **Score** function, which quantifies how similar the chosen motifs are to the consensus sequence.

The ``GreedyMotifSearch`` function is a simple yet effective algorithm to identify motifs (conserved patterns) in a set of DNA sequences. It starts with an initial guess and iteratively refines it to find the best possible motifs by incorporating probabilistic scoring and profiles.

🛠️ Program Code
```Python
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
```
Output of which looks something like this:
```
Best Motifs: ['CAG', 'CAG', 'CAA', 'CAA', 'CAA']
```

💡 Key Concept:
1. Input:
 - ``DNA`` : A list of DNA sequences in which motifs are to be found
 - ``k`` : The length of the motif
 - ``t`` : The number of sequences in the DNA list
2. Start with all possible k-mers from the first sequence
3. Use a Profile matrix to iteratively refine the motifs for subsequent sequences
4. Update the best motifs based on the Score function
5. The best motifs (one from each sequence) that minimize the score

### End of the episode with Greedy algorithm Search being the final algorithm which involves conglomeration of all the functions used and understood earlier.

### 🐍 Ep Extended - Finding Patterns with Mismatches
This script identifies occurrences of a specific pattern within a list of DNA strings, allowing for a defined number of mismatches. It is particularly useful for locating approximate matches in sequences where errors or variations might occur.

🛠️ Program Code
```Python
def find_pattern_with_mismatches(text, pattern, max_mismatches):
    """
    Find all occurrences of a pattern in a text with at most a given number of mismatches.

    Args:
        text (str): The main sequence to search.
        pattern (str): The pattern to find.
        max_mismatches (int): Maximum allowed mismatches.

    Returns:
        list: Starting positions of the pattern in the text with at most max_mismatches.
    """
    locations = []
    for i in range(len(text) - len(pattern) + 1):
        mismatches = 0
        for j in range(len(pattern)):
            if text[i + j] != pattern[j]:
                mismatches += 1
                if mismatches > max_mismatches:
                    break
        if mismatches <= max_mismatches:
            locations.append(i)
    return locations


# Example Input: Multiple DNA strings
sample_strings = [
    "XABACDXABACD",
    "CDEFCDEHABAC",
    "GHIJKLABACDa",
    "PQRSTUABACDX",
    "VWXYZABACDcA",
    "HIJKLMNOPQRST",
    "CDEFGHIJKLMN",
    "QRSTUVWXYZAB",
    "EFGHIJKLMNOP",
    "UVWXYZABCDEF"
]

pattern = "ABACD"  # Pattern to search for
max_mismatches = 2  # Allow up to 2 mismatches

# Process each string to find the pattern with mismatches
for text in sample_strings:
    locations = find_pattern_with_mismatches(text, pattern, max_mismatches)
    if locations:
        print(f"String: {text}, Pattern locations (with at most {max_mismatches} mismatches): {locations}")
    else:
        print(f"String: {text}, Pattern not found with at most {max_mismatches} mismatches.")
```
Output of which looks something like this:
```
String: XABACDXABACD, Pattern locations (with at most 2 mismatches): [1, 7]
String: CDEFCDEHABAC, Pattern not found with at most 2 mismatches.
String: GHIJKLABACDa, Pattern locations (with at most 2 mismatches): [6]
String: PQRSTUABACDX, Pattern locations (with at most 2 mismatches): [6]
String: VWXYZABACDcA, Pattern locations (with at most 2 mismatches): [5]
String: HIJKLMNOPQRST, Pattern not found with at most 2 mismatches.
String: CDEFGHIJKLMN, Pattern not found with at most 2 mismatches.
String: QRSTUVWXYZAB, Pattern not found with at most 2 mismatches.
String: EFGHIJKLMNOP, Pattern not found with at most 2 mismatches.
String: UVWXYZABCDEF, Pattern not found with at most 2 mismatches.
```
💡 Key concept:
1. Input parameters:
 - ``Text`` : DNA strings where we search for approximate matches
 - ``Pattern`` : The pattern we're trying to locate
 - ``max_mismatches`` : The tolerance for mismatches between the pattern and substrings
2. For each possible substring of ``text`` compare it to pattern character by character and Count mismatches and stop early if they exceed ``max_mismatches``
3. Add valid starting indices to the locations list
4. The starting indices of all substrings in ``text`` that match the pattern within the mismatch threshold

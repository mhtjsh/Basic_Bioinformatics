## 🧬 Season 3 - Adding Concrete slurry in the Skeletal Structure
Now season 3 and the upcoming last season involves a bit of Umbrealla type of Code where I will walk you through each and every type of functions seperately and then finally after using that small function codes we will make something very meaningful out of it. 

Therefore, I encourage you to jus stick with me and just enjoy the miracle or magic that you are about to experince at the end durning this journey.

Now to undertstand the further episodes (Ep) I would urge you to take a look at _NumPy_ an intresting and very efficient Python library to sort and visualize the Array is been used. I would urge you to go through the [_NumPy_ Documentation](https://numpy.org/doc/) a bit and this will help you to understand all the further episodes with a breeze.

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
 {'A': [2, 1, 0, 1, 0, 2],
  'C': [2, 3, 4, 1, 2, 0],
  'G': [0, 1, 1, 2, 1, 2],
  'T': [1, 0, 0, 3, 4, 1]}
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
 {'A': [0.4, 0.2, 0.0, 0.2, 0.0, 0.4],
  'C': [0.4, 0.6, 0.8, 0.2, 0.4, 0.0],
  'G': [0.0, 0.2, 0.2, 0.4, 0.2, 0.4],
  'T': [0.2, 0.0, 0.0, 0.6, 0.8, 0.2]}
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

💡 Key Concept of ``ConsenusMotif``
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

### 🐍 Ep 5 - 


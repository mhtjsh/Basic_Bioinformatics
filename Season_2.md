## 🧬 Season 2 - Adding steel rods as skeleton for the building
A bit more complex and tough concepts are being touched upon which needed a lot more efforts to be understood are taken into account here, and I am giving my level best to make it most simpler interpretation of that concept to understand and trying to address each and every question which I had and possibly you could have during your learning journey

### 🐍 Ep 1 - Symbol array in DNA sequence
🧠 What exactly is symbol array?
  - The Symbol Array calculates how frequently a specific nucleotide (e.g., "A", "C", "G", or "T") appears in a half-length sliding window across a given DNA sequence.
  This approach is particularly useful when analyzing circular genomes or trying to identify regions that are enriched with a particular symbol.

🧠 How does it actually work?
  - Suppose, input Genome ``ATGATAGTCCGAAA`` has a length of _n_ = 14
  - The function extends this genome to simulate a circular genome by appending half of the genome (_n/2_) to the end
    - Extended Genome ``ATGATAGTCCGAAA`` → ``ATGATAGTCCGAAAATGATAG``
  - For each position _i_ in the genome, it calculates the count of the symbol "A" within a half-length window (n/2=7) starting at _i_

🧠 Why Half ``n/2`` the Genome length?
- Because of the circularized Genome nature i.e. By extending the genome (``ExtendedGenome``), you ensure that sliding windows near the end of the sequence "wrap around" and include the beginning of the sequence. This avoids missing patterns that cross the genome's end. Also, dividing the genome into segments of half its length is a common approach to detect local enrichments of a nucleotide.

  This primarily involves 2 major programs
  1. ``SymbolArray`` : Extends the genome and calculates the symbol counts.
  2. ``PatternCount``: Counts how many times a symbol appears in a window.

🛠️ Program Code
```Python
def SymbolArray(Genome, symbol):
    array = {}  # Stores counts for each starting position
    n = len(Genome)
    ExtendedGenome = Genome + Genome[0:n//2]  # Extend the genome to handle circular overlap
    for i in range(n):
        array[i] = PatternCount(symbol, ExtendedGenome[i:i+(n//2)])  # Count occurrences in each half-window
    return array

def PatternCount(Pattern, Text):
    count = 0
    for i in range(len(Text) - len(Pattern) + 1):
        if Text[i:i+len(Pattern)] == Pattern:  # Check for matches
            count += 1
    return count

# Example usage
print(SymbolArray("ATGATAGTCCGAAA", "A"))
```

Output of which something looks like this:
```
{0: 3, 1: 2, 2: 2, 3: 2, 4: 1, 5: 2, 6: 2, 7: 3, 8: 4, 9: 4, 10: 4, 11: 5, 12: 4, 13: 4}
```
which basically explains us:
1. Each key (e.g., ``0``, ``1``, ``2``...) represents a starting position in the genome.
2. Each value (e.g., ``3``, ``2``) indicates how many times the symbol "A" appears in a window of size 7 starting from that position.

For example :
- At position ``0``, the window is ``ATGATA`` → "A" appears 3 times.
- At position ``3``, the window is ``ATAGTC`` → "A" appears 2 times.
- At position ``10``, the window is ``AAAATG`` → "A" appears 4 times.

💡 In a nutshell,
1. Regions with high counts indicate **"A"-rich regions** (or symbol-rich areas).
2. By sliding the window across the genome, it allows us to detect where the symbol of interest is concentrated.

### 🐍 Ep 1 **Extended** - Symbol array in DNA but faster
The ``FasterSymbolArray`` function is an optimized version of the ``SymbolArray`` function. Instead of recalculating the PatternCount for every sliding window (as done previously), it reduces redundant calculations by incrementally updating the count.

🧠 Why is it Faster?
 - Compared to SymbolArray, which recalculates the PatternCount from scratch for every position ``FasterSymbolArray`` avoids redundant computations by only adjusting the count when the window slides.

    Now I know this doesn't makes a much of difference at a smaller level but as we upscale the production it increases the time required to compile the program, for CS enthusiasts it is known as ``Time Complexity``, so here by using ``FasterSymbolArray`` it reduces the time complexity significantly by linearizing it to _O(n)_ as compared to the normal approach.

🛠️ Program Code
```Python
def PatternCount(Pattern, Text):
    # Counts occurrences of a Pattern in Text
    count = 0
    for i in range(len(Text) - len(Pattern) + 1):
        if Text[i:i+len(Pattern)] == Pattern:
            count += 1
    return count

def FasterSymbolArray(Genome, symbol):
    # Optimized version of SymbolArray using sliding window technique
    array = {}
    n = len(Genome)
    ExtendedGenome = Genome + Genome[0:n//2]

    # Initialize with the count of the symbol in the first half of the genome
    array[0] = PatternCount(symbol, Genome[0:n//2])

    for i in range(1, n):
        # Start with previous array value
        array[i] = array[i-1]

        # Adjust count for the sliding window
        if ExtendedGenome[i-1:i-1+len(symbol)] == symbol:
            array[i] -= 1  # Subtract 1 if symbol exits the window
        if ExtendedGenome[i + (n//2) - 1:i + (n//2) - 1 + len(symbol)] == symbol:
            array[i] += 1  # Add 1 if symbol enters the window

    return array

# Example execution
print(FasterSymbolArray("ATGCATATGACTACTAGATACTGATACTGATACATA", "AT"))
```

Output of which looks something like this:
```
{0: 3, 1: 3, 2: 3, 3: 3, 4: 3, 5: 3, 6: 3, 7: 3, 8: 3, 9: 3, 10: 3, 11: 3, 12: 3, 13: 3, 14: 3,
15: 3, 16: 3, 17: 3, 18: 3, 19: 3, 20: 3, 21: 3, 22: 3, 23: 3, 24: 3, 25: 3, 26: 3, 27: 3, 28: 3,
29: 3, 30: 3, 31: 3, 32: 3, 33: 3, 34: 3, 35: 3}
```

💡 Key concept of ``FasterSymbolArray`` is:
1. ``PatternCount`` function initialization which Counts how many times a pattern (``symbol``) occurs in a given text.
2. Calculation of ``ExtendedGenome`` to simulate a **circular genome** and ensure patterns near the end of the genome can "wrap around."
   ```Python
   ExtendedGenome = Genome + Genome[0:n//2]
   ```
3. Initialization of ``FasterSymbolArray`` which start by calculating the count of the symbol in the first half of the genome.
   ```Python
   array[0] = PatternCount(symbol, Genome[0:n//2])
   ```
4. Sliding the window in ``FasterSymbolArray``
 - The array is initialized with the count of the symbol in the first half of the genome.
 - As the window slides by one position, the count is adjusted efficiently.
    - Subtract **1** if the symbol exits the window. ``(i-1)`` The symbol which leaves the window (leftmost position).
    - ``(i+n//2-1)`` : The symbol enters the 
    - Add **1** if the symbol enters the window.
5. Finally returning the result of array of counts for each starting position in the Genome.

### 🐍 Ep 2 Skew Array in DNA
🧠 What exactly is Skew an array means?

Skew function calculates the skewness of the data set. 
- skewness = 0 : normally distributed.
- skewness > 0 : more weight in the left tail of the distribution.
- skewness < 0 : more weight in the right tail of the distribution.

This is how the graph looks like:

![Skew Graphs](https://github.com/user-attachments/assets/fe9f9120-67a5-46a3-b5ce-f5db04abdc6c)

🧠 How is it used in BioInformatics?

The ``SkewArray`` function calculates the cumulative skew of a DNA sequence, which is the difference between the number of G and C nucleotides encountered at each position. This is particularly useful for identifying replication origins where the imbalance between G and C is at its minimum, which we will also discuss which will be ``MinimumSkew``
- **G** increases the skew value by +1.
- **C** decreases the skew value by -1.
- **A** and **T** do not affect the skew.

🛠️ Program Code ↔ By listing method (Method 1)
```Python
def SkewArray(Genome):
    # Initialize skew array starting at 0
    skew = [0]
    for i in range(len(Genome)):
        if Genome[i] == 'A' or Genome[i] == 'T':
            skew.append(skew[-1])  # No change for A or T
        elif Genome[i] == 'G':
            skew.append(skew[-1] + 1)  # Increment skew for G
        elif Genome[i] == 'C':
            skew.append(skew[-1] - 1)  # Decrement skew for C
    return skew

# Example Execution
print(SkewArray("CATGGGCATCGGCCATACGCC"))
```
or 

🛠️ Program Code ↔ By Dictionary method (Method 2)
```Python
def SkewArray(Genome):
    skew = {0: 0}  # Initialize dictionary with skew at index 0 as 0
    for i in range(len(Genome)):
        if Genome[i] == 'A' or Genome[i] == 'T':
            skew[i + 1] = skew[i]  # Same skew if A or T
        elif Genome[i] == 'G':
            skew[i + 1] = skew[i] + 1  # Increment skew for G
        elif Genome[i] == 'C':
            skew[i + 1] = skew[i] - 1  # Decrement skew for C
    return skew

print(SkewArray("GATACACTTCCCGAGTAGGTACTG"))
```

Output of both looks something like this:
``[0, 0, -1, 0, 1, 2, 3, 2, 1, 1, 1, 2, 3, 2, 1, 0, -1, -1, -2, -1, -2]``

💡 Key concept of ``SkewArray``
1. Skew Initialisation : The skew starts at 0, representing a balance between G and C at the beginning of the genome. ``skew = [0]``
2. Skew Update Rule: For each nucleotide:
   - G - +1 to current Skew
   - C - -1 to current Skew
```Python
if Genome[i] == 'G':
    skew.append(skew[-1] + 1)
elif Genome[i] == 'C':
    skew.append(skew[-1] - 1)
else:
    skew.append(skew[-1])
```
3. Cumulative calculation: The skew is calculated cumulatively along the genome, building a sequence of skew values.
4. Returning the output in a list where each element represents the skew at that position in the genome.
   
### 🐍 Ep 3 Minimum Skew in a DNA Sequence
🧠 Why exactly use ``MinimumSkew`` function?

The positions of minimum skew typically represent regions of high C content relative to G, marking potential **"replication origin"** sites.

🛠️ Program Code
```Python
def SkewArray(Genome):
    # Initialize skew array starting at 0
    skew = {0: 0}
    for i in range(len(Genome)):
        if Genome[i] == 'A' or Genome[i] == 'T':
            skew[i + 1] = skew[i]  # No change for A or T
        elif Genome[i] == 'G':
            skew[i + 1] = skew[i] + 1  # Increment skew for G
        elif Genome[i] == 'C':
            skew[i + 1] = skew[i] - 1  # Decrement skew for C
    return skew

def MinimumSkew(Genome):
    # Identify positions with the minimum skew
    positions = []
    skew_raw = SkewArray(Genome)
    min_value = min(skew_raw.values())  # Find minimum skew value
    for key, value in skew_raw.items():
        if value == min_value:
            positions.append(key)  # Collect positions with minimum skew
    return positions

# Example Execution
print(MinimumSkew("GATACACTTCCCGAGTAGGTACTG"))
```

Output of which looks something like this ``[10, 23]``

💡 Key Concept of ``MinimumSkew`` is:
1. SkewArray: Calculates the cumulative skew for each position in the genome and gives a dictionary where keys represent genome positions and values represent cumulative skew values as output which looks like
```Python
SkewArray("CATGGGCATCGGCC")
# Output: {0: 0, 1: 0, 2: -1, 3: 0, 4: 1, 5: 2, 6: 3, 7: 2, ...}
```
2. Finding the Minimum Skew:
 - Compute the skew array using the SkewArray function
 - Identify the minimum skew value in the array.
 - Collect all positions where this minimum value occurs.
```Python
min_value = min(skew_raw.values())
for key, value in skew_raw.items():
    if value == min_value:
        positions.append(key)
```
3. Output: The minimum skew value occurs at positions 10 and 23 (here), which signifies that these positions represent regions where the cumulative G-C imbalance is at its lowest. Also, these are likely candidates for replication origin sites in the genome. 

### 🐍 Ep 4 Calculating Hamming Distance between DNA sequence
🧠 What is Hamming Distance?

The Hamming distance between two strings of equal length is the number of positions at which the corresponding elements (characters or nucleotides) are different.

🧠 Where exactly  is ``HammingDistance`` used?

It is primarirly used in **Mutation Detection** and **Sequence Comparison**

🛠️ Program Code
```Python
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
```

Output of which is ``37``

💡 Key Concept of ``HammingFunction`` is:
1. It iterate over the position of both the sequences
2. Compare the corresponding characters at each position.
3. If the characters differ, increment the counter.
4. Output ie. ``37`` shows that the two sequences differs at 37 positions, which shows that there are 37 mutations in the between the DNA sequence.

### 🐍 Ep 5 Approximate Pattern Matching in DNA Sequences

🧠 What is Approximate Pattern Matching?

Approximate Pattern Matching identifies locations where a pattern almost matches a part of the sequence, allowing for a limited number of mismatches.


🧠 Where is ``ApproximatePatternMatching`` function used?

Primarily it is used in **Mutation analysis** and **Sequence Search**.

🛠️ Program Code
```Python
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
```
Output of which looks something like this: ``[6, 17, 23]``

💡 Key concept of ``ApproximatePatternMatching``
1. We are giving 3 input variables:
  - **Text**: DNA sequence to search through
  - **Pattern**: The motif or sequence to search for
  - **d**: Maximum allowed mismatches (Hamming distance)
2. Logic:
  - Extract every substring of ``Text`` with the same length as ``Pattern``
  - Calculate the **Hamming Distance** between the substring and the pattern
  - If the Hamming Distance is ≤ ``d``, record the position
3. Output tells that the pattern "ATAG" (here) appears approximately at positions 6, 17, and 23 in the sequence, where each match has exact 1 mismatch i.e. ``d=1``

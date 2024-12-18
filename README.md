# 🦠 Bioinformatics Algorithm 
## An  understanding of how Biopython and its algorithm works generally on Python
I am sharing my code repository to understand the basic working and principle of various BioInformatics code in Python; I am writting this github repository with 2 main intensions of mine first and foremost is to build myself a code base to look back onto how I learnt the Bioinformatics and understood how the skeleton of code works/interacts with my subject of interest i.e. Biology and second is to make you people understand how the BioPython a popluar Biological library in Python which almost feel like a cheat sheet to solve complex algorithms in Biologial genetic data.

Hence, Let's dive into the journey with me where I start from basic to the complex concepts which is actually necessary to solve the real world calculations in Bioinformatics.

## 🧬 Season 1 - Diggin down deep
A very basic intial Strucutral concept explained which primarirly involves this concepts:
1. Calcuate pattern
2. Frequecney map of DNA sequence
3. Most frequent K-mers
4. Complementary DNA string
5. Pattern Matching

This all seems very intimaidating at first but trust me all of 'em have been explained in a very intresting manner, I urge you to take a look at it and you'll surely thank me later considering that.

## 🧬 Season 2 - Adding steel rods as skeleton for the building
A bit more complex and tough concepts are being touched upon which needed a lot more efforts to be understood are taken into account here, and I am giving my level best to make it most simpler interpreation of that concept to understand and trying to address each and every question which I had and possibly you could have durning your learning journey

### 🐍 Ep 1 - Symbol array in DNA sequence
🧠 What exactly is symbol array?
  - The Symbol Array calculates how frequently a specific nucleotide (e.g., "A", "C", "G", or "T") appears in a half-length sliding window across a given DNA sequence.
  This approach is particularly useful when analyzing circular genomes or trying to identify regions that are enriched with a particular symbol.

🧠 How does it actually work?
  - Suppose, input Genome ``ATGATAGTCCGAAA`` have length of _n_ = 13
  - The function extends this genome to simulate a circular genome by appending half of the genome (_n/2_) to the end
    - Extended Genome ``ATGATAGTCCGAAA`` → ``ATGATAGTCCGAAAATGA``
  - For each position _i_ in the genome, it calculates the count of the symbol "A" within a half-length window (n/2=6) starting at _i_

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
{0: 4, 1: 4, 2: 4, 3: 3, 4: 3, 5: 3, 6: 3, 7: 3, 8: 3, 9: 4, 10: 4, 11: 4, 12: 4}
```
which basically explains us:
1. Each key (e.g., ``0``, ``1``, ``2``...) represents a starting position in the genome.
2. Each value (e.g., ``4``, ``3``) indicates how many times the symbol "A" appears in a window of size 6 starting from that position.

For example :
- At position ``0``, the window is ``ATGATA`` → "A" appears 4 times.
- At position ``3``, the window is ``ATAGTC`` → "A" appears 3 times.
- At position ``10``, the window is ``AAAATG`` → "A" appears 4 times.

💡 In a nutshell,
1. Regions with high counts indicate **"A"-rich regions** (or symbol-rich areas).
2. By sliding the window across the genome, it allows us to detect where the symbol of interest is concentrated.

### 🐍 Ep 1 **Extended** - Symbol array in DNA but faster
The ``FasterSymbolArray`` function is an optimized version of the ``SymbolArray`` function. Instead of recalculating the PatternCount for every sliding window (as done previously), it reduces redundant calculations by incrementally updating the count.

🧠 Why is it Faster?
 - Compared to SymbolArray, which recalculates the PatternCount from scratch for every position ``FasterSymbolArray`` avoids redundant computations by only adjusting the count when the window slides.

    Now I know this dosen't makes a much of difference at a smaller level but as we upscale the production it increases the time required to compile the program, for CS enthusiasts it is known as ``Time Complexity``, so here by using ``FasterSymbolArray`` it reduces the time complexity significantly by linearizing it to _O(n)_ as compared to the normal appraoch.

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
1. ``PatternCount`` fucntion intialization which Counts how many times a pattern (``symbol``) occurs in a given text.
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

The ``SkewArray`` function calculates the cumulative skew of a DNA sequence, which is the difference between the number of G and C nucleotides encountered at each position. This is particularly useful for identifying replication origins where the imbalance between G and C is at its minimum, which we will also disucss which will be ``MinimumSkew``
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
1. Skew Initialsation : The skew starts at 0, representing a balance between G and C at the beginning of the genome. ``skew = [0]``
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
3. Cummulative calculation: The skew is calculated cumulatively along the genome, building a sequence of skew values.
4. Returning the output in a list where each element represents the skew at that position in the genome.
   
### 🐍 Ep 3 Minimum Skew in a DNA Sequence

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

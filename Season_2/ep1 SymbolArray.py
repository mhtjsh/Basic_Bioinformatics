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

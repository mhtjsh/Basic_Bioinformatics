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

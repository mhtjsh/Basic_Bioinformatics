def PatternMatching(Pattern, Genome):
    Positions = []  # Output list to store positions
    for i in range(len(Genome) - len(Pattern) + 1):
        if Genome[i:i+len(Pattern)] == Pattern:  # Check for the match
            Positions.append(i)  # Append position if match found
    return Positions

# Example usage
text = 'CGCGATACGTTACATACATGATAGACCGCGCGCGATCATATCGCGATTATC'
pattern = 'CGCG'
print(PatternMatching(pattern, text))

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

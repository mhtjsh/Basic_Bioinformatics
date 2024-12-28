# Method 1 By List method:

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

print(SkewArray("CATGGGCATCGGCCATACGCC"))



# Method 2 By Dictionary Method:

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

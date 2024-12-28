def FrequencyMap(Text, k):
    freq = {}
    n = len(Text)
    for i in range(n - k + 1):
        Pattern = Text[i:i + k]
        if Pattern in freq:
            freq[Pattern] += 1
        else:
            freq[Pattern] = 1
    return freq

# Example usage
text = 'ACGTTGCATGTCGCATGAGCATGAGAGCT'
k = 3
print(FrequencyMap(text, k))

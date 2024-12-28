def Complement(Pattern):
    compl = ''
    for i in Pattern:
        if i == 'A':
            compl = compl + 'T'
        elif i == 'T':
            compl = compl + 'A'
        elif i == 'C':
            compl = compl + 'G'
        elif i == 'G':
            compl = compl + 'C'
    return compl

# Example usage
print(Complement('ACGTTGCATGTCGCATGAGCATGAGAGCT'))

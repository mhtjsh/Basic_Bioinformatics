# 🦠 Basic Bioinformatics - An  understanding of how Biopython works generally on Python
I am sharing my code repository to understand the basic working and principle of various BioInformatics code in Python; I am writting this github repository with 2 main intensions of mine first and foremost is to build myself a code base to look back onto how I learnt the Bioinformatics and understood how the skeleton of code works/internacts with my subject of intrest i.e. Biology and second is to make you people undertand how the BioPython a popluar Biological library in Python which almost feel like a cheat sheet to solve complex algorithms in Biologial genetic data.

Hence, Let's dive into the journey with me where I start from basic to the complex concepts which is actually necessary to solve the real world calculations in Bioinformatics.

## 🧬 Season 1 - Diggin down deep
### Ep 1 - Let's calcualte the basic Pattern given in the series of Genetic sequence code
```Python
text = 'ATGCATATGACTACTAGATACTGATACTGATACATA'
pattern = 'ATA'
print(len(pattern))
print(len(text))
print(text[0:3])
```

Let's suppose this is the basic dataset given to us, in which we are asked to find the pattern "ATA" from the Sequence i.e. how much time this particlaur motif occurs or is present in the Parent sequence.

```Python
text = 'ATGCATATGACTACTAGATACTGATACTGATACATA'
pattern = 'ATA'
def PatternCount(Text, Pattern):
    count = 0
    for i in range(len(Text)-len(Pattern)+1):
        if Text[i:i+len(Pattern)] == Pattern:
            count = count+1 
    return count 
print(PatternCount(text, pattern))
```

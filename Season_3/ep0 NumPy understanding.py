import numpy as np

# Define a motif matrix as a list of lists
motif = [
    ['A', 'T', 'G', 'C'],
    ['A', 'G', 'T', 'G'],
    ['T', 'C', 'C', 'A'],
    ['G', 'T', 'A', 'G'],
    ['C', 'G', 'A', 'T']
]

# Convert the list of lists into a NumPy array
motif_array = np.array(motif)

# Extract the third column (index 2) using NumPy slicing
print("Extracted Column (using NumPy):", motif_array[:, 2])

# Alternative method: Extract the third column using a loop
column = []  # Initialize an empty list
for row in motif_array:
    column.append(row[2])  # Append the element at index 2 in each row

# Print the manually extracted column
print("Extracted Column (using Loop):", column)

# Display the original motif array
print("Motif Array:\n", motif_array)

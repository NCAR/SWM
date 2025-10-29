import numpy as np
import sys

# Function to read a matrix from a text file
def read_matrix(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()
    matrix = []
    for line in lines:
        if line.strip():
            matrix.append([float(x) for x in line.strip().split()])
    return np.array(matrix)

# Read in two output txt files for comparison
if len(sys.argv) != 3:
    print("Usage: python comparison_results.py <file1> <file2>")
    sys.exit(1)

file1 = sys.argv[1]
file2 = sys.argv[2]

# Read matrices
matrix1 = read_matrix(file1)
matrix2 = read_matrix(file2)

debug = False
# Compare matrices
if matrix1.shape != matrix2.shape:
    print("Matrices have different dimensions.")
else:
    diff = matrix1 - matrix2
    if (debug):
        print("Difference matrix:")
        print(diff)
    print("Are matrices equal?", np.array_equal(matrix1, matrix2))
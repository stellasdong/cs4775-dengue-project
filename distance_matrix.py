from Bio import AlignIO
from Bio.Phylo.TreeConstruction import DistanceCalculator

# Load the aligned sequences in FASTA format
alignment = AlignIO.read("sequence.fasta", "fasta")

# Create a DistanceCalculator object (using p-distance or another model)
calculator = DistanceCalculator('identity')  # 'identity' is p-distance
distance_matrix = calculator.get_distance(alignment)

# Convert the distance matrix to a formatted string with tabs
def distance_matrix_to_string_with_tabs(distance_matrix):
    labels = distance_matrix.names
    result = []

    # Add the top row of labels, joined by tabs
    result.append("\t" + "\t".join(labels))

    # Add each row of the matrix
    for i, label in enumerate(labels):
        row = [f"{distance_matrix[label, other]:.6f}" for other in labels]
        result.append(f"{label}\t" + "\t".join(row))

    # Join the result list into a single string with newline characters
    return "\n".join(result)

# Generate the string representation of the distance matrix
matrix_string = distance_matrix_to_string_with_tabs(distance_matrix)

# Write the string to a text file
output_file = "distance_matrix.txt"
with open(output_file, 'w') as f:
    f.write(matrix_string)

print(f"Distance matrix written to {output_file}")

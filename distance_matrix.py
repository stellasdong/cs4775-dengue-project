from Bio import AlignIO
from Bio.Phylo.TreeConstruction import DistanceCalculator

alignment = AlignIO.read("sequence.fasta", "fasta")
calculator = DistanceCalculator('identity')
distance_matrix = calculator.get_distance(alignment)

def distance_matrix_create(distance_matrix):
    labels = distance_matrix.names
    result = []

    result.append("\t" + "\t".join(labels))

    for i, label in enumerate(labels):
        row = [f"{distance_matrix[label, other]:.6f}" for other in labels]
        result.append(f"{label}\t" + "\t".join(row))

    return "\n".join(result)

matrix_string = distance_matrix_create(distance_matrix)

output_file = "distance_matrix.txt"
with open(output_file, 'w') as f:
    f.write(matrix_string)

print(f"Distance matrix written to {output_file}")

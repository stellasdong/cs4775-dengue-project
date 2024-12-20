import random
from Bio import SeqIO
from Bio.Phylo import PhyloXML
from Bio.Phylo import BaseTree
import math
import itertools

# Read sequences from FASTA
def read_fasta(fasta_file):
    sequences = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequences[record.id] = str(record.seq)
    return sequences

# Generate initial tree using maximum parsimony
def generate_initial_tree(sequences):
    alignments = list(sequences.values())
    alignment = [''.join(seq) for seq in zip(*alignments)]
    
    constructor = ParsimonyTreeConstructor()
    tree = constructor.build_tree(alignment)
    return tree

# Placeholder function to calculate tree likelihood using JC69
def calculate_jukes_cantor_likelihood(tree, sequences):
    """
    Calculate the likelihood of a phylogenetic tree under the Jukes-Cantor model.
    """
    alpha = 1.0  # Substitution rate parameter for JC69
    likelihood = 0.0
    
    # Compute pairwise likelihoods for each pair of sequences
    for seq_id1, seq_id2 in itertools.combinations(sequences.keys(), 2):
        seq1 = sequences[seq_id1]
        seq2 = sequences[seq_id2]
        
        # Calculate the pairwise distance using JC69
        pairwise_distance = 0.0
        for nucleotide1, nucleotide2 in zip(seq1, seq2):
            if nucleotide1 != nucleotide2:
                pairwise_distance += 1
        
        # JC69 model: P = e^(-4 * alpha * d), where d is the pairwise distance
        # Calculate the likelihood contribution for this pair
        if pairwise_distance > 0:
            # Compute the probability using JC69 formula
            prob = (1/4) * (1 + 3 * math.exp(-4 * alpha * pairwise_distance))
            likelihood += math.log(prob)
        else:
            # If the sequences are identical, there's no change, likelihood is 1
            likelihood += 0
    
    # Return the overall likelihood of the tree
    return likelihood

# Apply an NNI move (simplified)
def apply_nni(tree):
    """
    Apply a simplified NNI move to the tree. This operation selects an internal node
    and swaps its two subtrees to generate a new topology.
    """
    # We randomly pick two internal nodes (clades) to swap
    internal_clades = [clade for clade in tree.get_nonterminals() if isinstance(clade, BaseTree.Clade)]

    # If there are less than two internal nodes, we can't perform NNI
    if len(internal_clades) < 2:
        return tree

    # Select two random internal clades
    clade1, clade2 = random.sample(internal_clades, 2)

    # Swap the subtrees
    tree.root.clades = [clade2, clade1]

    return tree

# Hill climbing NNI optimization
def hill_climbing_nni(tree, sequences):
    best_likelihood = calculate_jukes_cantor_likelihood(tree, sequences)
    improved = True
    while improved:
        improved = False
        nni_moves = []  # Generate possible NNI moves
        for move in nni_moves:
            apply_nni(tree)
            new_likelihood = calculate_jukes_cantor_likelihood(tree, sequences)
            if new_likelihood > best_likelihood:
                best_likelihood = new_likelihood
                improved = True
        if improved:
            tree = best_likelihood_tree
    return tree

# Stochastic NNI step to escape local optima
def stochastic_nni(tree, sequences):
    n = len(sequences)
    perturbations = int(0.5 * (n - 3))
    best_tree = tree
    best_likelihood = calculate_jukes_cantor_likelihood(tree, sequences)

    for _ in range(perturbations):
        apply_nni(tree)
        perturbed_tree = hill_climbing_nni(tree, sequences)
        new_likelihood = calculate_jukes_cantor_likelihood(perturbed_tree, sequences)

        if new_likelihood > best_likelihood:
            best_tree = perturbed_tree
            best_likelihood = new_likelihood

    return best_tree

# Save tree to Newick file
def save_tree_to_file(tree, filename):
    with open(filename, 'w') as f:
        Phylo.write(tree, f, 'newick')

# Main function to reconstruct tree
def reconstruct_tree(fasta_file, output_file):
    sequences = read_fasta(fasta_file)
    initial_tree = generate_initial_tree(sequences)
    optimized_tree = hill_climbing_nni(initial_tree, sequences)
    final_tree = stochastic_nni(optimized_tree, sequences)
    save_tree_to_file(final_tree, output_file)

# Run the tree reconstruction
reconstruct_tree("test.fasta", "output_tree.newick")

import numpy as np
from itertools import combinations
from ete3 import Tree

# Example: Define a simple substitution model (JC69)
def jc69_likelihood(tree, alignment, branch_lengths):
    """
    Compute the likelihood of the given tree under the JC69 model using
    Felsenstein's pruning algorithm.
    """
    num_sites = len(alignment[0])
    states = ['A', 'T', 'C', 'G']
    
    def jc_transition(i, j, branch_length):
        """Compute the transition probability under JC69."""
        if i == j:
            return 0.25 + 0.75 * np.exp(-4 * branch_length / 3)
        else:
            return 0.25 - 0.25 * np.exp(-4 * branch_length / 3)
    
    def conditional_prob(node, site):
        if node['is_leaf']:
            # node['id'] is an integer index here
            return [1.0 if state == alignment[node['id']][site] else 0.0 for state in states]
        else:
            left = node['left']
            right = node['right']
            left_probs = conditional_prob(left, site)
            right_probs = conditional_prob(right, site)
            probs = []
            for i, state in enumerate(states):
                p_left = sum(left_probs[j] * jc_transition(state, states[j], branch_lengths[left['id']]) for j in range(4))
                p_right = sum(right_probs[j] * jc_transition(state, states[j], branch_lengths[right['id']]) for j in range(4))
                probs.append(p_left * p_right)
            return probs

    likelihood = 1.0
    for site in range(num_sites):
        root_probs = conditional_prob(tree['root'], site)
        likelihood *= sum(root_probs) / 4  # Averaging over all possible root states
    
    return likelihood

def generate_random_tree(sequences):
    """Generate a random bifurcating tree."""
    n = len(sequences)
    tree = {i: {'id': i, 'is_leaf': True, 'left': None, 'right': None} for i in range(n)}
    
    next_node_id = n
    while len(tree) > 1:
        # Randomly pair two nodes to create a new internal node
        left, right = np.random.choice(list(tree.keys()), size=2, replace=False)
        new_node = {'id': next_node_id, 'is_leaf': False, 'left': tree.pop(left), 'right': tree.pop(right)}
        tree[next_node_id] = new_node
        next_node_id += 1
    
    return {'root': list(tree.values())[0]}

# Apply a random tree modification (NNI)
def modify_tree(tree):
    """Randomly modify the tree using NNI."""
    def swap_subtrees(node):
        if node['is_leaf']:
            return
        if np.random.rand() > 0.5:
            node['left'], node['right'] = node['right'], node['left']
        swap_subtrees(node['left'])
        swap_subtrees(node['right'])
    
    swap_subtrees(tree['root'])
    return tree

def replace_ids_with_names(tree, id_to_name_map):
    """Replace integer IDs in leaf nodes with sequence names."""
    def traverse(node):
        if node['is_leaf']:
            node['id'] = id_to_name_map[node['id']]
        else:
            traverse(node['left'])
            traverse(node['right'])
    traverse(tree['root'])
    return tree

def stochastic_search(alignment, max_iterations=1000):
    """Run the stochastic search algorithm."""
    tree = generate_random_tree(alignment)
    branch_lengths = {i: 0.1 for i in range(len(alignment) * 2 - 1)}  # Initialize branch lengths
    best_tree = tree
    best_likelihood = -np.inf
    
    for iteration in range(max_iterations):
        # Propose a new tree
        new_tree = modify_tree(tree)
        
        # Calculate likelihood
        likelihood = jc69_likelihood(new_tree, alignment, branch_lengths)
        
        # Accept or reject the new tree
        if likelihood > best_likelihood:
            best_tree = new_tree
            best_likelihood = likelihood
        else:
            # Optionally accept with a small probability
            if np.random.rand() < 0.1:
                best_tree = new_tree
                best_likelihood = likelihood

    return best_tree

def tree_to_newick(tree):
    def traverse(node):
        if node['is_leaf']:
            return f"{node['id']}"
        else:
            left = traverse(node['left'])
            right = traverse(node['right'])
            return f"({left},{right})"
    return traverse(tree['root']) + ";"

###################################
# Example usage
###################################
# Suppose you have a dictionary like this:
alignment_dict = {
    "seq1": "ATGCGATGCGTGGGAATAGGCAACAGAGACTTCGTTGAAGGACTGTCAGGAGCAACATGGGTGGATGTGGTACTGGAGCATGGAAGCTGCGTCACCACCATGGCAAAAAATAAACCAACATTGGACATTGAACTCTTGAAGACGGAGGTCACGAACCCTGCCGTCTTGCGCAAACTGTGCATTGAAGCTAAAATATCAAACACCACTACCGATTCAAGATGTCCAACACAAGGAGAAGCTACACTGGTGGAAGAACAAGACGCAAACTTTGTGTGTCGACGAACATTCGTGGACAGAGGCTGGGGTAATGGTTGTGGACTATTCGGGAAGGGAAGCTTACTAACGTGTGCTAAGTTTAAGTGTGTGACAAAACTTGAAGGAAAGATAGTTCAATATGAAAACTTAAAATATTCGGTGATAGTCACTGTCCACACTGGGGACCAGCACCAGGTAGGAAATGAAACTACAGAACATGGAACAATTGCAACCATAACACCTCAAGCTCCCACGTCGGAAATACAGCTGACTGACTACGGAGCCCTTACATTGGATTGCTCACCTAGAACAGGGCTGGACTTTAATGAGATGGTGCTGTTGACAATGAAAGAGAAATCATGGCTTGTCCACAAACAATGGTTTCTAGACTTACCATTGCCCTGGACCT",
    "seq2": "ATGCGATGCGTGGGAATAGGCAACAGAGACTTCGTTGAAGGACTGTCAGGAGCAACATGGGTGGATGTGGTACTGGAGCATGGAAGCTGCGTCACCACCATGGCAAAAAATAAACCAACATTGGACATTGAACTCTTGAAGACGGAGGTCACGAACCCTGCCGTCTTGCGCAAACTGTGCATTGAAGCTAAAATATCAAACACCACTACCGATTCAAGATGTCCAACACAAGGAGAAGCTACACTGGTGGAAGAACAAGACGCAAACTTTGTGTGTCGACGAACATTCGTGGACAGAGGCTGGGGTAATGGTTGTGGACTATTCGGGAAGGGAAGCTTACTAACGTGTGCTAAGTTTAAGTGTGTGACAAAACTTGAAGGAAAGATAGTTCAATATGAAAACTTAAAATATTCGGTGATAGTCACTGTCCACACTGGGGACCAGCACCAGGTAGGAAATGAAACTACAGAACATGGAACAATTGCAACCATAACACCTCAAGCTCCCACGTCGGAAATACAGCTGACTGACTACGGAGCCCTTACATTGGATTGCTCACCTAGAACAGGGCTGGACTTTAATGAGATGGTGCTGTTGACAATGAAAGAGAAATCATGGCTTGTCCACAAACAATGGTTTCTAGACTTACCATTGCCCTGGACCT",
    "seq3": "ATGCGATGCGTGGGAATAGGCAACAGAGACTTCGTTGAAGGACTGTCAGGAGCAACATGGGTGGATGTGGTACTGGAGCATGGAAGCTGCGTCACCACCATGGCAAAAAATAAACCAACATTGGACATTGAACTCTTGAAGACGGAGGTCACGAACCCTGCCGTCTTGCGCAAACTGTGCATTGAAGCTAAAATATCAAACACCACTACCGATTCAAGATGTCCAACACAAGGAGAAGCTACACTGGTGGAAGAACAAGACGCAAACTTTGTGTGTCGACGAACATTCGTGGACAGAGGCTGGGGTAATGGTTGTGGACTATTCGGGAAGGGAAGCTTACTAACGTGTGCTAAGTTTAAGTGTGTGACAAAACTTGAAGGAAAGATAGTTCAATATGAAAACTTAAAATATTCGGTGATAGTCACTGTCCACACTGGGGACCAGCACCAGGTAGGAAATGAAACTACAGAACATGGAACAATTGCAACCATAACACCTCAAGCTCCCACGTCGGAAATACAGCTGACTGACTACGGAGCCCTTACATTGGATTGCTCACCTAGAACAGGGCTGGACTTTAATGAGATGGTGCTGTTGACAATGAAAGAGAAATCATGGCTTGTCCACAAACAATGGTTTCTAGACTTACCATTGCCCTGGACCT"
}

# Create an alignment list from the values
alignment = list(alignment_dict.values())

# Create a map from integer ID to sequence name based on order
id_to_name_map = {i: name for i, name in enumerate(alignment_dict.keys())}

# Run the stochastic search
best_tree = stochastic_search(alignment)

# Replace IDs with names
best_tree = replace_ids_with_names(best_tree, id_to_name_map)

print("Best tree with sequence IDs:", best_tree)

# Convert the tree to Newick format
newick_tree = tree_to_newick(best_tree)

with open("tree.nwk", "w") as f:
    f.write(newick_tree)

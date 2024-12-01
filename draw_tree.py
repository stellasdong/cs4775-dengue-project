import make_tree
import argparse
from Bio import Phylo
import matplotlib.pyplot as plt


def main():
    parser = argparse.ArgumentParser(
    description='Neighbor-joining algorithm on a set of n sequences')
    parser.add_argument('-f', action="store", dest="f", type=str, default='distance_matrix.txt')
    parser.add_argument('-nwk', action="store", dest="nwk", type=str, default='tree.nwk')
    parser.add_argument('-o', action="store", dest="o", type=str, default='OM920075.1')
    args = parser.parse_args()
    distances_file = args.f
    og_ = args.o
    nwk_ = args.nwk

    D, mapping = make_tree.read_data(distances_file)
    og = dict(map(reversed, mapping.items()))[og_]

    E, uD, fake_root = make_tree.neighbor_join(D, og) # TIP: Please change the arguments to pass and return here if you don't want to follow the provided structure
    mapping[fake_root] = "root"
    tree_map = make_tree.assemble_tree(fake_root, E) # TIP: Please change the arguments to pass and return here if you don't want to follow the provided structure
    nwk_str = make_tree.generate_newick(fake_root, tree_map, uD, mapping) # TIP: Please change the arguments to pass and return here if you don't want to follow the provided structure
    
    tree = Phylo.read('tree.nwk', 'newick')

    # Root the tree using Green_monkey
    tree.root_with_outgroup("OM920075.1")

    # Set the branch lengths to display with 6 decimal places
    for clade in tree.find_clades():
        if clade.branch_length:
            clade.branch_length = round(clade.branch_length, 6)

    # Plot the tree with labeled tips and branch lengths
    Phylo.draw(tree, branch_labels=lambda c: '%.6f' % c.branch_length if c.branch_length else '')

# Plot the tree and save it as an image (e.g., PNG)
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(1, 1, 1)
    Phylo.draw(tree, do_show=False, branch_labels=lambda c: '%.6f' % c.branch_length if c.branch_length else '', axes=ax)
    plt.savefig('primate_tree.png')

main()
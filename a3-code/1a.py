#!/usr/bin/env python3

''' Generates a Newick formatted tree using the Neighbor-Joining algorithm, with an outgroup species specified.

Arguments:
    -f: The path to the distance matrix file (a symmetric matrix with 0 on the diagonal).
        (Default: dist10.txt)
    -nwk: The path to the output file to store the nwk format 
        (Default: tree.nwk)
    -o: Name of the outgroup species to root the tree.

Outputs:
    A Newick formatted tree.

Example usage:
    python 1a.py -f dist10.txt -o Green_monkey -nwk tree.nwk


***************** IMPORTANT:   ***************************************************************************************
********    If any, please remove any print statements except the one provided in the main function before    ********
********    submission, as the autograder evaluation based solely on the print output of this script.         ********
**********************************************************************************************************************

'''

import argparse
from copy import deepcopy


''' Reads the input distance file between species, encode the distance matrix into a dictionary
    with the species names encoded as integer indices.

Arguments:
    distances_file: Path to the file containing the distance matrix between species.
Returns:
    D: A dictionary of dictionaries, defining distances between all species, every key is a species index,
        and the corresponding value is a dictionary containing all species indexes as keys. The values of
        these keys are the distance between species. For example {1: {1: 0.0, 2: 1.0}, 2: {1: 1.0, 2: 0.0}}
        defines the distance between two species, 1 and 2.
    mapping: A dictionary mapping indices to species names. For example, {1: 'Chimp', 2: 'Bonobo'}.
'''
def read_data(distances_file):
    with open(distances_file, "r") as f:
        lines = [l.strip().split() for l in f.readlines()]
        mapping = {i: s for i, s in enumerate(lines[0])}
        lines = [l[1:] for l in lines[1:]]
        D = {i: {} for i in range(len(lines))}
        for i, l in enumerate(lines):
            for j, sval in enumerate(l):
                D[i][j] = float(sval)
    return D, mapping


''' Performs the neighbor joining algorithm on the distance matrix and the index of the outgroup species.

Arguments:
    D: A dictionary of dictionaries, defining distances between all species, every key is a species index,
        and the corresponding value is a dictionary containing all species indexes as keys. (See the description
        in `read_data` for details).
    og: outgroup index, defining which species serves as an outgroup.
        A fake root should be inserted in the middle of the pendant edge
        leading to this outgroup node.

Returns:
    E : A list storing the edges chosen from the NJ algorithm in the form of tuples: (index, index). 
        For example [(3,1),(3,2)] represents an unrooted NJ tree of two edges, 
        3<-->1 and 3<-->2, where 1 & 2 are indexes of leaf nodes in the tree,
        and 3 is the index of the internal node you added.
    uD: A dictionary of dictionary, defining distances between all nodes (leaves and internal nodes),
        it's of the same format as D, storing all edge lengths of the NJ tree whose topology is specified by E.
        For example, {1: {1: 0.0, 2: 1.0, 3: 1.5}, 2: {1: 1.0, 2: 0.0, 3: 2.0}, 3: {1: 1.5, 2: 2.0, 3: 0.0}}
        will fully specify the edge lengths for the tree represented by the E example ([(3,1),(3,2)]):
        Length(3<-->1) = 1.5, Length(3<-->2) = 2.0.
    fake_root: which node index represent the root

        *************************************************************************************
        ***** Please note that this return value format is just for guiding purpose,    *****
        ***** and we are not grading based on the output of this function (neighbor_join).***
        ***** Please feel free to use return structure / indexing system you like.      *****
        *************************************************************************************
'''
def neighbor_join(D, og):
    ''' Complete this function. '''
    uD = deepcopy(D)
    temp_D = deepcopy(D)
    E = []
    removed = []

    for z in range(len(D)+1, len(D)+(len(D)-1)):
        n = len(temp_D)
        d_sum = column_sums(temp_D)
        if n > 3:
            x, y = find_Q(temp_D, d_sum)
            E += [(z, x), (z, y)]
            removed += [x, y]
            new_row = {}
            new_row[x] = (0.5*uD[x][y])+(d_sum[x]-d_sum[y])/(2*(n-2))
            new_row[y] = (0.5*uD[x][y])+(d_sum[y]-d_sum[x])/(2*(n-2))
            # print(removed)
            for k in uD:
                if k not in removed:
                    # print(x, y, k)
                    # print(uD)

                    new_row[k] = 0.5*(uD[x][k]+uD[y][k]-uD[x][y])
            
            add_row(uD, new_row, z)
            add_row(temp_D, new_row, z)
            # print(x, y)
            # print(uD)
            temp_D = remove_rows(temp_D, x, y)
        if n == 3:
            # print("run")
            x = list(temp_D.keys())[0]
            y = list(temp_D.keys())[1]
            k = list(temp_D.keys())[2]
            # print(x, y, k)
            # print(z)
            E += [(x, z), (y, z), (k, z)]
            new_row = {}
            new_row[x] = (0.5*uD[x][y])+(d_sum[x]-d_sum[y])/(2*(n-2))
            new_row[y] = uD[x][y]-new_row[x]
            new_row[k] = 0.5*(uD[x][k]+uD[y][k]-uD[x][y])
            # print(new_row)
            add_row(uD, new_row, z)
            # print(uD)
            add_row(temp_D, new_row, z)
    fake_root = len(uD)+1
    node = 0
    dist = 0
    # print(E)
    # print(uD)
    for x, y in E:
        if y == og:
            node = x
            dist = uD[x][y]
            E.remove((x, y))
        elif x == og:
            node = y
            dist = uD[x][y]
            E.remove((x, y))

    E += [(fake_root, og), (node, fake_root)]
    uD[fake_root] = {og: dist/2.0, node: dist/2.0}
        
    return E, uD, fake_root

def remove_rows(uD, x, y):
    temp_D = deepcopy(uD)
    del temp_D[x]
    del temp_D[y]
    for row in temp_D:
        if x in temp_D[row]:
            del temp_D[row][x]
        if y in temp_D[row]:
            del temp_D[row][y]
    return temp_D

def add_row(uD, new_row, z):
    uD[z] = new_row
    for row in new_row:
        if row != z:
            uD[row][z] = new_row[row]
    uD[z][z] = 0
    # print(uD)

def column_sums(D):
    d_sum = {}
    for row in D:
        for col in D[row]:
            if col not in d_sum:
                d_sum[col] = 0
            d_sum[col] += D[row][col]
    return d_sum

def find_Q(D, d_sum):
    Q = {}
    m = None
    m_index = (0, 0)
    for row in D:
        dict = {}
        for col in D[row]:
            if row > col:
                dict[col] = (len(D) - 2)*D[row][col] - d_sum[row] - d_sum[col]
                if m is None or dict[col] < m:
                    m = dict[col]
                    m_index = (col, row)
        if row != 1:
            Q[row] = dict
        if len(D) == 3:
            m_index = (list(D.keys())[0], list(D.keys())[2])
    return m_index[0], m_index[1]


''' Helper function for defining a tree data structure.
    First finds the root node and find its children, and then generates 
    the whole binary tree based on .

Arguments:
    E：A list storing the edges chosen from the NJ algorithm in the form of tuples: (index, index). 
         (See the description in `neighbor_join` for details).
    fake_root: which node index represent the root
Returns:
    tree_map：A dictionary storing the topology of the tree, where each key is a node index and each value is a list of
              node indexs of the direct children nodes of the key, of at most length 2. For example, {3: [1, 2], 2: [], 1: []}
              represents the tree with 1 internal node (3) and two leaves (1, 2). the [] value shows the leaf node status.

    *************************************************************************************
    ***** Please note that this return value format is just for guiding purpose,    *****
    ***** and we are not grading based on the output of this function (assemble_tree).***
    ***** Please feel free to use return structure / indexing system you like.      *****
    *************************************************************************************
'''
def assemble_tree(fake_root, E):
    ''' Complete this function. '''
    tree_map = {fake_root: []}
    build_tree(fake_root, E, tree_map)

    return tree_map

def build_tree(parent, E, tree_map):
    for node1, node2 in E[:]:

        if node1 == parent:
            tree_map[parent] += [node2]
            tree_map[node2] = []
            if (node1, node2) in E:
                E.remove((node1, node2))
                build_tree(node2, E, tree_map)
        if node2 == parent:
            tree_map[parent] += [node1]
            tree_map[node1] = []
            if (node1, node2) in E:
                E.remove((node1, node2)) 
                build_tree(node1, E, tree_map)


''' Returns a string of the Newick tree format for the tree, rooted at a pre-defined node (fake_root).

Arguments:
    fake_root: which node index represent the root
    tree_map：A dictionary storing the topology of the tree (See the description in `assemble_tree` for details).
    uD: A dictionary of dictionary, defining distances between all nodes (leaves and internal nodes)
        (See the description in `neighbor_join` for details)
    mapping: A dictionary mapping indices to species names. (See the description in `read_data` for details)
Returns:
    output: rooted tree in Newick tree format (string). The branch lengths should be in 6 decimal digits.
            For example, you could do "string = '%s:%.6f' % (name, length)"

    *********************************************************************************************
    ***** We will grade on the newick string output by this function (generate_newick).     *****
    *********************************************************************************************
'''
def generate_newick(fake_root, tree_map, uD, mapping = None):
    # print(tree_map)
    # print(mapping)
    # print(uD)
    ''' Complete this function. '''
    def display(fake_root):
        if len(tree_map[fake_root]) == 0:
            if mapping == None: return str(fake_root)
            else: return mapping[fake_root]
        
        [left_child, right_child] = tree_map[fake_root]
        return '(%s:%.6f, %s:%.6f)' % (display(left_child), uD[fake_root][left_child],
                                       display(right_child), uD[fake_root][right_child])
    
    return display(fake_root) + ';'

# def generate_newick(fake_root, tree_map, uD, mapping=None):
#     print(fake_root)
#     print(tree_map)
#     print(mapping)
#     print(uD)
#     ''' Complete this function. '''
#     visited = set()  
#     print(tree_map)


#     def display(node):
#         visited.add(node)
#         if len(tree_map[node]) == 0:
#             if mapping is None: return str(node)
#             else: return mapping.get(node, str(node))
#         subtrees = []
#         for child in tree_map[node]:
#             print(node)
#             print(child)
#             subtrees.append('%s:%.6f' % (display(child), uD[node][child]))
#         return '(%s)' % ', '.join(subtrees)

#     newick_string = display(fake_root)
#     for node in tree_map:
#         if node not in visited:
#             newick_string += ', ' + display(node)
    
#     return newick_string + ';'


def main():
    parser = argparse.ArgumentParser(
        description='Neighbor-joining algorithm on a set of n sequences')
    parser.add_argument('-f', action="store", dest="f", type=str, default='dist10.txt')
    parser.add_argument('-nwk', action="store", dest="nwk", type=str, default='tree.nwk')
    parser.add_argument('-o', action="store", dest="o", type=str, default='Green_monkey')
    args = parser.parse_args()
    distances_file = args.f
    og_ = args.o
    nwk_ = args.nwk

    D, mapping = read_data(distances_file)
    og = dict(map(reversed, mapping.items()))[og_]

    E, uD, fake_root = neighbor_join(D, og) # TIP: Please change the arguments to pass and return here if you don't want to follow the provided structure
    mapping[fake_root] = "root"
    tree_map = assemble_tree(fake_root, E) # TIP: Please change the arguments to pass and return here if you don't want to follow the provided structure
    nwk_str = generate_newick(fake_root, tree_map, uD, mapping) # TIP: Please change the arguments to pass and return here if you don't want to follow the provided structure
    
    # Print and save the Newick string.
    ''' 
    ****************************************************************************************
    ***** Please note that we will grade on this print statement, so please make sure  *****
    ***** to delete any print statement except for this given one before submission!!! *****
    ****************************************************************************************
    '''
    print(nwk_str)
    with open(nwk_, "w") as nwk_file:
        print(nwk_str, file=nwk_file)


if __name__ == "__main__":
    main()

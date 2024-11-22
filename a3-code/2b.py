#!/usr/bin/env python3

''' Find the maximum a posteriori distribution for each base based on the
    Jukes-Cantor model using Felsenstein's algorithm.

Arguments:
    -f: fasta file containing the multiple alignment
        (default is apoe.fasta)
    -t: topology index (1-3), corresponding to the maximum likelihood topology
        (from 2a)
Outputs:
    alignment of input sequences and maximum likelihood root sequence for the
    ML topology from Figure 1

Example usage:
    python 2b.py -f test.fasta -t 2
'''

import argparse
import numpy as np
import math


class Node():
    ''' Initializes a node with given parameters.

    Arguments:
        name: name of node (only relevant for leaves)
        left: left child (Node)
        right: right child (Node)
        branch_length: length of branch that leads to this node (float)
        branch_id: id of branch that leads to this node (int)
        probs: probability of observed bases beneath this node
                [list of 4 probs for 'ACGT'] (initialized to None]
    '''
    def __init__(self, name, left, right, branch_length, branch_id):
        self.name = name
        self.left = left
        self.right = right
        self.branch_length = branch_length
        self.branch_id = branch_id
        self.probs = [None for _ in range(4)]


''' Reads data from ```filename``` in fasta format.

Arguments:
    filename: name of fasta file to read
Returns:
    sequences: dictionary of outputs (string (sequence id) -> sequence (string))
    size: length of each sequence
'''
def read_data(filename):
    with open(filename, "r") as f:
        lines = f.readlines()
        sequences = {}
        output = ""
        size = 0
        curr = ""
        flag = False
        for l in lines:
            if l[0] == ">": 
                if (len(output) != 0):
                    sequences[curr] = output
                    size = len(output)
                    output = ""
                curr = l[2:].strip()
            else:
                output += l.strip()
        sequences[curr] = output
    return sequences, size


''' Evaluates P(b|a, t) under the Jukes-Cantor model

Arguments:
    b: descendant base (string)
    a: ancestral base (string)
    t: branch length (float)
    u: mutation rate (float, defaults to 1)
Returns:
    prob: float probability P(b|a, t)
'''
def jcm(b, a, t, u = 1.0):
    ''' Complete this function. '''
    if b == a:
        return 1/4 + (3/4 * math.exp(-4 * u * t / 3))
    else:
        return 1/4 - (1/4 * math.exp(-4 * u * t / 3))    
    pass


''' Constructs the ordering of the post-order traversal of ```index```
    topology from the pset.
Arguments:
    index: which topology to use
Returns:
    list of Nodes corresponding to post-order traversal of the topology
    branch_probs: 6x4x4 matrices, indexed as:
                  branch_probs[branch_id][a][b] = P(b | a, t_branch_id)
'''
def initialize_topology(index):
    bases = 'ACGT'

    branch_lengths = np.array(
        [[0.07517, 0.03059, 0.03161, 0.11761, 0.14289],
        [0.20843, 0.03397, 0.03497, 0.24952, 0.00000],
        [0.20843, 0.03397, 0.03497, 0.24952, 0.00000]], dtype = float)

    names = ['human', 'mouse', 'rat', 'dog']
    branches = [0, 1, 2, 3]
    leaves = [Node(s, None, None, bl, i) for (s, i, bl) in 
                zip(names, branches, branch_lengths[index, :])]
    ordering = None
    branch_probs = [np.zeros((4,4), dtype = float) for _ in range(6)]
    # Note that branch 5 (or 6 in 1-index) is the branch of 0-length
    if (index == 0):
        hum_dog = Node(None, leaves[0], leaves[3], 0, 5)
        mouse_rat = Node(None, leaves[1], leaves[2], branch_lengths[index,4], 4)
        root = Node('root', hum_dog, mouse_rat, None, None)
        ordering = [leaves[0], leaves[3], hum_dog, leaves[1], leaves[2], \
                    mouse_rat, root]
    elif (index == 1):
        hum_mouse = Node(None, leaves[0], leaves[1], 0, 5)
        rat_dog = Node(None, leaves[2], leaves[3], branch_lengths[index, 4], 4)
        root = Node('root', hum_mouse, rat_dog, None, None)
        ordering = [leaves[0], leaves[1], hum_mouse, leaves[2], leaves[3], \
                    rat_dog, root]
    else:
        mouse_dog = Node(None, leaves[1], leaves[3], 0, 5)
        hum_rat = Node(None, leaves[0], leaves[2], branch_lengths[index, 4], 4)
        root = Node('root', mouse_dog, hum_rat, None, None)
        ordering = [leaves[1], leaves[3], mouse_dog, leaves[0], leaves[2], \
                    hum_rat, root]

    ''' Assign 6x4x4 branch_probs values: branch_probs[branch_id][ancestor_base][descendant_base] '''
    for id in range(6):
        for a in range(4):
            for b in range(4):
                if id == 5:
                    length = 0
                else:
                    length = branch_lengths[index][id]
                branch_probs[id][a][b] = jcm(b, a, length)    
    return ordering, branch_probs


''' Computes maximum posterior distribution of bases at the root of the tree
    given the topology specified by ordering

Arguments:
    data: sequence data (dict: name of sequence owner -> sequence)
    seqlen: length of sequences
    ordering: postorder traversal of our topology
    bp: branch probabilities for the given branches: 6x4x4 matrix indexed as
        branch_probs[branch_id][a][b] = P(b | a, t_branch_id)
Returns:
    output: maximum a posteriori (MAP) estimate of root sequence
'''
def map_estimate(data, seqlen, ordering, bp):
    ''' Complete this function. '''
    bases = ['A', 'C', 'G', 'T']
    for node in ordering:
        probs = [[0 for _ in range(4)] for _ in range(seqlen)]
        if node.left is None and node.right is None:
            seq = data[node.name]
            for s in range(seqlen):
                for i in range(4):
                    if bases[i] == seq[s]:
                        probs[s][i]= 1.0
                    else:
                        probs[s][i]= 0.0
        else:
            for s in range(seqlen):
                for i in range(4):
                    left = 0
                    right = 0
                    for j in range(4):
                        left += (bp[node.left.branch_id][i][j]*node.left.probs[s][j])
                        right += (bp[node.right.branch_id][i][j]*node.right.probs[s][j])
                    probs[s][i] = left * right        
        node.probs = probs

    root = ordering[-1] 
    total_log_prob = 0.0
    posterior_probs = [[0 for _ in range(4)] for _ in range(seqlen)]        

    for s in range(seqlen):
        root_prob = 0
        for i in range(4):
            root_prob += 0.25 * root.probs[s][i]        
        total_log_prob += np.log(root_prob)
        
        for i in range(4):
            posterior_probs[s][i] = (0.25 * root.probs[s][i]) / root_prob

    root_sequence = ''
    for s in range(seqlen):
        root_sequence += bases[np.argmax(posterior_probs[s])]
    return root_sequence
    

        
''' Outputs the MAP estimate and the data in a way that can easily be examined
    for comparison.

Arguments:
    data: dictionary of sequences (name -> sequence)
    map_output: a string containing maximum a posteriori sequence at root
'''
def output_alignment(data, map_output, per_line = 70):
    iter = len(map_output) // per_line
    names = ['human: ', 'mouse: ', 'rat:   ', 'dog:   ']
    if (len(map_output) % per_line != 0): iter += 1

    for i in range(iter):
        if (i == iter - 1):
            print('root:  ' + map_output[i * per_line:])
            for name in names: print(name + data[name.strip()[:-1]][i * per_line:])
        else:
            print('root:  ' + map_output[i * per_line: (i+1) * per_line])
            for name in names: print(name + data[name.strip()[:-1]][i * per_line: (i+1) * per_line])
            print('\n')


def main():
    parser = argparse.ArgumentParser(
        description='Reconstruct MAP sequence and output this aligned with human, mouse, rat, and dog sequences.')
    parser.add_argument('-f', action="store", dest="f", type=str, default='apoe.fasta')
    parser.add_argument('-t', action="store", dest="t", type=int, required=True)
    args = parser.parse_args()
    fasta_file = args.f
    topology_index = args.t

    assert topology_index in [0, 1, 2], "Invalid topology index."
    
    data, seqlen = read_data(fasta_file)
    ordering, probs = initialize_topology(topology_index)
    map_output = map_estimate(data, seqlen, ordering, probs)
    output_alignment(data, map_output)


if __name__ == "__main__":
    main()

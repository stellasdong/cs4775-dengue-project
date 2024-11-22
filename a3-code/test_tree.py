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
            E.remove((node1, node2))
            build_tree(node2, E, tree_map)
        if node2 == parent:
            tree_map[parent] += [node1]
            tree_map[node1] = []
            E.remove((node1, node2)) 
            build_tree(node1, E, tree_map)            


E = [(9, 1), (9, 2), (10, 5), (10, 6), (11, 7), (12, 10), (12, 11), (13, 3), (13, 9), (4, 14), (12, 14), (13, 14), (15, 8), (11, 15)]
print(assemble_tree(15, E))
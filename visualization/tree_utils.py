from Bio import Phylo


# It would be more efficient to stick this inside the "set_node_x_values"
# method of the gene tree visuals module. But I'm leaving it here
# for now for modularity.
def print_node_distances(tree, out_file):

    X = Phylo.to_networkx(tree)
    nodes = list(X.nodes)
    col_headers=["node_name","dist","leaves"]
    data_lines = []

    for node in nodes:
            name="no_name"
            if node.name:
                name = node.name
            leaves = node.get_terminals()
            leaf_names = " ".join([leaf.name for leaf in leaves])
            dist_string=str(tree.distance(node))
            data=[name ,dist_string,leaf_names]
            data_lines.append(",".join(data) + "\n")

    with open(out_file, 'w') as f:
        f.writelines(",".join(col_headers) +"\n")
        f.writelines(data_lines)



#TODO, put a test case around this
def get_parent(tree, child_clade):
        node_path = tree.get_path(child_clade)
        if len(node_path)>2:
            return node_path[-2]
        else:
            return False
def get_sibs(tree, child_clade):
       parent = get_parent(tree, child_clade)
       if not parent:
           return []
       sibs=parent.clades
       return sibs
def birth_order(tree, child_clade):
    sibs = get_sibs(tree, child_clade)
    num_sibs=len(sibs)
    if num_sibs < 2:
        return (0,1)
    birth_order=sibs.index(child_clade)
    return (birth_order,num_sibs)
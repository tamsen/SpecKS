
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
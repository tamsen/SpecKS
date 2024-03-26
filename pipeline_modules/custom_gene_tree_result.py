from io import StringIO
from Bio import Phylo

class custom_gene_tree_result():

    original_newick=""
    simple_newick=""
    terminal_leaf_names=0
    gene_tree_name=""
    #gene_tree_file_name=""
    #leaf_map_file_name=""
    leaves_by_species={}
    tree=False

    def __init__(self, gene_tree_name,original_newick):

        self.original_newick=original_newick
        self.simple_newick=original_newick
        self.tree= Phylo.read(StringIO(self.simple_newick), "newick")
        self.terminal_leaf_names = [t.name for t in self.tree.get_terminals()]
        self.num_terminal_leaves  = len(self.terminal_leaf_names)
        self.gene_tree_name = gene_tree_name

        leaves_by_species={}
        genomes=["O","P1","P2"]
        for g in genomes:
            leaves_by_species[g]=[]
            for leaf in self.terminal_leaf_names:
                if g in leaf:
                    leaves_by_species[g].append(leaf)
        self.leaves_by_species=leaves_by_species

    def get_tree_length_as_in_PAML(self):
        # As per the PAML manual, “The tree length is the expected number of substitutions
        # per site along all branches in the phylogeny, calculated as the sum of the branch lengths”.

        total_branch_length = self.tree.total_branch_length()
        return total_branch_length

    def get_named_modes(self):
        X1 = Phylo.to_networkx(self.tree)
        nodes = list(X1.nodes)
        named_nodes=[]
        for n in nodes:
            if n.name:
                named_nodes.append(n.name)
        return named_nodes
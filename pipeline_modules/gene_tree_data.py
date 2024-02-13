import os
from Bio import Phylo
from io import StringIO
class gene_tree_result():

    original_newick=""
    simple_newick=""
    terminal_leaf_names=0
    gene_tree_name=""
    gene_tree_file_name=""
    leaf_map_file_name=""
    leaves_by_species={}
    tree=False

    def __init__(self, gene_tree_name,gene_tree_file_name, leaf_map_file):

        self.original_newick=get_newick_from_tree_file(gene_tree_file_name)
        self.simple_newick=clean_newick(self.original_newick)
        self.tree= Phylo.read(StringIO(self.original_newick), "newick")
        self.terminal_leaf_names = [t.name for t in self.tree.get_terminals()]
        self.num_terminal_leaves  = len(self.terminal_leaf_names)
        self.gene_tree_name = gene_tree_name
        self.gene_tree_file_name = gene_tree_file_name
        self.leaf_map_file_name=leaf_map_file

        self.set_leafmap_data(leaf_map_file)

    def add_back_outgroup(self, leg_distance):
        original_newick = self.simple_newick
        terminal_leaf_names = self.terminal_leaf_names

        if (not "O" in terminal_leaf_names) and (not "O1" in terminal_leaf_names):
            out_group_leaf="O1"
            origin_node_name="O"
            new_newick, new_tree,new_leaves = unprune_outgroup(self.tree,
                                                               leg_distance, out_group_leaf, origin_node_name)

            self.simple_newick = new_newick
            self.original_newick = original_newick
            self.tree= new_tree

            for new_leaf in new_leaves:
                if not "O" in self.leaves_by_species:
                    self.leaves_by_species["O"]=[]
                self.leaves_by_species["O"].append(new_leaf)
                self.terminal_leaf_names.append(new_leaf)

        self.num_terminal_leaves = len(self.terminal_leaf_names)
        return self

    def update_tree(self, new_tree):

        handle = StringIO()
        Phylo.write(new_tree, handle , "newick")
        new_newick = handle.getvalue()

        self.tree=new_tree
        self.original_newick = self.simple_newick
        self.simple_newick=clean_newick(new_newick)
        self.terminal_leaf_names = [t.name for t in self.tree.get_terminals()]
        self.num_terminal_leaves  = len(self.terminal_leaf_names)

        return self
    def set_leafmap_data(self, leafmap_file):

        leaves_by_species,num_extant_leaves= read_leaf_map_data(leafmap_file)
        self.leaves_by_species=leaves_by_species
        self.num_extant_leaves=num_extant_leaves
        return self
def unprune_outgroup(old_tree, leg_distance, out_group_leaf, origin_node_name):

        outgroup_clade = Phylo.BaseTree.Clade(branch_length=leg_distance, name=out_group_leaf)
        new_clade = Phylo.BaseTree.Clade(branch_length=leg_distance, name=origin_node_name,
                                         clades=[outgroup_clade, old_tree.clade])
        new_tree = Phylo.BaseTree.Tree.from_clade(new_clade)
        handle = StringIO()
        Phylo.write(new_tree, handle, "newick")
        new_newick = handle.getvalue()
        return new_newick, new_tree, [out_group_leaf]
def clean_newick(taggy_newick):
        open_brackets = [i for i in range(0, len(taggy_newick)) if taggy_newick[i] == "["]
        if len(open_brackets) == 0:
            return taggy_newick

        close_brackets_plus_zero = [0] + [i + 1 for i in range(0, len(taggy_newick)) if taggy_newick[i] == "]"]
        stuff_to_keep = [taggy_newick[close_brackets_plus_zero[i]:open_brackets[i]] for i in range(0, len(open_brackets))]
        clean_tree = "".join(stuff_to_keep) + ";"
        return clean_tree
def read_gene_tree_result_from_tree_and_leaf_map_files(tree_file, leaf_map_file):

    gene_tree_name = os.path.basename(tree_file).replace(".pruned.tree", "")
    result = gene_tree_result(gene_tree_name,tree_file,leaf_map_file)
    return result

def get_newick_from_tree_file(tree_file):

    with open(tree_file, "r") as f:
        lines = f.readlines()
        newick=lines[0]

    return newick


def read_leaf_map_data(leafmap_file):

    leaves_by_species={}
    num_extant_leaves=0
    with open(leafmap_file, "r") as f:
            lines = f.readlines()
            for line in lines:
                splat=line.split()
                leaf=splat[0]
                species=splat[1]
                if species not in leaves_by_species.keys():
                    leaves_by_species[species]=[]
                leaves_by_species[species].append(leaf)
                num_extant_leaves=num_extant_leaves+1

    return leaves_by_species,num_extant_leaves

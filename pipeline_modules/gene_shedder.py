import os
from io import StringIO

import numpy as np
from Bio import Phylo


def shed_gene_trees(polyploid, relaxed_gene_tree_results):

    config = polyploid.general_sim_config
    subfolder=os.path.join(polyploid.species_subfolder, str(polyploid.analysis_step_num) + "_post_WGD_gene_pruning")
    if not os.path.exists(subfolder):
        os.makedirs(subfolder)

    post_WGD_shedding_rate_per_MY = 0.10
    for gene_tree, gene_tree_results in relaxed_gene_tree_results.items():
        random_seed=random_seed+1
        tree = Phylo.read(StringIO(gene_tree_results.simple_newick), "newick")

    print("TODO: after WGD event, remove x% of genes per MY")
    polyploid.analysis_step_num=polyploid.analysis_step_num+1


def build_node_matrix(tree, max_sim_time):
    X = Phylo.to_networkx(tree)
    nodes = list(X.nodes)
    edges = list(X.edges)
    col_headers = ["node_name", "dist", "leaves"]
    data_lines = []
    bin_size = 1 #1 MY
    bins_post_WGD=  np.arange(0, max_sim_time, bin_size)


    #for node in nodes, figure out if the node is still extant
    #at time slice t
    for i in range(0,len(nodes)):

        node=nodes[i]
        name = "no_name"
        if node.name:
            name = node.name
        leaves = node.get_terminals()
        leaf_names = " ".join([leaf.name for leaf in leaves])
        dist = tree.distance(node)


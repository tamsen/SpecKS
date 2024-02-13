import os
from random import sample
from io import StringIO
import numpy as np
from Bio import Phylo

from visualization import tree_visuals_by_phylo


def shed_genes(polyploid, relaxed_gene_tree_results):


    if len(polyploid.subtree_subfolder) > 0:
        subfolder = os.path.join(polyploid.species_subfolder,
                                 str(polyploid.analysis_step_num) + "_gene_shedder_" + polyploid.subtree_subfolder)
    else:
        subfolder = os.path.join(polyploid.species_subfolder, str(polyploid.analysis_step_num) + "_gene_shedder")

    print(subfolder)
    if not os.path.exists(subfolder):
        os.makedirs(subfolder)

    max_tree_distance = polyploid.FULL_time_MYA
    tree_distance_defore_WGD = polyploid.FULL_time_MYA - polyploid.WGD_time_MYA
    bin_size = 1
    time_slices = np.arange(tree_distance_defore_WGD, max_tree_distance, bin_size)
    gene_trees_after_gene_shedding_by_gt=relaxed_gene_tree_results.copy()

    # At each time-slice post WGD, the number of extant genes are counted.
    # The total number of genes to be shed is calculated as a percent of the
    # total on the time slice. Each gene will come from a different gene tree.

    #TODO, this needs to be calculated based on some parameters.
    num_genes_to_shed = 2
    all_gene_trees=list(relaxed_gene_tree_results.keys())
    gene_trees_to_loose_a_gene = sample(all_gene_trees, num_genes_to_shed)

    # For each gene to be shed, a GT is randomly selected without replacement,
    # and a vertex from that GT which crosses the time-slice is removed.
    #   This process is repeated until the desired number of genes are shed.
    unprunable_leaves = ['O', 'O1', 'O2']
    time_slice = time_slices[-10]
    num_genes_to_remove_per_gene_tree = 1
    print("time slice:\t" + str(time_slice))

    for gt_name in gene_trees_to_loose_a_gene:

            gt_data_to_edit= gene_trees_after_gene_shedding_by_gt[gt_name]
            gt_tree = gt_data_to_edit.tree
            nodes_on_edges_that_cross_this_time = (
                get_nodes_on_edges_that_cross_this_time(gt_tree, time_slice))

            #just to visualize what is going on beforehand
            out_file=os.path.join(subfolder ,gt_name + "_newick_pre_gene_shedding.txt")
            Phylo.write(gt_tree, out_file, "newick")
            tree_visuals_by_phylo.save_tree_plot(gt_data_to_edit.simple_newick,
                                                 out_file.replace("txt","png"))

            list_of_terminal_leaves_to_remove = chose_leaves_to_remove(
                nodes_on_edges_that_cross_this_time, num_genes_to_remove_per_gene_tree, unprunable_leaves)

            for leaf in list_of_terminal_leaves_to_remove:
                gt_tree.prune(leaf)

            gene_trees_after_gene_shedding_by_gt[gt_name].update_tree(gt_tree)

            # just to visualize what is going on after shedding..
            out_file = os.path.join(subfolder, gt_name + "_newick_post_gene_shedding.txt")
            Phylo.write(gt_tree, out_file, "newick")
            tree_visuals_by_phylo.save_tree_plot(gt_data_to_edit.simple_newick,
                                                 out_file.replace("txt", "png"))

    polyploid.analysis_step_num=polyploid.analysis_step_num+1
    return gene_trees_after_gene_shedding_by_gt



def chose_leaves_to_remove(nodes_on_edges_that_cross_this_time,num_to_remove, unprunable_leaves):

        leaves_to_remove_by_node = {}
        nodes_allowed_to_prune = []
        for node in nodes_on_edges_that_cross_this_time:
            node_OK_to_prune = True
            leaves = node.get_terminals()
            for leaf in leaves:
                if leaf.name in unprunable_leaves:
                    node_OK_to_prune = False
                    continue
            if node_OK_to_prune:
                nodes_allowed_to_prune.append(node)
                leaves_to_remove_by_node[node] = leaves.copy()

        #Normally, you'd expect num_to_remove to be < nodes_allowed_to_prune
        #But *just in case* num_to_remove >= nodes_allowed_to_prune,
        #in that case, we just have to remove *all* the nodes.

        if (len(nodes_allowed_to_prune) >= num_to_remove):
            nodes_to_remove_list = sample(nodes_allowed_to_prune, num_to_remove)
        else: # just have to remove everything!
            nodes_to_remove_list = nodes_allowed_to_prune

        leaves_to_remove_list = []
        for node in nodes_to_remove_list:
            leaves_to_remove_list = leaves_to_remove_list + leaves_to_remove_by_node[node]
        return leaves_to_remove_list


def get_distance_intervals_for_edge(edges, tree):
    distance_intervals_by_edge_idx = {}
    nodes_by_edge_idx = {}
    for e_idx in range(0, len(edges)):
        edge = edges[e_idx]
        dist1 = tree.distance(edge[0])
        dist2 = tree.distance(edge[1])
        if dist1 < dist2:
            interval = (dist1, dist2)
            nodes = (edge[0], edge[1])
        else:
            interval = (dist2, dist1)
            nodes = (edge[1], edge[0])
        distance_intervals_by_edge_idx[e_idx] = interval
        nodes_by_edge_idx[e_idx] = nodes

    return distance_intervals_by_edge_idx, nodes_by_edge_idx


def get_nodes_on_edges_that_cross_this_time(tree,
                                            time_slice):
    X = Phylo.to_networkx(tree)
    edges = list(X.edges)
    distance_intervals_by_edge_idx, nodes_by_edge_idx = get_distance_intervals_for_edge(edges, tree)
    nodes_on_edges_that_cross_this_time = []
    for e_idx in range(0, len(edges)):
        distance_interval = distance_intervals_by_edge_idx[e_idx]

        # does edge cross the slice of time?
        if (distance_interval[0]) < time_slice and (distance_interval[1] >= time_slice):
            outer_node = nodes_by_edge_idx[e_idx][1]
            nodes_on_edges_that_cross_this_time.append(outer_node)

    return nodes_on_edges_that_cross_this_time
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


import os
from random import sample
from io import StringIO
import numpy as np
from Bio import Phylo
def execute_gene_shedding():

        max_sim_time = 500
        WGD_time = 200
        bin_size = 1
        time_slices = np.arange(WGD_time, max_sim_time, bin_size)

        # At each time-slice post WGD, the number of extant genes are counted.
        # The total number of genes to be shed is calculated as a percent of the
        # #total on the time slice.
        num_genes_to_shed = 2
        n1 = "(O:500,(P1:300,P2:300):200);"
        n2 = "(O:500,(P1:300,P2:300):200);"
        n3 = "(O:500,(P1:300,P2:300):200);"
        newicks = {1: n1, 2: n2, 3: n3, 4: n1, 5: n2, 6: n3}
        gt_to_loose_a_gene = sample(list(newicks.keys()), num_genes_to_shed)
        self.assertEqual(len(gt_to_loose_a_gene), 2)

        # For each gene to be shed, a GT is randomly selected without replacement,
        # and a vertex from that GT which crosses the time-slice is removed.
        #   This process is repeated until the desired number of genes are shed.
        unprunable_leaves = ['O', 'O1', 'O2']
        time_slice = time_slices[-10]
        num_genes_to_remove_per_gene_tree = 1
        print("time slice:\t" + str(time_slice))
        for gt in gt_to_loose_a_gene:
            tree_1 = Phylo.read(StringIO(newicks[gt]), "newick")
            nodes_on_edges_that_cross_this_time = (
                get_edges_that_cross_this_time(tree_1, time_slice))

            list_of_terminal_leaves_to_remove = self.chose_leaves_to_remove(
                nodes_on_edges_that_cross_this_time, num_genes_to_remove_per_gene_tree, unprunable_leaves)

            print("tree before:")
            terminals_after = [t.name for t in tree_1.get_terminals()]
            self.assertEqual(len(terminals_after), 3)
            self.assertEqual(terminals_after, ['O', 'P1', 'P2'])
            Phylo.draw_ascii(tree_1)

            for leaf in list_of_terminal_leaves_to_remove:
                print("pruning leaf " + str(leaf))
                tree_1.prune(leaf)

            print("tree after:")
            if len(tree_1.clade.clades) > 0:
                Phylo.draw_ascii(tree_1)
            else:
                print("No branches remaining")

            terminals_after = [t.name for t in tree_1.get_terminals()]
            self.assertEqual(len(terminals_after), 2)
            self.assertEqual(terminals_after[0], 'O')
            self.assertTrue((terminals_after[1] == 'P1') or (terminals_after[1] == 'P2'))

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

        nodes_to_remove_list = sample(nodes_allowed_to_prune, num_to_remove)
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


def get_edges_that_cross_this_time(tree,
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


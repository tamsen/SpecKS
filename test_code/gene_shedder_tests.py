import unittest
from random import sample

import numpy as np
from Bio import Phylo
from pathlib import Path
from io import StringIO

class GeneShedderTestCases(unittest.TestCase):
    def test_gene_sheddder(self):


        newick_1="(O:500,(P1:300,P2:300):200);"
        tree_1 = Phylo.read(StringIO(newick_1), "newick")
        Phylo.draw_ascii(tree_1)
        max_sim_time=500
        WGD_time = 200

        X = Phylo.to_networkx(tree_1 )
        edges = list(X.edges)

        #check we calculate the distance intervals right
        distance_intervals=get_distance_intervals_for_edge(edges, tree_1)
        print(distance_intervals)
        self.assertEqual(len(distance_intervals.keys()), 4)
        self.assertEqual(distance_intervals[0], (0,500))
        self.assertEqual(distance_intervals[1], (0,200))
        #{0: (0, 500.0), 1: (0, 200.0), 2: (200.0, 500.0), 3: (200.0, 500.0)}

        #check we know when a gene-tree branch crosses a time interval
        edges_that_cross_this_time= get_edges_that_cross_this_time(tree_1, 5)
        print(edges_that_cross_this_time)
        self.assertEqual(len(edges_that_cross_this_time), 2)
        self.assertEqual(edges_that_cross_this_time[0], 0)
        self.assertEqual(edges_that_cross_this_time[1], 1)


    def test_execute_gene_shedding_1(self):

        max_sim_time=500
        WGD_time = 200
        bin_size = 1
        time_slices = np.arange(WGD_time, max_sim_time, bin_size)

        #At each time-slice post WGD, the number of extant genes are counted.
        #The total number of genes to be shed is calculated as a percent of the
        # #total on the time slice.
        num_genes_to_shed=2
        n1="(O:500,(P1:300,P2:300):200);"
        n2="(O:500,(P1:300,P2:300):200);"
        n3="(O:500,(P1:300,P2:300):200);"
        newicks={1:n1,2:n2,3:n3}
        gt_to_loose_a_gene= sample(list(newicks.keys()), num_genes_to_shed)
        self.assertEqual(len(gt_to_loose_a_gene), 2)

        # For each gene to be shed, a GT is randomly selected without replacement,
        # and a vertex from that GT which crosses the time-slice is removed.
        #   This process is repeated until the desired number of genes are shed.

        time_slices= time_slices[0]
        for gt in gt_to_loose_a_gene:
            tree_1 = Phylo.read(StringIO(newicks[gt]), "newick")
            X = Phylo.to_networkx(tree_1)

            # test removing a tree node by edge...
            edges = list(X.edges)
            first_edge=edges[0]
            second_node_in_edge=first_edge[1]
            print("nodes to remove:" + str(second_node_in_edge))
            print("tree before:")
            Phylo.draw_ascii(tree_1)
            tree_1.prune(second_node_in_edge)
            print("tree after:")
            Phylo.draw_ascii(tree_1)
            #

            edges_that_cross_this_time, nodes_on_edges_that_cross_this_time = get_edges_that_cross_this_time(tree_1, 5)
            nodes_to_remove = sample(nodes_on_edges_that_cross_this_time , 1)
            print("nodes to remove:" + str(nodes_to_remove))
            print("tree before:")
            Phylo.draw_ascii(tree_1)
            for edge_idx in edges_that_cross_this_time:
                print("pruning " + str(edge))
                print("pruning " + str(edge_idx ))
                tree_1.prune(edge_idx )
            print("tree after:")
            Phylo.draw_ascii(tree_1)
        print()


def get_distance_intervals_for_edge(edges, tree):
    distance_intervals_by_edge_idx = {}
    nodes_by_edge_idx = {}
    for e_idx in range(0, len(edges)):
        edge = edges[e_idx]
        dist1 = tree.distance(edge[0])
        dist2 = tree.distance(edge[1])
        if dist1 < dist2:
            interval = (dist1, dist2)
            nodes = (edge[0],edge[1])
        else:
            interval = (dist2, dist1)
            nodes = (edge[1],edge[0])
        distance_intervals_by_edge_idx[e_idx] = interval
        nodes_by_edge_idx[e_idx] = nodes

    return distance_intervals_by_edge_idx,nodes_by_edge_idx


def get_edges_that_cross_this_time(tree,
                                   time_slice):
    X = Phylo.to_networkx(tree)
    edges = list(X.edges)
    distance_intervals_by_edge_idx,nodes_by_edge_idx=get_distance_intervals_for_edge(edges, tree)

    tree.prune(edges[0])
    list_of_edges_that_cross_this_time=[]
    nodes_on_edges_that_cross_this_time=[]
    for e_idx in range(0, len(edges)):
        distance_interval = distance_intervals_by_edge_idx[e_idx]
        nodes = nodes_by_edge_idx[e_idx]
        # does edge cross the slice of time?
        if (distance_interval[0]) < time_slice and (distance_interval[1] >= time_slice):
            list_of_edges_that_cross_this_time.append(e_idx)
            outer_node=  nodes_by_edge_idx[e_idx][1]
            nodes_on_edges_that_cross_this_time.append(outer_node)

    return list_of_edges_that_cross_this_time, nodes_on_edges_that_cross_this_time

#https://stackoverflow.com/questions/46718281/remove-percentage-of-items-in-a-list-randomly
def pick_edges_to_remove(my_list, fraction_to_remove):
    return sample(my_list, int(len(my_list) * fraction_to_remove))

if __name__ == '__main__':
    unittest.main()

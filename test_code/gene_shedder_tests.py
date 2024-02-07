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

        #TODO: check we remove the right num edges
        execute_gene_shedding(tree_1, max_sim_time, WGD_time)


def execute_gene_shedding(tree, max_sim_time, WGD_time):
    X = Phylo.to_networkx(tree)
    nodes = list(X.nodes)
    edges = list(X.edges)
    col_headers = ["node_name", "dist", "leaves"]
    data_lines = []
    bin_size = 1 #1 MY
    bins_post_WGD=  np.arange(WGD_time, max_sim_time, bin_size)
    print("bins_post_WGD:\t" + str(bins_post_WGD[0:5]) + "..." + str(bins_post_WGD[-5:]))
    fraction_to_remove=0.1
    edges_that_cross_time_list_by_time={}
    for time_slice in bins_post_WGD:

        edges_that_cross_this_time=get_edges_that_cross_this_time(tree,time_slice)

        #TODO - worry about rounding errors
        # If the number of edges is so very low, that the number to remove
        # ends up always being less that 0.5, ie, rounding to zero,
        # then we will not get any genes removed..
        # should probably change the algorithm not to pick n at random,
        # put to give each gene a chance of being removed and see if it executes.
        edges_to_remove=pick_edges_to_remove(edges_that_cross_this_time, fraction_to_remove)

        #for edge in edges_to_remove:
        print("edges_to_remove:\t" + str(edges))

def get_distance_intervals_for_edge(edges, tree):
    distance_intervals_by_edge_idx = {}
    for e_idx in range(0, len(edges)):
        edge = edges[e_idx]
        dist1 = tree.distance(edge[0])
        dist2 = tree.distance(edge[1])
        interval = (min(dist1, dist2), max(dist1, dist2))
        distance_intervals_by_edge_idx[e_idx] = interval
    return distance_intervals_by_edge_idx


def get_edges_that_cross_this_time(tree,
                                   time_slice):
    X = Phylo.to_networkx(tree)
    edges = list(X.edges)
    distance_intervals_by_edge_idx=get_distance_intervals_for_edge(edges, tree)
    list_of_edges_that_cross_this_time=[]
    for e_idx in range(0, len(edges)):
        distance_interval = distance_intervals_by_edge_idx[e_idx]

        # does edge cross the slice of time?
        if (distance_interval[0]) < time_slice and (distance_interval[1] >= time_slice):
            list_of_edges_that_cross_this_time.append(e_idx)

    return list_of_edges_that_cross_this_time

#https://stackoverflow.com/questions/46718281/remove-percentage-of-items-in-a-list-randomly
def pick_edges_to_remove(my_list, fraction_to_remove):
    return sample(my_list, int(len(my_list) * fraction_to_remove))

if __name__ == '__main__':
    unittest.main()

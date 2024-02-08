import os
import unittest
from random import sample
import numpy as np
from Bio import Phylo
from io import StringIO
from pipeline_modules import gene_shedder

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
        distance_intervals=gene_shedder.get_distance_intervals_for_edge(edges, tree_1)
        print(distance_intervals)
        self.assertEqual(len(distance_intervals.keys()), 4)
        self.assertEqual(distance_intervals[0], (0,500))
        self.assertEqual(distance_intervals[1], (0,200))
        #{0: (0, 500.0), 1: (0, 200.0), 2: (200.0, 500.0), 3: (200.0, 500.0)}

        #check we know when a gene-tree branch crosses a time interval
        edges_that_cross_this_time= gene_shedder.get_edges_that_cross_this_time(tree_1, 5)
        print(edges_that_cross_this_time)
        self.assertEqual(len(edges_that_cross_this_time), 2)
        self.assertEqual(edges_that_cross_this_time[0], 0)
        self.assertEqual(edges_that_cross_this_time[1], 1)


    def test_execute_gene_shedding_1(self):

        test_out = "test_out"
        if not os.path.exists(test_out):
            os.makedirs(test_out)

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
        newicks={1:n1,2:n2,3:n3,4:n1,5:n2,6:n3}
        gt_to_loose_a_gene= sample(list(newicks.keys()), num_genes_to_shed)
        self.assertEqual(len(gt_to_loose_a_gene), 2)

        # For each gene to be shed, a GT is randomly selected without replacement,
        # and a vertex from that GT which crosses the time-slice is removed.
        #   This process is repeated until the desired number of genes are shed.
        unprunable_leaves = ['O','O1','O2']
        time_slice = time_slices[-10]
        num_genes_to_remove_per_gene_tree=1
        print("time slice:\t" + str(time_slice))
        for gt in gt_to_loose_a_gene:
            tree_1 = Phylo.read(StringIO(newicks[gt]), "newick")
            nodes_on_edges_that_cross_this_time = (
                gene_shedder.get_edges_that_cross_this_time(tree_1, time_slice))


            list_of_terminal_leaves_to_remove = gene_shedder.chose_leaves_to_remove(
                nodes_on_edges_that_cross_this_time, num_genes_to_remove_per_gene_tree, unprunable_leaves)

            print("tree before:")
            terminals_after= [t.name for t in tree_1.get_terminals()]
            self.assertEqual(len(terminals_after), 3)
            self.assertEqual(terminals_after, ['O','P1','P2'])
            Phylo.draw_ascii(tree_1)

            for leaf in list_of_terminal_leaves_to_remove:
                    print("pruning leaf " + str(leaf))
                    tree_1.prune(leaf)

            print("tree after:")
            if len(tree_1.clade.clades) > 0:
                Phylo.draw_ascii(tree_1)
            else:
                print("No branches remaining")


            terminals_after= [t.name for t in tree_1.get_terminals()]
            self.assertEqual(len(terminals_after), 2)
            self.assertEqual(terminals_after[0],'O')
            self.assertTrue((terminals_after[1]=='P1') or (terminals_after[1]=='P2'))

            handle = StringIO()
            out_file=os.path.join(test_out,"text_newick_generation.txt")
            Phylo.write(tree_1, out_file, "newick")
            Phylo.write(tree_1, handle , "newick")
            print(handle )
            contents = handle.getvalue()
            print(contents)

if __name__ == '__main__':
    unittest.main()

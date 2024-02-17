import os
import unittest
from io import StringIO
from Bio import Phylo

from pipeline_modules import gene_tree_maker, gene_tree_data
#from pipeline_modules.gene_tree_maker import unprune_outgroup_2
from pipeline_modules.gene_tree_data import unprune_outgroup

class GeneTreeTests(unittest.TestCase):
    def test_newick_simplify(self):

        test_data_folder="gene_tree_maker_test_data"
        results_dict= gene_tree_maker.read_pruned_trees(test_data_folder, 500)

        files=list(results_dict.keys())
        self.assertEqual(len(files), 1)

        taggy_tree=results_dict[files[0]].original_newick
        expected_clean_tree="((((G0_1:43.39704,(G1_1:18.169443,G2_1:18.169443)G3_1:25.227596824464555)"+\
            "G4_1:10.494192025367472,G5_1:53.891232)G6_1:136.6172393732139,G7_1:190.50847)"+\
            "G8_1:9.491528753855258,G9_2:200.0)G10_3:300.0;"

        clean_tree = gene_tree_data.clean_newick(taggy_tree)
        print(taggy_tree)
        print(clean_tree)
        self.assertEqual(clean_tree, expected_clean_tree)

    def test_read_leaf_maps(self):

        test_data_file="GeneTree0.pruned.leafmap"
        test_data_folder="gene_tree_maker_test_data"
        full_path=os.path.join(test_data_folder,test_data_file)
        leaves_by_species,num_extant_leaves = gene_tree_data.read_leaf_map_data(full_path)
        print(leaves_by_species)
        self.assertEqual(len(leaves_by_species), 3)

    def test_get_distances_from_newick(self):

        gt_leaf_map="GeneTree4.pruned.leafmap"
        gt_file="GeneTree4.pruned.tree"
        test_data_folder="gene_tree_maker_test_data"
        leaf_map_path=os.path.join(test_data_folder,gt_leaf_map)
        tree_path=os.path.join(test_data_folder,gt_file)
        leaves_by_species,num_extant_leaves = gene_tree_data.read_leaf_map_data(leaf_map_path)
        print(leaves_by_species)
        self.assertEqual(len(leaves_by_species), 1)
        gene_tree_result = gene_tree_data.gene_tree_result("gt_name",tree_path, leaf_map_path)
        tree = gene_tree_result.tree
        X = Phylo.to_networkx(tree)
        nodes = list(X.nodes)
        sum_of_the_branch_lengths = 0
        root_distance=tree.clade.branch_length
        clades=tree.find_clades()
        for node in nodes:
            distance = tree.distance(node)
            sum_of_the_branch_lengths = sum_of_the_branch_lengths + distance
            name = "no name"
            if node.name:
                name = node.name
            print(name + " distance from first node:\t" + str(distance))
            print(name + " distance from root:\t" + str(distance+root_distance))
        total_branch_length = tree.total_branch_length()
        expected_total_branch_length = 399.86767+399.86767+100.13233

        print("total_branch_length:\t" + str(total_branch_length))
        print("expected_branch_length:\t" + str(expected_total_branch_length))
        self.assertAlmostEqual(total_branch_length, expected_total_branch_length,4)


#https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-13-209/tables/1
# https: // biopython.org / docs / 1.75 / api / Bio.Phylo.BaseTree.html?highlight = split  # Bio.Phylo.BaseTree.TreeMixin.split
    def test_unprune_outgroup(self):

        newick1="(G0_1:292.9058231549251,(G1_2:167.19426864617353,G2_2:166.32932582899295)G3_2:128.15959745471997)G4_3:198.68395473225863;"
        out_group_leaf="O1"
        origin_node_name="O"
        tree1= Phylo.read(StringIO(newick1), "newick")

        new_newick, tree4, new_leaves = gene_tree_data.unprune_outgroup(tree1, 500, out_group_leaf, origin_node_name)
        print(new_newick)
        terminal_leaves = tree4.get_terminals()
        print("tree4")
        Phylo.draw_ascii(tree4)

        for node in terminal_leaves:
            distance=tree4.distance(node)
            if node.name:
                print("node:\t" + node.name + "\tdistance:\t" + str(distance))
            else:
                print("node:\t" + "no name"+ "\tdistance:\t" + str(distance))

            self.assertTrue(distance>490)
            self.assertTrue(distance<510)

        self.assertEqual(len(terminal_leaves), 4)

if __name__ == '__main__':
    unittest.main()

import unittest
from io import StringIO
from Bio import Phylo

from pipeline_modules import gene_tree_maker
from pipeline_modules.gene_tree_maker import unprune_outgroup_2


class GeneTreeTests(unittest.TestCase):
    def test_newick_simplify(self):

        test_data_folder="gene_tree_maker_test_data"
        results_dict= gene_tree_maker.read_pruned_trees(test_data_folder)

        files=list(results_dict.keys())
        self.assertEqual(len(files), 1)

        taggy_tree=results_dict[files[0]].newick_with_tags
        expected_clean_tree="((((G0_1:43.39704,(G1_1:18.169443,G2_1:18.169443)G3_1:25.227596824464555)"+\
            "G4_1:10.494192025367472,G5_1:53.891232)G6_1:136.6172393732139,G7_1:190.50847)"+\
            "G8_1:9.491528753855258,G9_2:200.0)G10_3:300.0;"

        clean_tree = gene_tree_maker.clean_newick(taggy_tree)
        print(taggy_tree)
        print(clean_tree)
        self.assertEqual(clean_tree, expected_clean_tree)

    def test_read_leaf_maps(self):

        test_data_file="GeneTree0.pruned.leafmap"
        result= gene_tree_maker.gene_tree_result()
        result= gene_tree_maker.read_leaf_map(test_data_file, result)

        leaf_dict_by_species=result.leaves_by_species
        print(leaf_dict_by_species)
        self.assertEqual(len(leaf_dict_by_species), 3)

    def test_read_leaf_maps(self):

        test_data_file="GeneTree1.pruned.info"
        result= gene_tree_maker.gene_tree_result()
        result= gene_tree_maker.read_gene_tree_info(test_data_file, result)

        info_dict_by_species=result.info_dict
        print(info_dict_by_species)
        self.assertEqual(len(info_dict_by_species), 17)
        self.assertEqual(info_dict_by_species["No. of vertices"], "7")

#https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-13-209/tables/1
# https: // biopython.org / docs / 1.75 / api / Bio.Phylo.BaseTree.html?highlight = split  # Bio.Phylo.BaseTree.TreeMixin.split
    def test_unprune_outgroup(self):

        newick1="(G0_1:292.9058231549251,(G1_2:167.19426864617353,G2_2:166.32932582899295)G3_2:128.15959745471997)G4_3:198.68395473225863;"
        out_group_leaf="O1"
        origin_node_name="O"
        new_newick, tree4 = unprune_outgroup_2(newick1, 500, out_group_leaf, origin_node_name)
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

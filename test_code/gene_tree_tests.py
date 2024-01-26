import unittest
import gene_tree_maker

class GeneTreeTests(unittest.TestCase):
    def test_newick_simplify(self):

        test_data_folder="."
        results_dict=gene_tree_maker.read_pruned_trees(test_data_folder)

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
        result=gene_tree_maker.gene_tree_result()
        result=gene_tree_maker.read_leaf_map(test_data_file,result)

        leaf_dict_by_species=result.leaves_by_species
        print(leaf_dict_by_species)
        self.assertEqual(len(leaf_dict_by_species), 3)

    def test_read_leaf_maps(self):

        test_data_file="GeneTree1.pruned.info"
        result=gene_tree_maker.gene_tree_result()
        result=gene_tree_maker.read_gene_tree_info(test_data_file,result)

        info_dict_by_species=result.info_dict
        print(info_dict_by_species)
        self.assertEqual(len(info_dict_by_species), 17)
        self.assertEqual(info_dict_by_species["No. of vertices"], "7")

if __name__ == '__main__':
    unittest.main()

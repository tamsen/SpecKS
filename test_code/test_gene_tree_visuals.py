import os
import unittest

import networkx as nx
from Bio import Phylo
from matplotlib import pyplot as plt
import numpy as np
from pathlib import Path
from pipeline_modules import gene_tree_maker
from visualization import gene_tree_visuals, tree_visuals_by_phylo
from io import StringIO


class VisualizationTests2(unittest.TestCase):

    def test_newick_to_ascii(self):

        #newick_string = "(O:500,(P1:300,P2:300):200);"
        newick_string =get_gt1()
        tree = Phylo.read(StringIO(newick_string), "newick")
        Phylo.draw_ascii(tree)
        self.assertEqual(True, True)


    def test_gene_tree_plotting(self):

        species_filter = ["P1", "P2"]
        test_out = "test_out"
        test_in = "gene_tree_visuals_test_data"
        if not os.path.exists(test_out):
            os.makedirs(test_out)

        newick_strings=[get_gt1(),get_gt1(),get_gt2(),get_gt3()]
        plot_names=["gt_test0.png","gt_test1.png","gt_test2.png","gt_test3.png"]
        leaf_data_files=[os.path.join(test_in ,f) for f in
            ["GeneTree0.test.leafmap","GeneTree0.test.leafmap",
                         "GeneTree495.pruned.leafmap", "GeneTree521.pruned.leafmap"]]

        for i in range(0,len(newick_strings)):
            tree = Phylo.read(StringIO(newick_strings[i]), "newick")
            Phylo.draw_ascii(tree)

            leaf_map=get_leaf_map(leaf_data_files[i])
            test_file_to_save = os.path.join("./test_out", plot_names[i])
            expected_file_to_save = os.path.join("./test_out",
                                                 plot_names[i].replace(".png",".phylo.png"))
            tree_visuals_by_phylo.save_tree_plot(newick_strings[i], expected_file_to_save )
            gene_tree_visuals.plot_gene_tree_from_phylo_tree(test_file_to_save,
                species_filter,leaf_map, tree)
            made_a_file=os.path.exists(test_file_to_save)
            self.assertTrue(made_a_file)

    def test_plot_newick24(self):

        species_filter = ["P1", "P2"]
        leaf_data_file=os.path.join("gene_tree_visuals_test_data","GeneTree0.test.leafmap")
        leaf_map = get_leaf_map(leaf_data_file)
        newick_string =get_gt1()
        tree = Phylo.read(StringIO(newick_string), "newick")
        #leaves_to_prune = ["G0_0", "G1_0", "G2_0"]
        #for leaf in leaves_to_prune:
        #    tree.prune(leaf)
        Phylo.draw_ascii(tree)
        file_to_save = os.path.join("test_out", "GeneTree42.png")
        gene_tree_visuals.plot_gene_tree_from_phylo_tree(
            file_to_save,  species_filter,leaf_map,tree)

def get_leaf_map(leaf_map_file_path):
    leaf_map = gene_tree_maker.read_leaf_map(
        leaf_map_file_path, gene_tree_maker.gene_tree_result()).leaves_by_species
    return leaf_map

def get_gt1():
    #from GeneTree0.test.pruned.tree
    gt1="((G0_0:385.45802[ID=0 HOST=0 GUEST=null VERTEXTYPE=Leaf],(G1_0:303.7403[ID=1 HOST=0 GUEST=null VERTEXTYPE=Leaf],G2_0:303.7403[ID=2 HOST=0 GUEST=null VERTEXTYPE=Leaf])G3_0:81.7177139711581[ID=3 HOST=0 GUEST=null VERTEXTYPE=Duplication DISCPT=(1,1)])G4_0:114.54198363837152[ID=4 HOST=0 GUEST=null VERTEXTYPE=Duplication DISCPT=(1,1)],((G5_1:243.39763[ID=5 HOST=1 GUEST=null VERTEXTYPE=Leaf],G6_1:243.39763[ID=6 HOST=1 GUEST=null VERTEXTYPE=Leaf])G7_1:56.60236874941369[ID=7 HOST=1 GUEST=null VERTEXTYPE=Duplication DISCPT=(0,2)],(G8_2:290.54914[ID=8 HOST=2 GUEST=null VERTEXTYPE=Leaf],G9_2:290.54914[ID=9 HOST=2 GUEST=null VERTEXTYPE=Leaf])G10_2:9.45085528985734[ID=10 HOST=2 GUEST=null VERTEXTYPE=Duplication DISCPT=(0,2)])G11_3:200.0[ID=11 HOST=3 GUEST=null VERTEXTYPE=Speciation DISCPT=(1,0)])G12_4:0.0[ID=12 HOST=4 GUEST=null VERTEXTYPE=Speciation DISCPT=(2,0)][NAME=PrunedTree];"
    return gt1

def get_gt2():
    #GeneTree495.pruned.tree
    gt2="((G0_0:2.320191[ID=0 HOST=0 GUEST=null VERTEXTYPE=Leaf],G1_0:2.320191[ID=1 HOST=0 GUEST=null VERTEXTYPE=Leaf])G2_0:497.67980899130964[ID=2 HOST=0 GUEST=null VERTEXTYPE=Duplication DISCPT=(0,1)],((G3_2:30.389777[ID=3 HOST=2 GUEST=null VERTEXTYPE=Leaf],G4_2:30.389777[ID=4 HOST=2 GUEST=null VERTEXTYPE=Leaf])G5_2:159.5032415504242[ID=5 HOST=2 GUEST=null VERTEXTYPE=Duplication DISCPT=(0,1)],G6_2:189.89302[ID=6 HOST=2 GUEST=null VERTEXTYPE=Leaf])G7_2:310.10698[ID=7 HOST=2 GUEST=null VERTEXTYPE=Duplication DISCPT=(0,2)])G8_4:0.0[ID=8 HOST=4 GUEST=null VERTEXTYPE=Speciation DISCPT=(2,0)][NAME=PrunedTree];"
    return gt2

def get_gt3():
    #GeneTree521.pruned.tree
    gt3="(G0_0:500.0[ID=0 HOST=0 GUEST=null VERTEXTYPE=Leaf],((G1_1:300.0[ID=1 HOST=1 GUEST=null VERTEXTYPE=Leaf],G2_2:300.0[ID=2 HOST=2 GUEST=null VERTEXTYPE=Leaf])G3_3:181.08452[ID=3 HOST=3 GUEST=null VERTEXTYPE=Speciation DISCPT=(1,0)],((G4_2:166.45892[ID=4 HOST=2 GUEST=null VERTEXTYPE=Leaf],G5_2:166.45892[ID=5 HOST=2 GUEST=null VERTEXTYPE=Leaf])G6_2:106.52013404221599[ID=6 HOST=2 GUEST=null VERTEXTYPE=Duplication DISCPT=(0,2)],G7_2:272.97905[ID=7 HOST=2 GUEST=null VERTEXTYPE=Leaf])G8_2:208.10547[ID=8 HOST=2 GUEST=null VERTEXTYPE=Duplication DISCPT=(0,2)])G9_3:18.915483193535646[ID=9 HOST=3 GUEST=null VERTEXTYPE=Duplication DISCPT=(1,2)])G10_4:0.0[ID=10 HOST=4 GUEST=null VERTEXTYPE=Speciation DISCPT=(2,0)][NAME=PrunedTree];"
    return gt3
def get_config_file():
    par_dir = Path(__file__).parent.parent
    config_file = os.path.join(par_dir, "sample_configs","cos-bi-config.xml")
    return config_file


if __name__ == '__main__':
    unittest.main()

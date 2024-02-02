import os
import unittest

import networkx as nx
from Bio import Phylo
from matplotlib import pyplot as plt
import numpy as np
from pathlib import Path
from pipeline_modules import gene_tree_maker
from visualization import gene_tree_visuals
from io import StringIO


class VisualizationTests2(unittest.TestCase):

    def test_newick_to_ascii(self):

        #newick_string = "(O:500,(P1:300,P2:300):200);"
        newick_string =get_gt1()
        tree = Phylo.read(StringIO(newick_string), "newick")
        Phylo.draw_ascii(tree)
        self.assertEqual(True, True)


    def test_plot_newick24(self):
        #newick_string = "(O:500,(P1:300,P2:300):200);"

        #SHOULD HAVE 13 NODES AND 12 EDGES
        newick_string =get_gt1()
        tree = Phylo.read(StringIO(newick_string), "newick")
        #leaves_to_prune = ["G0_0", "G1_0", "G2_0"]
        #for leaf in leaves_to_prune:
        #    tree.prune(leaf)
        Phylo.draw_ascii(tree)
        file_to_save = os.path.join("./test_out", "GeneTree42.png")
        gene_tree_visuals.plot_gene_tree_from_phylo_tree(file_to_save, tree)


def get_gt1():
    gt1="((G0_0:385.45802[ID=0 HOST=0 GUEST=null VERTEXTYPE=Leaf],(G1_0:303.7403[ID=1 HOST=0 GUEST=null VERTEXTYPE=Leaf],G2_0:303.7403[ID=2 HOST=0 GUEST=null VERTEXTYPE=Leaf])G3_0:81.7177139711581[ID=3 HOST=0 GUEST=null VERTEXTYPE=Duplication DISCPT=(1,1)])G4_0:114.54198363837152[ID=4 HOST=0 GUEST=null VERTEXTYPE=Duplication DISCPT=(1,1)],((G5_1:243.39763[ID=5 HOST=1 GUEST=null VERTEXTYPE=Leaf],G6_1:243.39763[ID=6 HOST=1 GUEST=null VERTEXTYPE=Leaf])G7_1:56.60236874941369[ID=7 HOST=1 GUEST=null VERTEXTYPE=Duplication DISCPT=(0,2)],(G8_2:290.54914[ID=8 HOST=2 GUEST=null VERTEXTYPE=Leaf],G9_2:290.54914[ID=9 HOST=2 GUEST=null VERTEXTYPE=Leaf])G10_2:9.45085528985734[ID=10 HOST=2 GUEST=null VERTEXTYPE=Duplication DISCPT=(0,2)])G11_3:200.0[ID=11 HOST=3 GUEST=null VERTEXTYPE=Speciation DISCPT=(1,0)])G12_4:0.0[ID=12 HOST=4 GUEST=null VERTEXTYPE=Speciation DISCPT=(2,0)][NAME=PrunedTree];"
    return gt1

def get_gt2():
    gt2="(G0_0:500.0[ID=0 HOST=0 GUEST=null VERTEXTYPE=Leaf],G1_1:500.0[ID=1 HOST=1 GUEST=null VERTEXTYPE=Leaf])G2_4:0.0[ID=2 HOST=4 GUEST=null VERTEXTYPE=Speciation DISCPT=(2,0)][NAME=PrunedTree];"
    return gt2

def get_gt3():
    gt3="((((G0_0:62.147061[ID=0 HOST=0 GUEST=null VERTEXTYPE=Leaf],G1_0:62.147061[ID=1 HOST=0 GUEST=null VERTEXTYPE=Leaf])G2_0:120.57833966761288[ID=2 HOST=0 GUEST=null VERTEXTYPE=Duplication DISCPT=(0,1)],G3_0:182.7254[ID=3 HOST=0 GUEST=null VERTEXTYPE=Leaf])G4_0:211.12567[ID=4 HOST=0 GUEST=null VERTEXTYPE=Duplication DISCPT=(0,2)],((G5_0:200.72794[ID=5 HOST=0 GUEST=null VERTEXTYPE=Leaf],G6_0:200.72794[ID=6 HOST=0 GUEST=null VERTEXTYPE=Leaf])G7_0:173.17475[ID=7 HOST=0 GUEST=null VERTEXTYPE=Duplication DISCPT=(0,2)],(((G8_0:58.077819[ID=8 HOST=0 GUEST=null VERTEXTYPE=Leaf],G9_0:58.077819[ID=9 HOST=0 GUEST=null VERTEXTYPE=Leaf])G10_0:94.49669998269559[ID=10 HOST=0 GUEST=null VERTEXTYPE=Duplication DISCPT=(0,1)],G11_0:152.57452[ID=11 HOST=0 GUEST=null VERTEXTYPE=Leaf])G12_0:29.200733628117394[ID=12 HOST=0 GUEST=null VERTEXTYPE=Duplication DISCPT=(0,2)],G13_0:181.77525[ID=13 HOST=0 GUEST=null VERTEXTYPE=Leaf])G14_0:192.12743578315605[ID=14 HOST=0 GUEST=null VERTEXTYPE=Duplication DISCPT=(0,2)])G15_0:19.948383788675812[ID=15 HOST=0 GUEST=null VERTEXTYPE=Duplication DISCPT=(1,1)])G16_0:106.14893[ID=16 HOST=0 GUEST=null VERTEXTYPE=Duplication DISCPT=(1,1)],(G17_2:118.08417[ID=17 HOST=2 GUEST=null VERTEXTYPE=Leaf],G18_2:118.08417[ID=18 HOST=2 GUEST=null VERTEXTYPE=Leaf])G19_2:381.91583[ID=19 HOST=2 GUEST=null VERTEXTYPE=Duplication DISCPT=(0,1)])G20_4:0.0[ID=20 HOST=4 GUEST=null VERTEXTYPE=Speciation DISCPT=(2,0)][NAME=PrunedTree];"
    return gt3
def get_config_file():
    par_dir = Path(__file__).parent.parent
    config_file = os.path.join(par_dir, "sample_configs","cos-bi-config.xml")
    return config_file


if __name__ == '__main__':
    unittest.main()

import os
import unittest
from Bio import Phylo
import tree_visuals_by_phylo
from pathlib import Path
import polyploid_setup
import config
from pipeline_modules import species_tree_maker
from io import StringIO
import tree_visuals_as_network


class VisualizationTests(unittest.TestCase):
    def test_plot_newick(self):

        #plot a simple newick
        test_out = "./test_out"
        if not os.path.exists(test_out):
            os.makedirs(test_out)

        newick_string1="(O:500,(P1:300,P2:300):200);"
        newick_string2="((G0_0:345.38716[ID=0 HOST=0 GUEST=null VERTEXTYPE=Leaf],G1_0:345.38716[ID=1 HOST=0 GUEST=null VERTEXTYPE=Leaf])G2_0:154.61283965381836[ID=2 HOST=0 GUEST=null VERTEXTYPE=Duplication DISCPT=(1,1)],((G3_1:179.15393[ID=3 HOST=1 GUEST=null VERTEXTYPE=Leaf],G4_1:179.15393[ID=4 HOST=1 GUEST=null VERTEXTYPE=Leaf])G5_1:120.84607035635047[ID=5 HOST=1 GUEST=null VERTEXTYPE=Duplication DISCPT=(0,2)],G6_2:300.0[ID=6 HOST=2 GUEST=null VERTEXTYPE=Leaf])G7_3:200.0[ID=7 HOST=3 GUEST=null VERTEXTYPE=Speciation DISCPT=(1,0)])G8_4:0.0[ID=8 HOST=4 GUEST=null VERTEXTYPE=Speciation DISCPT=(2,0)][NAME=PrunedTree];"

        out_file_name=os.path.join(test_out,"newick1.png")
        tree_visuals_by_phylo.save_tree_plot(newick_string1, out_file_name)

        self.assertEqual(True, os.path.exists(out_file_name))  # add assertion here

        #plot the species tree outline using networkx
        config_file=get_config_file()
        general_sim_config=config.SpecKS_config(config_file)
        polyploid=polyploid_setup.polyploid_data(test_out,"poly1",
                                 300,200,general_sim_config)
        species_tree_out_file_name = os.path.join(test_out, "species1.png")
        species_tree_viz_data = species_tree_maker.plot_species_tree(
            species_tree_out_file_name, polyploid)


        #combine the two..?
        gene_tree_viz_data = tree_visuals_as_network.gene_tree_newick_to_tree_viz_data(newick_string1, "gt_name1")
        out_file_name = os.path.join(test_out, "gene_tree.png")
        #tree_viz_data=[species_tree_viz_data ,gene_tree_viz_data]
        tree_visuals_as_network.plot_gene_tree(out_file_name,gene_tree_viz_data)


def get_config_file():
    par_dir = Path(__file__).parent.parent
    config_file = os.path.join(par_dir, "sample_configs","cos-bi-config.xml")
    return config_file

if __name__ == '__main__':
    unittest.main()

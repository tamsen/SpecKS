import os
import unittest
from Bio import Phylo
from matplotlib import pyplot as plt

from visualization import tree_visuals_by_phylo
from pathlib import Path
import polyploid_setup
import config
from io import StringIO
import tree_visuals_as_network
import networkx as nx


class VisualizationTests(unittest.TestCase):

    def test_plot_newick2(self):
        newick_string = "(O:500,(P1:300,P2:300):200);"
        tree = Phylo.read(StringIO(newick_string), "newick")
        Phylo.draw_ascii(tree)

        file_to_save = os.path.join("../test_code/test_out", "species42.png")
        fig, ax = plt.subplots()
        #X = nx.Graph()
        X = Phylo.to_networkx(tree)
        print(list(X.nodes))
        print(list(X.edges))
        for node in X.nodes:

            print("Node:" + str(node))
            if node.branch_length:
                print(node.branch_length)

        nx.draw(X)
        ax.tick_params(left=True, bottom=True, labelleft=True, labelbottom=False)
        plt.xlim(-50, 50)
        plt.ylim(0, 500 + 50)
        y_ticks = plt.yticks()[0]
        new_ticks = [500 - y_tick for y_tick in y_ticks]
        num_ticks = len(y_ticks)
        plt.title('polyploid.species_name')
        # reverse, since back in time
        ax.set_yticks(y_ticks[0:num_ticks - 1])
        ax.set_yticklabels(new_ticks[0:num_ticks - 1])
        plt.legend()
        plt.ylabel("MYA")
        plt.savefig(file_to_save)
        plt.cla()
        plt.close()
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
        #species_tree_viz_data = species_tree_maker.plot_species_tree(
        #    species_tree_out_file_name, polyploid)


        #combine the two..?

        #gt_viz_data_0 = tree_visuals_as_network.gene_tree_newick_to_tree_viz_data(newick_string1, "gt_name1")
        gt_viz_data_1 = tree_visuals_as_network.gene_tree_newick_to_tree_viz_data(get_gt1(), "gt_name1")
        #gt_viz_data_2 = tree_visuals_as_network.gene_tree_newick_to_tree_viz_data(get_gt2(), "gt_name2")
        #gt_viz_data_3 = tree_visuals_as_network.gene_tree_newick_to_tree_viz_data(get_gt3(), "gt_name3")
        out_file_name = os.path.join(test_out, "gene_tree.png")
        #list_of_tree_viz_data=[species_tree_viz_data ,gt_viz_data_0]
        tree_visuals_as_network.plot_gene_tree(out_file_name, [gt_viz_data_1])

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

import os
import unittest
from Bio import Phylo
from pathlib import Path
from io import StringIO

import polyploid_setup
from pipeline_modules import gene_tree_maker, species_tree_maker, gene_tree_data
from visualization import gene_tree_visuals, tree_visuals_by_phylo
import config
from visualization.combined_tree_view import plot_combined_tree_view
from visualization.gene_tree_visuals import plot_polyploid_gene_tree_alone


class GeneTreeVisualsTests(unittest.TestCase):

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
            tree_visuals_by_phylo.save_tree_plot_from_newick(newick_strings[i], expected_file_to_save)
            sl= polyploid_setup.sim_time_interval_forward_in_time(0, 300,['P1','P2'])
            gene_tree_visuals.plot_polyploid_gene_tree_alone(
                sl,leaf_map, newick_strings[i], "gt_name",
                300, "species_name" ,test_file_to_save)
            made_a_file=os.path.exists(test_file_to_save)
            self.assertTrue(made_a_file)

    def test_plotting_a_species_tree_under_a_gene_tree(self):

        test_in = "gene_tree_visuals_test_data"
        newick_strings=[get_gt1(),get_gt1(),get_gt2(),get_gt3()]
        leaf_data_files=[os.path.join(test_in ,f) for f in
            ["GeneTree0.test.leafmap","GeneTree0.test.leafmap",
                         "GeneTree495.pruned.leafmap", "GeneTree521.pruned.leafmap"]]
        #plot a simple newick
        test_out = "test_out"
        if not os.path.exists(test_out):
            os.makedirs(test_out)

        newick_string1="(O:500,(P1:300,P2:300):200);"
        out_file_name=os.path.join(test_out,"species1_by_phylo.png")
        tree_visuals_by_phylo.save_tree_plot_from_newick(newick_string1, out_file_name)

        self.assertEqual(True, os.path.exists(out_file_name))  # add assertion here

        #plot the species tree outline using networkx
        config_file=get_config_file()
        general_sim_config=config.SpecKS_config(config_file)
        polyploid=polyploid_setup.polyploid_data(test_out,"poly1",
                                 300,200,general_sim_config)
        species_tree_out_file_name = os.path.join(test_out, "species1_by_specks.png")
        species_tree_viz_data = species_tree_maker.plot_species_tree(
                species_tree_out_file_name, polyploid)

        gt_tree_out_file_name = os.path.join(test_out, "gt1_by_specks.png")
        species_filter = ["P1", "P2"]

        gt_tree_viz_data_by_gt_name= {}
        time_of_WGD_MYA=200
        time_of_SPEC_MYA=300

        for i in range(1,4):
            leaf_map = get_leaf_map(leaf_data_files[i])
            sl= polyploid_setup.sim_time_interval_forward_in_time(0, 300,species_filter)
            gt_tree_viz_data =plot_polyploid_gene_tree_alone(sl, leaf_map, newick_strings[i],'gt_name',
                                                             300,'species_name',
                                                             gt_tree_out_file_name)
            gt_tree_viz_data_by_gt_name["gt"+str(i)]=gt_tree_viz_data


        full_sim_time=500
        s_and_gt_tree_out_file_name = os.path.join(test_out, "species_and_gt_by_specks.png")
        plot_combined_tree_view(species_tree_viz_data, gt_tree_viz_data_by_gt_name,
                                time_of_WGD_MYA, time_of_SPEC_MYA,full_sim_time,"species_name",
                            s_and_gt_tree_out_file_name)

    def test_short_gts(self):

        tree_files=["GeneTree0.pruned.tree","GeneTree1.pruned.tree",
                    "GeneTree2.pruned.tree","GeneTree3.pruned.tree",
                    "GeneTree4.pruned.tree"]

        test_out = "test_out"
        test_in = "gene_tree_visuals_test_data"
        if not os.path.exists(test_out):
            os.makedirs(test_out)

        newick_strings=[]
        leaf_maps=[]
        tree_data_files=[os.path.join(test_in,"short_gene_trees" ,f) for f in tree_files]

        for tree_file in tree_data_files:
            newick,leaf_map = get_gt_from_file(tree_file)
            newick_strings.append(newick)
            leaf_maps.append(leaf_map)

        #plot the species tree outline using networkx
        config_file=get_config_file()
        general_sim_config=config.SpecKS_config(config_file)
        polyploid=polyploid_setup.polyploid_data(test_out,"poly1",
                                 300,200,general_sim_config)
        species_tree_out_file_name = os.path.join(test_out, "species1_by_specks.png")
        species_tree_viz_data = species_tree_maker.plot_species_tree(
                species_tree_out_file_name, polyploid)

        gt_tree_out_file_name = os.path.join(test_out, "gt1_by_specks.png")
        species_filter = ["P1", "P2"]

        gt_tree_viz_data_by_gt_name= {}
        time_of_WGD_MYA=200
        time_of_SPEC_MYA=300

        for i in range(0,5):

            sl= polyploid_setup.sim_time_interval_forward_in_time(0, 300,species_filter)
            gt_tree_viz_data =plot_polyploid_gene_tree_alone(sl, leaf_maps[i], newick_strings[i],
                                                   "gt" + str(i),300,'species_name',
                                                             gt_tree_out_file_name)
            gt_tree_viz_data_by_gt_name["gt"+str(i)]=gt_tree_viz_data


        full_sim_time=500
        s_and_gt_tree_out_file_name = os.path.join(test_out, "species_and_short_gt_by_specks.png")
        plot_combined_tree_view(species_tree_viz_data, gt_tree_viz_data_by_gt_name,
                                time_of_WGD_MYA, time_of_SPEC_MYA,full_sim_time,'species_name',
                            s_and_gt_tree_out_file_name)

    def test_add_back_root(self):


        tree_files=["GeneTree0.pruned.tree","GeneTree1.pruned.tree",
                    "GeneTree2.pruned.tree","GeneTree3.pruned.tree",
                    "GeneTree4.pruned.tree"]

        test_in = "gene_tree_visuals_test_data"
        newick_strings=[]
        tree_data_files=[os.path.join(test_in,"short_gene_trees" ,f) for f in tree_files]

        for tree_file in tree_data_files:
            newick,leaf_map = get_gt_from_file(tree_file)
            newick_strings.append(newick)


        # we need some smarts to add back in the outgroup, if it got lost..
        for newick_1 in newick_strings:

            tree_1 = Phylo.read(StringIO(newick_1), "newick")
            Phylo.draw_ascii(tree_1)
            terminals = tree_1.get_terminals()

            for t in terminals:
                dist = tree_1.distance(t)
                print(dist)

            out_group_leaf="O1"
            origin_node_name="O"

            adjusted_newick, new_tree, [out_group_leaf] = gene_tree_data.unprune_outgroup(tree_1,500,
                                                              out_group_leaf, origin_node_name)
            tree_1_fixed = Phylo.read(StringIO(adjusted_newick), "newick")
            Phylo.draw_ascii(tree_1_fixed)
            terminals = tree_1_fixed.get_terminals()

            tolerance=0.005
            for t in terminals:
                dist = tree_1_fixed.distance(t)
                print(dist)
                self.assertEqual(True, dist>=(500-tolerance))

    def test_plot_newick24(self):

        species_filter = ["P1", "P2"]
        leaf_data_file=os.path.join("gene_tree_visuals_test_data","GeneTree0.test.leafmap")
        leaf_map = get_leaf_map(leaf_data_file)
        newick_string =get_gt1()
        tree = Phylo.read(StringIO(newick_string), "newick")
        #leaves_to_prune = ["G0_0", "G1_0", "G2_0"]
        #for leaf in leaves_to_prune:
        #    tree.prune(leaf)

        #new_clade = Phylo.BaseTree.Clade(branch_length=full_sim_time, name=origin_node_name,
        #                                 clades=[outgroup_clade, old_tree.clade])
        #new_tree = Phylo.BaseTree.Tree.from_clade(new_clade)

        my_clade1=Phylo.BaseTree.Clade( name="foo1",branch_length=500)
        my_clade2=Phylo.BaseTree.Clade( name="foo2",branch_length=500)
        my_clade3=Phylo.BaseTree.Clade( name="foo3",branch_length=500, clades=[my_clade1,my_clade2])
        new_tree = Phylo.BaseTree.Tree.from_clade(my_clade3)
        Phylo.draw_ascii(new_tree)
        file_to_save = os.path.join("test_out", "GeneTree42.png")
        sl = polyploid_setup.sim_time_interval_forward_in_time(0, 300, species_filter)
        gene_tree_visuals.plot_polyploid_gene_tree_alone(
              sl,leaf_map,newick_string, "gt_name",300, "species_name",file_to_save)

def get_leaf_map(leaf_map_file_path):
    leaf_map,num_extant_leaves = gene_tree_data.read_leaf_map_data(leaf_map_file_path)
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

def get_gt_from_file(tree_file):
    leafmap_file = tree_file.replace(".tree", ".leafmap")
    result= gene_tree_data.read_gene_tree_result_from_tree_and_leaf_map_files(tree_file,leafmap_file)
    result= result.add_back_outgroup(500)
    return result.simple_newick,result.leaves_by_species
def get_config_file():
    par_dir = Path(__file__).parent.parent
    config_file = os.path.join(par_dir, "sample_configs","cos-bi-config.xml")
    return config_file


if __name__ == '__main__':
    unittest.main()

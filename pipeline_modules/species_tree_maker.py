import os

import log
from visualization import tree_visuals_by_phylo
from visualization.species_tree_visuals import plot_species_tree


def get_example_autopolyploid_trees(polyploid):

    time_before_WGD=polyploid.FULL_time_MYA - polyploid.WGD_time_MYA
    before_WGD_nstring="(O:{0},P:{1});".format(time_before_WGD,time_before_WGD )
    after_WGD_nstrings=[
        "(O1:{0},P1:{1});".format(polyploid.WGD_time_MYA,polyploid.WGD_time_MYA ),
        "(O2:{0},P2:{1});".format(polyploid.WGD_time_MYA,polyploid.WGD_time_MYA)]

    return [before_WGD_nstring] + after_WGD_nstrings

def get_example_allopolyploid_tree(time_span, speciation_MYA):

    #Example tree 1: (O:500, (P1:200, P2:200):300);

    #Example tree 2: (O:500, (P1:100, P2:100)G2_0: 400);
    #speciation happened 100MYA. We'd expect Ks ~2. to have build up between P1 and P2.
    time_before_speciation = time_span - speciation_MYA
    nstring="(O:{0}, (P1:{1}, P2:{1}):{2});".format(time_span,speciation_MYA, time_before_speciation)
    log.write_to_log("tree generated:\t" + nstring)


    return [nstring]


def make_species_trees(polyploid):

    time_span=polyploid.FULL_time_MYA
    subfolder=os.path.join(polyploid.species_subfolder, str(polyploid.analysis_step_num) + "_species_trees")
    include_visualizations = polyploid.general_sim_config.include_visualizations
    log.write_to_log("making species tree for " + polyploid.species_name)

    if not os.path.exists(subfolder):
        os.makedirs(subfolder)

    if polyploid.is_allo():
        species_trees = get_example_allopolyploid_tree(time_span,polyploid.SPC_time_MYA)
    else:
        species_trees = get_example_autopolyploid_trees(polyploid)

    out_file_png = polyploid.species_name + "_tree_WGD{0}of{1}MY.png".format(polyploid.WGD_time_MYA, time_span)
    out_file_tree = out_file_png.replace(".png", ".tree")

    if include_visualizations:
        plot_species_tree(os.path.join(subfolder, out_file_png), polyploid)


    #save species tree netwick to file, and also plot the tree (if vis are desired)
    i=0
    tree_path=os.path.join(subfolder,out_file_tree)
    with open(tree_path, 'w') as f:

        for tree_as_newick in species_trees:
            f.writelines(tree_as_newick + "\n")

            if include_visualizations:
                out_file_png_by_phylo = out_file_png.replace(".png", "_"+str(i)+"_by_phylo.png")
                tree_visuals_by_phylo.save_tree_plot_from_newick(tree_as_newick, os.path.join(subfolder, out_file_png_by_phylo))

            i=i+1

    polyploid.analysis_step_num=polyploid.analysis_step_num+1
    return species_trees

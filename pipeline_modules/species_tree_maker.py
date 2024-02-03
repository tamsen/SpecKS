import os
from visualization.species_tree_visuals import plot_species_tree


def get_example_autopolyploid_trees(polyploid):

    time_before_WGD=polyploid.FULL_time_MYA - polyploid.WGD_time_MYA
    before_WGD_nstring="(O:{0},P:{1});".format(time_before_WGD,time_before_WGD )
    after_WGD_nstrings=[
        "(O1:{0},P1:{1});".format(polyploid.WGD_time_MYA,polyploid.WGD_time_MYA ),
        "(O2:{0},P2:{1});".format(polyploid.WGD_time_MYA,polyploid.WGD_time_MYA)]

    return [before_WGD_nstring] + after_WGD_nstrings

def get_example_allopolyploid_tree(time_span, time_before_WGD):

    #Example tree (O:500, (P1:200, P2:200):300);
    time_since_WGD=time_span - time_before_WGD
    nstring="(O:{0}, (P1:{1}, P2:{1}):{2});".format(time_span, time_since_WGD,time_before_WGD )
    print("species_tree:\t" + nstring)

    return [nstring]


def make_species_trees(polyploid):

    time_before_WGD=polyploid.WGD_time_MYA
    time_span=polyploid.FULL_time_MYA
    subfolder=os.path.join(polyploid.species_subfolder, str(polyploid.analysis_step_num) + "_species_trees")
    print(subfolder)

    if not os.path.exists(subfolder):
        os.makedirs(subfolder)

    if polyploid.is_allo():
        species_trees = get_example_allopolyploid_tree(time_span,time_before_WGD)
    else:
        species_trees = get_example_autopolyploid_trees(polyploid)

    out_file_png = polyploid.species_name + "_tree_WGD{0}of{1}MY.png".format(time_before_WGD, time_span)
    out_file_tree = out_file_png.replace(".png", ".tree")
    plot_species_tree(os.path.join(subfolder, out_file_png), polyploid)

    #save species tree netwick to file
    tree_path=os.path.join(subfolder,out_file_tree)
    with open(tree_path, 'w') as f:
        for tree in species_trees:
            f.writelines(tree + "\n")

    polyploid.analysis_step_num=polyploid.analysis_step_num+1
    return species_trees

import os
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors


def get_example_autopolyploid_trees(polyploid):

    time_before_WGD=polyploid.FULL_time_MYA - polyploid.WGD_time_MYA
    before_WGD_nstring="(O:{0},P:{1});".format(time_before_WGD,time_before_WGD )
    after_WGD_nstrings=[
        "(O1:{0},P1:{1});".format(time_before_WGD,time_before_WGD ),
        "(O2:{0},P2:{1});".format(polyploid.WGD_time_MYA,polyploid.WGD_time_MYA)]

    return [before_WGD_nstring] + after_WGD_nstrings

def get_example_allopolyploid_tree(time_span, time_before_WGD):

    #Example tree (O:500, (P1:200, P2:200):300);
    time_since_WGD=time_span - time_before_WGD
    nstring="(O:{0}, (P1:{1}, P2:{1}):{2});".format(time_span, time_since_WGD,time_before_WGD )
    print("species_tree:\t" + nstring)

    return [nstring]

def plot_species_tree(file_to_save, polyploid):

    time_before_WGD = polyploid.FULL_time_MYA - polyploid.WGD_time_MYA
    time_before_SPEC = polyploid.FULL_time_MYA - polyploid.SPC_time_MYA
    time_span = polyploid.FULL_time_MYA

    fig, ax = plt.subplots()
    X = nx.Graph()

    if polyploid.is_allo():
        offset_for_width=10
        pos = {0: (0, 0), 1: (0, time_before_SPEC),
               2: (20, time_span),
               3: (-20, time_span),
               4: (0, time_before_SPEC+offset_for_width),
               5: (0, time_before_SPEC - 10*offset_for_width),
               6: (0, time_before_SPEC - 10*offset_for_width)}

        edge_before_WGD = [(1, 0)]
        #edge_after_WGD = [(1, 2), (1, 3)]
        edge_after_WGD = [(5, 2), (6, 3)]
    else:
        pos = {0: (0, 0), 1: (0, time_before_WGD),
               2: (20, time_before_WGD), 3: (-20, time_before_WGD), 4: (20, time_span), 5: (-20, time_span)}

        edge_before_WGD = [(1, 0)]
        edge_after_WGD = [(2, 4), (3, 5)]

    X.add_edges_from(edge_before_WGD + edge_after_WGD)

    #nx.draw_networkx_nodes(X,pos,node_size=3,nodelist=[0,1,5,6],node_color='r',
    #                  ax=ax)

    nx.draw_networkx_edges(X, pos, edge_before_WGD, arrows=False,
                      edge_color=mcolors.CSS4_COLORS['lightgrey'],width=50)
    nx.draw_networkx_edges(X, pos, edge_after_WGD, arrows=False,
                      edge_color=mcolors.CSS4_COLORS['lightgrey'],width=40)

    plt.axhline(y=time_before_WGD, color='r', linestyle='--', label="WGD")
    plt.axhline(y=time_before_SPEC, color='b', linestyle=':', label="SPEC")
    ax.tick_params(left=True, bottom=True, labelleft=True, labelbottom=False)
    plt.xlim(-50,50)
    plt.ylim(0,time_span+50)
    y_ticks=plt.yticks()[0]
    new_ticks=[500-y_tick for y_tick in y_ticks]
    num_ticks=len(y_ticks)
    plt.title(polyploid.species_name)
    #reverse, since back in time
    ax.set_yticks(y_ticks[0:num_ticks-1])
    ax.set_yticklabels(new_ticks[0:num_ticks-1])
    plt.legend()
    plt.ylabel("MYA")
    plt.savefig(file_to_save)
    plt.cla()
    plt.close()

def make_species_trees(polyploid):

    time_before_WGD=polyploid.WGD_time_MYA
    time_span=polyploid.FULL_time_MYA
    subfolder=os.path.join(polyploid.species_subfolder, str(polyploid.analysis_step_num) + "_species_trees")
    print(subfolder)

    if not os.path.exists(subfolder):
        os.makedirs(subfolder)

    if polyploid.is_allo():
        species_trees = get_example_allopolyploid_tree(time_span,time_before_WGD)
        out_file_png = "allopolyploid_species_tree_WGD{0}of{1}MY.png".format(time_before_WGD, time_span)
    else:
        species_trees = get_example_autopolyploid_trees(polyploid)
        out_file_png = "autopolyploid_species_tree_WGD{0}of{1}MY.png".format(time_before_WGD, time_span)

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

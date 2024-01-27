import os

import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors


def get_example_allopolyploid_tree(time_span, time_before_WGD):


    #Time since WGD: 5,10, 15,50,100,200 MYA. Total tree length 500 MY. Make allo and autopoly examples.
    #Example tree (O:500, (P1:200, P2:200):300);

    time_since_WGD=time_span - time_before_WGD
    nstring="(O:{0}, (P1:{1}, P2:{1}):{2});".format(time_span, time_since_WGD,time_before_WGD )
    print("species_tree:\t" + nstring)

    return nstring

def plot_allopolyploid_species_tree(file_to_save,time_span, time_before_WGD):

    fig, ax = plt.subplots()
    pos = {0: (0, 0), 1: (0, time_before_WGD), 2: (20, time_span), 3: (-20, time_span)}
    X=nx.Graph()

    edge_before_WGD=[(1,0)]
    edge_after_WGD=[(1,2), (1,3)]
    X.add_edges_from(edge_before_WGD+edge_after_WGD)

    nx.draw_networkx_nodes(X,pos,node_size=3,nodelist=[0,1,2,3],node_color='r',
                      ax=ax)
    nx.draw_networkx_edges(X, pos, edge_before_WGD, arrows=False,
                      edge_color=mcolors.CSS4_COLORS['lightgrey'],width=50)
    nx.draw_networkx_edges(X, pos, edge_after_WGD, arrows=False,
                      edge_color=mcolors.CSS4_COLORS['lightgrey'],width=40)


    ax.tick_params(left=True, bottom=True, labelleft=True, labelbottom=False)
    plt.xlim(-50,50)
    plt.ylim(0,time_span+50)
    y_ticks=plt.yticks()[0]
    new_ticks=[500-y_tick for y_tick in y_ticks]
    num_ticks=len(y_ticks)
    #reverse, since back in time
    ax.set_yticks(y_ticks[0:num_ticks-1])
    ax.set_yticklabels(new_ticks[0:num_ticks-1])
    plt.ylabel("MYA")
    plt.savefig(file_to_save)
    plt.cla()
    plt.close()

def make_species_trees(config):

    species_results=[]
    time_before_WGD=300
    time_span=500
    out_dir = config.output_folder
    print(out_dir )
    subfolder=os.path.join(out_dir, str(config.sim_step_num) + "_species_trees")
    print(subfolder)

    if not os.path.exists(subfolder):
        os.makedirs(subfolder)

    out_file_png = "allopolyploid_species_tree_WGD{0}of{1}MY.png".format(time_before_WGD, time_span)
    out_file_tree= out_file_png .replace(".png",".tree")

    species_tree = get_example_allopolyploid_tree(time_span,time_before_WGD)
    plot_allopolyploid_species_tree(os.path.join(subfolder,out_file_png),
                                    time_span,time_before_WGD)

    #save species tree netwick to file
    tree_path=os.path.join(subfolder,out_file_tree)
    with open(tree_path, 'w') as f:
        f.writelines(species_tree + "\n")

    config.sim_step_num=config.sim_step_num+1
    return species_tree

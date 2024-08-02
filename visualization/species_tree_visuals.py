import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from visualization.combined_tree_view import tree_viz_data


def plot_species_tree(file_to_save, polyploid):

    time_before_WGD = polyploid.FULL_time_MYA - polyploid.WGD_time_MYA
    time_before_SPEC = polyploid.FULL_time_MYA - polyploid.DIV_time_MYA
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

    edges=edge_before_WGD + edge_after_WGD
    X.add_edges_from(edges)

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

    #TODO - clean up, this is a messy way to do it...
    species_tree_viz_data = save_tree_vis_data(edges, polyploid, pos)

    return species_tree_viz_data


def save_tree_vis_data(edges, polyploid, pos):
    species_tree_viz_data = tree_viz_data()
    species_tree_viz_data.verts = edges
    species_tree_viz_data.points = pos
    species_tree_viz_data.color = mcolors.CSS4_COLORS['lightgrey']
    species_tree_viz_data.name = polyploid.species_name
    species_tree_viz_data.width = 50
    species_tree_viz_data.alpha = 1.0
    return species_tree_viz_data

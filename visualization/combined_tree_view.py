import networkx as nx
from matplotlib import pyplot as plt

def plot_combined_tree_view(tree_viz_data_list,time_of_WGD_MYA, time_of_SPEC_MYA,
                            file_to_save):

    time_before_WGD=500-time_of_WGD_MYA
    time_before_SPEC=500-time_of_SPEC_MYA
    fig, ax = plt.subplots()
    X = nx.Graph()

    for data_to_plot in tree_viz_data_list:
        X.add_edges_from(data_to_plot.verts)
        nx.draw_networkx_edges(X, data_to_plot.points, data_to_plot.verts, arrows=False,
                           edge_color=data_to_plot.color,alpha=data_to_plot.alpha , width=data_to_plot.width)

        #if data_to_plot.labels:
        #    plot_node_labels(data_to_plot.points, data_to_plot.labels, plt)

    ax.tick_params(left=True, bottom=True, labelleft=True, labelbottom=True)
    plt.axhline(y=time_before_WGD, color='r', linestyle='--', label="WGD")
    plt.axhline(y=time_before_SPEC, color='b', linestyle=':', label="SPEC")
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

def plot_node_labels(points, labels, plt):
    for i, pos in points.items():
        plt.text(pos[0], pos[1], labels[i])
class tree_viz_data():
    points={}
    verts=[]
    color='k'
    name=""
    width=1
    alpha=0.5
    labels={}
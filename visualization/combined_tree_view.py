import networkx as nx
import random
from matplotlib import pyplot as plt
import matplotlib.colors as mcolors


def color_list():

    return [mcolors.CSS4_COLORS['black'], mcolors.CSS4_COLORS['purple'], mcolors.CSS4_COLORS['gold'],
            mcolors.CSS4_COLORS['hotpink'], mcolors.CSS4_COLORS['darkgreen'], mcolors.CSS4_COLORS['brown'],
            mcolors.CSS4_COLORS['green'],mcolors.CSS4_COLORS['red'],mcolors.CSS4_COLORS['cyan'],
            mcolors.CSS4_COLORS['blue'], mcolors.CSS4_COLORS['yellow'], mcolors.CSS4_COLORS['orange']
            ]

def plot_combined_tree_view(species_tree_viz_data,gt_tree_viz_data_by_name,
                            time_of_WGD_MYA, time_of_SPEC_MYA, full_sim_time,
                            species_name, file_to_save, limit_gt_number=False):

    #perturb the gene tree values so they dont sit right on top of each other for plotting
    random.seed(10)
    colors=color_list()
    gt_names=list(gt_tree_viz_data_by_name.keys())
    tree_viz_data_list=[species_tree_viz_data]
    for i in range(0,len(gt_names)):

        if limit_gt_number:
            if i > limit_gt_number:
                break

        name=gt_names[i]
        gt_viz_data = gt_tree_viz_data_by_name[name]
        gt_viz_data.color = colors[i]
        tree_viz_data_list.append(gt_viz_data)
        num_to_add=20*(random.random()-0.5) # random float between 0 and 1
        for id, point in gt_viz_data.points.items():
            new_x=point[0]+num_to_add
            old_y=point[1]
            gt_viz_data.points[id]=(new_x,old_y)

    time_before_WGD=full_sim_time-time_of_WGD_MYA
    time_before_SPEC=full_sim_time-time_of_SPEC_MYA
    fig, ax = plt.subplots()
    X = nx.Graph()

    for data_to_plot in tree_viz_data_list:
        label=None
        if data_to_plot.width<10:
            label=data_to_plot.name

        X.add_edges_from(data_to_plot.verts)
        nx.draw_networkx_edges(X, data_to_plot.points, data_to_plot.verts, arrows=False,
                           edge_color=data_to_plot.color,alpha=data_to_plot.alpha , width=data_to_plot.width,
                               label=label)

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
    plt.title(species_name)
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
    node_distances={}
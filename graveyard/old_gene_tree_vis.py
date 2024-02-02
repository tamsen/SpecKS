import networkx as nx
import numpy as np
from Bio import Phylo
from matplotlib import pyplot as plt
from six import StringIO
import matplotlib.colors as mcolors



def newick_to_points_and_verts2(newick_string):
    newick_string = "(O:500,(P1:300,P2:300):200);"
    tree = Phylo.read(StringIO(newick_string), "newick")
    Phylo.draw_ascii(tree)
    net = Phylo.to_networkx(tree)

def newick_to_points_and_verts(newick_string):

    tree = Phylo.read(StringIO(newick_string), "newick")
    Phylo.draw_ascii(tree)
    #do any pruning
    #ie, we dont want the "outgroup" nodes to be plotted
    leaves_to_prune=["G0_0","G1_0","G2_0"]
    for leaf in leaves_to_prune:
        tree.prune(leaf)
    Phylo.draw_ascii(tree)
    net = Phylo.to_networkx(tree)
    print("nodes:" + str(net.nodes))
    #G0_0	O
    #G1_0	O
    #G2_0
    i = 0
    origin_point = (0, 0)
    points = {i: origin_point}
    labels = {i : tree.clade[0].name}
    connections = []
    recursive_tree_data_to_points(i, origin_point, 0,tree.root,
                                  points, connections, labels)
    return points, connections, labels

def gene_tree_newick_to_tree_viz_data(newick_string, gt_name):

    pos, edges, labels = newick_to_points_and_verts(newick_string)
    gene_tree_viz_data = tree_viz_data()
    gene_tree_viz_data.verts = edges
    gene_tree_viz_data.points = pos
    gene_tree_viz_data.labels = labels
    gene_tree_viz_data.name = gt_name
    gene_tree_viz_data.color = mcolors.CSS4_COLORS['black']
    gene_tree_viz_data.alpha = 0.3
    gene_tree_viz_data.width = 1.0

    return gene_tree_viz_data
def recursive_tree_data_to_points(i, origin_point, recursion_depth, my_root,
                                  points, connections, labels):

    origin_name = i
    origin_point_x = origin_point[0]
    origin_point_y = origin_point[1]

    min = -10
    max = 10
    num_needed = len(my_root.clades) - 1
    step = round((max - min) / num_needed, 2)
    xs = np.arange(min, max + 1, step)
    print(xs)
    xs_index = 0
    recursion_depth = recursion_depth+1
    for clade in my_root.clades:

        i = i + 1
        print("clade " + str(i))
        if clade.name:
            print("clade " + str(i) + " has name " + clade.name)
            print("clade " + str(i) + ", length=" + str(clade.branch_length))
            point = (xs[xs_index] + origin_point_x, clade.branch_length + origin_point_y)
            node_name = "clade " + str(i) + " " + clade.name
            # we dont want to plot the outgroup data
            if clade.name == "O":
                continue
        else:
            node_name = "clade " + str(i)
            print(node_name + " has no name")
            print(node_name + " has length " + str(clade.branch_length))
            if recursion_depth > 1:
                point = (xs[xs_index] + origin_point_x, clade.branch_length + origin_point_y)
            else:
                point = (origin_point_x, clade.branch_length + origin_point_y)

        xs_index = xs_index + 1
        points[i] = point
        labels[i] = node_name
        connections.append((origin_name, i))

        if len(clade.clades) > 0:
            print("has children!")
            recursive_tree_data_to_points(i, point, recursion_depth, clade, points, connections, labels)

    print("connections\t" + str(connections))
    print("points\t" + str(points))
    return points, connections, labels


def plot_gene_tree(file_to_save, tree_viz_data_list):

    fig, ax = plt.subplots()
    X = nx.Graph()

    for data_to_plot in tree_viz_data_list:
        X.add_edges_from(data_to_plot.verts)
        nx.draw_networkx_edges(X, data_to_plot.points, data_to_plot.verts, arrows=False,
                           edge_color=data_to_plot.color,alpha=data_to_plot.alpha , width=data_to_plot.width)
        #if data_to_plot.width < 5:
        #nx.draw_networkx_labels(X, data_to_plot.points, font_size=8, font_family="sans-serif")

        for i in data_to_plot.points.keys():
            plt.text(*data_to_plot.points[i],data_to_plot.labels[i])
        #for p in data_to_plot.points.items:
        #    G.add_node(1, weight=2)
    #my_nodes=dict(X.nodes(data="weight", default=1))


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

class tree_viz_data():
    points={}
    verts=[]
    color='k'
    name=""
    width=1
    alpha=0.5
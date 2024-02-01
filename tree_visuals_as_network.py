import networkx as nx
import numpy as np
from Bio import Phylo
from matplotlib import pyplot as plt
from six import StringIO
import matplotlib.colors as mcolors

def newick_to_points_and_verts(newick_string):

    tree = Phylo.read(StringIO(newick_string), "newick")
    i = 0
    origin_point = (0, 0)
    points = {i: origin_point}
    connections = []
    recursive_tree_data_to_points(i, origin_point, tree.root, points, connections)
    return points, connections

def gene_tree_newick_to_tree_viz_data(newick_string, gt_name):

    pos, edges = newick_to_points_and_verts(newick_string)
    gene_tree_viz_data = tree_viz_data
    gene_tree_viz_data.verts = edges
    gene_tree_viz_data.points = pos
    gene_tree_viz_data.name = gt_name
    gene_tree_viz_data.color = mcolors.CSS4_COLORS['lightgrey']

    return gene_tree_viz_data
def recursive_tree_data_to_points(i, origin_point, my_root, points, connections):

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
    for clade in my_root.clades:

        i = i + 1
        print("clade " + str(i))
        if clade.name:
            print("clade " + str(i) + ", length=" + str(clade.branch_length))
            point = (xs[xs_index] + origin_point_x, clade.branch_length + origin_point_y)

            # we dont want to plot the outgoup data
            if clade.name == "O":
                continue
        else:
            print("clade " + str(i) + " has no name")
            print("clade " + str(i) + " has length " + str(clade.branch_length))
            point = (xs[xs_index] + origin_point_x, clade.branch_length + origin_point_y)

        xs_index = xs_index + 1
        points[i] = point
        connections.append((origin_name, i))

        if len(clade.clades) > 0:
            print("has children!")
            recursive_tree_data_to_points(i, point, clade, points, connections)

    print("connections\t" + str(connections))
    print("points\t" + str(points))
    return points, connections


def plot_gene_tree(file_to_save, tree_viz_data):

    points=tree_viz_data.points
    cxns=tree_viz_data.verts

    fig, ax = plt.subplots()
    X = nx.Graph()

    #for data_to_plot in tree_viz_data_list:

    X.add_edges_from(cxns)
    nx.draw_networkx_edges(X, points, cxns, arrows=False,
                           edge_color="k", width=1)


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
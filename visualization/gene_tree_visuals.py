import os

import networkx as nx
from io import StringIO
from Bio import Phylo
from matplotlib import pyplot as plt
import matplotlib.colors as mcolors
from pipeline_modules import gene_tree_maker, species_tree_maker
from visualization import tree_utils, tree_visuals_by_phylo, combined_tree_view
from visualization.combined_tree_view import tree_viz_data, plot_combined_tree_view


def plot_gene_trees_on_top_of_species_trees(polyploid, config,
                                            gt_tree_viz_data_by_name, out_folder):


    # plot the species tree outline using networkx
    species_tree_out_file_name = os.path.join(out_folder, "species1_by_specks.png")
    species_tree_viz_data = species_tree_maker.plot_species_tree(
        species_tree_out_file_name, polyploid)

    # plot the gene trees over the species tree
    s_and_gt_tree_out_file_name = os.path.join(out_folder, "species_and_gt_by_specks.png")
    combined_tree_view.plot_combined_tree_view(gt_tree_viz_data_by_name,
                                               species_tree_viz_data,
                            polyploid.WGD_time_MYA, polyploid.SPC_time_MYA,
                            polyploid.FULL_time_MYA, s_and_gt_tree_out_file_name)


def plot_gene_tree_alone(
        species_filter, leaf_map, tree_as_newick, gt_name, file_to_save):

    tree = Phylo.read(StringIO(tree_as_newick), "newick")
    time_since_speciation=300.0
    leaf_aim=12.0
    width=5.0 #speciation_tree_width, the envelope

    slope=leaf_aim/time_since_speciation  #this is where inside the parent tree the leaves should end up
    X = Phylo.to_networkx(tree)
    nodes = list(X.nodes)
    edges = list(X.edges)
    node_coordinates_by_i = {i:node_coordinate() for i in range(0,len(nodes))}
    species_by_leaf_dict = get_species_by_leaf_dict_from_leaf_map(leaf_map)

    #calculate x and y values for each node to graph
    node_i_by_name, node_names_by_i = set_node_y_values(node_coordinates_by_i, nodes, tree)

    #Only keep vertexes that touch the species of interest.
    #We dont want to clutter the diagram with the outgroup.
    nodes_to_visualize = set_node_x_values(node_coordinates_by_i, species_by_leaf_dict,
                                          species_filter, slope, width,nodes, tree)
    #print(nodes_to_visualize)

    #Translate so its indexed by "i", not by clade.
    edge_list_in_i_coords = get_edge_list_in_i_coords(edges, node_i_by_name,nodes_to_visualize)
    pos_by_i = pos_dict_in_i_coords(node_coordinates_by_i,nodes_to_visualize)

    #plot it
    plot_gene_tree(X, edge_list_in_i_coords,node_names_by_i, pos_by_i, file_to_save)

    gt_vis_data=save_tree_vis_data(edge_list_in_i_coords, pos_by_i,
                                   node_names_by_i,  gt_name)
    return gt_vis_data

def save_tree_vis_data(edges, pos, labels, name):
    gt_vis_data = tree_viz_data()
    gt_vis_data.verts = edges
    gt_vis_data.points = pos
    gt_vis_data.color = mcolors.CSS4_COLORS['green']
    gt_vis_data.name = name
    gt_vis_data.width = 3
    gt_vis_data.alpha = 0.3
    gt_vis_data.labels = labels
    return gt_vis_data
def plot_gene_tree(X, edge_list_in_i_coords,
                    node_names_by_i, pos_by_i, file_to_save):

    fig, ax = plt.subplots()
    X.add_edges_from(edge_list_in_i_coords)
    nx.draw_networkx_edges(X, pos_by_i, edge_list_in_i_coords)  # nodes_to_visulize)
    plot_node_labels(pos_by_i, node_names_by_i, plt)
    ax.tick_params(left=True, bottom=True, labelleft=True, labelbottom=False)
    # plt.xlim(-50, 50)
    plt.ylim(0, 500 + 50)
    y_ticks = plt.yticks()[0]
    new_ticks = [500 - y_tick for y_tick in y_ticks]
    num_ticks = len(y_ticks)
    plt.title('polyploid.species_name')
    # reverse, since back in time
    ax.set_yticks(y_ticks[0:num_ticks - 1])
    ax.set_yticklabels(new_ticks[0:num_ticks - 1])
    # plt.legend()
    plt.ylabel("MYA")
    plt.savefig(file_to_save)
    plt.cla()
    plt.close()


def pos_dict_in_i_coords(node_coordinates_by_i, nodes_to_visualize):
    pos = {}
    for i, nc in node_coordinates_by_i.items():

        if i in nodes_to_visualize:
            pos[i] = (nc.x, nc.y)

    return pos

def plot_node_labels(pos_by_i, node_names_by_i, plt):

    for i, pos in pos_by_i.items():
        name = "no name" + "_i" + str(i)
        if i in node_names_by_i:
            name = node_names_by_i[i] + "_i" + str(i)

        #print("{0}:{1}\t({2},{3})".format(
        #    i, name, nc.x, nc.y))
        #plt.plot(nc.x, nc.y)
        plt.text(pos[0], pos[1], name)

def get_edge_list_in_i_coords(edges, node_i_by_name,nodes_to_visualize):

    edges_to_visualize = []
    for edge in edges:
        n1 = edge[0]
        n2 = edge[1]
        index_1 = node_i_by_name[n1.name]
        index_2 = node_i_by_name[n2.name]

        if index_1 in nodes_to_visualize and index_2 in nodes_to_visualize:
            edges_to_visualize.append((index_1, index_2))

    return edges_to_visualize


def set_node_x_values(node_coordinates_by_i, species_by_leaf_dict,
                      species_filter, slope, width, nodes, tree):

    nodes_to_visulize={}
    for i in range(0, len(nodes)):
        leafs = nodes[i].get_terminals()
        names = [leaf.name for leaf in leafs]
        species = [species_by_leaf_dict[name] for name in names]

        if (species_filter[0] in species) and (species_filter[1] in species):#ie, does a leaf land in P1?
            f = .0
            nodes_to_visulize[i]=names.copy()
        elif species_filter[0] in species: #ie, does a leaf land in P1?
            f = slope
            nodes_to_visulize[i]=names.copy()
        elif species_filter[1] in species: #ie, does a leaf land in P2/
            f = -1.0 * slope
            nodes_to_visulize[i]=names.copy()
        else: #never conected with our species of interest
            f = 0
            #print(species)
        (b, num_sibs) = tree_utils.birth_order(tree, nodes[i])
        delta = b * width / num_sibs
        #print(names)
        node_coordinates_by_i[i].x = tree.distance(nodes[i]) * f + delta

    return nodes_to_visulize.keys() #we dont really need a full dict. that is just for troubleshooting

def set_node_y_values(node_coordinates_by_i, nodes, tree):
    node_names_by_i = {}
    node_i_by_name = {}
    for i in range(0, len(nodes)):
        node_coordinates_by_i[i].y = tree.distance(nodes[i])
        if nodes[i].name:
            node_names_by_i[i] = nodes[i].name
            node_i_by_name[nodes[i].name] = i
    return node_i_by_name, node_names_by_i


def get_species_by_leaf_dict():
    leaf_map = gene_tree_maker.read_leaf_map(
        "GeneTree0.test.leafmap", gene_tree_maker.gene_tree_result()).leaves_by_species
    leaf_map_by_leaf = {}
    for species in leaf_map:
        for leaf in leaf_map[species]:
            leaf_map_by_leaf[leaf] = species
    return leaf_map_by_leaf

def get_species_by_leaf_dict_from_leaf_map(leaf_map):
    leaf_map_by_leaf = {}
    for species in leaf_map:
        for leaf in leaf_map[species]:
            leaf_map_by_leaf[leaf] = species
    return leaf_map_by_leaf

class node_coordinate():
    y=0
    x=0
    name=""
class tree_data_for_nx():
    points={}
    verts=[]
    color='k'
    name=""
    width=1
    alpha=0.5
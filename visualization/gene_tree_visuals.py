import os

import networkx as nx
from io import StringIO
from Bio import Phylo
from matplotlib import pyplot as plt
import matplotlib.colors as mcolors
from pipeline_modules import sagephy_GBD_model, species_tree_maker
from visualization import tree_utils, tree_visuals_by_phylo, combined_tree_view
from visualization.combined_tree_view import tree_viz_data, plot_combined_tree_view


def histogram_node_distances(polyploid, gt_tree_viz_data_by_gene_tree,
                                               file_prefix,  out_folder):
    dist_histogram_out_file_name = os.path.join(out_folder,
                                              file_prefix + "_" + polyploid.species_name
                                                + "_node_dist_histogram.png")

    MYA_list=[]
    for gene_tree, gt_tree_viz_data in gt_tree_viz_data_by_gene_tree.items():

        #filter out the origin (dist==0) and reverse time, so its MYA and not distance
        MYA= [polyploid.FULL_time_MYA-d for d in gt_tree_viz_data.node_distances if d > 1.0]
        MYA_list= MYA + MYA_list

    #https: // matplotlib.org / stable / gallery / statistics / hist.html

    fig = plt.figure(figsize=(10, 10), dpi=100)
    n_bins=50
    n, bins, patches = plt.hist(MYA_list, bins=n_bins, facecolor='b', alpha=0.25, label='histogram data')
    plt.axvline(x=polyploid.WGD_time_MYA, color='r', linestyle='--', label="WGD")
    plt.axvline(x=polyploid.DIV_time_MYA, color='b', linestyle='--', label="SPEC")
    plt.legend()
    plt.xlim(0,polyploid.FULL_time_MYA + 50)
    plt.xlabel("Node dist (MYA)")
    plt.ylabel("Num nodes")
    plt.title("Node-distance histogram for " + polyploid.species_name  + ", " + file_prefix)
    plt.savefig(dist_histogram_out_file_name)
    plt.clf()
    plt.close()

    return dist_histogram_out_file_name
def plot_gene_trees_on_top_of_species_trees(polyploid,
                                            gt_tree_viz_data_by_name, file_prefix, out_folder):


    # plot the species tree outline using networkx
    species_tree_out_file_name = os.path.join(out_folder,
                                              file_prefix + "_" + polyploid.species_name + "_by_specks.png")
    species_tree_viz_data = species_tree_maker.plot_species_tree(
        species_tree_out_file_name, polyploid)

    # plot the gene trees over the species tree
    max_num_gt_to_visualize=5
    s_and_gt_tree_out_file_name = os.path.join(out_folder, file_prefix + "_" +
                                               polyploid.species_name + "_species_and_5gt_by_specks.png")
    combined_tree_view.plot_combined_tree_view(species_tree_viz_data, gt_tree_viz_data_by_name,
                                               polyploid.WGD_time_MYA, polyploid.DIV_time_MYA,
                                               polyploid.FULL_time_MYA, polyploid.species_name, s_and_gt_tree_out_file_name,
                                               max_num_gt_to_visualize)

    s_and_gt_tree_out_file_name = os.path.join(out_folder, file_prefix + "_" +
                                               polyploid.species_name + "_species_and_all_gt_by_specks.png")
    combined_tree_view.plot_combined_tree_view(species_tree_viz_data, gt_tree_viz_data_by_name,
                                               polyploid.WGD_time_MYA, polyploid.DIV_time_MYA,
                                               polyploid.FULL_time_MYA, polyploid.species_name, s_and_gt_tree_out_file_name)


def plot_polyploid_gene_tree_alone(simulation_leg, leaf_map, tree_as_newick,
                                   gt_name, time_since_speciation,
                                   species_name, file_to_save=False):

    tree = Phylo.read(StringIO(tree_as_newick), "newick")
    tree.clade.name = "Origin"
    leaf_aim=12.0
    width=5.0 #speciation_tree_width, the envelope

    slope=leaf_aim/time_since_speciation  #this is where inside the parent tree the leaves should end up
    X = Phylo.to_networkx(tree)
    nodes = list(X.nodes)
    edges = list(X.edges)
    node_coordinates_by_i = {i:node_coordinate() for i in range(0,len(nodes))}
    species_by_leaf_dict = get_species_by_leaf_dict_from_leaf_map(leaf_map)

    #calculate x and y values for each node to graph
    node_i_by_name, node_names_by_i = set_node_y_values(node_coordinates_by_i,
                                                        nodes,  tree, simulation_leg)

    #Only keep vertexes that touch the species of interest.
    #We dont want to clutter the diagram with the outgroup.
    nodes_to_visualize = set_node_x_values(node_coordinates_by_i, species_by_leaf_dict,
                                           slope, width, nodes, tree, simulation_leg)
    #print(nodes_to_visualize)

    #Translate so its indexed by "i", not by clade.
    edge_list_in_i_coords = get_edge_list_in_i_coords(edges, node_i_by_name,nodes_to_visualize)
    pos_by_i = pos_dict_in_i_coords(node_coordinates_by_i,nodes_to_visualize)
    node_distances=[node_coordinates_by_i[i].distance for i in range(0,len(node_coordinates_by_i)) ]
    #plot it & print it
    if file_to_save:
        plot_gene_tree(X, edge_list_in_i_coords,node_names_by_i, pos_by_i, species_name, file_to_save)
        css_to_save=file_to_save.replace(".png",".dist.csv")
        tree_utils.print_node_distances(tree, css_to_save)

    gt_vis_data=save_tree_vis_data(edge_list_in_i_coords, pos_by_i,
                                   node_names_by_i, node_distances, gt_name)
    return gt_vis_data

def save_tree_vis_data(edges, pos, labels, node_distances, name):
    gt_vis_data = tree_viz_data()
    gt_vis_data.verts = edges
    gt_vis_data.points = pos
    gt_vis_data.color = mcolors.CSS4_COLORS['green']
    gt_vis_data.name = name
    gt_vis_data.width = 3
    gt_vis_data.alpha = 0.5
    gt_vis_data.labels = labels
    gt_vis_data.node_distances = node_distances
    return gt_vis_data
def plot_gene_tree(X, edge_list_in_i_coords,
                    node_names_by_i, pos_by_i, polyploid_species_name, file_to_save):

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
    plt.title(polyploid_species_name)
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

        if (n1.name and n2.name):
            index_1 = node_i_by_name[n1.name]
            index_2 = node_i_by_name[n2.name]

            if index_1 in nodes_to_visualize and index_2 in nodes_to_visualize:
                edges_to_visualize.append((index_1, index_2))
        else:
            print("what happened here?")

    return edges_to_visualize

def set_node_x_values(node_coordinates_by_i, species_by_leaf_dict, slope,
                      width, nodes, tree, simulation_leg):

    species_filter=simulation_leg.subgenomes_during_this_interval
    nodes_to_visulize={}
    for i in range(0, len(nodes)):
        leafs = nodes[i].get_terminals()
        names = [leaf.name for leaf in leafs]
        species = [species_by_leaf_dict[name] for name in names]
        distance =tree.distance(nodes[i]) + simulation_leg.interval_start_time_MY

        if len(species_filter)==1: #ie, "P" is the only subgenome we are looking at...
            if (species_filter[0] in species):  # ie, does a leaf land in "P" ?
                f = .0
                nodes_to_visulize[i] = names.copy()
            else:
                f = 0.0
                distance = -1
        else:
            if (species_filter[0] in species) and (species_filter[1] in species):#ie, does a leaf land in P1 AND P2?
                f = .0
                nodes_to_visulize[i]=names.copy()
            elif species_filter[0] in species: #ie, does a leaf land in P1?
                f = -1.0*slope
                nodes_to_visulize[i]=names.copy()
            elif species_filter[1] in species: #ie, does a leaf land in P2/
                f = 1.0 * slope
                nodes_to_visulize[i]=names.copy()
            else: #never conected with our species of interest
                f = 0.0
                distance = -1

        (b, num_sibs) = tree_utils.birth_order(tree, nodes[i])
        delta = b * width / num_sibs
        #print(names)
        node_coordinates_by_i[i].x = distance * f + delta
        node_coordinates_by_i[i].distance=distance

    return nodes_to_visulize.keys() #we dont really need a full dict. that is just for troubleshooting

def set_node_y_values(node_coordinates_by_i, nodes, tree, simulation_leg):
    node_names_by_i = {}
    node_i_by_name = {}
    for i in range(0, len(nodes)):
        node_coordinates_by_i[i].y = tree.distance(nodes[i]) + simulation_leg.interval_start_time_MY
        if nodes[i].name:
            node_names_by_i[i] = nodes[i].name
            node_i_by_name[nodes[i].name] = i
        else:
            if i==0:
                origin_name="0"
                node_names_by_i[i] = origin_name
                node_i_by_name[origin_name] = i

    return node_i_by_name, node_names_by_i


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
    distance=0
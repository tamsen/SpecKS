from collections import OrderedDict
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
#from file_io import microsat_file_utils
#import distance_metrics
import os


# from pygraphviz import *

# https://stackoverflow.com/questions/46784028/edge-length-in-networkx
# https://networkx.org/documentation/stable/auto_examples/geospatial/plot_points.html#sphx-glr-auto-examples-geospatial-plot-points-py
# https://stackoverflow.com/questions/48065567/legend-based-on-edge-color-in-networkx

def get_vertex_data(input_files,
                    species_filter=None, locality_filter=None, loci_filter=None):

    raw_microsat_data, raw_locality_data = microsat_file_utils.read_microsat_files(input_files,
                                                                          allow_nulls=False,
                                                                          species_filter=species_filter,
                                                                          locality_filter=locality_filter)

    microsat_data, locality_data = distance_metrics.consolidate_all_clones(raw_microsat_data, raw_locality_data )

    vertex_data_by_extraction, locality_data = \
        distance_metrics.calculate_all_similarity_scores(microsat_data,locality_data,
                                                         loci_filter)

    locality_data_by_vertex = OrderedDict()
    for vertex in vertex_data_by_extraction:
        locality_data_by_vertex[vertex] = locality_data[vertex]

    # print histogram of link strengths
    return vertex_data_by_extraction, locality_data_by_vertex


def plot_a_microsat_csv_as_a_network(input_files, out_dir, species=None, sites=None,
                                     loci_filter=None,
                                     plot_title=None, fig_size_tuple=None, force_node_color_by_site=None,
                                     output_file=None, legend_placement=None,

                                     random_seed=7):

    file_base_name = os.path.basename(input_files[0]).replace(".csv", "")

    similarity_scores_dict_by_plantA_by_plant_B, locality_dict = get_vertex_data(input_files,
                                                            species_filter=species,
                                                             loci_filter=loci_filter)

    extractions = similarity_scores_dict_by_plantA_by_plant_B.keys()
    locality_list = list(set(locality_dict.values()))
    extractions_by_locality_by_num_clones = {locality: {} for locality in locality_list}

    colors_by_locality = OrderedDict()
    for locality in locality_list:
        colors_by_locality[locality] = np.array([0.5, 0.5, 0.5])
    if force_node_color_by_site:
        for locality in force_node_color_by_site:
            colors_by_locality[locality] = force_node_color_by_site[locality]


    plt.clf()
    plt.close()
    fig = plt.figure(figsize=(8, 10),dpi=100)
    if fig_size_tuple:
        fig = plt.figure(figsize=fig_size_tuple, dpi=100)

    #fig = plt.figure(figsize=(20, 30), dpi=100)

    line_fatness = 10000
    AS_threshold2 = 0.8
    AS_threshold1 = 0.5
    node_size=700
    clones_by_extraction={}

    G = nx.Graph()
    for extractionA in extractions:

        #print("extractionA:" + extraction)
        clones_by_extraction[extractionA]=[]
        vertex_data_by_extractionB = similarity_scores_dict_by_plantA_by_plant_B[extractionA]
        ordered_edges_list = [[k, v] for k, v in
                              sorted(vertex_data_by_extractionB.items(), reverse=True, key=lambda item: item[1][0])]
        nearest_3 = [edge[0] for edge in ordered_edges_list[0:3]]

        # for extractionB in edges_by_extractionB:
        for extractionB in nearest_3:

            strength = round(vertex_data_by_extractionB[extractionB][0], 2)
            different_alleles_between_them=vertex_data_by_extractionB[extractionB][1]
            diff_alleles_str = ",".join(different_alleles_between_them)

            if strength > AS_threshold2:  # > 1.0, black
                G.add_edge(extractionA, extractionB, label=diff_alleles_str, AS=strength,
                           weight=line_fatness * 3)
            elif strength > AS_threshold1:  # between 0.5 and 1.0, dashed
                G.add_edge(extractionA, extractionB, label=diff_alleles_str, AS=strength,
                           weight=line_fatness * 2)
            else:  #less that 0.5
                G.add_edge(extractionA, extractionB, label=diff_alleles_str, AS=strength,
                           weight=line_fatness * 1)

    elarge = [(u, v) for (u, v, d) in G.edges(data=True) if d["AS"] > AS_threshold2]
    emedium = [(u, v) for (u, v, d) in G.edges(data=True) if d["AS"] > AS_threshold1]
    esmall = [(u, v) for (u, v, d) in G.edges(data=True) if d["AS"] <= AS_threshold1]

    pos = nx.spring_layout(G, seed=random_seed)  # positions for all nodes - seed for reproducibility



    for extraction in locality_dict:
        if extraction in pos:
            locality = locality_dict[extraction]
            num_clones = extraction.count(",")
            data_so_far=extractions_by_locality_by_num_clones[locality]
            if not (num_clones in data_so_far):
                extractions_by_locality_by_num_clones[locality][num_clones]=[]

            extractions_by_locality_by_num_clones[locality][num_clones].append(extraction)

    have_BMW_data=False
    for locality in extractions_by_locality_by_num_clones:
       if ("Original" in locality) or ("BMW" in locality):
           have_BMW_data=True

    print("drawing nodes")
    print("extractions_by_locality_by_num_clones:")
    print(str(extractions_by_locality_by_num_clones))
    BMW_background_color="black"
    TDJ_background_color="c"

    if have_BMW_data:
        nx.draw_networkx_nodes(G, pos, node_size=node_size, node_color=BMW_background_color,
                           nodelist=[extraction], label="BMW collections", node_shape="o")

    nx.draw_networkx_nodes(G, pos, node_size=node_size, node_color=TDJ_background_color,
                           nodelist=[extraction], label="TJD collections", node_shape="o")

    for locality in extractions_by_locality_by_num_clones:

        background_color=TDJ_background_color
        if ("Original" in locality) or ("BMW" in locality):
            background_color = BMW_background_color

        for num_clones in extractions_by_locality_by_num_clones[locality]:
            print("num clones: " + str(num_clones))
            print("locality: " + locality)
            nodes = extractions_by_locality_by_num_clones[locality][num_clones]

            print("nodes: " + str(nodes))
            color = colors_by_locality[locality]
            print("color: " + str(color))
            colors_repeated= [color] * len(nodes)

            node_sizes = []
            background_sizes=[]
            for node in nodes:
                node_sizes.append( node_size*(1+num_clones))
                background_sizes.append( node_size*(1+num_clones)*1.2)

            if num_clones == 0:

                #draw backgound
                nx.draw_networkx_nodes(G, pos, node_size=background_sizes, node_color=background_color,
                                           nodelist=nodes, label=None, node_shape="o")

                #draw real thing
                nx.draw_networkx_nodes(G, pos, node_size=node_sizes, node_color=colors_repeated,
                               nodelist=nodes, label=locality,node_shape = "p")


            else:

                #draw backgound
                nx.draw_networkx_nodes(G, pos, node_size=background_sizes, node_color=background_color,
                                           nodelist=nodes, label=None, node_shape="o")

                if not(0 in extractions_by_locality_by_num_clones[locality]):
                    nx.draw_networkx_nodes(G, pos, node_size=node_size, node_color=colors_repeated,
                               nodelist=nodes, label=locality, node_shape = "p")

                #draw real thing
                nx.draw_networkx_nodes(G, pos, node_size=node_sizes, node_color=colors_repeated,
                               nodelist=nodes,node_shape = "p")

    # edges
    edge_collection1 = nx.draw_networkx_edges(G, pos, edgelist=elarge, edge_color="k", width=4,
                                              label='AS > ' + str(AS_threshold2))
    edge_collection2 = nx.draw_networkx_edges(
        G, pos, edgelist=emedium, width=2, alpha=0.5, edge_color="k", style="dashed",
        label='AS > ' + str(AS_threshold1))

    edge_collection2 = nx.draw_networkx_edges(
        G, pos, edgelist=esmall, width=1, alpha=0.2, edge_color="k", style="dashed",
        label='AS <= ' + str(AS_threshold1))

    # node labels
    nx.draw_networkx_labels(G, pos, font_size=8, font_family="sans-serif")
    # edge weight labels
    # edge_labels = nx.get_edge_attributes(G, "weight")
    edge_labels = nx.get_edge_attributes(G, "label")
    nx.draw_networkx_edge_labels(G, pos, edge_labels, font_size=8, alpha=1.0)

    ax = plt.gca()
    ax.margins(0.08)
    plt.axis("off")
    # plt.tight_layout()

    if species:
        file_base_name = file_base_name + "_" + species
    if locality:
        file_base_name = file_base_name + "_" + locality

    file_name = os.path.join(out_dir, file_base_name + "_network.png")
    if output_file:
        file_name = os.path.join(out_dir, output_file)

    if os.path.exists(file_name):
        os.remove(file_name)

    # props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    # ax.text(0.05, 0.05, "textstr", transform=ax.transAxes, fontsize=14,
    #        verticalalignment='top', bbox=props)

    plt.title(file_base_name + " Genotype Network")
    if plot_title:
        plt.title(plot_title)

    if legend_placement:
        plt.legend(loc=legend_placement, prop={'size': 12})
    else:
        plt.legend(loc="upper right", prop={'size': 12})
    plt.savefig(file_name)
    # plt.show()

    found_clones=False
    for extraction in clones_by_extraction:
        clones_for_extraction=clones_by_extraction[extraction]
        if len(clones_for_extraction) > 0:
            print(extraction + " clones: " + str(clones_for_extraction))
            found_clones= True

    if not found_clones:
        print("No clones found in this dataset")
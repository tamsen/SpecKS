import os
import unittest

import networkx as nx
from Bio import Phylo
from matplotlib import pyplot as plt
import numpy as np
from pathlib import Path
from pipeline_modules import gene_tree_maker
from io import StringIO


class VisualizationTests2(unittest.TestCase):


    def test_plot_newick24(self):
        #newick_string = "(O:500,(P1:300,P2:300):200);"

        #SHOULD HAVE 13 NODES AND 12 EDGES
        newick_string =get_gt1()
        tree = Phylo.read(StringIO(newick_string), "newick")
        #leaves_to_prune = ["G0_0", "G1_0", "G2_0"]
        #for leaf in leaves_to_prune:
        #    tree.prune(leaf)
        Phylo.draw_ascii(tree)

        file_to_save = os.path.join("./test_out", "species42.png")
        fig, ax = plt.subplots()

        X = Phylo.to_networkx(tree)
        print(list(X.nodes))
        print(list(X.edges))

        nodes=list(X.nodes)
        edges=list(X.edges)
        #self.assertEqual(len(nodes), 13) #now,7
        #self.assertEqual(len(edges), 12)


        node_distances_x={}
        node_distances_y={}
        node_names_by_i={}
        node_i_by_name={}
        for i in range(0,len(nodes)):
            node_distances_y[i]=tree.distance(nodes[i])
            if nodes[i].name:
                node_names_by_i[i] = nodes[i].name
                node_i_by_name[nodes[i].name] = i

        leaf_map= gene_tree_maker.read_leaf_map(
            "GeneTree0.test.leafmap", gene_tree_maker.gene_tree_result()).leaves_by_species

        leaf_map_by_leaf={}
        for species in leaf_map:
            for leaf in leaf_map[species]:
                leaf_map_by_leaf[leaf]=species

        nodes_to_visulize=[]
        for i in range(0, len(nodes)):
            leafs = nodes[i].get_terminals()
            names=[leaf.name for leaf in leafs]
            species=[leaf_map_by_leaf[name] for name in names]
            print("species" + str(species))
            if "O" in species:
                f=0
            elif "P1" in species:
                f=.3
                nodes_to_visulize.append(i)
            elif "P2" in species:
                f = -.3
                nodes_to_visulize.append(i)
            else:
                print("bad")
                print(species)
            (b, num_sibs) =birth_order(tree, nodes[i])
            delta= b*100.0/num_sibs
            print(names)
            node_distances_x[i] = tree.distance(nodes[i])*f+delta

        new_edge_list=[]
        for edge in edges:
            n1=edge[0]
            n2=edge[1]
            index_1=node_i_by_name[n1.name]
            index_2=node_i_by_name[n2.name]
            new_edge_list.append((index_1,index_2))
            #TODO - check which species it is in before adding it to the edge list..

        #self.assertEqual(len( node_distances_y), 13)
        #self.assertEqual(len( node_distances_x), 13)

        pos={}
        for i in node_distances_y.keys():

            name="no name"+ "_i" +str(i)
            if i in node_names_by_i:
                name=node_names_by_i[i] + "_i" +str(i)

            print("{0}:{1}\t({2},{3})".format(
                i,name,node_distances_x[i],node_distances_y[i]))
            plt.plot(node_distances_x[i],node_distances_y[i])

            plt.text(node_distances_x[i],node_distances_y[i],name)
            pos[i]=(node_distances_x[i],node_distances_y[i])

        X.add_edges_from(new_edge_list)

        #todo - dont visualize the outgroup. last time this worked with nodelist= [1,2,3]
        #short_pos_list={i:pos[i] for i in [1,2,4]}
        for edge in new_edge_list:
            nx.draw_networkx_edges(X, pos, new_edge_list)
        #    nx.draw_networkx_edges(X, pos, new_edge_list) #nodes_to_visulize)

        ax.tick_params(left=True, bottom=True, labelleft=True, labelbottom=False)
        #plt.xlim(-50, 50)
        plt.ylim(0, 500 + 50)
        y_ticks = plt.yticks()[0]
        new_ticks = [500 - y_tick for y_tick in y_ticks]
        num_ticks = len(y_ticks)
        plt.title('polyploid.species_name')
        # reverse, since back in time
        ax.set_yticks(y_ticks[0:num_ticks - 1])
        ax.set_yticklabels(new_ticks[0:num_ticks - 1])
        #plt.legend()
        plt.ylabel("MYA")
        plt.savefig(file_to_save)
        plt.cla()
        plt.close()

    def get_node_distances(self, i, node_points, cxns, clade_to_id,
                           recursion_depth, root, tree):

        start_i=i
        recursion_depth = recursion_depth + 1
        index_at_level = 0
        parent_x=node_points[start_i][0]

        min = -10
        max = 10
        num_needed = len(root.clades) - 1
        step = round((max - min) / num_needed, 2)
        xs = np.arange(min, max + 1, step)
        print(xs)
        xs_index = 0

        for node in root.clades:  # X.nodes:
            i=i+1
            print("Node:" + str(node))
            if node.branch_length:
                print(node.branch_length)
                print("distance:" + str(tree.distance(node)))
                node_points[i] = [parent_x+xs[index_at_level],tree.distance(node)]
                clade_to_id[i] = node
                if start_i > 0:
                    cxns.append((start_i,i))
                for child in node.clades:
                    self.get_node_distances(i, node_points, cxns, clade_to_id,
                                            recursion_depth, child, tree)

                index_at_level = index_at_level + 1

        return node_points,cxns, clade_to_id


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
            #if clade.name == "O":
            #    continue
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

def get_gt1():
    gt1="((G0_0:385.45802[ID=0 HOST=0 GUEST=null VERTEXTYPE=Leaf],(G1_0:303.7403[ID=1 HOST=0 GUEST=null VERTEXTYPE=Leaf],G2_0:303.7403[ID=2 HOST=0 GUEST=null VERTEXTYPE=Leaf])G3_0:81.7177139711581[ID=3 HOST=0 GUEST=null VERTEXTYPE=Duplication DISCPT=(1,1)])G4_0:114.54198363837152[ID=4 HOST=0 GUEST=null VERTEXTYPE=Duplication DISCPT=(1,1)],((G5_1:243.39763[ID=5 HOST=1 GUEST=null VERTEXTYPE=Leaf],G6_1:243.39763[ID=6 HOST=1 GUEST=null VERTEXTYPE=Leaf])G7_1:56.60236874941369[ID=7 HOST=1 GUEST=null VERTEXTYPE=Duplication DISCPT=(0,2)],(G8_2:290.54914[ID=8 HOST=2 GUEST=null VERTEXTYPE=Leaf],G9_2:290.54914[ID=9 HOST=2 GUEST=null VERTEXTYPE=Leaf])G10_2:9.45085528985734[ID=10 HOST=2 GUEST=null VERTEXTYPE=Duplication DISCPT=(0,2)])G11_3:200.0[ID=11 HOST=3 GUEST=null VERTEXTYPE=Speciation DISCPT=(1,0)])G12_4:0.0[ID=12 HOST=4 GUEST=null VERTEXTYPE=Speciation DISCPT=(2,0)][NAME=PrunedTree];"
    return gt1

def get_gt2():
    gt2="(G0_0:500.0[ID=0 HOST=0 GUEST=null VERTEXTYPE=Leaf],G1_1:500.0[ID=1 HOST=1 GUEST=null VERTEXTYPE=Leaf])G2_4:0.0[ID=2 HOST=4 GUEST=null VERTEXTYPE=Speciation DISCPT=(2,0)][NAME=PrunedTree];"
    return gt2

def get_gt3():
    gt3="((((G0_0:62.147061[ID=0 HOST=0 GUEST=null VERTEXTYPE=Leaf],G1_0:62.147061[ID=1 HOST=0 GUEST=null VERTEXTYPE=Leaf])G2_0:120.57833966761288[ID=2 HOST=0 GUEST=null VERTEXTYPE=Duplication DISCPT=(0,1)],G3_0:182.7254[ID=3 HOST=0 GUEST=null VERTEXTYPE=Leaf])G4_0:211.12567[ID=4 HOST=0 GUEST=null VERTEXTYPE=Duplication DISCPT=(0,2)],((G5_0:200.72794[ID=5 HOST=0 GUEST=null VERTEXTYPE=Leaf],G6_0:200.72794[ID=6 HOST=0 GUEST=null VERTEXTYPE=Leaf])G7_0:173.17475[ID=7 HOST=0 GUEST=null VERTEXTYPE=Duplication DISCPT=(0,2)],(((G8_0:58.077819[ID=8 HOST=0 GUEST=null VERTEXTYPE=Leaf],G9_0:58.077819[ID=9 HOST=0 GUEST=null VERTEXTYPE=Leaf])G10_0:94.49669998269559[ID=10 HOST=0 GUEST=null VERTEXTYPE=Duplication DISCPT=(0,1)],G11_0:152.57452[ID=11 HOST=0 GUEST=null VERTEXTYPE=Leaf])G12_0:29.200733628117394[ID=12 HOST=0 GUEST=null VERTEXTYPE=Duplication DISCPT=(0,2)],G13_0:181.77525[ID=13 HOST=0 GUEST=null VERTEXTYPE=Leaf])G14_0:192.12743578315605[ID=14 HOST=0 GUEST=null VERTEXTYPE=Duplication DISCPT=(0,2)])G15_0:19.948383788675812[ID=15 HOST=0 GUEST=null VERTEXTYPE=Duplication DISCPT=(1,1)])G16_0:106.14893[ID=16 HOST=0 GUEST=null VERTEXTYPE=Duplication DISCPT=(1,1)],(G17_2:118.08417[ID=17 HOST=2 GUEST=null VERTEXTYPE=Leaf],G18_2:118.08417[ID=18 HOST=2 GUEST=null VERTEXTYPE=Leaf])G19_2:381.91583[ID=19 HOST=2 GUEST=null VERTEXTYPE=Duplication DISCPT=(0,1)])G20_4:0.0[ID=20 HOST=4 GUEST=null VERTEXTYPE=Speciation DISCPT=(2,0)][NAME=PrunedTree];"
    return gt3
def get_config_file():
    par_dir = Path(__file__).parent.parent
    config_file = os.path.join(par_dir, "sample_configs","cos-bi-config.xml")
    return config_file

def get_parent(tree, child_clade):
        node_path = tree.get_path(child_clade)
        if len(node_path)>2:
            return node_path[-2]
        else:
            return False
def get_sibs(tree, child_clade):
       parent = get_parent(tree, child_clade)
       if not parent:
           return []
       sibs=parent.clades
       return sibs
def birth_order(tree, child_clade):
    sibs = get_sibs(tree, child_clade)
    num_sibs=len(sibs)
    if num_sibs < 2:
        return (0,1)
    birth_order=sibs.index(child_clade)
    return (birth_order,num_sibs)

if __name__ == '__main__':
    unittest.main()

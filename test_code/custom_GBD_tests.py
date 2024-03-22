import unittest
from io import StringIO
from Bio import Phylo

from pipeline_modules import custom_GBD_maker


class TestCustomGBD(unittest.TestCase):
    def test_custom_GBD_model(self):
        base_gene_tree_newick = '(O:100, (P1:73.36, P2:73.36): 26.64);'
        # custom_GBD_maker.add_GBD_to_GT(base_gene_tree_newick, 1)

        tree = Phylo.read(StringIO(base_gene_tree_newick), "newick")
        print("start tree")
        Phylo.draw_ascii(tree)
        print("newick")
        print(base_gene_tree_newick)

        nodes_to_add_by_branch_name = {}
        internal_node_idx = 0
        duplicate_idx = 0
        self.iterate_though_clades(tree.clade, nodes_to_add_by_branch_name, internal_node_idx, duplicate_idx)


        new_branch_inx=0
        #tree.prune("P1")
        duplicate_idx=0

        for branch_name, branches_to_add in nodes_to_add_by_branch_name.items():

            if len(branches_to_add)==0:
                continue

            for new_branch_data in branches_to_add:
                self.split_branch_with_this_name(branch_name, internal_node_idx,
                                                 new_branch_data, tree.clade.clades)

            #new_n=self.write_new_newick(tree)
            #print("branch:\t" +  branch_name)
            #print("new_newick:\t" + new_n)

        print("new tree")
        Phylo.draw_ascii(tree)
        self.assertEqual(True, False)  # add assertion here

    def split_branch_with_this_name(self, branch_name, internal_node_idx, new_branch_data, clades):

        for c in clades:

            if c.name == branch_name:

                print("found "+ branch_name)
                new_branch_length = new_branch_data.relative_start_time
                orignal_branch_length=c.branch_length

                #https://biopython.org/docs/1.75/api/Bio.Phylo.BaseTree.html
                c.split(n=2, branch_length=orignal_branch_length-new_branch_data.relative_start_time)
                c.branch_length=new_branch_length
                c.name="internal_branch_" + str(internal_node_idx)
                c.clades[0].name=branch_name
                c.clades[1].name = new_branch_data.new_branch_name
                print("at distance " + str(c.branch_length))
                print("adding "+ c.clades[1].name)
                print("preserving " + c.clades[0].name)
                internal_node_idx=internal_node_idx+1
                return

            self.split_branch_with_this_name(branch_name, internal_node_idx,  new_branch_data, c.clades)



    # if death rate = 0.001359 gene per MY, thats a death rate of 1 per 700 MY
    #
    def get_nodes_to_add_to_branch(self, parent_clade, duplicate_idx):
        length = parent_clade.branch_length
        parent_branch_name = parent_clade.name
        spaces_between_nodes = 50
        first_space = spaces_between_nodes / 2.0
        node_life_span = 30
        num_nodes_to_add = int(length / spaces_between_nodes)
        distances_between_nodes_to_add = [spaces_between_nodes for i in range(0, num_nodes_to_add)]
        distances_between_nodes_to_add = [first_space] + distances_between_nodes_to_add

        list_of_nodes_to_add = []
        distance_traveled_along_branch = 0
        for i in range(0, num_nodes_to_add):
            # start_pos=distance_traveled_along_branch+ distances_between_nodes_to_add[i]
            # end_pos=start_pos + node_life_span
            relative_start_pos = distances_between_nodes_to_add[i]
            relative_end_pos = relative_start_pos + node_life_span
            new_branch_name=parent_branch_name+ "_duplicate_" + str(duplicate_idx)
            new_node = node_to_add(parent_branch_name, new_branch_name, relative_start_pos, relative_end_pos)
            list_of_nodes_to_add.append(new_node)
            distance_traveled_along_branch = distance_traveled_along_branch + relative_start_pos
            duplicate_idx= duplicate_idx+1

        return list_of_nodes_to_add

    def write_new_newick(self, tree):
        handle = StringIO()
        Phylo.write(tree, handle, "newick")
        new_newick = handle.getvalue()
        return new_newick

    def iterate_though_clades(self, clade, nodes_to_add_by_branch_name, internal_node_idx, duplicate_idx):

        for c in clade.clades:

            if not c.name:
                internal_node_name = "internal_node_" + str(internal_node_idx)
                c.name = internal_node_name
                internal_node_idx = internal_node_idx + 1


            print(c.name + ":\t" + str(c.branch_length))
            nodes_to_add_to_branch = self.get_nodes_to_add_to_branch(c,duplicate_idx)
            print("nodes_to_add_to_branch:" + str(nodes_to_add_to_branch))
            nodes_to_add_by_branch_name[c.name] = nodes_to_add_to_branch
            self.iterate_though_clades(c, nodes_to_add_by_branch_name, internal_node_idx,duplicate_idx)


if __name__ == '__main__':
    unittest.main()


class node_to_add():
    parent_branch_name = ""
    new_branch_name = ""
    relative_start_time = 0
    relative_end_time = 0

    def __init__(self, parent_branch_name, new_branch_name, relative_start_time, relative_end_time):
        self.parent_branch_name = parent_branch_name
        self.relative_start_time = relative_start_time
        self.relative_end_time = relative_end_time
        self.new_branch_name = new_branch_name

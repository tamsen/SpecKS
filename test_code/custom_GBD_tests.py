import random
import unittest
from io import StringIO
from Bio import Phylo
from scipy.stats import lognorm, expon, poisson
from pipeline_modules import custom_GBD_maker


class TestCustomGBD(unittest.TestCase):
    def test_custom_GBD_model(self):
        base_gene_tree_newick = '(O:100, (P1:73.36, P2:73.36): 26.64);'
        
        #Gene family evolution in green plants with emphasis on the origination and evolution of Arabidopsis thaliana genes
        #gene_birth_rate=0.001359

        gene_birth_rate = 0.1359
        mean_time_between_gene_births= int( 1.0/gene_birth_rate )


        #https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4672517/
        #https: // getd.libs.uga.edu / pdfs / yan_zhen_201605_phd.pdf
        #The Exponential Distribution is the time between events in a Poisson process.
        # If the number of occurrences follows a Poisson distribution,
        # the lapse of time between these events is distributed exponentially.
        distribution_fxn=expon

        # custom_GBD_maker.add_GBD_to_GT(base_gene_tree_newick, 1)

        tree = Phylo.read(StringIO(base_gene_tree_newick), "newick")
        print("start tree")
        Phylo.draw_ascii(tree)
        print("newick")
        print(base_gene_tree_newick)

        total_tree_length=tree.total_branch_length()
        max_number_gene_births=int(total_tree_length*gene_birth_rate*1.5)+10 # the one x1.5 is just to make sure

        #TODO, print out these distributions and make sure they really do what I want

        #SSD_life_span_halflife = 1 #MY
        # => 0.5A=Ae^kt, at t=1MY
        # kt=ln(1/2), set t=1
        #=k = ln(0.5)t = -0.6931
        #so, expon shape_factor= 0.6731
        #SSD_life_spans=expon.rvs(0.6731,size=max_number_gene_births,random_state=4)
        SSD_life_spans = expon.rvs(10, size=max_number_gene_births, random_state=4)
        #SSD_time_between_gene_birth_events=expon.rvs(mean_time_between_gene_births,
        #                                             size=max_number_gene_births,random_state=4)

        SSD_life_spans=[30 for r in range(0,max_number_gene_births)]
        SSD_time_between_gene_birth_events=[10 for r in range(0,max_number_gene_births)]

        nodes_to_add_by_branch_name = {}
        internal_node_idx = 0
        duplicate_idx = 0
        randomness_idx = 0
        end_of_sim=100
        self.recursively_get_new_branches_to_add(tree,tree.clade,
                                                 nodes_to_add_by_branch_name,SSD_time_between_gene_birth_events,
                                                 SSD_life_spans, internal_node_idx, duplicate_idx,randomness_idx)

        print("")
        print("Genes born:")
        print_dict_of_node_to_add(nodes_to_add_by_branch_name)

        retained_nodes_by_branch_name,deleted_nodes_by_branch_name=(
            self.prune_any_branches_that_would_be_dead_before_the_end_of_the_sim(nodes_to_add_by_branch_name, end_of_sim))
        print("")
        print("Genes died:")
        print_dict_of_node_to_add(deleted_nodes_by_branch_name)
        print("")

        print("adding the new genes.. ")
        for branch_name, branches_to_add in retained_nodes_by_branch_name.items():

            if branch_name=='O':
                branches_to_add= []
            else:
                branches_to_add= [branches_to_add[0],branches_to_add[1]]
            print("\nFor branch "  + branch_name)
            for new_branch_data in branches_to_add:
                self.recursively_split_branch_with_this_name(tree, branch_name, internal_node_idx,
                                                             new_branch_data, tree.clade.clades)

        print("new tree")
        Phylo.draw_ascii(tree)
        #self.assertEqual(True, False)  # add assertion here

        clades=tree.clade.clades
        for c in clades:
                if c.name == "O":
                    print("length of O = " + str(tree.distance(c)))

                    self.assertEqual(100, tree.distance(c))  # add assertion here
    def recursively_split_branch_with_this_name(self, full_tree, branch_name, internal_node_idx, new_branch_data, clades):

        for c in clades:

            if c.name == branch_name:

                print("\n\nfound "+ branch_name)
                print("new branch data: ")
                new_branch_data.print_data()
                original_branch_length=c.branch_length
                original_parent_branch_end=full_tree.distance(c)
                original_parent_branch_start=full_tree.distance(c)-original_branch_length

                new_branch_starts_here=new_branch_data.absolute_start_time
                relative_new_branch_start=new_branch_starts_here-original_parent_branch_start


                print("original parent branch length " + str(original_branch_length))
                print("abs original parent branch start " + str(original_parent_branch_start)) #shoudl be ~26 for P2
                print("relative start to new child branch " + str(new_branch_data.relative_start_time))
                print("abs new child branch start (calculated)" + str(original_parent_branch_start+new_branch_data.relative_start_time))
                print("abs new child branch start (saved)" + str(new_branch_data.absolute_start_time))

                child_branch_length=original_parent_branch_end-new_branch_data.absolute_start_time
                print("child branch length (calculated from abs) " + str(child_branch_length))
                print("abs new child branch start (saved)" + str(new_branch_data.absolute_start_time))

                #https://biopython.org/docs/1.75/api/Bio.Phylo.BaseTree.html
                #print("where_to_place_split (abs coords) " + str(new_branch_starts_here))


                #c.split:
                #New clades have the given branch_length and the same name as this clade’s
                # root plus an integer suffix (counting from 0). For example, splitting
                # a clade named “A” produces sub-clades named “A0” and “A1”.
                # If the clade has no name, the prefix “n” is used for child nodes, e.g. “n0” and “n1”.

                c.split(n=2, branch_length=child_branch_length) #new_branch_data.relative_start_time)
                c.branch_length=original_branch_length-child_branch_length

                c.name="internal_branch_" + str(internal_node_idx)
                c.clades[0].name=branch_name
                c.clades[1].name = new_branch_data.new_branch_name

                print("child0 branch length after splitting " + str(c[0].branch_length))
                print("child1 branch length after splitting " + str(c[1].branch_length))
                print("adding "+ c.clades[1].name)
                print("preserving " + c.clades[0].name + " but its start pos has changed")
                internal_node_idx=internal_node_idx+1

                #add distave since last duplication gene to new_branch_data
                return internal_node_idx

            self.recursively_split_branch_with_this_name(full_tree, branch_name, internal_node_idx, new_branch_data, c.clades)



    # if death rate = 0.001359 gene per MY, thats a death rate of 1 per 700 MY
    #

    def get_nodes_to_add_to_branch(self, parents_distance_from_root, child_clade,
                                   SSD_time_between_gene_birth_events,
                                   SSD_life_spans, duplicate_idx,randomness_idx):

        child_branch_length = child_clade.branch_length
        child_branch_name = child_clade.name
        jump_to_next_event=0
        list_of_nodes_to_add = []

        #how long since the last gene birth for this gene family
        #distance_traveled_along_branch = random.randint(0, int(SSD_time_between_gene_birth_events[randomness_idx]))
        distance_traveled_along_branch = 15
        print("random start:\t" + str( distance_traveled_along_branch))

        while distance_traveled_along_branch < child_branch_length:

            relative_start_pos = distance_traveled_along_branch
            relative_end_pos = relative_start_pos + SSD_life_spans[randomness_idx]
            absolute_end_time = parents_distance_from_root + relative_end_pos
            absolute_start_time = parents_distance_from_root + relative_start_pos
            new_branch_name=child_branch_name+ "_duplicate_" + str(duplicate_idx)
            new_node = node_to_add(child_branch_name, new_branch_name,
                                   relative_start_pos, absolute_end_time, absolute_start_time, jump_to_next_event)
            list_of_nodes_to_add.append(new_node)
            new_node.print_data()
            jump_to_next_event=SSD_time_between_gene_birth_events[randomness_idx]
            distance_traveled_along_branch = distance_traveled_along_branch + jump_to_next_event
            duplicate_idx= duplicate_idx+1
            randomness_idx = randomness_idx +1

        return list_of_nodes_to_add, duplicate_idx, randomness_idx

    def write_new_newick(self, tree):
        handle = StringIO()
        Phylo.write(tree, handle, "newick")
        new_newick = handle.getvalue()
        return new_newick

    def prune_any_branches_that_would_be_dead_before_the_end_of_the_sim(self,
                                                                        nodes_to_add_by_branch_name,
                                                                        end_of_sim):

        print("pruning duplicates")
        print("\t")
        deleted_nodes_by_branch_name={}
        retained_nodes_by_branch_name={}
        branches=nodes_to_add_by_branch_name.keys()
        for parent_branch, list_of_nodes_to_add_by_branch in nodes_to_add_by_branch_name.items():

            for duplicate_branch in list_of_nodes_to_add_by_branch:
                print("duplicate_branch_name:\t" + duplicate_branch.new_branch_name)
                print("relative_start_time:\t" + str(duplicate_branch.relative_start_time))
                print("absolute_end_time:\t" + str(duplicate_branch.absolute_end_time))
                print("sim_end_time:\t" + str(end_of_sim))


                if end_of_sim > duplicate_branch.absolute_end_time:
                    print("killing branch.\t")
                    print("\t")
                    if parent_branch not in deleted_nodes_by_branch_name:
                        deleted_nodes_by_branch_name[parent_branch]=[]
                    deleted_nodes_by_branch_name[parent_branch].append(duplicate_branch)
                else:
                    print("keeping branch.\t")
                    print("\t")
                    if parent_branch not in retained_nodes_by_branch_name:
                        retained_nodes_by_branch_name[parent_branch]=[]
                    retained_nodes_by_branch_name[parent_branch].append(duplicate_branch)

        return retained_nodes_by_branch_name,deleted_nodes_by_branch_name

    def recursively_get_new_branches_to_add(self, tree, parent_clade, nodes_to_add_by_branch_name,
                                            SSD_time_between_gene_birth_events, SSD_life_spans,
                                            internal_node_idx, duplicate_idx, randomness_idx):

        parents_distance_from_root = tree.distance(parent_clade)
        for child_clade in parent_clade.clades:

            if not child_clade.name:
                internal_node_name = "internal_node_" + str(internal_node_idx)
                child_clade.name = internal_node_name
                internal_node_idx = internal_node_idx + 1


            print("branch " + child_clade.name + ":\t length of " + str(child_clade.branch_length))

            nodes_to_add_to_branch, duplicate_idx, randomness_idx = self.get_nodes_to_add_to_branch(
                parents_distance_from_root, child_clade,
                SSD_time_between_gene_birth_events,SSD_life_spans,
                duplicate_idx,randomness_idx)
            print("nodes_to_add_to_branch:" + str(nodes_to_add_to_branch))
            nodes_to_add_by_branch_name[child_clade.name] = nodes_to_add_to_branch
            self.recursively_get_new_branches_to_add(tree,child_clade, nodes_to_add_by_branch_name,
                                                     SSD_time_between_gene_birth_events,
                                                     SSD_life_spans,internal_node_idx, duplicate_idx,randomness_idx)


if __name__ == '__main__':
    unittest.main()


class node_to_add():
    parent_branch_name = ""
    new_branch_name = ""
    relative_start_time = 0
    absolute_end_time = 0
    time_between_splits = 0
    absolute_start_time = 0
    def __init__(self, parent_branch_name, new_branch_name,
                 relative_start_time, absolute_end_time,
                 absolute_start_time,time_between_splits):
        self.parent_branch_name = parent_branch_name
        self.relative_start_time = relative_start_time
        self.absolute_end_time = absolute_end_time
        self.absolute_start_time = absolute_start_time
        self.new_branch_name = new_branch_name
        self.time_between_splits = time_between_splits


    def print_data(self):
        print("~ data for new node " + self.new_branch_name + " ~")
        print("rel start pos:\t" + str(self.relative_start_time))
        print("abs end time:\t" + str(self.absolute_end_time))
        print("jump:\t" + str(self.time_between_splits))

def print_dict_of_node_to_add(dict_of_node_to_add):

    for branch_name, nodes_to_add in dict_of_node_to_add.items():
        print("")
        print(branch_name+":")

        for node in nodes_to_add:
            node.print_data()
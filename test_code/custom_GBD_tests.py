import unittest
from io import StringIO
from Bio import Phylo
from scipy.stats import lognorm, expon, poisson
from pipeline_modules import custom_GBD_maker


class TestCustomGBD(unittest.TestCase):
    def test_custom_GBD_model(self):
        base_gene_tree_newick = '(O:100, (P1:73.36, P2:73.36): 26.64);'
        gene_birth_rate=0.001359
        gene_death_rate=0.001359
        mean_time_between_gene_births=1.0/gene_birth_rate
        mean_time_between_gene_deaths=1.0/gene_death_rate

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
        num_gene_birth_intervals_needed=int(total_tree_length*gene_birth_rate*1.5)+10 # the one x1.5 is just to make sure
        num_gene_death_intervals_needed=int(total_tree_length*gene_death_rate*1.5)+10# we always have extra by a good margin
        random_start_times =random.sample(range(0, mean_time_between_gene_births), num_gene_birth_intervals_needed)
        #scale = 1 mode.
        random_intervals_between_gene_births=distribution_fxn.rvs(mean_time_between_gene_births,
                                                                               size=num_gene_birth_intervals_needed)
        random_intervals_between_gene_deaths=distribution_fxn.rvs(mean_time_between_gene_deaths,
                                                                               size=num_gene_death_intervals_needed)

        reservoirs_of_randomness=[random_intervals_between_gene_births,random_intervals_between_gene_deaths, random_start_times]
        nodes_to_add_by_branch_name = {}
        internal_node_idx = 0
        duplicate_idx = 0

        end_of_sim=100
        self.recursively_get_new_branches_to_add(tree.clade,
                                                 nodes_to_add_by_branch_name,
                                                 reservoirs_of_randomness, internal_node_idx, duplicate_idx)
        print("Genes born:" + str(nodes_to_add_by_branch_name))

        retained_nodes_by_branch_name,deleted_nodes_by_branch_name=(
            self.prune_any_branches_that_would_be_dead_before_the_end_of_the_sim(nodes_to_add_by_branch_name, end_of_sim))
        print("Genes died:" + str(deleted_nodes_by_branch_name))

        for branch_name, branches_to_add in retained_nodes_by_branch_name.items():

            if len(branches_to_add)==0:
                continue

            for new_branch_data in branches_to_add:
                self.split_branch_with_this_name(branch_name, internal_node_idx,
                                                 new_branch_data, tree.clade.clades)

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
    def get_nodes_to_add_to_branch(self, parent_clade, parent_distance_from_root, reservoirs_of_randomness, duplicate_idx):

        randomness_idx=5
        [gb_intervals,gd_intervals,start_times]=reservoirs_of_randomness
        length = parent_clade.branch_length
        parent_branch_name = parent_clade.name

        #distances_between_nodes_to_add = [spaces_between_nodes for i in range(0, num_nodes_to_add)]
        #distances_between_nodes_to_add = [first_space] + distances_between_nodes_to_add

        list_of_nodes_to_add = []
        distance_traveled_along_branch = 0
        last_start_pos=start_times[randomness_idx]
        #for i in range(0, num_nodes_to_add):
        while distance_traveled_along_branch < length:
            # start_pos=distance_traveled_along_branch+ distances_between_nodes_to_add[i]
            # end_pos=start_pos + node_life_span
            relative_start_pos = last_start_pos + gb_intervals[randomness_idx]]
            last_start_pos = relative_start_pos
            relative_end_pos = relative_start_pos + gd_intervals[randomness_idx]
            absolute_end_time = parent_distance_from_root + relative_end_pos
            new_branch_name=parent_branch_name+ "_duplicate_" + str(duplicate_idx)
            new_node = node_to_add(parent_branch_name, new_branch_name, relative_start_pos, absolute_end_time)
            list_of_nodes_to_add.append(new_node)
            distance_traveled_along_branch = distance_traveled_along_branch + relative_start_pos
            duplicate_idx= duplicate_idx+1

        return list_of_nodes_to_add, duplicate_idx

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

    def recursively_get_new_branches_to_add(self, clade, nodes_to_add_by_branch_name,
                                            reservoirs_of_randomness, internal_node_idx, duplicate_idx):

        for c in clade.clades:

            if not c.name:
                internal_node_name = "internal_node_" + str(internal_node_idx)
                c.name = internal_node_name
                internal_node_idx = internal_node_idx + 1


            print(c.name + ":\t" + str(c.branch_length))
            nodes_to_add_to_branch, duplicate_idx = self.get_nodes_to_add_to_branch(c,
                                                                                    reservoirs_of_randomness,
                                                                                    duplicate_idx)
            print("nodes_to_add_to_branch:" + str(nodes_to_add_to_branch))
            nodes_to_add_by_branch_name[c.name] = nodes_to_add_to_branch
            self.recursively_get_new_branches_to_add(c, nodes_to_add_by_branch_name,
                                                     reservoirs_of_randomness,internal_node_idx, duplicate_idx)


if __name__ == '__main__':
    unittest.main()


class node_to_add():
    parent_branch_name = ""
    new_branch_name = ""
    relative_start_time = 0
    absolute_end_time = 0

    def __init__(self, parent_branch_name, new_branch_name, relative_start_time, absolute_end_time):
        self.parent_branch_name = parent_branch_name
        self.relative_start_time = relative_start_time
        self.absolute_end_time = absolute_end_time
        self.new_branch_name = new_branch_name

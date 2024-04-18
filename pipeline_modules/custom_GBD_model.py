
import os
import random
import numpy as np
from matplotlib import pyplot as plt
from scipy.stats import lognorm, expon, poisson
from Bio import Phylo
from io import StringIO
import log
from pipeline_modules import gene_tree_maker
from pipeline_modules.gene_tree_info import custom_gene_tree_result

def run_custom_GBD_model(polyploid, all_genomes_of_interest, simulation_leg, base_gene_tree_newicks_by_tree_name):

    config = polyploid.general_sim_config
    num_gene_trees_needed = config.num_gene_trees_per_species_tree
    include_visualizations = config.include_visualizations
    length_of_leg = simulation_leg.interval_end_time_MY-simulation_leg.interval_start_time_MY
    gt_index_formatter = gene_tree_maker.get_gt_index_format(num_gene_trees_needed)
    subfolder = os.path.join(polyploid.species_subfolder, str(polyploid.analysis_step_num) + "_custom_GBD")

    if not os.path.exists(subfolder):
        os.makedirs(subfolder)


    num_values_needed=config.num_gene_trees_per_species_tree * config.num_replicates_per_gene_tree
    SSD_life_spans, SSD_time_between_gene_birth_events, skip_GBD_model = set_up(config,subfolder,
                                                                                num_values_needed)
    gene_tree_data_by_tree_name = {}
    total_num_dup_added=0
    total_num_dup_shed=0
    randomness_idx=idx_by_reference(0)
    for gt_idx in range(0, num_gene_trees_needed):


        gene_tree_name = "GeneTree" + gt_index_formatter.format(gt_idx)
        base_gene_tree_newick = base_gene_tree_newicks_by_tree_name[gene_tree_name]
        if skip_GBD_model:
            gene_tree_newick_with_GBD = base_gene_tree_newick
        else:
            gene_tree_newick_with_GBD, num_dup_added, num_dup_shed, randomness_idx = add_GBD_to_newick(base_gene_tree_newick,
                                                      gene_tree_name,
                                                      gt_idx,randomness_idx,
                                                      SSD_life_spans,SSD_time_between_gene_birth_events,
                                                      length_of_leg,subfolder,include_visualizations)
            total_num_dup_added=total_num_dup_added+num_dup_added
            total_num_dup_shed=total_num_dup_shed+num_dup_shed
        gene_tree_data=custom_gene_tree_result(gene_tree_name,gene_tree_newick_with_GBD,all_genomes_of_interest)
        gene_tree_data_by_tree_name[gene_tree_name] = gene_tree_data

    log.write_to_log("Overall, " +
                     str(total_num_dup_added + total_num_dup_shed) +
                     " total num SSD duplicates occurred over time ")
    log.write_to_log(str(total_num_dup_shed) + " were shed ")
    log.write_to_log(str(total_num_dup_added) + " were retained (not yet shed) ")

    polyploid.analysis_step_num = polyploid.analysis_step_num + 1
    return gene_tree_data_by_tree_name


def set_up(config, subfolder,num_gene_trees_needed):

    include_visualizations = config.include_visualizations
    mean_gene_birth_rate  = config.mean_gene_birth_rate #0.001359
    mean_SSD_life_span = config.mean_SSD_life_span #1 #
    skip_GBD_model=False

    if (not mean_gene_birth_rate) or (not mean_SSD_life_span):
        skip_GBD_model = True
        log.write_to_log("skipping GBD model")
        return [], [], skip_GBD_model
    else:
        log.write_to_log("running GBD model")
        max_possible_base_tree_length = config.full_sim_time * 2
        multiplier=max(1,mean_gene_birth_rate)
        max_number_gene_births = int(
            max_possible_base_tree_length * multiplier * 2 * num_gene_trees_needed) + 10  # the one x1.5 is just to make sure
        mean_time_between_gene_births = int(1.0 / mean_gene_birth_rate)

        # loc=mode of data (location of peak of distribution), scale=location of average of data (cm of distribution)
        SSD_life_spans = expon.rvs(loc=0, scale=mean_SSD_life_span, size=max_number_gene_births)
        SSD_time_between_gene_birth_events = poisson.rvs(mean_time_between_gene_births, size=max_number_gene_births)

        if include_visualizations:
            visualize_distribution(SSD_life_spans, "SSD_life_spans", mean_SSD_life_span, subfolder)
            visualize_distribution(SSD_time_between_gene_birth_events,
                                   "SSD_time_between_gene_birth_events", mean_time_between_gene_births, subfolder)
    return SSD_life_spans, SSD_time_between_gene_birth_events, skip_GBD_model


def add_GBD_to_newick(base_gene_tree_newick, gene_tree_name,
                      gt_idx, randomness_idx,
                      SSD_life_spans, SSD_time_between_gene_birth_events,
                      end_of_leg, outfolder, include_visuals):


        tree = Phylo.read(StringIO(base_gene_tree_newick), "newick")
        nodes_to_add_by_branch_name = {}
        internal_node_idx = 0
        duplicate_idx = 0
        tree.clade.name="root_" + str(gt_idx)

        if include_visuals:
            print(base_gene_tree_newick)
            save_ascii_tree(base_gene_tree_newick, gene_tree_name, outfolder, "_before_GBD.txt", tree)

        recursively_get_new_branches_to_add(tree, tree.clade,
                                                 nodes_to_add_by_branch_name, SSD_time_between_gene_birth_events,
                                                 SSD_life_spans, internal_node_idx, duplicate_idx, randomness_idx)

        retained_nodes_by_branch_name, deleted_nodes_by_branch_name = (
            prune_any_branches_that_would_be_dead_before_the_end_of_the_sim(nodes_to_add_by_branch_name,
                                                                            end_of_leg))


        num_dup_shed=num_genes_in_dict(deleted_nodes_by_branch_name)
        num_dup_retained=num_genes_in_dict(retained_nodes_by_branch_name)
        log.write_to_log(str(num_dup_retained)+ " were retained (not yet shed) ")

        for branch_name, branches_to_add in retained_nodes_by_branch_name.items():

            for new_branch_data in branches_to_add:
                recursively_split_branch_with_this_name(tree, branch_name, internal_node_idx,
                                                             new_branch_data, tree.clade.clades)


        gene_tree_newick_with_GBD=tree_to_newick(tree)

        if include_visuals:
            save_ascii_tree(gene_tree_newick_with_GBD, gene_tree_name, outfolder, "_after_GBD.txt", tree)

        return gene_tree_newick_with_GBD, num_dup_retained, num_dup_shed,randomness_idx


def save_ascii_tree(base_gene_tree_newick, gene_tree_name, outfolder, suffix, tree):
    tree_file_path = os.path.join(outfolder, gene_tree_name + suffix)
    with open(tree_file_path, 'w') as f:

        f.writelines("newick:" + str(base_gene_tree_newick) + "\n")

        if len(tree.get_terminals()) < 2:
            f.writelines("The tree has less than two branches, so nothing to draw.")
        else:
            Phylo.draw_ascii(tree, f)


def recursively_split_branch_with_this_name(full_tree, branch_name, internal_node_idx, new_branch_data,
                                                clades):

        for c in clades:

            if c.name == branch_name:

                original_branch_length = c.branch_length
                original_parent_branch_end = full_tree.distance(c)

                child_branch_length = original_parent_branch_end - new_branch_data.absolute_start_time

                # https://biopython.org/docs/1.75/api/Bio.Phylo.BaseTree.html
                # c.split:
                # New clades have the given branch_length and the same name as this clade’s
                # root plus an integer suffix (counting from 0). For example, splitting
                # a clade named “A” produces sub-clades named “A0” and “A1”.
                # If the clade has no name, the prefix “n” is used for child nodes, e.g. “n0” and “n1”.

                c.split(n=2, branch_length=child_branch_length)  # new_branch_data.relative_start_time)
                c.branch_length = original_branch_length - child_branch_length
                c.name = "internal_branch_" + str(internal_node_idx)
                c.name = None
                c.clades[0].name = branch_name
                c.clades[1].name = new_branch_data.new_branch_name
                internal_node_idx = internal_node_idx + 1

                return internal_node_idx

            recursively_split_branch_with_this_name(full_tree, branch_name, internal_node_idx, new_branch_data,
                                                         c.clades)

# if death rate = 0.001359 gene per MY, thats a death rate of 1 per 700 MY
#

def get_nodes_to_add_to_branch(parents_distance_from_root, child_clade,
                                   SSD_time_between_gene_birth_events,
                                   SSD_life_spans, duplicate_idx, randomness_idx):

    child_branch_length = child_clade.branch_length
    child_branch_name = child_clade.name
    list_of_nodes_to_add = []

    # how long since the last gene birth for this gene family
    distance_traveled_along_branch = random.randint(0, int(SSD_time_between_gene_birth_events[randomness_idx.value]))

    while distance_traveled_along_branch < child_branch_length:


        absolute_start_time = parents_distance_from_root + distance_traveled_along_branch
        absolute_end_time = absolute_start_time + SSD_life_spans[randomness_idx.value]


        new_branch_name = child_branch_name + "_duplicate_" + str(duplicate_idx)
        new_node = node_to_add(child_branch_name, new_branch_name, absolute_end_time, absolute_start_time)
        list_of_nodes_to_add.append(new_node)

        jump_to_next_event = SSD_time_between_gene_birth_events[randomness_idx.value]
        distance_traveled_along_branch = distance_traveled_along_branch + jump_to_next_event
        duplicate_idx = duplicate_idx + 1
        randomness_idx.increment()

    return list_of_nodes_to_add, duplicate_idx, randomness_idx

def write_new_newick(tree):
    handle = StringIO()
    Phylo.write(tree, handle, "newick")
    new_newick = handle.getvalue()
    return new_newick

def prune_any_branches_that_would_be_dead_before_the_end_of_the_sim(nodes_to_add_by_branch_name,
                                                                    end_of_sim):

    deleted_nodes_by_branch_name = {}
    retained_nodes_by_branch_name = {}
    for parent_branch, list_of_nodes_to_add_by_branch in nodes_to_add_by_branch_name.items():

        for duplicate_branch in list_of_nodes_to_add_by_branch:


            if end_of_sim > duplicate_branch.absolute_end_time:
                if parent_branch not in deleted_nodes_by_branch_name:
                    deleted_nodes_by_branch_name[parent_branch] = []
                deleted_nodes_by_branch_name[parent_branch].append(duplicate_branch)
            else:
                if parent_branch not in retained_nodes_by_branch_name:
                    retained_nodes_by_branch_name[parent_branch] = []
                retained_nodes_by_branch_name[parent_branch].append(duplicate_branch)

    return retained_nodes_by_branch_name, deleted_nodes_by_branch_name

def recursively_get_new_branches_to_add(tree, parent_clade, nodes_to_add_by_branch_name,
                                        SSD_time_between_gene_birth_events, SSD_life_spans,
                                        internal_node_idx, duplicate_idx, randomness_idx):

    parents_distance_from_root = tree.distance(parent_clade)
    for child_clade in parent_clade.clades:

        if not child_clade.name:
            internal_node_name = "internal_node_" + str(internal_node_idx)
            child_clade.name = internal_node_name
            internal_node_idx = internal_node_idx + 1

        nodes_to_add_to_branch, duplicate_idx, randomness_idx = get_nodes_to_add_to_branch(
            parents_distance_from_root, child_clade,
            SSD_time_between_gene_birth_events, SSD_life_spans,
            duplicate_idx, randomness_idx)

        nodes_to_add_by_branch_name[child_clade.name] = nodes_to_add_to_branch

        recursively_get_new_branches_to_add(tree, child_clade, nodes_to_add_by_branch_name,
                                                 SSD_time_between_gene_birth_events,
                                                 SSD_life_spans, internal_node_idx, duplicate_idx, randomness_idx)



def tree_to_newick(tree):
    handle = StringIO()
    Phylo.write(tree, handle, "newick")
    newick = handle.getvalue()
    return newick

def log_dict_of_node_to_add(dict_of_node_to_add):
    for branch_name, nodes_to_add in dict_of_node_to_add.items():

        log.write_to_log(branch_name + ":")
        for node in nodes_to_add:
            log.write_to_log(node.data_to_string())

def num_genes_in_dict(dict_of_node_to_add):

    total_num_genes=0
    for branch_name, nodes_to_add in dict_of_node_to_add.items():
        total_num_genes = total_num_genes+len(nodes_to_add)
    return total_num_genes

def visualize_distribution(data, title, theoretical_mean, out_folder):

    data_sum=sum(data)
    mean=data_sum/len(data)
    raw_data_plot_title=title +"_raw_data"
    hist_data_plot_title=title +"_hist_data"
    xs = [i for i in range(0, len(data))]

    fig, ax = plt.subplots(1, 1, figsize=(10, 10))
    plt.scatter(xs, data, label=raw_data_plot_title, marker="o", color='r')
    plt.axhline(y=mean, color='b', linestyle='-', label="sim avg at " + str(mean))
    plt.axhline(y=mean, color='y', linestyle='--', label="expected avg at " + str(theoretical_mean))
    plt.title(raw_data_plot_title)
    plt.legend()
    ax.set(ylabel="MY")
    ax.set(xlabel="datapoints")
    file_to_save = os.path.join(out_folder, raw_data_plot_title+".png")
    plt.savefig(file_to_save)
    plt.cla()
    plt.close()

    fig, ax = plt.subplots(1, 1, figsize=(10, 10))
    bin_size=1
    bins = np.arange(0, max(data), bin_size)
    ax.hist(data, density=True, bins=bins, alpha=0.2, label=hist_data_plot_title)
    ax.set(ylabel="density")
    ax.set(xlabel="MY (bins)")
    plt.axvline(x=mean, color='b', linestyle='-', label="sim avg at " + str(mean))
    plt.axvline(x=mean, color='y', linestyle='--', label="expected avg at " + str(theoretical_mean))
    plt.title(raw_data_plot_title)
    plt.legend()
    file_to_save = os.path.join(out_folder,hist_data_plot_title+".png")
    plt.savefig(file_to_save)
    plt.cla()
    plt.close()

class idx_by_reference():
    value=0
    def __init__(self,value):
        self.value=value

    def increment(self):
        self.value = self.value+1

    def __str__(self):
        return str(self.value)

class node_to_add():
    parent_branch_name = ""
    new_branch_name = ""
    absolute_end_time = 0
    absolute_start_time = 0

    def __init__(self, parent_branch_name, new_branch_name, absolute_end_time,
                 absolute_start_time):
        self.parent_branch_name = parent_branch_name
        self.absolute_end_time = absolute_end_time
        self.absolute_start_time = absolute_start_time
        self.new_branch_name = new_branch_name

    def print_data(self):
        print("~ data for new node " + self.new_branch_name + " ~")
        print("abs end time:\t" + str(self.absolute_end_time))

    def data_to_string(self):
        s1="~ data for new node " + self.new_branch_name + " ~"
        s3="abs end time:\t" + str(self.absolute_end_time)
        s="\n".join([s1,s3])
        return s

import os
import random
import numpy as np
from matplotlib import pyplot as plt
from scipy.stats import lognorm, expon, poisson
from Bio import Phylo
from io import StringIO
import log
from pipeline_modules import sagephy_GBD_model
from pipeline_modules.gene_tree_info import custom_gene_tree_result


def run_custom_GBD_model(polyploid, simulation_leg, base_gene_tree_newicks_by_tree_name):

    config = polyploid.general_sim_config
    include_visualizations = config.include_visualizations
    num_gene_trees_needed = config.num_gene_trees_per_species_tree
    mean_gene_birth_rate  = config.mean_gene_birth_rate #0.001359
    mean_SSD_life_span = config.mean_SSD_life_span #1 #MY
    end_of_sim = config.full_sim_time
    gt_index_formatter = sagephy_GBD_model.get_gt_index_format(num_gene_trees_needed)
    skip_GBD_model=False
    subfolder = os.path.join(polyploid.species_subfolder, str(polyploid.analysis_step_num) + "_custom_GBD")

    if not os.path.exists(subfolder):
        os.makedirs(subfolder)

    if (not mean_gene_birth_rate) or (not mean_SSD_life_span):
        skip_GBD_model=True
        log.write_to_log("skipping GBD model")
    else:
        log.write_to_log("running GBD model")
        max_possible_base_tree_length=config.full_sim_time*2
        max_number_gene_births=int(max_possible_base_tree_length*1.5*num_gene_trees_needed)+10 # the one x1.5 is just to make sure
        mean_time_between_gene_births = int(1.0 / mean_gene_birth_rate)

        #loc=mode of data (location of peak of distribution), scale=location of average of data (cm of distribution)
        SSD_life_spans = expon.rvs(loc=0, scale=mean_SSD_life_span,size=max_number_gene_births)
        SSD_time_between_gene_birth_events = poisson.rvs(mean_time_between_gene_births, size=max_number_gene_births)

        if include_visualizations:
            visualize_distribution(SSD_life_spans,"SSD_life_spans", mean_SSD_life_span, subfolder)
            visualize_distribution(SSD_time_between_gene_birth_events,
                               "SSD_time_between_gene_birth_events", mean_time_between_gene_births,subfolder)
    gene_tree_data_by_tree_name = {}
    total_num_dup_added=0
    randomness_idx=0
    for gt_idx in range(0, num_gene_trees_needed):


        gene_tree_name = "GeneTree" + gt_index_formatter.format(gt_idx)
        base_gene_tree_newick = base_gene_tree_newicks_by_tree_name[gene_tree_name]
        if skip_GBD_model:
            gene_tree_newick_with_GBD = base_gene_tree_newick
        else:
            gene_tree_newick_with_GBD, num_dup_added, randomness_idx = add_GBD_to_newick(base_gene_tree_newick,
                                                      gene_tree_name,
                                                      gt_idx,randomness_idx,
                                                      SSD_life_spans,SSD_time_between_gene_birth_events,
                                                      end_of_sim,subfolder)
            total_num_dup_added=total_num_dup_added+num_dup_added
        gene_tree_data=custom_gene_tree_result(gene_tree_name,gene_tree_newick_with_GBD)
        gene_tree_data_by_tree_name[gene_tree_name] = gene_tree_data

    new_gene_trees = gene_tree_data_by_tree_name.keys()
    #log.write_to_log("pruned_gtrees:\t" + str(new_gene_trees))
    log.write_to_log("total num duplicates added:\t" + str(total_num_dup_added))

    polyploid.analysis_step_num = polyploid.analysis_step_num + 1
    return gene_tree_data_by_tree_name



def add_GBD_to_newick(base_gene_tree_newick, gene_tree_name,
                      gt_idx, randomness_idx,
                      SSD_life_spans, SSD_time_between_gene_birth_events,
                      end_of_sim, outfolder):


        tree = Phylo.read(StringIO(base_gene_tree_newick), "newick")
        save_ascii_tree(base_gene_tree_newick, gene_tree_name, outfolder, "_before_GBD.txt", tree)

        nodes_to_add_by_branch_name = {}
        internal_node_idx = 0
        duplicate_idx = 0
        tree.clade.name="root_" + str(gt_idx)
        recursively_get_new_branches_to_add(tree, tree.clade,
                                                 nodes_to_add_by_branch_name, SSD_time_between_gene_birth_events,
                                                 SSD_life_spans, internal_node_idx, duplicate_idx, randomness_idx)

        #log.write_to_log("Genes born:")
        #log_dict_of_node_to_add(nodes_to_add_by_branch_name)

        retained_nodes_by_branch_name, deleted_nodes_by_branch_name = (
            prune_any_branches_that_would_be_dead_before_the_end_of_the_sim(nodes_to_add_by_branch_name,
                                                                                 end_of_sim))
        #log.write_to_log("Genes died:")
        #log_dict_of_node_to_add(deleted_nodes_by_branch_name)

        num_dup_to_add=num_genes_in_dict(retained_nodes_by_branch_name)
        log.write_to_log("adding " + str(num_dup_to_add)+ " new genes.. ")

        for branch_name, branches_to_add in retained_nodes_by_branch_name.items():

            #log.write_to_log("\nFor branch " + branch_name)
            for new_branch_data in branches_to_add:
                recursively_split_branch_with_this_name(tree, branch_name, internal_node_idx,
                                                             new_branch_data, tree.clade.clades)


        gene_tree_newick_with_GBD=tree_to_newick(tree)
        save_ascii_tree(gene_tree_newick_with_GBD, gene_tree_name, outfolder, "_after_GBD.txt", tree)
        return gene_tree_newick_with_GBD, num_dup_to_add, randomness_idx


def save_ascii_tree(base_gene_tree_newick, gene_tree_name, outfolder, suffix, tree):
    tree_file_path = os.path.join(outfolder, gene_tree_name + suffix)
    with open(tree_file_path, 'w') as f:
        f.writelines("newick:" + str(base_gene_tree_newick) + "\n")
        Phylo.draw_ascii(tree, f)


def recursively_split_branch_with_this_name(full_tree, branch_name, internal_node_idx, new_branch_data,
                                                clades):

        for c in clades:

            if c.name == branch_name:
                #print("\n\nfound " + branch_name)
                #print("new branch data: ")
                new_branch_data.print_data()
                original_branch_length = c.branch_length
                original_parent_branch_end = full_tree.distance(c)
                original_parent_branch_start = full_tree.distance(c) - original_branch_length

                new_branch_starts_here = new_branch_data.absolute_start_time

                #print("original parent branch length " + str(original_branch_length))
                #print("abs original parent branch start " + str(original_parent_branch_start))
                #print("relative start to new child branch " + str(new_branch_data.relative_start_time))
                #print("abs new child branch start (calculated)" + str(
                #    original_parent_branch_start + new_branch_data.relative_start_time))
                #print("abs new child branch start (saved)" + str(new_branch_data.absolute_start_time))

                child_branch_length = original_parent_branch_end - new_branch_data.absolute_start_time
                #print("child branch length (calculated from abs) " + str(child_branch_length))
                #print("abs new child branch start (saved)" + str(new_branch_data.absolute_start_time))

                # https://biopython.org/docs/1.75/api/Bio.Phylo.BaseTree.html
                # print("where_to_place_split (abs coords) " + str(new_branch_starts_here))

                # c.split:
                # New clades have the given branch_length and the same name as this clade’s
                # root plus an integer suffix (counting from 0). For example, splitting
                # a clade named “A” produces sub-clades named “A0” and “A1”.
                # If the clade has no name, the prefix “n” is used for child nodes, e.g. “n0” and “n1”.

                c.split(n=2, branch_length=child_branch_length)  # new_branch_data.relative_start_time)
                c.branch_length = original_branch_length - child_branch_length

                c.name = "internal_branch_" + str(internal_node_idx)
                c.clades[0].name = branch_name
                c.clades[1].name = new_branch_data.new_branch_name

                #print("child0 branch length after splitting " + str(c[0].branch_length))
                #print("child1 branch length after splitting " + str(c[1].branch_length))
                #print("adding " + c.clades[1].name)
                #print("preserving " + c.clades[0].name + " but its start pos has changed")
                internal_node_idx = internal_node_idx + 1

                # add distave since last duplication gene to new_branch_data
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
    jump_to_next_event = 0
    list_of_nodes_to_add = []

    # how long since the last gene birth for this gene family
    distance_traveled_along_branch = random.randint(0, int(SSD_time_between_gene_birth_events[randomness_idx]))
    #distance_traveled_along_branch = 15
    #print("random start:\t" + str(distance_traveled_along_branch))

    while distance_traveled_along_branch < child_branch_length:
        relative_start_pos = distance_traveled_along_branch
        relative_end_pos = relative_start_pos + SSD_life_spans[randomness_idx]
        absolute_end_time = parents_distance_from_root + relative_end_pos
        absolute_start_time = parents_distance_from_root + relative_start_pos
        new_branch_name = child_branch_name + "_duplicate_" + str(duplicate_idx)
        new_node = node_to_add(child_branch_name, new_branch_name,
                               relative_start_pos, absolute_end_time, absolute_start_time, jump_to_next_event)
        list_of_nodes_to_add.append(new_node)
        new_node.print_data()
        jump_to_next_event = SSD_time_between_gene_birth_events[randomness_idx]
        distance_traveled_along_branch = distance_traveled_along_branch + jump_to_next_event
        duplicate_idx = duplicate_idx + 1
        randomness_idx = randomness_idx + 1

    return list_of_nodes_to_add, duplicate_idx, randomness_idx

def write_new_newick(tree):
    handle = StringIO()
    Phylo.write(tree, handle, "newick")
    new_newick = handle.getvalue()
    return new_newick

def prune_any_branches_that_would_be_dead_before_the_end_of_the_sim(nodes_to_add_by_branch_name,
                                                                    end_of_sim):

    #print("pruning duplicates")
    #print("\t")
    deleted_nodes_by_branch_name = {}
    retained_nodes_by_branch_name = {}
    for parent_branch, list_of_nodes_to_add_by_branch in nodes_to_add_by_branch_name.items():

        for duplicate_branch in list_of_nodes_to_add_by_branch:
            #print("duplicate_branch_name:\t" + duplicate_branch.new_branch_name)
            #print("relative_start_time:\t" + str(duplicate_branch.relative_start_time))
            #print("absolute_end_time:\t" + str(duplicate_branch.absolute_end_time))
            #print("sim_end_time:\t" + str(end_of_sim))

            if end_of_sim > duplicate_branch.absolute_end_time:
                #print("killing branch.\t")
                #print("\t")
                if parent_branch not in deleted_nodes_by_branch_name:
                    deleted_nodes_by_branch_name[parent_branch] = []
                deleted_nodes_by_branch_name[parent_branch].append(duplicate_branch)
            else:
                #print("keeping branch.\t")
                #print("\t")
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

        #print("branch " + child_clade.name + ":\t length of " + str(child_clade.branch_length))

        nodes_to_add_to_branch, duplicate_idx, randomness_idx = get_nodes_to_add_to_branch(
            parents_distance_from_root, child_clade,
            SSD_time_between_gene_birth_events, SSD_life_spans,
            duplicate_idx, randomness_idx)
        #print("nodes_to_add_to_branch:" + str(nodes_to_add_to_branch))
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

class node_to_add():
    parent_branch_name = ""
    new_branch_name = ""
    relative_start_time = 0
    absolute_end_time = 0
    time_between_splits = 0
    absolute_start_time = 0

    def __init__(self, parent_branch_name, new_branch_name,
                 relative_start_time, absolute_end_time,
                 absolute_start_time, time_between_splits):
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
    def data_to_string(self):
        s1="~ data for new node " + self.new_branch_name + " ~"
        s2="rel start pos:\t" + str(self.relative_start_time)
        s3="abs end time:\t" + str(self.absolute_end_time)
        s4="jump:\t" + str(self.time_between_splits)
        s="\n".join([s1,s2,s3,s4])
        return s
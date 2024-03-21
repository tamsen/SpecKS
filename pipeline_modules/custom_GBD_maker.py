import glob
import os
import matplotlib.pyplot as plt
from scipy.stats import beta
from Bio import Phylo
from io import StringIO
import log
import process_wrapper
from pipeline_modules import gene_tree_data
from pipeline_modules import gene_tree_GBD_maker
from visualization import tree_visuals_by_phylo
from visualization import gene_tree_visuals

def run_custom_GBD_model(polyploid, simulation_leg, base_gene_tree_newicks_by_tree_name):

    config = polyploid.general_sim_config
    include_visualizations = config.include_visualizations
    num_gene_trees_needed = config.num_gene_trees_per_species_tree
    dup_rate_parameters = config.dup_rate_parameters
    loss_rate_parameters = config.loss_rate_parameters
    gt_index_formatter = gene_tree_GBD_maker.get_gt_index_format(num_gene_trees_needed)
    subfolder = os.path.join(polyploid.species_subfolder, str(polyploid.analysis_step_num) + "_custom_GBD")

    if not os.path.exists(subfolder):
        os.makedirs(subfolder)

    dup_values, loss_values = gene_tree_GBD_maker.get_randomized_dup_and_loss_rates(
        dup_rate_parameters, loss_rate_parameters, num_gene_trees_needed)

    if include_visualizations:
        gene_tree_GBD_maker.visualize_dup_and_loss_rates(dup_values, loss_values, subfolder)

    gene_tree_data_by_tree_name = {}
    for i in range(0, num_gene_trees_needed):

        gene_tree_name = "GeneTree" + gt_index_formatter.format(i)
        base_gene_tree_newick = base_gene_tree_newicks_by_tree_name[gene_tree_name]
        gene_tree_newick_with_GBD = add_GBD_to_GT(base_gene_tree_newick, i)
        gene_tree_data_by_tree_name[gene_tree_name] = gene_tree_newick_with_GBD

    pruned_tree_names = gene_tree_data_by_tree_name.keys()
    log.write_to_log("pruned_tree_files:\t" + str(pruned_tree_names))

    if include_visualizations:
        print("do nothing yet")

    polyploid.analysis_step_num = polyploid.analysis_step_num + 1
    return gene_tree_data_by_tree_name


def add_GBD_to_GT(base_gene_tree_newick, i):
    # Idea
    # for each species_tree_newick,
    # figure out how long each branch is
    # figure how many new branches to make
    # recursion would add new branches along new branches,
    # (if they were long enough & lasted long enough to duplicate before they died)
    # give each branch a start time and an end time
    # figure out which branches would still be alive by the end of the simulation
    # and add *only* those ones to the final new GT newick.
    random_seed = i  # to make the results repeatable, run to run, but different between GT
    tree = Phylo.read(StringIO(base_gene_tree_newick), "newick")
    branches=tree.get_branches()
    for b in branches:
        print(b)
    print(tree)
    #self.terminal_leaf_names = [t.name for t in self.tree.get_terminals()]
    gene_tree_newick_with_GBD = base_gene_tree_newick
    return gene_tree_newick_with_GBD

class custom_gene_tree_result():

    original_newick=""
    simple_newick=""
    terminal_leaf_names=0
    gene_tree_name=""
    gene_tree_file_name=""
    leaf_map_file_name=""
    leaves_by_species={}
    tree=False

    def __init__(self, gene_tree_name,original_newick):

        self.original_newick=original_newick
        self.simple_newick=original_newick
        self.tree= Phylo.read(StringIO(self.simple_newick), "newick")
        self.terminal_leaf_names = [t.name for t in self.tree.get_terminals()]
        self.num_terminal_leaves  = len(self.terminal_leaf_names)
        self.gene_tree_name = gene_tree_name

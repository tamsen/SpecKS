import os
from random import sample
from io import StringIO
from matplotlib import pyplot as plt
from scipy.stats import expon
import numpy as np
from Bio import Phylo

import config
import log
from pipeline_modules import custom_gene_tree_result, custom_GBD_model
from pipeline_modules.gene_tree_maker import plot_distribution

def shed_genes(polyploid, gene_data_by_gt_name, leaf_to_prune, genomes):


    if len(polyploid.subtree_subfolder) > 0:
        subfolder = os.path.join(polyploid.species_subfolder,
                                 str(polyploid.analysis_step_num) + "_WGD_gene_shedder_" + polyploid.subtree_subfolder)
    else:
        subfolder = os.path.join(polyploid.species_subfolder, str(polyploid.analysis_step_num) + "_gene_shedder")

    if not os.path.exists(subfolder):
        os.makedirs(subfolder)

    #Genes dupicated by WGD live hundreds of MY. 100 MY
    #so, to get the number remaining after t MY, we need an exponential decay fxn, with avg at 500MY
    time_since_WGD=polyploid.WGD_time_MYA
    avg_WGD_gene_lifespan=500 #my
    xscale=avg_WGD_gene_lifespan
    # scale = 1 / lambda, the decay parameter.
    # the mean of exp is 1/lambda
    fraction_WGD_genes_remaining_at_time_since_WGD = xscale*expon.pdf(time_since_WGD, loc=0, scale=xscale)
    bin_size = max(time_since_WGD/100.0,0.0001)
    xs = np.arange(0, xscale*3, bin_size)
    ys = [xscale*y for y in  [expon.pdf(x, loc=0, scale=avg_WGD_gene_lifespan) for x in xs]]

    plot_distribution(avg_WGD_gene_lifespan,time_since_WGD, "exponential", subfolder, 0,
                      config.SpecKS_config.full_sim_time, xs, ys, " decay for genes duplicated by WGD")

    all_gt_newicks=list(gene_data_by_gt_name.keys())
    num_original_gt=float(len(all_gt_newicks))
    num_genes_to_shed = int(num_original_gt*(1.0-fraction_WGD_genes_remaining_at_time_since_WGD))
    gene_trees_to_loose_a_duplicate_gene = sample(all_gt_newicks, num_genes_to_shed)

    # For each gene to be shed, a GT is randomly selected without replacement,
    # and a vertex from that GT which crosses the time-slice is removed.
    #   This process is repeated until the desired number of genes are shed.

    gt_after_everyone_that_needed_pruning_is_pruned={}
    pruned_gt=[]
    for gt in all_gt_newicks:
            if gt in gene_trees_to_loose_a_duplicate_gene:
                original_gt_newick=gene_data_by_gt_name[gt]
                print("pruning time!")
                new_gt_data =remove_a_duplicate(original_gt_newick, leaf_to_prune, genomes)
                gt_after_everyone_that_needed_pruning_is_pruned[gt] = new_gt_data
                pruned_gt.append(new_gt_data.gene_tree_name)
            else:
               new_gt_data =gene_data_by_gt_name[gt]

            gt_after_everyone_that_needed_pruning_is_pruned[gt] = new_gt_data
            custom_GBD_model.save_ascii_tree(new_gt_data.original_newick, new_gt_data.gene_tree_name, subfolder,
                            "_after_WGD_shedding.txt", new_gt_data.tree)

    log.write_to_log("total num WGD duplicates shedded:\t" + str(len(pruned_gt)))
    log.write_to_log("GT loosing WGD duplicates:\t" + str(pruned_gt))
    polyploid.analysis_step_num=polyploid.analysis_step_num+1
    return gt_after_everyone_that_needed_pruning_is_pruned



def remove_a_duplicate(gt_data, name_of_branch_to_prune, genomes_to_have_in_leaf_dict):

        #if this is the alloployploid, there are two choices: P1 or P2.
        #since left and right are totaly random, we will prune left,
        # and the result is still random
        new_gt_data= new_tree_with_a_branch_renamed(
            gt_data.original_newick, gt_data.gene_tree_name, name_of_branch_to_prune,genomes_to_have_in_leaf_dict)
        return new_gt_data

def new_tree_with_a_branch_renamed(original_newick, original_name,
                                   name_of_branch_to_prune, genomes):

    tree_copy = Phylo.read(StringIO(original_newick), "newick")
    terminal_clades = tree_copy.get_terminals()
    for c in terminal_clades:

        #this removes "P" and all duplicated with "P" in the name
        if name_of_branch_to_prune in c.name:
            tree_copy.prune(c)

    handle = StringIO()
    Phylo.write(tree_copy, handle, "newick")
    new_newick = handle.getvalue()
    new_gt_result=custom_gene_tree_result.custom_gene_tree_result(original_name,new_newick,genomes)
    return new_gt_result



def plot_distribution(theoretical_avg, time_since_WGD, distribution_name, out_folder, start,
                      xaxis_limit, xs, ys, title):
    fig, ax = plt.subplots(1, 1, figsize=(10, 10))

    plt.plot(xs, ys, label='underlying distribution')
    plt.axvline(x=theoretical_avg, color='b', linestyle='-', label="mean WGD duplicate life span ")
    plt.axvline(x=time_since_WGD, color='g', linestyle='-', label="Present time (time since WGD=" +
                                                                  str(time_since_WGD) + "MY)")
    out_file_name = os.path.join(out_folder, "Raw data " + distribution_name +title+ ".png")
                                # " Distribution in bifurcation time of gene trees for orthologs.png")
    title = distribution_name + title #' Distribution in bifurcation time of gene trees for orthologs'
    fig.suptitle(title)
    ax.set(xlabel="MY since WGD")
    #ax.set(xlim=[start, xaxis_limit])
    ax.legend()
    plt.savefig(out_file_name)
    plt.close()


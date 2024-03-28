import os
from random import sample

from matplotlib import pyplot as plt
from scipy.stats import expon
import numpy as np
from Bio import Phylo

import config
import log
from pipeline_modules.gene_tree_maker import plot_distribution
from visualization import tree_visuals_by_phylo


def shed_genes(polyploid, gene_tree_results, leaf_to_prune):


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

    all_gene_trees=list(gene_tree_results.keys())
    num_genes_to_shed = len(all_gene_trees)*(1.0-fraction_WGD_genes_remaining_at_time_since_WGD)
    gene_trees_to_loose_a_duplicate_gene = sample(all_gene_trees, num_genes_to_shed)

    # For each gene to be shed, a GT is randomly selected without replacement,
    # and a vertex from that GT which crosses the time-slice is removed.
    #   This process is repeated until the desired number of genes are shed.

    unprunable_leaves = ['O', 'O1', 'O2'] #no pruning the outgroup
    gt_after_everyone_that_needed_pruning_is_pruned={}

    for gt in all_gene_trees:
            if gt in gene_trees_to_loose_a_duplicate_gene:
                original_gt=gene_tree_results[gt]
                print("pruning time!")
                new_gt_data=remove_a_duplicate(original_gt,leaf_to_prune)
            else:
               gt_after_everyone_that_needed_pruning_is_pruned[gt]=gene_tree_results[gt]

    polyploid.analysis_step_num=polyploid.analysis_step_num+1
    return gt_after_everyone_that_needed_pruning_is_pruned



def remove_a_duplicate(gt_data,leaf_to_prune):

        #if this is the alloployploid, there are two choices: P1 or P2.
        #since left and right are totaly random, we will prune left,
        # and the result is still
        return gt_data




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


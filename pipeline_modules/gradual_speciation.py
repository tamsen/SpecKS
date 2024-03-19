import os
from scipy.stats import norm,lognorm
import numpy as np
from matplotlib import pyplot as plt
from Bio import Phylo
from io import StringIO

from pipeline_modules import species_tree_maker
from pipeline_modules.gene_tree_maker import get_gt_index_format

def get_mode_and_cm(xs, pdfs):

    ymax = max(pdfs)
    xs_of_ymax = []
    weighted_mass = 0
    weights = 0

    for i in range(0, len(xs)):
        x = xs[i]
        y = pdfs[i]
        weighted_mass = weighted_mass + (x * y)
        weights = weights + y
        if y == ymax:
            xs_of_ymax.append(x)

    x_value_of_ymax = sum(xs_of_ymax) / len(xs_of_ymax)
    center_of_mass = weighted_mass / weights
    return center_of_mass, x_value_of_ymax

def get_per_gene_tree_variation_on_speciation_time(out_folder,
                                                   num_gt_needed, include_vis):
    # for 90% < 0.55MY
    # shape_parameter = 0.8
    # xscale = 0.2
    # start = -0.55
    
    #for 10 MY spread
    # for 90% < 10MY
    time_span_MY = 10 * 1.1
    shape_parameter=0.5
    xscale = 5.27
    #start = -4.1
    #loc=start
    start=0
    loc=0

    #https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.lognorm.html
    random_draws_from_distribution = lognorm.rvs(shape_parameter, size=num_gt_needed, scale=xscale, loc=loc)

    #if include_vis:
    time_span_MY = time_span_MY
    bin_size = 0.1
    xaxis_limit = time_span_MY + 0.1

    xs = np.arange(start, time_span_MY+ 0.01, bin_size)
    ys = [lognorm.pdf(x,shape_parameter, scale=xscale,loc=loc ) for x in xs]
    center_of_mass, x_value_of_ymax = get_mode_and_cm(xs, ys)

    fig, ax = plt.subplots(1, 1, figsize=(10, 10))
    bins = np.arange(start, xaxis_limit, bin_size)
    ax.hist(random_draws_from_distribution, density=True, bins=bins, alpha=0.2, label='distribution of simulated data')
    plt.plot(xs, ys, label='underlying distribution')
    ax.scatter(x_value_of_ymax, 0.01,
                            color='darkgreen', marker='^', label="mode="+str(round(x_value_of_ymax,2)), s=100)

    out_file_name = os.path.join(out_folder, "Distribution in bifurcation time of gene trees for orthologs.png")
    title='Distribution in bifurcation time of gene trees for orthologs'
    fig.suptitle(title)
    ax.set(xlabel="MYA")
    ax.set(xlim=[start, xaxis_limit])
    ax.legend()
    plt.savefig(out_file_name)
    plt.close()

    bifurcaton_variations=[ri-x_value_of_ymax for  ri in random_draws_from_distribution ]

    if include_vis:
        fig, ax = plt.subplots(1, 1, figsize=(10, 10))
        bins = np.arange(start-x_value_of_ymax, xaxis_limit-x_value_of_ymax, bin_size)
        ax.hist(bifurcaton_variations, density=True, bins=bins, alpha=0.2, label='re-centered distribution of simulated data')
        out_file_name = os.path.join(out_folder, "Re-centered Distribution in bifurcation time of gene trees for orthologs.png")
        title='Distribution in bifurcation time of gene trees for orthologs'
        fig.suptitle(title)
        ax.set(xlabel="MYA")
        ax.set(xlim=[start-x_value_of_ymax, xaxis_limit-x_value_of_ymax])
        ax.legend()
        plt.savefig(out_file_name)
        plt.close()

    return bifurcaton_variations

def make_randomized_gene_trees(polyploid, simulation_leg, species_tree_newick):

    config = polyploid.general_sim_config
    include_visualizations = config.include_visualizations
    num_gene_trees_needed = config.num_gene_trees_per_species_tree
    gt_index_formatter = get_gt_index_format(num_gene_trees_needed)

    subfolder=os.path.join(polyploid.species_subfolder, str(polyploid.analysis_step_num) + "_randomized_gene_trees")

    if not os.path.exists(subfolder):
        os.makedirs(subfolder)


    #if include_visualizations:

    base_speciation_time=polyploid.SPC_time_MYA
    time_span=polyploid.FULL_time_MYA
    r=get_per_gene_tree_variation_on_speciation_time(subfolder,num_gene_trees_needed,include_visualizations)

    gene_tree_newicks_by_tree_name={}
    for i in range(0, num_gene_trees_needed):

        gt_name = "GeneTree" +  gt_index_formatter.format(i)
        #perturb parent species tree
        variation=round(r[i],3)
        gene_bifurcation_time=base_speciation_time+variation
        [new_newick_string]= species_tree_maker.get_example_allopolyploid_tree(time_span,gene_bifurcation_time)
        gene_tree_newicks_by_tree_name[gt_name ]=new_newick_string

        delta_path = os.path.join(subfolder, "variations_in_bifurtaion_time.txt")
        with open(delta_path, 'a') as f:
            f.writelines(str(variation) + "\n")

        tree_path = os.path.join(subfolder, "gene_trees_with_lognorm_distributed_bifurcation.txt")
        with open(tree_path, 'a') as f:
            f.writelines(new_newick_string + "\n")


    polyploid.analysis_step_num = polyploid.analysis_step_num + 1
    return gene_tree_newicks_by_tree_name


import os
from scipy.stats import lognorm, expon
import numpy as np
from matplotlib import pyplot as plt
from Bio import Phylo
from io import StringIO

import log
from pipeline_modules import species_tree_maker
from pipeline_modules.sagephy_GBD_model import get_gt_index_format

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


def get_per_gene_tree_variation_on_speciation_time(out_folder,num_gt_needed,
                                                   distribution_parameters, include_vis):
    # for 90% < 0.55MY
    # shape_parameter = 0.8
    # xscale = 0.2
    # start = -0.55
    
    #for 10 MY spread
    # for 90% < 10MY
    time_span_MY = 10 * 1.1
    #shape_parameter=0.5
    #xscale = 5.27
    #start = -4.1
    #loc=start
    start=0

    time_span_MY = time_span_MY

    #xaxis_limit = time_span_MY + 0.1

    distribution_name=distribution_parameters[0].upper() #exponential, lognorm..
    shape_parameter = float(distribution_parameters[1])
    loc=float(distribution_parameters[1])
    xscale = float(distribution_parameters[2])
    bin_size = min(0.1,xscale*.1)
    xs = np.arange(start, bin_size*50+ 0.01, bin_size)
    xaxis_limit = max(xs)
    #for expon
    #loc = 0, scale = mean_SSD_life_span
    if distribution_name=="LOGNORM":
        distribution_fxn=lognorm
    elif (distribution_name=="EXPONENTIAL") or (distribution_name=="EXPON") :
        distribution_fxn=expon
    elif (distribution_name == "IMPULSE"):
        distribution_fxn = "IMPULSE"
    else:
        message="The distribution requested is not supported"
        log.write_to_log(message)
        raise ValueError(message)

    if distribution_fxn == "IMPULSE":
        random_draws_from_distribution = [0 for i in range(0, num_gt_needed)]
        ys = [1 if x == 0 else 0 for x in xs]
    else:
        #https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.lognorm.html
        random_draws_from_distribution = distribution_fxn.rvs(loc=loc, size=num_gt_needed, scale=xscale)
        ys = [distribution_fxn.pdf(x, loc=loc, scale=xscale) for x in xs]

    center_of_mass, x_value_of_ymax = get_mode_and_cm(xs, ys)

    fig, ax = plt.subplots(1, 1, figsize=(10, 10))
    bins = np.arange(start, xaxis_limit, bin_size)
    ax.hist(random_draws_from_distribution, density=True, bins=bins, alpha=0.2, label='distribution of simulated data')
    plt.plot(xs, ys, label='underlying distribution')
    ax.scatter(x_value_of_ymax, 0.01,
                            color='darkgreen', marker='^', label="mode="+str(round(x_value_of_ymax,2)), s=100)
    ax.scatter(center_of_mass, 0.01,
                            color='lightgreen', marker='^', label="mean="+str(round(center_of_mass,2)), s=100)

    out_file_name = os.path.join(out_folder, "Raw data " + distribution_name +
                                 " Distribution in bifurcation time of gene trees for orthologs.png")
    title= distribution_name + ' Distribution in bifurcation time of gene trees for orthologs'
    fig.suptitle(title)
    ax.set(xlabel="MYA")
    ax.set(xlim=[start,xaxis_limit])
    ax.legend()
    plt.savefig(out_file_name)
    plt.close()

    bifurcaton_variations=[ri-x_value_of_ymax for  ri in random_draws_from_distribution ]

    if include_vis:
        fig, ax = plt.subplots(1, 1, figsize=(10, 10))
        bins = np.arange(start-x_value_of_ymax, xaxis_limit-x_value_of_ymax, bin_size)
        ax.hist(bifurcaton_variations, density=True, bins=bins, alpha=0.2, label='distribution centered at mode')
        out_file_name = os.path.join(out_folder,
                                     "Shifted distribution " + distribution_name + " of bifurcation time of gene trees for orthologs.png")
        title= (distribution_name + ' Distribution in bifurcation time of gene trees for orthologs, subtracted '
                +str(x_value_of_ymax))
        fig.suptitle(title)
        ax.set(xlabel="MYA")
        ax.set(xlim=[start-x_value_of_ymax, xaxis_limit-x_value_of_ymax])
        #ax.set(xlim=[start - x_value_of_ymax, center_of_mass*2 - x_value_of_ymax])
        ax.legend()
        plt.savefig(out_file_name)
        plt.close()

    return bifurcaton_variations



def make_randomized_gene_trees(polyploid, species_tree_newick):

    config = polyploid.general_sim_config
    include_visualizations = config.include_visualizations
    num_gene_trees_needed = config.num_gene_trees_per_species_tree
    gene_divergence_distribution_parameters = config.divergence_distribution_parameters
    base_speciation_time=polyploid.SPC_time_MYA
    time_span=polyploid.FULL_time_MYA
    gt_index_formatter = get_gt_index_format(num_gene_trees_needed)
    subfolder=os.path.join(polyploid.species_subfolder, str(polyploid.analysis_step_num) + "_randomized_gene_trees")

    if not os.path.exists(subfolder):
        os.makedirs(subfolder)

    if polyploid.is_auto():
        log.write_to_log("All gene trees diverge at the same instant for an autopolyploid. Nothing to do here.")
        gene_divergence_distribution_parameters = ["impulse",1,1]
    elif not gene_divergence_distribution_parameters:
        log.write_to_log("Config specifies no variation in divergence time for orthologs.")
        gene_divergence_distribution_parameters = ["impulse",1,1]
    else:
        log.write_to_log("Calculating randomized divergence times for gene trees within allopolyploid species tree")

    variations_in_gene_div_time = get_per_gene_tree_variation_on_speciation_time(
            subfolder,num_gene_trees_needed,gene_divergence_distribution_parameters, include_visualizations)

    gene_tree_newicks_by_tree_name={}
    for i in range(0, num_gene_trees_needed):

        gt_name = "GeneTree" +  gt_index_formatter.format(i)
        variation=round(variations_in_gene_div_time[i],3)

        if gene_divergence_distribution_parameters[0] == "impulse":
            #In this case, the gene trees will be identical to the species tree.
            new_newick_string = species_tree_newick

        else:    #else, for an allopolyploid...we do something more fancy.
            #perturb parent species tree
            gene_div_time=base_speciation_time+variation
            [new_newick_string] = species_tree_maker.get_example_allopolyploid_tree(time_span,gene_div_time)

        gene_tree_newicks_by_tree_name[gt_name]=new_newick_string
        delta_path = os.path.join(subfolder, "variations_in_div_time.txt")
        with open(delta_path, 'a') as f:
            f.writelines(str(variation) + "\n")

        tree_path = os.path.join(subfolder, "gene_trees_with_lognorm_distributed_div_time.txt")
        with open(tree_path, 'a') as f:
            f.writelines(new_newick_string + "\n")


    polyploid.analysis_step_num = polyploid.analysis_step_num + 1
    return gene_tree_newicks_by_tree_name


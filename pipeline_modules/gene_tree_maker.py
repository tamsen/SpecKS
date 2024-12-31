# """gene_tree_maker.py
# Builds initial gene trees based on a given species tree.
# Gene trees are unique, each individually perturbed from the species
# tree, with the variation drawn from a distribution given by the
# gene_div_time_distribution_parameters stipulated in the config file
# """

import os
from scipy.stats import lognorm, expon
import numpy as np
from matplotlib import pyplot as plt
import log
from pipeline_modules import species_tree_maker

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
                                                   distribution_parameters, include_vis,
                                                   center_distribution_on_cm):

    #TODO consider moving distribution_parameters parsing / sanitization to the config class.
    start=0
    distribution_name=distribution_parameters[0].upper() #exponential, lognorm..
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
        xs = np.arange(start, xaxis_limit + bin_size, bin_size)
        ys = [1 if x == 0 else 0 for x in xs]
    elif distribution_fxn == expon:
        random_draws_from_distribution = distribution_fxn.rvs(loc=loc, size=num_gt_needed, scale=xscale)
        xs = np.arange(start, xaxis_limit + bin_size, bin_size)
        ys = [distribution_fxn.pdf(x, loc=loc, scale=xscale) for x in xs]

    else: #lognorm
        #https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.lognorm.html
        loc = 0
        xaxis_limit = 10
        bin_size = 0.1
        shape_parameter = float(distribution_parameters[1])
        xscale = float(distribution_parameters[2])
        random_draws_from_distribution = lognorm.rvs(shape_parameter, size=num_gt_needed, scale=xscale, loc=loc)

        xs = np.arange(start, xaxis_limit + bin_size, bin_size)
        ys = [distribution_fxn.pdf(x, shape_parameter, scale=xscale, loc=0) for x in xs]

    center_of_mass, x_value_of_ymax = get_mode_and_cm(xs, ys)

    plot_distribution(bin_size, center_of_mass, distribution_name, out_folder, random_draws_from_distribution, start,
                      x_value_of_ymax, xaxis_limit, xs, ys, "Distribution in bifurcation time of gene trees for orthologs.png")

    #Todo, make using x_value_of_ymax configurable.
    if center_distribution_on_cm:
        bifurcaton_variations=[ri-x_value_of_ymax for  ri in random_draws_from_distribution ]
    else:
        bifurcaton_variations=[ri for  ri in random_draws_from_distribution ]

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


def plot_distribution(bin_size, center_of_mass, distribution_name, out_folder, random_draws_from_distribution, start,
                      x_value_of_ymax, xaxis_limit, xs, ys, title):
    fig, ax = plt.subplots(1, 1, figsize=(10, 10))
    bins = np.arange(start, xaxis_limit, bin_size)
    ax.hist(random_draws_from_distribution, density=True, bins=bins, alpha=0.2, label='distribution of simulated data')
    plt.plot(xs, ys, label='underlying distribution')
    ax.scatter(x_value_of_ymax, 0.01,
               color='darkgreen', marker='^', label="mode=" + str(round(x_value_of_ymax, 2)), s=100)
    ax.scatter(center_of_mass, 0.01,
               color='lightgreen', marker='^', label="mean=" + str(round(center_of_mass, 2)), s=100)
    out_file_name = os.path.join(out_folder, "Raw data " + distribution_name +title+ ".png")
                                # " Distribution in bifurcation time of gene trees for orthologs.png")
    title = distribution_name + title #' Distribution in bifurcation time of gene trees for orthologs'
    fig.suptitle(title)
    ax.set(xlabel="MYA")
    ax.set(xlim=[start, xaxis_limit])
    ax.legend()
    plt.savefig(out_file_name)
    plt.close()


def make_randomized_gene_trees(polyploid, species_tree_newick):

    config = polyploid.general_sim_config
    include_visualizations = config.include_visualizations
    num_gene_trees_needed = config.num_gene_trees_per_species_tree
    center_distribution_on_cm = config.recenter_gene_divergence_distribution_on_cm
    gene_divergence_distribution_parameters_list = polyploid.list_of_divergence_distribution_parameters
    base_speciation_time=polyploid.DIV_time_MYA
    time_span=polyploid.FULL_time_MYA
    gt_index_formatter = get_gt_index_format(num_gene_trees_needed)
    subfolder=os.path.join(polyploid.species_subfolder, str(polyploid.analysis_step_num) + "_randomized_gene_trees")

    if not os.path.exists(subfolder):
        os.makedirs(subfolder)

    if (not gene_divergence_distribution_parameters_list) or (len(gene_divergence_distribution_parameters_list)==0):
        log.write_to_log("Config specifies no variation in divergence time for orthologs.")
        log.write_to_log("Defaulting to force all gene trees to diverge at the same instant.")
        gene_divergence_distribution_parameters_list = [["impulse",1,1,1]]
    else:
        log.write_to_log("Calculating randomized divergence times for gene trees within species tree")

    num_gt_needed_per_distribution=[]
    for gene_divergence_distribution_parameters in gene_divergence_distribution_parameters_list:

        if len(gene_divergence_distribution_parameters) < 4:
            percent_genes_to_follow_this_distribution = 1.0
        else:
            percent_genes_to_follow_this_distribution=float(gene_divergence_distribution_parameters[3])

        num_gene_trees_needed_for_this_distribution=(round(percent_genes_to_follow_this_distribution*
                                                          num_gene_trees_needed))
        num_gt_needed_per_distribution.append(num_gene_trees_needed_for_this_distribution)

    #don't accidentally loose a single gt due to integer roudning issues
    deficit= num_gene_trees_needed - sum(num_gt_needed_per_distribution)
    num_gt_needed_per_distribution[0]=num_gt_needed_per_distribution[0]+deficit

    #accept that the GT variations might come from different distributions
    # ie, segmental allo/autopolyploid

    variations_including_all_distributions=[]
    for i in range(0,len(gene_divergence_distribution_parameters_list)):
        gene_divergence_distribution_parameters= gene_divergence_distribution_parameters_list[i]
        variations_in_gene_div_time_for_dist = get_per_gene_tree_variation_on_speciation_time(
            subfolder,
            num_gt_needed_per_distribution[i],gene_divergence_distribution_parameters, include_visualizations,
            center_distribution_on_cm)
        variations_including_all_distributions=variations_including_all_distributions+variations_in_gene_div_time_for_dist

    gene_tree_newicks_by_tree_name={}
    for i in range(0, num_gene_trees_needed):

        gt_name = "GeneTree" +  gt_index_formatter.format(i)
        variation=round(variations_including_all_distributions[i],3)

        if gene_divergence_distribution_parameters[0] == "impulse":
            #In this case, the gene trees will be identical to the species tree.
            new_newick_string = species_tree_newick

        else:    #else, for an allopolyploid...we do something more fancy.
            #perturb parent species tree
            gene_div_time=base_speciation_time+variation

            if gene_div_time > time_span:
                message=["The coalescent for this gene tree reaches further back than the sim time",
                         "Please increase your sim time",
                         "gene tree divergence time, MYA: " + str(gene_div_time ),
                         "earliest sim time, MYA: " + str(time_span),]
                log.write_to_log(message[0])
                log.write_to_log(message[1])
                log.write_to_log(message[2])
                log.write_to_log(message[3])
                raise ValueError(" ".join(message))

            [new_newick_string] = species_tree_maker.get_polyploid_species_tree(time_span, gene_div_time)


        gene_tree_newicks_by_tree_name[gt_name]=new_newick_string
        delta_path = os.path.join(subfolder, "variations_in_div_time.txt")
        with open(delta_path, 'a') as f:
            f.writelines(str(variation) + "\n")

        tree_path = os.path.join(subfolder, "gene_trees_with_lognorm_distributed_div_time.txt")
        with open(tree_path, 'a') as f:
            f.writelines(new_newick_string + "\n")


    polyploid.analysis_step_num = polyploid.analysis_step_num + 1
    return gene_tree_newicks_by_tree_name

def get_gt_index_format(num_gene_trees_needed):
    decimals_needed = len(str(num_gene_trees_needed))
    formatter = "{:0" + str(decimals_needed) + "d}"
    return formatter

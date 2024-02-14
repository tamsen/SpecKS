import glob
import os

from Bio import Phylo
from io import StringIO
import common
import matplotlib.pyplot as plt

from pipeline_modules import gene_tree_data
from visualization import tree_visuals_by_phylo
from visualization import gene_tree_visuals
from scipy.stats import beta


def write_SaGePhy_GuestTreeGen_commands(config, species_tree_newick, dup_rate, loss_rate,
                                        random_seed, out_file_name):


    cmd = ["java", "-jar", config.path_to_sagephy,
         "GuestTreeGen", species_tree_newick,
         str(dup_rate), str(loss_rate), "0.0",out_file_name, "-s", str(random_seed)]

    #"nox" is for "no auxillary tags" as per the manual
    # If nox is ON you get the simple netwick that evolver wants
    # If nox is OFF you don't get the *.pruned.leafmap
    # files that help figure out how many sequsences you need
    # to propagate through the tree
    return cmd

def get_randomized_dup_and_loss_rates(dup_rate_parameters,loss_rate_parameters,num_values_needed):
    dup_values = beta.rvs(dup_rate_parameters[0],dup_rate_parameters[1], size=num_values_needed)
    loss_values = beta.rvs(loss_rate_parameters[0],loss_rate_parameters[1], size=num_values_needed)
    return dup_values,loss_values

def visualize_dup_and_loss_rates(dup_values,loss_values,out_folder):

    xs=[i for i in range(0,len(dup_values))]
    plt.scatter(xs,dup_values, label="duplication rate",marker="o",color='r')
    plt.scatter(xs,loss_values, label="loss rate",marker="o",color='b')

    if len(dup_values) > 0 and  len(loss_values) > 0:
        mean_dup_values=sum(dup_values)/len(dup_values)
        mean_loss_values=sum(loss_values)/len(loss_values)
        plt.axhline(y=mean_dup_values, color='r', linestyle='--', label="avg dup rate")
        plt.axhline(y=mean_loss_values, color='b', linestyle=':', label="avg loss rate")
        print("Mean Gene Duplication Rate:\t" + str(mean_dup_values))
        print("Mean Gene Loss Rate:\t" + str(mean_loss_values))

    plt.title("Gene duplication and loss rates for each gene tree")
    plt.legend()
    file_to_save=os.path.join(out_folder,"Gene duplication and loss rates")
    plt.savefig(file_to_save)
    plt.cla()
    plt.close()

def read_pruned_trees(subfolder, leg_distance):

    tree_files = glob.glob(subfolder + "/*.pruned.tree")
    results_by_tree_name = {}
    for tree_file in tree_files:

        leafmap_file=tree_file.replace(".tree",".leafmap")
        result = gene_tree_data.read_gene_tree_result_from_tree_and_leaf_map_files(tree_file, leafmap_file)
        result.add_back_outgroup(leg_distance)
        results_by_tree_name[result.gene_tree_name]=result

        new_tree_file=tree_file.replace(".tree",".updated.tree")
        with open(new_tree_file, 'w') as f:
            f.writelines(result.simple_newick + "\n")

    return results_by_tree_name

def run_sagephy(polyploid, simulation_leg, species_tree_newick):

    config = polyploid.general_sim_config
    num_gene_trees_needed = config.num_gene_trees_per_species_tree
    dup_rate_parameters = config.dup_rate_parameters
    loss_rate_parameters = config.loss_rate_parameters

    if len(polyploid.subtree_subfolder)>0:
        subfolder = os.path.join(polyploid.species_subfolder,
                                 str(polyploid.analysis_step_num) + "_gene_trees_" + polyploid.subtree_subfolder)
    else:
        subfolder=os.path.join(polyploid.species_subfolder, str(polyploid.analysis_step_num) + "_gene_trees")

    print(subfolder)
    if not os.path.exists(subfolder):
        os.makedirs(subfolder)

    dup_values, loss_values = get_randomized_dup_and_loss_rates(
        dup_rate_parameters, loss_rate_parameters, num_gene_trees_needed)

    visualize_dup_and_loss_rates(dup_values, loss_values, subfolder)

    for i in range(0, num_gene_trees_needed):
        out_file_name = "GeneTree" + str(i)
        random_seed = i #to make the results repeatable, run to run, but different between GT
        cmd = write_SaGePhy_GuestTreeGen_commands(config, species_tree_newick,
                                                  dup_values[i], loss_values[i],
                                                  random_seed,
                                                  os.path.join(subfolder,out_file_name))
        common.run_and_wait_on_process(cmd, subfolder)


    gene_tree_data_by_tree_name = read_pruned_trees(subfolder, simulation_leg.leg_length())
    pruned_tree_names=gene_tree_data_by_tree_name.keys()
    print("pruned_tree_files:\t" + str(pruned_tree_names))

    gt_tree_viz_data_by_name={}
    j=0
    max_num_gt_to_plot=10
    for gt_name in pruned_tree_names:
        if j> max_num_gt_to_plot:
            break
        plot_file_name_1= os.path.join(subfolder,gt_name +"_phylo.png")
        plot_file_name_2= os.path.join(subfolder,gt_name +"_specks.png")
        gt_newick=gene_tree_data_by_tree_name[gt_name].simple_newick
        leaf_map = gene_tree_data_by_tree_name[gt_name].leaves_by_species
        print("newick to plot:\t" +gt_newick)
        tree_visuals_by_phylo.save_tree_plot(gt_newick, plot_file_name_1)
        gt_tree_viz_data=gene_tree_visuals.plot_polyploid_gene_tree_alone(
            simulation_leg, leaf_map,gt_newick, gt_name,
            polyploid.SPC_time_MYA,polyploid.species_name,
            plot_file_name_2)
        gt_tree_viz_data_by_name[gt_name]=gt_tree_viz_data
        i=i+1

    gene_tree_visuals.histogram_node_distances(polyploid, gt_tree_viz_data_by_name,
                                               "raw_gene_tree",  subfolder)
    gene_tree_visuals.plot_gene_trees_on_top_of_species_trees(polyploid,
                                                              gt_tree_viz_data_by_name,
                                                              "raw_gene_tree", subfolder)
    polyploid.analysis_step_num=polyploid.analysis_step_num+1
    return gene_tree_data_by_tree_name

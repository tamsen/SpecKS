import glob
import os
import matplotlib.pyplot as plt
from scipy.stats import beta

import process_wrapper
from pipeline_modules import gene_tree_data
from visualization import tree_visuals_by_phylo
from visualization import gene_tree_visuals



def write_SaGePhy_GuestTreeGen_commands(config, species_tree_newick, dup_rate, loss_rate,
                                        random_seed, out_file_name):


    cmd = ["java", "-jar", config.path_to_sagephy,
         "GuestTreeGen", species_tree_newick,
         str(dup_rate), str(loss_rate), "0.0",out_file_name, "-s", str(random_seed)]

    #"nox" is for "no auxillary tags" as per the manual
    # If nox is ON you get the simple netwick that evolver wants
    # If nox is OFF you don't get the *.pruned.leafmap
    # files that help figure out how many sequences you need
    # to propagate through the tree
    return cmd

def get_randomized_dup_and_loss_rates(dup_rate_parameters,loss_rate_parameters,num_values_needed):

    if dup_rate_parameters:
        dup_values = beta.rvs(dup_rate_parameters[0],dup_rate_parameters[1], size=num_values_needed)
    else:
        dup_values = [0 for x in range(0,num_values_needed)]

    if loss_rate_parameters:
        loss_values = beta.rvs(loss_rate_parameters[0],loss_rate_parameters[1], size=num_values_needed)
    else:
        loss_values = [0 for x in range(0,num_values_needed)]

    return dup_values,loss_values

def get_gt_index_format(num_gene_trees_needed):
        decimals_needed = len(str(num_gene_trees_needed))
        formatter = "{:0" + str(decimals_needed) + "d}"
        return formatter

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
        #result.add_back_outgroup(leg_distance)
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
    gt_index_formatter = get_gt_index_format(num_gene_trees_needed)

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
        out_file_name = "GeneTree" +  gt_index_formatter.format(i)
        random_seed = i #to make the results repeatable, run to run, but different between GT
        cmd = write_SaGePhy_GuestTreeGen_commands(config, species_tree_newick,
                                                  dup_values[i], loss_values[i],
                                                  random_seed,
                                                  os.path.join(subfolder,out_file_name))
        process_wrapper.run_and_wait_on_process(cmd, subfolder)


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
        plot_file_name_3= os.path.join(subfolder,gt_name +"_tree_phylo.png")
        gt_tree=gene_tree_data_by_tree_name[gt_name].tree
        gt_newick=gene_tree_data_by_tree_name[gt_name].simple_newick
        leaf_map = gene_tree_data_by_tree_name[gt_name].leaves_by_species
        print("newick to plot:\t" +gt_newick)
        tree_visuals_by_phylo.save_tree_plot_from_newick(gt_newick, plot_file_name_1)
        tree_visuals_by_phylo.save_tree_plot_directly_from_tree(gt_tree, plot_file_name_3)
        gt_tree_viz_data=gene_tree_visuals.plot_polyploid_gene_tree_alone(
            simulation_leg, leaf_map,gt_newick, gt_name,
            polyploid.SPC_time_MYA,polyploid.species_name,
            plot_file_name_2)
        gt_tree_viz_data_by_name[gt_name]=gt_tree_viz_data
        j=j+1

    gene_tree_visuals.histogram_node_distances(polyploid, gt_tree_viz_data_by_name,
                                               "raw_gene_tree",  subfolder)
    gene_tree_visuals.plot_gene_trees_on_top_of_species_trees(polyploid,
                                                              gt_tree_viz_data_by_name,
                                                              "raw_gene_tree", subfolder)
    polyploid.analysis_step_num=polyploid.analysis_step_num+1
    return gene_tree_data_by_tree_name

def run_sagephy_with_root_seq(polyploid, simulation_leg, species_tree_newick,
                              root_seq_files_written_by_gene_tree_by_child_tree):

    config = polyploid.general_sim_config
    dup_rate_parameters = config.dup_rate_parameters
    loss_rate_parameters = config.loss_rate_parameters

    original_gene_trees=root_seq_files_written_by_gene_tree_by_child_tree.keys()
    num_children_per_GT=[len(root_seq_files_written_by_gene_tree_by_child_tree[gt]) for gt in original_gene_trees]
    num_rates_needed=sum(num_children_per_GT)


    if len(polyploid.subtree_subfolder)>0:
        subfolder = os.path.join(polyploid.species_subfolder,
                                 str(polyploid.analysis_step_num) + "_gene_trees_" + polyploid.subtree_subfolder)
    else:
        subfolder=os.path.join(polyploid.species_subfolder, str(polyploid.analysis_step_num) + "_gene_trees")

    print(subfolder)
    if not os.path.exists(subfolder):
        os.makedirs(subfolder)

    dup_values, loss_values = get_randomized_dup_and_loss_rates(
        dup_rate_parameters, loss_rate_parameters, num_rates_needed)

    visualize_dup_and_loss_rates(dup_values, loss_values, subfolder)

    idx_for_new_tree_generated=0
    for original_GT_name, sequence_data_file_by_child_for_leaves in root_seq_files_written_by_gene_tree_by_child_tree.items():

        for child_GT_name in sequence_data_file_by_child_for_leaves:

            #GeneTree05_childG6_1.rep1.txt
            out_file_name= original_GT_name + "_child" + child_GT_name
            random_seed = idx_for_new_tree_generated #to make the results repeatable, run to run, but different between GT
            cmd = write_SaGePhy_GuestTreeGen_commands(config, species_tree_newick,
                                                  dup_values[idx_for_new_tree_generated], loss_values[idx_for_new_tree_generated],
                                                  random_seed,
                                                  os.path.join(subfolder,out_file_name))
            process_wrapper.run_and_wait_on_process(cmd, subfolder)
            idx_for_new_tree_generated = idx_for_new_tree_generated + 1


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
        plot_file_name_3= os.path.join(subfolder,gt_name +"_tree_phylo.png")
        gt_tree=gene_tree_data_by_tree_name[gt_name].tree
        gt_newick=gene_tree_data_by_tree_name[gt_name].simple_newick
        leaf_map = gene_tree_data_by_tree_name[gt_name].leaves_by_species
        print("newick to plot:\t" +gt_newick)
        tree_visuals_by_phylo.save_tree_plot_from_newick(gt_newick, plot_file_name_1)
        tree_visuals_by_phylo.save_tree_plot_directly_from_tree(gt_tree, plot_file_name_3)
        gt_tree_viz_data=gene_tree_visuals.plot_polyploid_gene_tree_alone(
            simulation_leg, leaf_map,gt_newick, gt_name,
            polyploid.SPC_time_MYA,polyploid.species_name,
            plot_file_name_2)
        gt_tree_viz_data_by_name[gt_name]=gt_tree_viz_data
        j=j+1

    gene_tree_visuals.histogram_node_distances(polyploid, gt_tree_viz_data_by_name,
                                               "raw_gene_tree",  subfolder)
    gene_tree_visuals.plot_gene_trees_on_top_of_species_trees(polyploid,
                                                              gt_tree_viz_data_by_name,
                                                              "raw_gene_tree", subfolder)
    polyploid.analysis_step_num=polyploid.analysis_step_num+1
    return gene_tree_data_by_tree_name

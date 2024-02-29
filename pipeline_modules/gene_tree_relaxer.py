import os
import pipeline_modules.gene_tree_maker as gene_tree_maker
import process_wrapper
from pipeline_modules import gene_tree_data
from visualization import tree_visuals_by_phylo, gene_tree_visuals

def relax(polyploid, simulation_leg, gene_tree_results_by_tree_name):

    config = polyploid.general_sim_config

    if len(polyploid.subtree_subfolder) > 0:
        subfolder = os.path.join(polyploid.species_subfolder,
                                 str(polyploid.analysis_step_num) + "_relaxed_gene_trees_" + polyploid.subtree_subfolder)
    else:
        subfolder = os.path.join(polyploid.species_subfolder, str(polyploid.analysis_step_num) + "_relaxed_gene_trees")

    print(subfolder)
    if not os.path.exists(subfolder):
        os.makedirs(subfolder)

    relaxed_gene_tree_results_by_gene_tree={}
    gt_tree_viz_data_by_gene_tree={}
    random_seed=0
    for gene_tree, gene_tree_results in gene_tree_results_by_tree_name.items():
        random_seed=random_seed+1
        relaxed_tree_file_out=gene_tree +".relaxed.tree"
        relaxation_model=config.branch_relaxation_parameters

        if not relaxation_model:
            relaxed_gene_tree_results_by_gene_tree[gene_tree] = gene_tree_results
            continue

        in_file_name =gene_tree_results.gene_tree_file_name
        cmd= ["java","-jar",
             config.path_to_sagephy,"BranchRelaxer",
              "-x","-innms","-o" ,relaxed_tree_file_out, in_file_name,
              *relaxation_model, "-s", str(random_seed)]

        #-x for "include auxillary tags
        #-innms for "Keep interior vertex names in the output tree."
        # relaxation_model of "ACRY07 1 0.000001" is Autocorrelated lognormal rate in accordance with Rannala-
        # Yang ’07, with parameters as specified.

        process_wrapper.run_and_wait_on_process(cmd, subfolder)
        full_path_to_relaxed_tree_file=os.path.join(subfolder,relaxed_tree_file_out)
        relaxed_gene_tree_results = gene_tree_data.read_gene_tree_result_from_tree_and_leaf_map_files(
            full_path_to_relaxed_tree_file, gene_tree_results.leaf_map_file_name)
        #relaxed_gene_tree_results.add_back_outgroup(simulation_leg.leg_length())
        relaxed_gene_tree_results_by_gene_tree[gene_tree] =relaxed_gene_tree_results

        plot_file_name_1= full_path_to_relaxed_tree_file +"_phylo.png"
        plot_file_name_2= os.path.join(subfolder,gene_tree +"_specks.png")
        print("newick to plot:\t" +relaxed_gene_tree_results.simple_newick)
        tree_visuals_by_phylo.save_tree_plot_from_newick(relaxed_gene_tree_results.simple_newick, plot_file_name_1)

        new_tree_file=full_path_to_relaxed_tree_file.replace(".tree",".updated.tree")
        with open(new_tree_file, 'w') as f:
            f.writelines(relaxed_gene_tree_results.simple_newick + "\n")

        gt_newick=relaxed_gene_tree_results.simple_newick
        leaf_map = relaxed_gene_tree_results.leaves_by_species
        gt_tree_viz_data=gene_tree_visuals.plot_polyploid_gene_tree_alone(
            simulation_leg, leaf_map,gt_newick, gene_tree, polyploid.SPC_time_MYA,
            polyploid.species_name, plot_file_name_2)
        gt_tree_viz_data_by_gene_tree[gene_tree]=gt_tree_viz_data

    gene_tree_visuals.histogram_node_distances(polyploid, gt_tree_viz_data_by_gene_tree,
                                               "relaxed",  subfolder)
    gene_tree_visuals.plot_gene_trees_on_top_of_species_trees(polyploid, gt_tree_viz_data_by_gene_tree,
                                                              "relaxed",  subfolder)
    polyploid.analysis_step_num=polyploid.analysis_step_num+1
    return relaxed_gene_tree_results_by_gene_tree
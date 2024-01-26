import os
import tree_visuals
import gene_tree_maker
import common


def relax(config,gene_tree_results_by_tree_name):

    out_dir = config.output_folder
    subfolder = os.path.join(out_dir, str(config.sim_step_num) + "_relaxed_gene_trees")

    print(subfolder)
    if not os.path.exists(subfolder):
        os.makedirs(subfolder)

    relaxed_gene_tree_results_by_gene_tree={}
    for gene_tree, gene_tree_results in gene_tree_results_by_tree_name.items():

        out_file_name=gene_tree +".relaxed.tree"
        in_file_name =gene_tree_results.gene_tree_file_name
        cmd= ["java","-jar",
             config.path_to_sagephy,"BranchRelaxer",
              "-x","-innms","-o" ,out_file_name, in_file_name,
              "ACRY07","1","0.000001"]

        common.run_and_wait_on_process(cmd, subfolder)

        full_path_out_file=os.path.join(subfolder,out_file_name)
        relaxed_gene_tree_results=gene_tree_maker.read_tree_file(full_path_out_file)
        relaxed_gene_tree_results.info_dict=gene_tree_results.info_dict.copy()
        relaxed_gene_tree_results.leaves_by_species =gene_tree_results.leaves_by_species.copy()
        relaxed_gene_tree_results_by_gene_tree[gene_tree] =relaxed_gene_tree_results

        plot_file_name= full_path_out_file +".png"
        print("newick to plot:\t" +relaxed_gene_tree_results.simple_newick)
        tree_visuals.save_tree_plot(relaxed_gene_tree_results.simple_newick, plot_file_name)

    config.sim_step_num=config.sim_step_num+1
    return relaxed_gene_tree_results_by_gene_tree
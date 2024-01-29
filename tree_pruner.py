import os


def prune_gene_trees(polyploid):

    config = polyploid.general_sim_config
    subfolder=os.path.join(polyploid.species_subfolder, str(polyploid.analysis_step_num) + "_post_WGD_gene_pruning")
    if not os.path.exists(subfolder):
        os.makedirs(subfolder)

    print("TODO: after WGD event, remove x% of genes per MY")
    polyploid.analysis_step_num=polyploid.analysis_step_num+1
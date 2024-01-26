import os


def prune_gene_trees(config):

    out_dir = config.output_folder
    subfolder = os.path.join(out_dir, str(config.sim_step_num) +"_post_WGD_gene_pruning")
    if not os.path.exists(subfolder):
        os.makedirs(subfolder)

    print("TODO: after WGD event, remove x% of genes per MY")
    config.sim_step_num=config.sim_step_num+1
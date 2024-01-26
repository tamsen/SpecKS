import os


def prune_gene_trees(config,step_num):

    out_dir = config.output_folder
    subfolder = os.path.join(out_dir, str(step_num) +"_gene_pruning")
    if not os.path.exists(subfolder):
        os.makedirs(subfolder)

    print("TODO: after WGD event, remove x% of genes per MY")
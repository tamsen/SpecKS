import os


def prune_gene_trees(config):

    out_dir = config.output_folder
    subfolder = os.path.join(out_dir, "4_gene_pruning")
    if not os.path.exists(subfolder):
        os.makedirs(subfolder)

    print("TODO: after WGD event, remove x% of genes per MY")
import os


def relax(config,species_tree, num_gene_trees, step_num):

    out_dir = config.output_folder
    subfolder = os.path.join(out_dir, str(step_num) + "_relax_gene_trees")

    print(subfolder)
    if not os.path.exists(subfolder):
        os.makedirs(subfolder)
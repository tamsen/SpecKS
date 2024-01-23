import glob
import os
import subprocess
import matplotlib.pyplot as plt
import tree_visuals
from scipy.stats import beta


def write_SaGePhy_GuestTreeGen_commands(species_tree_newick, dup_rate, loss_rate, out_file_name):


    cmd = ["java", "-jar", "/home/tamsen/Apps/sagephy/sagephy-1.0.0.jar",
         "GuestTreeGen", species_tree_newick,
         str(dup_rate), str(loss_rate), "0.0",out_file_name, "-nox"]

    #"nox" is for "no auxillary tags" as per the manual

    print(cmd)
    return cmd

def get_randomized_dup_and_loss_rates(dup_rate_parameters,loss_rate_parameters,num_values_needed):
    dup_values = beta.rvs(dup_rate_parameters[0],dup_rate_parameters[1], size=num_values_needed)
    loss_values = beta.rvs(loss_rate_parameters[0],loss_rate_parameters[1], size=num_values_needed)
    return dup_values,loss_values

def visualize_dup_and_loss_rates(dup_values,loss_values,out_folder):

    xs=[i for i in range(0,len(dup_values))]
    plt.scatter(xs,dup_values, label="duplication rate",marker="o")
    plt.scatter(xs,loss_values, label="loss rate",marker="o")

    plt.title("Gene duplication and loss rates for each gene tree")
    plt.legend()
    file_to_save=os.path.join(out_folder,"Gene duplication and loss rates")
    plt.savefig(file_to_save)
    plt.cla()
    plt.close()

def read_pruned_trees(subfolder):

    tree_files = glob.glob(subfolder + "/*.pruned.tree")
    newicks_by_file={}

    for tree_file in tree_files:
        with open(tree_file, "r") as f:
            lines = f.readlines()
            newicks_by_file[tree_file]=lines[0]

    return newicks_by_file

def run_sagephy(config, species_tree_newick, num_gene_trees_needed):

    out_dir = config.output_folder
    dup_rate_parameters = config.dup_rate_parameters
    loss_rate_parameters = config.loss_rate_parameters
    subfolder = os.path.join(out_dir, "gene_trees")

    print(subfolder)
    if not os.path.exists(subfolder):
        os.makedirs(subfolder)

    dup_values, loss_values = get_randomized_dup_and_loss_rates(
        dup_rate_parameters, loss_rate_parameters, num_gene_trees_needed)

    visualize_dup_and_loss_rates(dup_values, loss_values, subfolder)

    for i in range(0, num_gene_trees_needed):
        out_file_name = "GeneTree" + str(i)
        cmd = write_SaGePhy_GuestTreeGen_commands(species_tree_newick,
                                                  dup_values[i], loss_values[i],
                                                  os.path.join(subfolder,out_file_name))

        result = subprocess.run(cmd, capture_output=True)
        print("result:\t" + str(result))

    newicks_by_file = read_pruned_trees(subfolder)
    pruned_tree_files=newicks_by_file.keys()
    print("pruned_tree_files:\t" + str(pruned_tree_files))

    for pruned_tree_file in pruned_tree_files:
        plot_file_name=pruned_tree_file.replace(".tree",".png")
        tree_visuals.save_tree_plot(newicks_by_file[pruned_tree_file], plot_file_name)

    return newicks_by_file

#https://docs.python.org/3/library/subprocess.html
#https://stackoverflow.com/questions/8953119/waiting-for-external-launched-process-finish
#https://stackoverflow.com/questions/14762048/subprocess-call-in-python-to-invoke-java-jar-files-with-java-opts
#java -jar /home/tamsen/Apps/sagephy/sagephy-1.0.0.jar BranchRelaxer -x -innms -o genetree.relaxed.tree genetree_test1.unpruned.tree ACRY07 1 0.000001
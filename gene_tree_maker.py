import glob
import os
import common
import matplotlib.pyplot as plt
import tree_visuals
from scipy.stats import beta


def write_SaGePhy_GuestTreeGen_commands(config, species_tree_newick, dup_rate, loss_rate, out_file_name):


    cmd = ["java", "-jar", config.path_to_sagephy,
         "GuestTreeGen", species_tree_newick,
         str(dup_rate), str(loss_rate), "0.0",out_file_name]

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
    #results_by_file={}
    results_by_tree_name = {}
    for tree_file in tree_files:

        #TODO this could really be cleaned up into a single constructor
        leafmap_file=tree_file.replace(".tree",".leafmap")
        info_file=tree_file.replace(".tree",".info")

        result = read_tree_file(tree_file)
        result = read_leaf_map(leafmap_file, result)
        result = read_gene_tree_info(info_file, result)
        results_by_tree_name[result.gene_tree_name]=result

    return results_by_tree_name


def read_tree_file(tree_file):
    with open(tree_file, "r") as f:
        lines = f.readlines()
        clean_tree_string = clean_newick(lines[0])
        full_tree_string = lines[0]

    gene_tree_name=os.path.basename(tree_file).replace(".pruned.tree","")
    result = gene_tree_result()
    result.gene_tree_file_name=tree_file
    result.simple_newick=clean_tree_string
    result.newick_with_tags=full_tree_string
    result.gene_tree_name=gene_tree_name
    return result


def read_leaf_map(leafmap_file, my_gene_tree_result):

    leaves_by_species={}
    num_extant_leaves=0
    with open(leafmap_file, "r") as f:
            lines = f.readlines()
            for line in lines:
                splat=line.split()
                leaf=splat[0]
                species=splat[1]
                if species not in leaves_by_species.keys():
                    leaves_by_species[species]=[]
                leaves_by_species[species].append(leaf)
                num_extant_leaves=num_extant_leaves+1

    my_gene_tree_result.leaves_by_species=leaves_by_species
    my_gene_tree_result.num_extant_leaves=num_extant_leaves
    return my_gene_tree_result

def read_gene_tree_info(gene_tree_info_file, my_gene_tree_result):

    info_dict={}
    with open(gene_tree_info_file, "r") as f:
            lines = f.readlines()
            for line in lines:
                if "#" in line:
                    continue
                if "Arguments" in line:
                    continue

                splat=line.split(":")
                if len(splat) < 2:
                    continue

                splat=line.split(":")
                info_dict[splat[0]]=splat[1].strip()

    my_gene_tree_result.info_dict=info_dict
    return my_gene_tree_result
def clean_newick(taggy_tree):

        open_brackets = [i for i in range(0, len(taggy_tree)) if taggy_tree[i] == "["]
        if len(open_brackets) ==0:
            return taggy_tree

        close_brackets_plus_zero = [0] + [i + 1 for i in range(0, len(taggy_tree)) if taggy_tree[i] == "]"]
        stuff_to_keep = [taggy_tree[close_brackets_plus_zero[i]:open_brackets[i]] for i in range(0, len(open_brackets))]
        clean_tree = "".join(stuff_to_keep) +";"
        return clean_tree
def run_sagephy(config, species_tree_newick, num_gene_trees_needed, step_num):

    out_dir = config.output_folder
    dup_rate_parameters = config.dup_rate_parameters
    loss_rate_parameters = config.loss_rate_parameters
    subfolder = os.path.join(out_dir, str(step_num) + "_gene_trees")

    print(subfolder)
    if not os.path.exists(subfolder):
        os.makedirs(subfolder)

    dup_values, loss_values = get_randomized_dup_and_loss_rates(
        dup_rate_parameters, loss_rate_parameters, num_gene_trees_needed)

    visualize_dup_and_loss_rates(dup_values, loss_values, subfolder)

    for i in range(0, num_gene_trees_needed):
        out_file_name = "GeneTree" + str(i)
        cmd = write_SaGePhy_GuestTreeGen_commands(config, species_tree_newick,
                                                  dup_values[i], loss_values[i],
                                                  os.path.join(subfolder,out_file_name))
        common.run_and_wait_on_process(cmd, subfolder)


    gene_tree_data_by_tree_name = read_pruned_trees(subfolder)
    pruned_tree_names=gene_tree_data_by_tree_name.keys()
    print("pruned_tree_files:\t" + str(pruned_tree_names))

    for pruned_tree_file in pruned_tree_names:
        plot_file_name= os.path.join(subfolder,pruned_tree_file +".png")
        print("newick to plot:\t" +gene_tree_data_by_tree_name[pruned_tree_file].simple_newick)
        tree_visuals.save_tree_plot(gene_tree_data_by_tree_name[pruned_tree_file].simple_newick, plot_file_name)

    return gene_tree_data_by_tree_name

class gene_tree_result():
    simple_newick=""
    newick_with_tags=""
    num_extant_leaves=0
    gene_tree_name=""
    gene_tree_file_name=""
    leaves_by_species={}
    info_dict={}
#https://docs.python.org/3/library/subprocess.html
#https://stackoverflow.com/questions/8953119/waiting-for-external-launched-process-finish
#https://stackoverflow.com/questions/14762048/subprocess-call-in-python-to-invoke-java-jar-files-with-java-opts
#java -jar /home/tamsen/Apps/sagephy/sagephy-1.0.0.jar BranchRelaxer -x -innms -o genetree.relaxed.tree genetree_test1.unpruned.tree ACRY07 1 0.000001
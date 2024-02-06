import glob
import os

from Bio import Phylo
from io import StringIO
import common
import matplotlib.pyplot as plt
from visualization import tree_visuals_by_phylo
from visualization import gene_tree_visuals
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

def read_pruned_trees(subfolder,full_sim_time):

    tree_files = glob.glob(subfolder + "/*.pruned.tree")
    results_by_tree_name = {}
    for tree_file in tree_files:

        #TODO this could really be cleaned up into a single constructor
        leafmap_file=tree_file.replace(".tree",".leafmap")
        info_file=tree_file.replace(".tree",".info")

        result = read_tree_file(tree_file)
        result = read_leaf_map(leafmap_file, result)
        result = read_gene_tree_info(info_file, result)
        result = add_back_outgroup(result, full_sim_time)
        results_by_tree_name[result.gene_tree_name]=result

    return results_by_tree_name

def unprune_outgroup(newick_1, full_sim_time):

        tree_1 = Phylo.read(StringIO(newick_1), "newick")
        terminals = tree_1.get_terminals()

        #this is the pruned tree, so all terminals should go the length of the tree
        if len(terminals)>0:
            dist = tree_1.distance(terminals[0])
        else:
            return "(O1:500,O2:500)",["O1","O2"] #I guess all the nodes died..

        missing_distance = full_sim_time - dist
        prefix = "(O1:500,"
        suffix = ":" + str(missing_distance) + ");"
        adjusted_newick = prefix + newick_1.replace(";", suffix)
        return adjusted_newick,["O1"]

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

def add_back_outgroup(result,full_sim_time):
    original_newick = result.simple_newick
    new_newick, new_leaves=unprune_outgroup(result.simple_newick,full_sim_time)
    result.simple_newick = new_newick
    result.original_newick = original_newick

    for new_leaf in new_leaves:
        if not "O" in result.leaves_by_species:
            result.leaves_by_species["O"]=[]
        result.leaves_by_species["O"].append(new_leaf)

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
def run_sagephy(polyploid, species_tree_newick):

    config = polyploid.general_sim_config
    num_gene_trees_needed = config.num_gene_trees_per_species_tree
    dup_rate_parameters = config.dup_rate_parameters
    loss_rate_parameters = config.loss_rate_parameters

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
        out_file_name = "GeneTree" + str(i)
        cmd = write_SaGePhy_GuestTreeGen_commands(config, species_tree_newick,
                                                  dup_values[i], loss_values[i],
                                                  os.path.join(subfolder,out_file_name))
        common.run_and_wait_on_process(cmd, subfolder)


    gene_tree_data_by_tree_name = read_pruned_trees(subfolder, polyploid.FULL_time_MYA)
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
        gt_newick=gene_tree_data_by_tree_name[gt_name].simple_newick
        leaf_map = gene_tree_data_by_tree_name[gt_name].leaves_by_species
        print("newick to plot:\t" +gt_newick)
        tree_visuals_by_phylo.save_tree_plot(gt_newick, plot_file_name_1)
        gt_tree_viz_data=gene_tree_visuals.plot_gene_tree_alone(
            polyploid.subgenome_names, leaf_map,gt_newick, gt_name, plot_file_name_2)
        gt_tree_viz_data_by_name[gt_name]=gt_tree_viz_data
        i=i+1

    gene_tree_visuals.plot_gene_trees_on_top_of_species_trees(polyploid, config,
                                                              gt_tree_viz_data_by_name, subfolder)
    polyploid.analysis_step_num=polyploid.analysis_step_num+1
    return gene_tree_data_by_tree_name

class gene_tree_result():
    original_newick=""
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
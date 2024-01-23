import os
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import glob
from scipy.stats import beta


# java -jar /home/tamsen/Apps/sagephy/sagephy-1.0.0.jar GuestTreeGen '(O:500,(P1:200,P2:200):300);' 0.00162 0.00194 0.0 genetree_test1

def write_SaGePhy_GuestTreeGen_commands(species_tree_newick, dup_rate, loss_rate, out_file_name):
    cmd = "java -jar /home/tamsen/Apps/sagephy/sagephy-1.0.0.jar GuestTreeGen " + \
          species_tree_newick + " " + str(dup_rate) + " " + str(loss_rate) + " 0.0 " + out_file_name

    print(cmd)
    return cmd

species_tree='(O:500,(P1:200,P2:200):300);'
dup_rate=0.00162
loss_rate=0.00194
out_file_name="foo"
write_SaGePhy_GuestTreeGen_commands(species_tree,dup_rate,loss_rate,out_file_name)

dup_rate_parameters=(4, 2460)
loss_rate_parameters=(4,2053)

def get_randomized_dup_and_loss_rates(dup_rate_parameters,loss_rate_parameters,num_values_needed):
    dup_values = beta.rvs(dup_rate_parameters[0],dup_rate_parameters[1], size=num_values_needed)
    loss_values = beta.rvs(loss_rate_parameters[0],loss_rate_parameters[1], size=num_values_needed)
    return dup_values,loss_values

    dup_values,loss_values=get_randomized_dup_and_loss_rates(
        dup_rate_parameters,loss_rate_parameters,10)

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



def write_gene_tree_commands(species_tree_newick, num_gene_trees_needed,out_folder):

    #species_tree_newick = '\'(O:500,(P1:200,P2:200):300);\''

    dup_values, loss_values = get_randomized_dup_and_loss_rates(
        dup_rate_parameters, loss_rate_parameters, num_gene_trees_needed)

    visualize_dup_and_loss_rates(dup_values, loss_values, out_folder)

    for i in range(0, num_gene_trees_needed):
        out_file_name = "GeneTree" + str(i)
        cmd = write_SaGePhy_GuestTreeGen_commands(species_tree_newick,
                                                  dup_values[i], loss_values[i], out_file_name)


#java -jar /home/tamsen/Apps/sagephy/sagephy-1.0.0.jar BranchRelaxer -x -innms -o genetree.relaxed.tree genetree_test1.unpruned.tree ACRY07 1 0.000001
import os
import numpy as np
import matplotlib.pyplot as plt
import glob
import shutil

def get_Ks_from_file(paml_out_file):

    KS_values = []
    with open(paml_out_file, "r") as f:
        lines = f.readlines()
        for l in lines:
            data = l.split()
            for d in data[1:len(data)]:
                KS_values.append(float(d))
    return KS_values


def get_Ks_from_folder(paml_out_folder, alg_name):

    res_files = glob.glob(paml_out_folder + "/*"+alg_name+".dS")
    KS_values = []

    with open(os.path.join(paml_out_folder,"Ks_by_GeneTree.csv"), 'w') as f:

        f.writelines("GeneTree,Ks\n")
        for paml_out_file in res_files:
            print(paml_out_file)
            Ks_for_og = get_Ks_from_file(paml_out_file)
            KS_values = KS_values + Ks_for_og
            Ks_string = [str(ks) for ks in KS_values]
            f.writelines(paml_out_file+","+ ",".join(Ks_string)+ "\n")

    return KS_values


def plot_Ks_histogram(PAML_hist_out_file, species_name, Ks_results, max_Ks, max_y, alg_name, color, bin_size):

    fig = plt.figure(figsize=(10, 10), dpi=100)
    bins = np.arange(0, max_Ks + 0.1, bin_size)
    x = Ks_results
    print(PAML_hist_out_file)

    n, bins, patches = plt.hist(x, bins=bins, facecolor=color, alpha=0.25, label='histogram data')

    plt.xlim([0, max_Ks * (1.1)])
    if max_y:
        plt.ylim([0, max_y])

    plt.xlabel("Ks")
    plt.ylabel("Count in Bin")
    plt.title("Ks histogram for " + species_name + ", last ~" + str(max_Ks * 100) + " MY\n" +
              "algorithm: PAML " + alg_name)
    plt.savefig(PAML_hist_out_file)
    plt.clf()
    plt.close()

def run_Ks_histogramer(polyploid,codeml_results_by_replicate_num,gene_tree_results_by_file_name):

    config = polyploid.general_sim_config
    subfolder = os.path.join(polyploid.species_subfolder, str(polyploid.analysis_step_num) + "_ks_histograms")

    if not os.path.exists(subfolder):
        os.makedirs(subfolder)

    replicates = list(codeml_results_by_replicate_num.keys())
    for replicate in replicates:
        rep_subfolder = os.path.join(subfolder, "replicate" + str(replicate))
        if not os.path.exists(rep_subfolder ):
            os.makedirs(rep_subfolder)

        ks_result_files=codeml_results_by_replicate_num[replicate]
        gene_trees=ks_result_files.keys()
        for gene_tree in gene_trees:

            codeml_result=ks_result_files[gene_tree]
            files=[codeml_result.ML_file,codeml_result.NG_file]
            for file in files:
                base=os.path.basename(file)
                dst=os.path.join(rep_subfolder, gene_tree + "_" +base)
                if os.path.exists(file):
                    shutil.copyfile(file, dst)
                else:
                    print("Warning: no codeml result for " + gene_tree + ", replicate " + str(replicate))
                    print("Maybe gene tree had no extant leaves? Codeml result file is here: "
                          + file)

        subgenomes_of_concern=["P1","P2"]
        species_str="between subgenomes " + "and".join(subgenomes_of_concern)
        print("making histograms " + species_str + ", replicate " + str(replicate))
        summarize_ks(rep_subfolder, subgenomes_of_concern,
                     config.max_ks_for_hist_plot, config.max_y_for_hist_plot,"pink", 0.1)
    polyploid.analysis_step_num=polyploid.analysis_step_num+1
def summarize_ks(paml_out_folder, subgenomes_of_concern,max_ks, max_y,color, step):

    for alg_name in ["ML","NG"]:

        print("getting results for PAML alg name " + alg_name)
        species_name="Polyploid " + "and".join(subgenomes_of_concern)
        ks_results = get_Ks_from_folder(paml_out_folder, alg_name)
        paml_hist_file = os.path.join(paml_out_folder, species_name + "_paml_hist_maxKS" + str(max_ks) +
                                  "_"+ alg_name + ".png")

        plot_Ks_histogram(paml_hist_file, species_name, ks_results, max_ks, max_y,
                                                alg_name, color, step)

    return





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
    res = glob.glob(paml_out_folder + "/*")
    og_list = [r for r in res if os.path.isdir(r)]
    # og_list = [og for og in os.listdir(paml_out_folder) if os.path.isdir(og)]
    KS_values = []
    for og in og_list:
        # print(og)
        # Ks_file="2NG.dS" or "2ML.dS"
        paml_out_file = os.path.join(paml_out_folder, og, "2" + alg_name + ".dS")
        print(paml_out_file)
        Ks_for_og = get_Ks_from_file(paml_out_file)
        KS_values = KS_values + Ks_for_og
    return KS_values


def plot_Ks_histogram(PAML_hist_out_file, species_name, Ks_results, max_Ks, max_y, alg_name, color, step):
    bins = np.arange(0, max_Ks + 0.1, step)
    x = Ks_results
    print(PAML_hist_out_file)

    n, bins, patches = plt.hist(x, bins=bins, facecolor=color, alpha=0.25, label='histogram data')

    plt.xlim([0, max_Ks * (1.1)])
    if not (max_y == False):
        plt.ylim([0, max_y])

    plt.xlabel("Ks")
    plt.ylabel("Count in Bin")
    plt.title("Ks histogram for " + species_name + ", last ~" + str(max_Ks * 100) + " MY\n" +
              "algorithm: PAML " + alg_name)
    plt.savefig(PAML_hist_out_file)
    return [n, bins, patches, plt]

def run_Ks_histogramer(config,codeml_results_by_replicate_num):
    out_dir = config.output_folder
    subfolder = os.path.join(out_dir, "6_ks_histograms")
    if not os.path.exists(subfolder):
        os.makedirs(subfolder)

    replicates = list(codeml_results_by_replicate_num.keys())
    for replicate in replicates:
        rep_subfolder = os.path.join(subfolder, "replicate" + str(replicate))
        if not os.path.exists(rep_subfolder ):
            os.makedirs(rep_subfolder )

        ks_result_files=codeml_results_by_replicate_num[replicate]
        gene_trees=ks_result_files.keys()
        for gene_tree in gene_trees:
            codeml_result=ks_result_files[gene_tree]
            files=[codeml_result.ML_file,codeml_result.NG_file]
            for file in files:
                base=os.path.basename(file)
                dst=os.path.join(rep_subfolder, gene_tree + "_" +base)
                shutil.copyfile(file, dst)

def summarize_Ks(paml_out_folder, species_name, max_ks, max_y, alg_name, color, step):
    print("making histograms")
    #ks_results = get_Ks_from_folder(paml_out_folder, alg_name)
    #paml_hist_file = os.path.join(paml_out_folder, species_name + "_paml_hist_maxKS" + str(max_ks) +
    #                              alg_name + ".png")

    #[n, bins, patches, plt] = plot_Ks_histogram(paml_hist_file, species_name, ks_results, max_ks, max_y,
    #                                            alg_name, color, step)

    #return [n, bins, patches, plt]

    #print("making histograms")
#my_paml_out_folder="/home/tamsen/Data/Ks_Genome_Simulator/Codeml"
#summarize_Ks(my_paml_out_folder,"my_species_name",5,False,"NG", 'c', 0.01)



import os
import numpy as np
import matplotlib.pyplot as plt
import glob
import shutil

from pipeline_modules import ks_calculator


def get_Ks_from_file(paml_out_file):

    ks_results =[]
    ordered_ortholog_list = []
    with open(paml_out_file, "r") as f:
        lines = f.readlines()

    row_index=0
    for l in lines[1:len(lines)]:
            data = l.split()
            if len(data)==0:
                continue

            #print("line: " + l)
            ordered_ortholog_list.append(data[0])
            for col_index in range(1,len(data)):
                d=data[col_index]
                new_ks_result=ks_data(float(d),row_index,col_index-1)
                ks_results.append(new_ks_result)
            row_index = row_index+1

    for ks_result in ks_results:
        ks_result.set_orthologs(ordered_ortholog_list)

    return ks_results


def get_Ks_from_folder(paml_out_folder, replicate, alg_name):

    res_files = glob.glob(paml_out_folder + "/*"+alg_name+".dS")
    csv_file_out=os.path.join(paml_out_folder,alg_name+ "_rep" +str(replicate) +
                                           "_Ks_by_GeneTree.csv")
    KS_values = []

    with open(csv_file_out, 'w') as f:

        f.writelines("GeneTree,Ks,full_path\n")
        for paml_out_file in res_files:
            print(paml_out_file)
            base_name=os.path.basename(paml_out_file)
            Ks_for_og = get_Ks_from_file(paml_out_file)
            for Ks_data in Ks_for_og:
                Ks_value = Ks_data.ks_between_orthologs
                ortholog_names_str = str(Ks_data.ortholog_pair).replace(","," ")
                f.writelines(base_name + "," +  ortholog_names_str + "," +
                             str(Ks_value)  +","+paml_out_file + "\n")
                KS_values.append(Ks_value)

    return KS_values,csv_file_out


def plot_Ks_histogram(PAML_hist_out_file, species_name, Ks_results, WGD_as_Ks, SPEC_as_Ks,
                      max_Ks, max_y, alg_name, color, bin_size):

    fig = plt.figure(figsize=(10, 10), dpi=100)
    #bins = np.arange(0, max_Ks + 0.1, bin_size)
    nBins=50
    x = Ks_results
    print(PAML_hist_out_file)

    n, bins, patches = plt.hist(x, bins=nBins, facecolor=color, alpha=0.25, label='histogram data')

    if max_Ks:
        plt.xlim([0, max_Ks * (1.1)])

    if max_y:
        plt.ylim([0, max_y])

    plt.axvline(x=WGD_as_Ks, color='r', linestyle='--', label="WGD time as Ks")
    plt.axvline(x=SPEC_as_Ks, color='b', linestyle='--', label="SPEC time as Ks")
    plt.legend()
    plt.xlabel("Ks")
    plt.ylabel("Count in Bin")
    plt.title("Ks histogram for " + species_name + ", last ~" + str(max_Ks * 100) + " MY\n" +
              "algorithm: PAML " + alg_name)
    plt.savefig(PAML_hist_out_file)
    plt.clf()
    plt.close()

def run_Ks_histogramer(polyploid,codeml_results_by_replicate_num):

    config = polyploid.general_sim_config
    subfolder = os.path.join(polyploid.species_subfolder, str(polyploid.analysis_step_num) + "_ks_histograms")
    results_files_by_replicate={}

    if not os.path.exists(subfolder):
        os.makedirs(subfolder)

    replicate_index_formatter = ks_calculator.get_replicate_index_format(config.num_replicates_per_gene_tree)
    replicates = list(codeml_results_by_replicate_num.keys())
    for replicate in replicates:

        replicate_string = replicate_index_formatter.format(replicate)
        rep_subfolder = os.path.join(subfolder, "replicate_" + replicate_string)
        if not os.path.exists(rep_subfolder ):
            os.makedirs(rep_subfolder)

        ks_result_files=codeml_results_by_replicate_num[replicate]
        gene_trees=ks_result_files.keys()
        for gene_tree in gene_trees:

            codeml_result=ks_result_files[gene_tree]
            files=[codeml_result.ML_file,codeml_result.NG_file]
            for file in files:
                base=os.path.basename(file)
                dst=os.path.join(rep_subfolder, gene_tree + "_" + base)
                if os.path.exists(file):
                    shutil.copyfile(file, dst)
                else:
                    print("Warning: no codeml result for " + gene_tree + ", replicate " + replicate_string )
                    print("Maybe gene tree had no extant leaves? Codeml result file should be here: "
                          + file)

        subgenomes_of_concern=["P1","P2"]
        species_str="between subgenomes " + "and".join(subgenomes_of_concern)
        print("making histograms " + species_str + ", replicate " + replicate_string )
        WGD_as_Ks=polyploid.WGD_time_as_ks()
        SPEC_as_Ks=polyploid.SPEC_time_as_ks()
        outfiles = summarize_ks(rep_subfolder, replicate_string, polyploid.species_name, WGD_as_Ks, SPEC_as_Ks,
                     config.max_ks_for_hist_plot, config.max_y_for_hist_plot,"pink", 0.1)
        results_files_by_replicate[replicate]=outfiles


    polyploid.analysis_step_num=polyploid.analysis_step_num+1
    return results_files_by_replicate
def summarize_ks(paml_out_folder, replicate, species_name, WGD_as_Ks, SPEC_as_Ks, max_ks, max_y,color, step):

    outfiles=[]
    for alg_name in ["ML","NG"]:

        print("getting results for PAML alg name " + alg_name)
        ks_results,csv_file_out = get_Ks_from_folder(paml_out_folder, replicate, alg_name)
        paml_hist_file = os.path.join(paml_out_folder, species_name + "_rep" + str(replicate) +
                                      "_paml_hist_maxKS" + str(max_ks) +
                                  "_"+ alg_name + ".png")

        plot_Ks_histogram(paml_hist_file, species_name, ks_results,WGD_as_Ks, SPEC_as_Ks, max_ks, max_y,
                                                alg_name, color, step)

        paml_hist_file = os.path.join(paml_out_folder, species_name + "_rep" + str(replicate) +
                                      "_paml_hist" +
                                  "_"+ alg_name + ".png")
        plot_Ks_histogram(paml_hist_file, species_name, ks_results,WGD_as_Ks, SPEC_as_Ks, False, max_y,
                                                alg_name, color, step)

        outfiles.append(csv_file_out)

    return outfiles



class ks_data():

    ks_between_orthologs=0
    ortholog_pair=[]
    row_index=-1
    col_index=-1
    def __init__(self, ks_data_result,row_index,col_index):
        self.ks_between_orthologs = ks_data_result
        self.ortholog_pair = []
        self.row_index = row_index
        self.col_index = col_index

    def set_orthologs(self, ordered_ortholog_list):

        ortholog1=ordered_ortholog_list[self.row_index]
        ortholog2=ordered_ortholog_list[self.col_index]
        self.ortholog_pair = [ortholog1,ortholog2]
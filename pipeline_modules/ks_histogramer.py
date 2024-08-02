import os
import numpy as np
import matplotlib.pyplot as plt
import glob
import shutil
import log
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


def get_Ks_from_folder(paml_out_folder, replicate, alg_name, version_string):

    res_files = glob.glob(paml_out_folder + "/*"+alg_name+".dS")
    csv_file_out=os.path.join(paml_out_folder,alg_name+ "_rep" +str(replicate) +
                                           "_Ks_by_GeneTree.csv")
    KS_values = extract_K_values(csv_file_out, res_files, version_string)

    return KS_values,csv_file_out
def get_Kn_from_folder(paml_out_folder, replicate, alg_name, version_string):

    res_files = glob.glob(paml_out_folder + "/*"+alg_name+".dN")
    csv_file_out=os.path.join(paml_out_folder,alg_name+ "_rep" +str(replicate) +
                                           "_Kn_by_GeneTree.csv")
    KS_values = extract_K_values(csv_file_out, res_files, version_string)

    return KS_values,csv_file_out

def extract_K_values(csv_file_out, res_files, version_string):
    KS_values = []
    with open(csv_file_out, 'w') as f:

        f.writelines("SpecKS " + version_string + "\n")  # version_info.to_string())
        f.writelines("GeneTree,leaf names,Ks,path to original output file\n")
        for paml_out_file in res_files:
            log.write_to_log(paml_out_file)
            base_name = os.path.basename(paml_out_file)
            Ks_for_og = get_Ks_from_file(paml_out_file)
            for Ks_data in Ks_for_og:
                Ks_value = Ks_data.ks_between_ortholog_and_LCA
                ortholog_names_str = str(Ks_data.ortholog_pair).replace(",", " ")
                f.writelines(base_name + "," + ortholog_names_str + "," +
                             str(Ks_value) + "," + paml_out_file + "\n")
                KS_values.append(Ks_value)
    return KS_values


def plot_Ks_histogram(PAML_hist_out_file, species_name, Ks_results, WGD_as_Ks, SPEC_as_Ks,
                      max_Ks, max_y, alg_name, color, bin_size):

    fig = plt.figure(figsize=(10, 10), dpi=100)
    x = Ks_results

    if max_Ks:
        bins = np.arange(0, max_Ks + 0.1, bin_size)
        log.write_to_log(PAML_hist_out_file)
        n, bins, patches = plt.hist(x, bins=bins, facecolor=color, alpha=0.25, label='histogram data')
        plt.xlim([0, max_Ks * (1.1)])
    else:
        nBins=50
        n, bins, patches = plt.hist(x, bins=nBins, facecolor=color, alpha=0.25, label='histogram data')
    if max_y:
        plt.ylim([0, max_y])

    plt.axvline(x=WGD_as_Ks, color='b', linestyle='-', label="WGD time as Ks")
    plt.axvline(x=SPEC_as_Ks, color='r', linestyle='--', label="SPEC time as Ks")
    plt.legend()
    plt.xlabel("Ks")
    plt.ylabel("Count in Bin")
    plt.title("Ks histogram for " + species_name + ", last ~" + str(max_Ks * 100) + " MY\n" +
              "algorithm: PAML " + alg_name)
    plt.savefig(PAML_hist_out_file)
    plt.clf()
    plt.close()

def run_Ks_histogramer(polyploid,genomes_of_interest_by_species,Ks_results_by_species_by_replicate_num):

    config = polyploid.general_sim_config
    specks_repo= config.version_info.repo_url
    specks_version_num= config.version_info.version_num
    specks_version = str(specks_version_num) +"," +specks_repo

    subfolder = os.path.join(polyploid.species_subfolder, str(polyploid.analysis_step_num) + "_ks_histograms")

    if not os.path.exists(subfolder):
        os.makedirs(subfolder)

    replicate_index_formatter = ks_calculator.get_replicate_index_format(config.num_replicates_per_gene_tree)

    results_files_by_species_by_replicate = {}
    for species, subgenomes in genomes_of_interest_by_species.items():

        codeml_results_by_replicate_num=Ks_results_by_species_by_replicate_num[species]
        replicates = list(codeml_results_by_replicate_num.keys())
        species_subfolder = os.path.join(subfolder,species)
        results_files_for_species_by_replicate = {}

        if not os.path.exists(species_subfolder):
            os.makedirs(species_subfolder)

        for replicate in replicates:

            replicate_string = replicate_index_formatter.format(replicate)
            rep_subfolder = os.path.join(species_subfolder, "replicate_" + replicate_string)
            if not os.path.exists(rep_subfolder ):
                os.makedirs(rep_subfolder)

            ks_result_files=codeml_results_by_replicate_num[replicate]
            gene_trees=ks_result_files.keys()
            for gene_tree in gene_trees:

                codeml_result=ks_result_files[gene_tree]
                files=[codeml_result.ML_dN_file, codeml_result.ML_dS_file]
                for file in files:
                    base=os.path.basename(file)
                    dst=os.path.join(rep_subfolder, gene_tree + "_" + base)
                    if os.path.exists(file):
                        shutil.copyfile(file, dst)
                    else:
                        log.write_to_log("Warning: no codeml result for " + gene_tree + ", replicate " + replicate_string )
                        log.write_to_log("Maybe gene tree had no extant leaves? Codeml result file should be here: "
                          + file)

            bin_size = 0.1
            subgenomes_of_concern=["P1","P2"]
            species_str="between subgenomes " + "and".join(subgenomes_of_concern)
            log.write_to_log("making histograms " + species_str + ", replicate " + replicate_string )
            WGD_as_Ks=polyploid.WGD_time_as_ks()
            SPEC_as_Ks=polyploid.DIV_time_as_ks()
            outfiles = summarize_ks(rep_subfolder, specks_version, replicate_string, polyploid.species_name, WGD_as_Ks, SPEC_as_Ks,
                     config.max_ks_for_hist_plot, config.max_y_for_hist_plot,"pink", bin_size)
            results_files_for_species_by_replicate[replicate]=outfiles

        results_files_by_species_by_replicate[species] =results_files_for_species_by_replicate

    polyploid.analysis_step_num=polyploid.analysis_step_num+1
    return results_files_by_species_by_replicate
def summarize_ks(paml_out_folder, specks_version, replicate, species_name,
                 WGD_as_Ks, SPEC_as_Ks, max_ks, max_y, color, bin_size):

    outfiles=[]
    #for alg_name in ["ML","NG"]:
    for alg_name in ["ML"]:

        log.write_to_log("getting results for PAML alg name " + alg_name)
        ks_results,ks_csv_file_out = get_Ks_from_folder(paml_out_folder, replicate, alg_name, specks_version)
        kn_results,kn_csv_file_out = get_Kn_from_folder(paml_out_folder, replicate, alg_name, specks_version)

        paml_hist_file = os.path.join(paml_out_folder, species_name + "_rep" + str(replicate) +
                                      "_paml_hist_maxKS" + str(max_ks) +
                                  "_"+ alg_name + ".png")

        plot_Ks_histogram(paml_hist_file, species_name, ks_results, WGD_as_Ks, SPEC_as_Ks, max_ks, max_y,
                          alg_name, color, bin_size)

        paml_hist_file = os.path.join(paml_out_folder, species_name + "_rep" + str(replicate) +
                                      "_paml_hist" +
                                  "_"+ alg_name + ".png")
        plot_Ks_histogram(paml_hist_file, species_name, ks_results, WGD_as_Ks, SPEC_as_Ks, False, max_y,
                          alg_name, color, bin_size)

        outfiles.append(ks_csv_file_out)
        outfiles.append(kn_csv_file_out)

    return outfiles



class ks_data():

    round_trip_ks_between_orthologs=0
    ks_between_ortholog_and_LCA=0
    ortholog_pair=[]
    row_index=-1
    col_index=-1
    def __init__(self, round_trip_ks_data_result, row_index, col_index):
        self.round_trip_ks_between_orthologs = round_trip_ks_data_result
        self.ks_between_ortholog_and_LCA = self.round_trip_ks_between_orthologs * 0.5
        self.ortholog_pair = []
        self.row_index = row_index
        self.col_index = col_index

    def set_orthologs(self, ordered_ortholog_list):

        ortholog1=ordered_ortholog_list[self.row_index]
        ortholog2=ordered_ortholog_list[self.col_index]
        self.ortholog_pair = [ortholog1,ortholog2]
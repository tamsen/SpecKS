import os
import unittest
import numpy as np
from matplotlib import pyplot as plt
import config

class BatchHistogrammer(unittest.TestCase):

    def test_make_histograms_for_batch(self):

        batch_name="sim40_5p0" ##"sim37_N20" #sim37_N0p1,sim37_N5,sim40_1p0,sim40_5p0,sim40_10p0"
        plot_title=("Simulation with custom GBD model, \n" +\
                    "with Ne-driven allopolyploid ortholog divergence ({0})".format(batch_name))

        input_folder=os.path.join( "/home/tamsen/Data/Specks_outout_from_mesx",batch_name)
        output_folder=os.path.join(input_folder,"analysis")
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)

        print("Reading csv files..")
        csvfiles_by_polyploid_by_species_rep_by_algorithm = get_ks_data_from_folders(input_folder)
        params_by_polyploid = get_truth_from_name_dict(csvfiles_by_polyploid_by_species_rep_by_algorithm)
        example_sim=list(csvfiles_by_polyploid_by_species_rep_by_algorithm.keys())[0]
        species=list(csvfiles_by_polyploid_by_species_rep_by_algorithm[example_sim].keys())
        replicates=list(csvfiles_by_polyploid_by_species_rep_by_algorithm[example_sim][species[0]].keys())
        alg="ML"

        print("Making plots..")
        bin_size = 0.001
        max_Ks = 0.1
        max_Y = False
        for spec in ['polyploid']:#species:
            print(spec)
            for replicate in replicates:
                get_histograms_for_runs_in_batch(output_folder, plot_title,
                                                     csvfiles_by_polyploid_by_species_rep_by_algorithm, spec,
                                                     replicate, alg, params_by_polyploid,
                                                 max_Ks, max_Y, bin_size)

    def test_Single_histogram(self):

        out_folder=("/home/tamsen/Data/SpecKS_output/" +
                    "SpecKS_m04d25y2024_h16m55s36/Allo_Maize/8_final_results")
        csv_file="Allo_Maize_ML_rep0_LCA_to_Ortholog_Ks_by_GeneTree.csv"
        out_png=os.path.join(out_folder,"out.png")
        full_path=os.path.join(out_folder,csv_file)
        Ks_results = read_Ks_csv(full_path)

        bin_size=0.01
        max_Ks=1
        color='blue'

        fig = plt.figure(figsize=(10, 10), dpi=100)
        x = Ks_results
        #print(PAML_hist_out_file)

        if max_Ks:
            bins = np.arange(0, max_Ks + 0.1, bin_size)
            n, bins, patches = plt.hist(x, bins=bins, facecolor=color, alpha=0.25, label='histogram data')
            plt.xlim([0, max_Ks * (1.1)])
        #plt.ylim([0, max_y])

        #plt.axvline(x=WGD_as_Ks, color='b', linestyle='-', label="WGD time as Ks")
        #plt.axvline(x=SPEC_as_Ks, color='r', linestyle='--', label="SPEC time as Ks")
        plt.legend()
        plt.xlabel("Ks")
        plt.ylabel("Count in Bin")
        #plt.title("Ks histogram for " + species_name + ", last ~" + str(max_Ks * 100) + " MY\n" +
        #          "algorithm: PAML " + alg_name)
        plt.savefig(out_png)
        plt.clf()
        plt.close()
def get_truth_from_name_dict(dict_by_name):
    names=dict_by_name.keys()
    return get_truth_from_name_list(names)


def get_truth_from_name_list(names):
    params_by_polyploid = {}
    for name in names:
        print(name)
        spec_time = int(name[7:10])
        wgd_time = int(name[11:14])
        params_by_polyploid[name] = config.PolyploidParams(spec_time, wgd_time, name)
    return params_by_polyploid


def get_ks_data_from_folders(output_folder):
    polyploid_data_folders = os.listdir(output_folder)
    csvfiles_by_polyploid_by_species_rep_by_algorthim = {}
    for polyploid_folder in polyploid_data_folders:
        print(polyploid_folder)

        if not ("Allo" in polyploid_folder) and not ("Auto" in polyploid_folder):
            continue

        csvfiles_by_species_rep_by_algorthim = {}
        full_polyploid_folder_path = os.path.join(output_folder, polyploid_folder)
        if (os.path.isdir(full_polyploid_folder_path)):
            files = os.listdir(full_polyploid_folder_path)
            csvfiles_by_replicate = {}
            for csv_file in files:

                #if "outgroup" in csv_file:
                #    continue

                if not ".csv" in csv_file:
                    continue

                if not "Ks" in csv_file:
                    continue

                # example file name: ML_rep0_Ks_by_GeneTree.csv
                splat = csv_file.split("_")
                len_splat=len(splat)
                if splat[0]== "outgroup":
                    species = "outgroup"
                    algorithm = splat[1]
                    replicate = splat[2]
                else:
                    species = "polyploid"
                    algorithm = splat[2]
                    replicate = splat[3]

                if not species in csvfiles_by_species_rep_by_algorthim:
                    csvfiles_by_species_rep_by_algorthim[species] = {}

                if not replicate in csvfiles_by_species_rep_by_algorthim[species]:

                    csvfiles_by_species_rep_by_algorthim[species][replicate] = {}

                # csvfiles_by_replicate[replicate][algorithm]=csv_file
                full_csv_path = os.path.join(full_polyploid_folder_path, csv_file)
                print("reading " + full_csv_path)
                ks_result_for_file = read_Ks_csv(full_csv_path)
                csvfiles_by_species_rep_by_algorthim[species][replicate][algorithm] = ks_result_for_file

            csvfiles_by_polyploid_by_species_rep_by_algorthim[polyploid_folder] = csvfiles_by_species_rep_by_algorthim
    return csvfiles_by_polyploid_by_species_rep_by_algorthim


def read_Ks_csv(csv_file):

    ks_results = []
    with open(csv_file, "r") as f:

        reading_header=True
        while True:
            line = f.readline()
            if "ersion" in line:
                continue
            if "Git" in line:
                continue
            if not line:
                break
            if len(line)==0:
                break
            if reading_header:
                reading_header=False
                continue
            data = line.split(",")
            #print(data)
            ks_value=float(data[2])
            ks_results.append(ks_value)

    return ks_results

def get_histograms_for_runs_in_batch(out_folder, sample_name, csvfiles_by_polyploid_by_rep_by_algorthim,
                                     spec, replicate, alg, params_by_polyploid, max_Ks_for_x_axis,
                                        max_Y_for_y_axis,bin_size):

    result_names=(list(csvfiles_by_polyploid_by_rep_by_algorthim.keys()))
    result_names.sort()
    ordered_allo1_results=[n for n in result_names if "Allo1" in n]
    ordered_auto_results=[n for n in result_names if "Auto" in n]
    metrics_by_result_names= {}

    # making subplots
    spec_times= [params_by_polyploid[auto_name].SPC_time_MYA for auto_name in ordered_auto_results]
    wgd_offsets= [params_by_polyploid[allo_name].SPC_time_MYA - params_by_polyploid[allo_name].WGD_time_MYA
                  for allo_name in ordered_allo1_results] + [0]
    num_spec_times=len(spec_times)

    num_autos=len(ordered_allo1_results)
    num_wgd_times = num_autos + 1

    fig, ax = plt.subplots(num_spec_times, num_wgd_times,figsize=(10, 10))
    fig.suptitle(sample_name + " for " + replicate +", "+ alg + " algorithm")
    ax[0, 0].set_title("Allopolyploid\n",fontsize=20)
    ax[0, num_autos].set_title("Autopolyploid\n",fontsize=20)

    for sim_idx in range(0, num_spec_times):

        allo_group="Allo"+str(sim_idx+1)
        ordered_allo_results=[n for n in result_names if allo_group in n]
        ordered_allo_results.sort()
        results_for_this_sim_index=ordered_allo_results + [ordered_auto_results[sim_idx]]
        num_omitted_sims= num_wgd_times - len(results_for_this_sim_index)
        for i in range(0, num_omitted_sims):
            results_for_this_sim_index = [False] + results_for_this_sim_index

        for result_idx in range(0,num_wgd_times):

            allo_result_name=results_for_this_sim_index[result_idx]
            if not allo_result_name:
                continue
            this_ax = ax[sim_idx, result_idx]
            csvs_for_allo_result= csvfiles_by_polyploid_by_rep_by_algorthim[allo_result_name]
            ks_for_allo_result= csvs_for_allo_result[spec][replicate][alg]
            params=params_by_polyploid[allo_result_name]
            hist_data=make_histogram_subplot(this_ax, spec, ks_for_allo_result,
                                        bin_size, params,
                                        max_Ks_for_x_axis, max_Y_for_y_axis,"blue")

            out_file_name = os.path.join(out_folder, allo_result_name + ".hist.csv")
            save_hist_to_csv(hist_data, out_file_name)


    for i in range(0,num_wgd_times):
        last_row=num_spec_times-1
        ax[last_row, i].set(xlabel="(wgd offset: " + str(wgd_offsets[i]) + "MYA\n<-- Ks -->")

    for i in range(0,num_spec_times):
        first_col=0
        ax[i, first_col].set(ylabel="<- density ->\nspec: "+ str(spec_times[i]) + "MYA")

    out_file_name=os.path.join(out_folder, "histogram" + "_plot_" + spec +
                               "_" + replicate + "_" + str(max_Ks_for_x_axis) + ".png")
    plt.savefig(out_file_name)
    plt.close()
    return metrics_by_result_names


def make_histogram_subplot(this_ax, spec, Ks_results, bin_size, params,
                           max_Ks, maxY, plot_color):

    WGD_as_Ks = params.WGD_time_MYA * 0.01 #/ 1.04 ..the peak max is about 96% off from where it should be
    SPEC_as_Ks = params.SPC_time_MYA * 0.01 #/ 1.04
    default_xaxis_limit =SPEC_as_Ks + 0.2
    x = Ks_results
    num_gene_pairs_str=str(len(Ks_results))

    bins = np.arange(0, max_Ks + 0.1, bin_size)
    n, bins, patches = this_ax.hist(x, bins=bins, facecolor=plot_color, alpha=0.25,
            label='ks hist ({0} pairs)'.format(num_gene_pairs_str))

    hist_result=[n, bins]
    hist_maximum=max(n)
    ymax_suggestion=hist_maximum*1.6


    this_ax.axvline(x=WGD_as_Ks, color='b', linestyle='-', label="WGD time, " + str(params.WGD_time_MYA) + " MYA")
    this_ax.axvline(x=SPEC_as_Ks, color='r', linestyle='--', label="SPEC time, " + str(params.SPC_time_MYA) + " MYA")


    if maxY:
        yaxis_limit=maxY
        ymax_suggestion=False
    else:
        yaxis_limit= ymax_suggestion

    #this_ax.legend()
    this_ax.set(ylim=[0, yaxis_limit])


    if max_Ks:
        this_ax.set(xlim=[0, max_Ks * 1.1])
        #if max_Ks < 0.51:
        #    this_ax.set(xlim=[0, 5])
    else:
        this_ax.set(xlim=[0, default_xaxis_limit])

    return hist_result


def save_hist_to_csv(hist_result, out_file_name):

    [n, bins]=hist_result
    with open(out_file_name, 'w') as f:
        data_headers= ['bin','# in bin']
        f.writelines(",".join(data_headers) +"\n")

        for i in range(0,len(n)):
            data_list=[str(bins[i]),str(n[i])]
            f.writelines(",".join(data_list) +"\n")

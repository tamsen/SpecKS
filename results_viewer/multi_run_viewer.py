import os
import unittest
from scipy.optimize import curve_fit
import numpy as np
from matplotlib import pyplot as plt

import config
import process_wrapper
from results_viewer import curve_fitting


#https://stackoverflow.com/questions/14770735/how-do-i-change-the-figure-size-with-subplots
class MulitRunViewerTests(unittest.TestCase):

    def test_single_run_viewer(self):

        #test_out_folder="/home/tamsen/Git/SpecKS/SpecKS/test_code/test_out/test_main"
        test_out_folder = "/home/tamsen/Data/Specks_outout_from_mesx/sim19_N1"
        csv_folder="/home/tamsen/Data/Specks_outout_from_mesx/sim19_N1/Auto5"
        csv_file_base0="Auto5_S010W010_ML_rep0_Ks_by_GeneTree.csv"
        csv_file_base1="Allo5_S010W005_ML_rep0_Ks_by_GeneTree.csv"
        csv_file_base2="outgroup_ML_rep0_Ks_by_GeneTree.csv"
        full_csv_path1=os.path.join(csv_folder,csv_file_base1)
        full_csv_path2=os.path.join(csv_folder,csv_file_base2)
        #csv_folder="/home/tamsen/Data/SpecKS_output/SpecKS_m03d15y2024_h15m35s38/Allo1_S150W100/8_final_results"
        #csv_folder = "/home/tamsen/Data/SpecKS_output/SpecKS_m03d15y2024_h16m27s56/Allo1_S070W065/8_final_results"

        full_csv_path1=os.path.join(csv_folder,csv_file_base0)
        full_csv_path2=os.path.join(csv_folder,csv_file_base0)

        bin_size = 0.01
        WGD_time_MYA=5
        SPC_time_MYA=10
        max_Ks_for_x_axis = 1
        csv_files=[full_csv_path1,full_csv_path2]
        plot_titles=['polyploid_with_gene birth and death','outgroup_with_gene birth and death']
        specs=['polyploid','outgroup']
        for i in range(0,len(csv_files)):
            csv_file=csv_files[i]
            plot_title=plot_titles[i]
            self.histogram_a_single_csv_file(specs[i], SPC_time_MYA, WGD_time_MYA, bin_size, max_Ks_for_x_axis,
                                             csv_file, plot_title)

    def histogram_a_single_csv_file(self,spec, SPC_time_MYA, WGD_time_MYA, bin_size, max_Ks_for_x_axis, full_csv_path, plot_title):
        ks_result_for_file = read_Ks_csv(full_csv_path)
        # making subplots
        fig, ax = plt.subplots(1, 1, figsize=(10, 10))
        fig.suptitle(plot_title)
        # ax[0, 0].set_title("Allopolyploid\n", fontsize=20)
        #max_Ks_for_x_axis = 8
        this_ax = ax
        this_ax, ymax_suggestion, metrics = make_subplot(this_ax,spec, ks_result_for_file, bin_size, WGD_time_MYA,
                                                SPC_time_MYA,
                                                max_Ks_for_x_axis, False, "blue", True)
        this_ax.set(xlabel="Ks")
        out_file_name = full_csv_path.replace(".csv", ".hist.png")
        if max_Ks_for_x_axis:
            out_file_name = out_file_name.replace(".png", "_maxKs" + str(max_Ks_for_x_axis) + ".png")
        plt.savefig(out_file_name)
        plt.close()

    def test_download_mesx_results(self):

        batch_folder="sim20_log"
        runs=range(0,7)
        me_at_remote_URL='tdunn@mesx.sdsu.edu'
        local_output_folder="/home/tamsen/Data/Specks_outout_from_mesx"
        remote_output_folder="/usr/scratch2/userdata2/tdunn/SpecKS_Output"
        local_batch_folder=os.path.join(local_output_folder, batch_folder)
        remote_batch_folder=os.path.join(remote_output_folder, batch_folder)
        os.makedirs(local_batch_folder)
        for run in runs:

            allo_run="Allo" + str(run)
            auto_run="Auto" + str(run)
            local_allo_folder = os.path.join(local_batch_folder, allo_run)
            local_auto_folder = os.path.join(local_batch_folder,auto_run)
            os.makedirs(local_allo_folder)
            os.makedirs(local_auto_folder)

            remote_allo_folder=os.path.join(remote_batch_folder,"specks_" + allo_run.upper()+"*")
            remote_to_match = remote_allo_folder + "/A*/*final*/*.csv"
            cmd2 = ['scp', '-r', me_at_remote_URL + ':' + remote_to_match, local_allo_folder ]
            print(" ".join(cmd2))
            out_string, error_string = process_wrapper.run_and_wait_on_process(cmd2, local_batch_folder)

            remote_auto_folder=os.path.join(remote_batch_folder,"specks_" + auto_run.upper()+"*")
            remote_to_match = remote_auto_folder + "/A*/*final*/*.csv"
            cmd2 = ['scp', '-r', me_at_remote_URL + ':' + remote_to_match, local_auto_folder ]
            print(" ".join(cmd2))
            out_string, error_string = process_wrapper.run_and_wait_on_process(cmd2, local_batch_folder)

        cmd2=['scp','-r',me_at_remote_URL+':' +remote_batch_folder +'/specks*/*.xml','.']
        out_string, error_string = process_wrapper.run_and_wait_on_process(cmd2, local_batch_folder)

        cmd2=['scp','-r',me_at_remote_URL+':' +remote_batch_folder +'/specks*/*log*','.']
        out_string, error_string = process_wrapper.run_and_wait_on_process(cmd2, local_batch_folder)
    def test_multi_run_viewer(self):

        plot_title='Simulation with custom GBD model, \nwith Ne-driven allopolyploid ortholog divergence'
        #suppose you have lots of results (cvs files) with all the KS results from many specks runs,
        #and you want to see them all together on one plot.

        output_folder="/home/tamsen/Data/Specks_outout_from_mesx/sim20_log"
        #params_by_polyploid = self.get_truth_for_1MY_sim() #self.get_truth_for_5MY_sim()
        params_by_polyploid = self.get_truth_for_Fig1_sim()
        print("Reading csv files..")
        csvfiles_by_polyploid_by_species_rep_by_algorithm = self.get_ks_data_from_folders(output_folder)
        example_sim='Allo1' #list(csvfiles_by_polyploid_by_species_rep_by_algorithm.keys())[0]
        species=list(csvfiles_by_polyploid_by_species_rep_by_algorithm[example_sim].keys())
        replicates=list(csvfiles_by_polyploid_by_species_rep_by_algorithm[example_sim][species[0]].keys())
        #algs=list(csvfiles_by_polyploid_by_rep_by_algorthim[example_sim][replicates[0]].keys())
        algs=["ML"]

        print("Making plots..")
        metric_by_sim_name=[]


        bin_size = 0.01
        for spec in species:#['outgroup']:#species:
            print(spec)
            for replicate in replicates:
                for alg in algs:

                    max_Ks = 1.0
                    metrics_by_result_names = plot_histograms_for_the_sim_runs(output_folder, plot_title,
                                         csvfiles_by_polyploid_by_species_rep_by_algorithm,spec,
                                         replicate, alg, params_by_polyploid, max_Ks, bin_size, True )

                    max_Ks = False
                    plot_histograms_for_the_sim_runs(output_folder, plot_title,
                                         csvfiles_by_polyploid_by_species_rep_by_algorithm,spec,
                                         replicate, alg, params_by_polyploid, max_Ks, bin_size, True)

                    max_Ks = 0.1
                    plot_histograms_for_the_sim_runs(output_folder, plot_title,
                                         csvfiles_by_polyploid_by_species_rep_by_algorithm,spec,
                                         replicate, alg, params_by_polyploid,max_Ks, 0.001, False )

                    plot_histograms_for_the_sim_runs(output_folder, plot_title,
                                         csvfiles_by_polyploid_by_species_rep_by_algorithm,spec,
                                         replicate, alg, params_by_polyploid,0.5, 0.001 , False)

                    print(metrics_by_result_names)

        self.plot_and_save_metrics(metrics_by_result_names, output_folder, params_by_polyploid)

        self.assertEqual(True,True)

    def plot_and_save_metrics(self, metrics_by_result_names, output_folder, params_by_polyploid):

        allo_results={}
        auto_results={}
        for sim_name, metric in metrics_by_result_names.items():
            if "Allo" in sim_name:
                allo_results[sim_name]=metric
            else:
                auto_results[sim_name] = metric

        self.save_metrics_to_csv(allo_results, auto_results, output_folder, params_by_polyploid)

        allo_xs_spc_times = []  # times
        auto_xs_spc_times = []  # times
        allo_ys_diffs = []  # metrics
        auto_ys_diffs = []  # metrics

        for sim_name, metric in metrics_by_result_names.items():

                spec_time = params_by_polyploid[sim_name].SPC_time_MYA
                mode = metric[0]
                cm = metric[1]
                num_paralogs=metric[2]
                diff = abs(mode - cm)
                if "Allo" in sim_name:
                    allo_xs_spc_times.append(spec_time)
                    allo_ys_diffs.append(diff)
                else:
                    auto_xs_spc_times.append(spec_time)
                    auto_ys_diffs.append(diff)

        plot_allo_vs_auto_metrics(output_folder, allo_xs_spc_times, allo_ys_diffs,
                                  auto_xs_spc_times, auto_ys_diffs, "Allo vs Auto Discrimination plot")

    def save_metrics_to_csv(self, allo_results, auto_results, output_folder, params_by_polyploid):
        out_csv = "metrics.csv"
        out_metrics_path = os.path.join(output_folder, out_csv)
        with open(out_metrics_path, 'w') as f:
            f.writelines("sim_name,spc_time,wgd_time,mode,cm,num_paralogs\n")

            for sim_name, metric in allo_results.items():
                spec_time = params_by_polyploid[sim_name].SPC_time_MYA
                wgd_time = params_by_polyploid[sim_name].WGD_time_MYA
                data = [sim_name, str(spec_time), str(wgd_time)] + [str(m) for m in metric]
                f.writelines(",".join(data) + "\n")

            for sim_name, metric in auto_results.items():
                spec_time = params_by_polyploid[sim_name].SPC_time_MYA
                wgd_time = params_by_polyploid[sim_name].WGD_time_MYA
                data = [sim_name, str(spec_time), str(wgd_time)] + [str(m) for m in metric]
                f.writelines(",".join(data) + "\n")

    def get_truth_for_Fig1_sim(self):
        params_by_polyploid = {}
        params_by_polyploid["Allo1"] = config.PolyploidParams(80, 75, "Allo1")
        params_by_polyploid["Allo2"] = config.PolyploidParams(60, 55, "Allo2")
        params_by_polyploid["Allo3"] = config.PolyploidParams(40, 35, "Allo3")
        params_by_polyploid["Allo4"] = config.PolyploidParams(20, 14, "Allo4")
        params_by_polyploid["Allo5"] = config.PolyploidParams(10, 5, "Allo5")
        params_by_polyploid["Allo6"] = config.PolyploidParams(5, 1, "Allo6")
        params_by_polyploid["Auto1"] = config.PolyploidParams(80, 80, "Auto1")
        params_by_polyploid["Auto2"] = config.PolyploidParams(60, 60, "Auto2")
        params_by_polyploid["Auto3"] = config.PolyploidParams(40, 40, "Auto3")
        params_by_polyploid["Auto4"] = config.PolyploidParams(20, 20, "Auto4")
        params_by_polyploid["Auto5"] = config.PolyploidParams(10, 10, "Auto5")
        params_by_polyploid["Auto6"] = config.PolyploidParams(5, 5, "Auto6")
        return params_by_polyploid
    def get_truth_for_1MY_sim(self):
        params_by_polyploid = {}
        params_by_polyploid["Allo0"] = config.PolyploidParams(20, 15, "Allo0")
        params_by_polyploid["Allo1"] = config.PolyploidParams(40, 30, "Allo1")
        params_by_polyploid["Allo2"] = config.PolyploidParams(60, 50, "Allo2")
        params_by_polyploid["Allo3"] = config.PolyploidParams(80, 70, "Allo3")
        params_by_polyploid["Auto0"] = config.PolyploidParams(20, 20, "Auto0")
        params_by_polyploid["Auto1"] = config.PolyploidParams(40, 40, "Auto1")
        params_by_polyploid["Auto2"] = config.PolyploidParams(60, 60, "Auto2")
        params_by_polyploid["Auto3"] = config.PolyploidParams(80, 80, "Auto3")
        return params_by_polyploid
    def get_truth_for_5MY_sim(self):
        params_by_polyploid = {}
        params_by_polyploid["Allo0"] = config.PolyploidParams(200, 150, "Allo0")
        params_by_polyploid["Allo1"] = config.PolyploidParams(150, 100, "Allo1")
        params_by_polyploid["Allo2"] = config.PolyploidParams(50, 25, "Allo2")
        params_by_polyploid["Allo3"] = config.PolyploidParams(25, 20, "Allo3")
        params_by_polyploid["Auto0"] = config.PolyploidParams(200, 200, "Auto0")
        params_by_polyploid["Auto1"] = config.PolyploidParams(150, 150, "Auto1")
        params_by_polyploid["Auto2"] = config.PolyploidParams(50, 50, "Auto2")
        params_by_polyploid["Auto3"] = config.PolyploidParams(25, 25, "Auto3")
        return params_by_polyploid

    def get_ks_data_from_folders(self, output_folder):
        polyploid_data_folders = os.listdir(output_folder)
        csvfiles_by_polyploid_by_species_rep_by_algorthim = {}
        for polyploid_folder in polyploid_data_folders:
            print(polyploid_folder)
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

                    # example file name: ML_rep0_Ks_by_GeneTree.csv
                    splat = csv_file.split("_")
                    len_splat=len(splat)
                    if splat[0]== "outgroup":
                        species = "outgroup"
                    else:
                        species = "polyploid"
                    algorithm = splat[len_splat-5]
                    replicate = splat[len_splat-4]

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

def plot_histograms_for_the_sim_runs(run_folder, sample_name, csvfiles_by_polyploid_by_rep_by_algorthim,
                                     spec,replicate, alg, params_by_polyploid, max_Ks_for_x_axis,
                                     bin_size, do_curve_fit):

    result_names=(list(csvfiles_by_polyploid_by_rep_by_algorthim.keys()))
    result_names.sort()
    ordered_allo_results=[n for n in result_names if "Allo" in n]
    ordered_auto_results=[n for n in result_names if "Auto" in n]
    metrics_by_result_names= {}

    # making subplots
    num_sims=len(ordered_allo_results)
    fig, ax = plt.subplots(num_sims, 2,figsize=(10, 10))
    fig.suptitle(sample_name + " for " + replicate +", "+ alg + " algorithm")
    ax[0, 0].set_title("Allopolyploid\n",fontsize=20)
    ax[0, 1].set_title("Autopolyploid\n",fontsize=20)


    for sim_idx in range(0, num_sims):

        #if len(ordered_allo_results) >0:
        allo_result_name=ordered_allo_results[sim_idx]
        csvs_for_allo_result= csvfiles_by_polyploid_by_rep_by_algorthim[allo_result_name]

        if spec in csvs_for_allo_result:
            ks_for_allo_result= csvs_for_allo_result[spec][replicate][alg]

            #plot allo result
            params=params_by_polyploid[allo_result_name]
            this_ax = ax[sim_idx, 0]
            this_ax, ymax_suggestion, metrics = make_subplot(this_ax, spec, ks_for_allo_result, bin_size,params.WGD_time_MYA, params.SPC_time_MYA,
                max_Ks_for_x_axis, False,"blue", do_curve_fit)

            if metrics:
                metrics_by_result_names[allo_result_name]=metrics

            this_ax.set(ylabel="simulation #" + str(sim_idx))
            if (sim_idx==3):
                this_ax.set(xlabel="Ks")

        #plot auto result

        #if len(ordered_auto_results) >0:
        auto_result_name=ordered_auto_results[sim_idx]
        csvs_for_auto_result= csvfiles_by_polyploid_by_rep_by_algorthim[auto_result_name]

        if spec in csvs_for_auto_result:
            ks_for_auto_result= csvs_for_auto_result[spec][replicate][alg]
            params=params_by_polyploid[auto_result_name]
            this_ax = ax[sim_idx, 1]
            this_ax, ymax_suggestion, metrics= make_subplot(this_ax,spec,
                ks_for_auto_result, bin_size, params.WGD_time_MYA, params.SPC_time_MYA,
                max_Ks_for_x_axis, ymax_suggestion,"blue", do_curve_fit)

            if metrics:
                metrics_by_result_names[auto_result_name]=metrics

            text_string="SPC: {0} MYA\nWGD: {1} MYA".format(params.SPC_time_MYA,params.WGD_time_MYA)
            #this_ax.text(0.8, 0.8, text_string,
            #             horizontalalignment='center', verticalalignment='center',
            #             transform=this_ax.transAxes)
            if (sim_idx==(num_sims-1)):
                this_ax.set(xlabel="Ks")

    out_file_name=os.path.join(run_folder,"histogram" + "_plot_" + spec +"_" + replicate +"_"+ alg+".png")
    if max_Ks_for_x_axis:
        out_file_name = out_file_name.replace(".png","_maxKs" +str(max_Ks_for_x_axis)+ ".png")

    plt.savefig(out_file_name)

    plt.close()
    return metrics_by_result_names


def make_subplot(this_ax, spec, Ks_results, bin_size,WGD_time_MYA, SPC_time_MYA,
                 max_Ks, maxY, plot_color, do_curve_fit):

    WGD_as_Ks = WGD_time_MYA * 0.01 #/ 1.04 ..the peak max is about 96% off from where it should be
    SPEC_as_Ks =SPC_time_MYA * 0.01 #/ 1.04
    default_xaxis_limit =SPEC_as_Ks + 0.2
    x = Ks_results
    num_gene_pairs_str=str(len(Ks_results))
    if max_Ks:
        bins = np.arange(0, max_Ks + 0.1, bin_size)
        n, bins, patches = this_ax.hist(x, bins=bins, facecolor=plot_color, alpha=0.25,
            label='ks hist ({0} pairs)'.format(num_gene_pairs_str))
    else:
        bins = np.arange(0, default_xaxis_limit, bin_size)
        n, bins, patches = this_ax.hist(x, bins=bins, facecolor=plot_color, alpha=0.25,
            label='ks hist ({0} pairs)'.format(num_gene_pairs_str))

    hist_maximum=max(n)
    ymax_suggestion=hist_maximum*1.6

    #do_curve_fitting=False
    #do_curve_fitting = (spec != "outgroup" ) and (max_Ks and (max_Ks > 0.1))
    curve_fit_done=False

    metrics=False
    if do_curve_fit and (spec != "outgroup" ):

        this_ax.axvline(x=WGD_as_Ks, color='b', linestyle='-', label="WGD time, "+ str(WGD_time_MYA)+ " MYA")
        this_ax.axvline(x=SPEC_as_Ks, color='r', linestyle='--', label="SPEC time, "+ str(SPC_time_MYA)+ " MYA")

        fit_curve_ys1, xs_for_wgd, mode,cm = curve_fitting.fit_curve_to_hist(bins, n)
        print("hist_maximum " + str(hist_maximum))
        print("mode " + str(mode))
        print("cm " + str(cm))
        curve_fit_done=True
        metrics=[mode,cm,len(Ks_results)]
        if fit_curve_ys1 and (hist_maximum>0):
            this_ax.plot(xs_for_wgd,fit_curve_ys1,
                 color='green', linestyle=':', label="lognorm fit")

        if SPC_time_MYA == 20:
            ymax_suggestion=400


    if maxY:
        yaxis_limit=maxY
        ymax_suggestion=False
    else:
        yaxis_limit= ymax_suggestion

    if curve_fit_done:
        if fit_curve_ys1 and (hist_maximum>0):
            this_ax.scatter(cm,0.05*yaxis_limit,
                 color='darkgreen', marker='o', s=100)# label="cm",)

            this_ax.scatter(mode,0.05*yaxis_limit,
                 color='cyan', marker='^', s=80)# label="mode")



    this_ax.legend()

    this_ax.set(ylim=[0, yaxis_limit])

    if max_Ks:
        this_ax.set(xlim=[0, max_Ks * 1.1])
        if max_Ks < 1.0:
            this_ax.set(ylim=[0, 5])
    else:
        this_ax.set(xlim=[0, default_xaxis_limit])

    return this_ax,ymax_suggestion, metrics

def plot_allo_vs_auto_metrics(out_folder, allo_xs, allo_ys, auto_xs, auto_ys, title):


    fig, ax = plt.subplots(1, 1, figsize=(10, 10))
    plt.plot(allo_xs, allo_ys, label='Allo',c='r')
    plt.plot(auto_xs, auto_ys, label='Auto',c='b')
    out_file_name = os.path.join(out_folder, title+ ".png")
    fig.suptitle(title)
    ax.set(xlabel="MYA")
    ax.set(ylabel="Metric (mode-cm)")

    ax.legend()
    plt.savefig(out_file_name)
    plt.close()
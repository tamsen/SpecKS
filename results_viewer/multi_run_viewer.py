import os
import unittest
from scipy.optimize import curve_fit
import numpy as np
from matplotlib import pyplot as plt

import config
from results_viewer import curve_fitting


#https://stackoverflow.com/questions/14770735/how-do-i-change-the-figure-size-with-subplots
class MulitRunViewerTests(unittest.TestCase):
    def test_multi_run_viewer(self):

        plot_title='Simulation with NO gene birth/death'
        #suppose you have lots of results (cvs files) with all the KS results from many specks runs,
        #and you want to see them all together on one plot.

        #output_folder="/home/tamsen/Data/SpecKS_mesx_data/mesx_sim1_no_genebirth_or_death"
        output_folder="/home/tamsen/Data/Specks_outout_from_mesx/sim1_redo"
        #output_folder="/home/tamsen/Data/Specks_outout_from_mesx/mesx_sim2_genebirth_and_death"
        #output_folder = "/home/tamsen/Data/SpecKS_mesx_data/mesx_sim2_genebirth_and_death"

        print("Reading csv files..")
        csvfiles_by_polyploid_by_rep_by_algorthim = self.get_ks_data_from_folders(output_folder)
        example_sim=list(csvfiles_by_polyploid_by_rep_by_algorthim.keys())[0]
        replicates=list(csvfiles_by_polyploid_by_rep_by_algorthim[example_sim].keys())
        #algs=list(csvfiles_by_polyploid_by_rep_by_algorthim[example_sim][replicates[0]].keys())
        algs=["ML"]

        print("Making plots..")
        params_by_polyploid={}
        params_by_polyploid["Allo0"]= config.PolyploidParams(200,150,"Allo0")
        params_by_polyploid["Allo1"]= config.PolyploidParams(150,100,"Allo1")
        params_by_polyploid["Allo2"]= config.PolyploidParams(50, 25,"Allo2")
        params_by_polyploid["Allo3"]= config.PolyploidParams( 25, 20,"Allo3")

        params_by_polyploid["Auto0"]= config.PolyploidParams(200,200,"Auto0")
        params_by_polyploid["Auto1"]= config.PolyploidParams(150,150,"Auto1")
        params_by_polyploid["Auto2"]= config.PolyploidParams( 50, 50,"Auto2")
        params_by_polyploid["Auto3"]= config.PolyploidParams( 25, 25,"Auto3")

        max_Ks = 8
        for replicate in replicates:
            for alg in algs:
                plot_histograms_for_the_sim_runs(output_folder, plot_title,
                                         csvfiles_by_polyploid_by_rep_by_algorthim,
                                         replicate, alg, params_by_polyploid, max_Ks )

                plot_histograms_for_the_sim_runs(output_folder, plot_title,
                                         csvfiles_by_polyploid_by_rep_by_algorthim,
                                         replicate, alg, params_by_polyploid, False )

        self.assertEqual(True,True)

    def get_ks_data_from_folders(self, output_folder):
        polyploid_data_folders = os.listdir(output_folder)
        csvfiles_by_polyploid_by_rep_by_algorthim = {}
        for polyploid_folder in polyploid_data_folders:
            print(polyploid_folder)
            full_polyploid_folder_path = os.path.join(output_folder, polyploid_folder)
            if (os.path.isdir(full_polyploid_folder_path)):
                files = os.listdir(full_polyploid_folder_path)
                csvfiles_by_replicate = {}
                for csv_file in files:

                    if not ".csv" in csv_file:
                        continue

                    # example file name: ML_rep0_Ks_by_GeneTree.csv
                    splat = csv_file.split("_")
                    algorithm = splat[0]
                    replicate = splat[1]
                    if not replicate in csvfiles_by_replicate:
                        csvfiles_by_alg = {}
                        csvfiles_by_replicate[replicate] = csvfiles_by_alg

                    # csvfiles_by_replicate[replicate][algorithm]=csv_file
                    full_csv_path = os.path.join(full_polyploid_folder_path, csv_file)
                    print("reading " + full_csv_path)
                    ks_result_for_file = read_Ks_csv(full_csv_path)
                    csvfiles_by_replicate[replicate][algorithm] = ks_result_for_file

                csvfiles_by_polyploid_by_rep_by_algorthim[polyploid_folder] = csvfiles_by_replicate
        return csvfiles_by_polyploid_by_rep_by_algorthim


def read_Ks_csv(csv_file):

    ks_results = []
    with open(csv_file, "r") as f:

        reading_header=True
        while True:
            line = f.readline()
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
                                     replicate, alg, params_by_polyploid, max_Ks_for_x_axis):

    result_names=(list(csvfiles_by_polyploid_by_rep_by_algorthim.keys()))
    result_names.sort()
    ordered_allo_results=[n for n in result_names if "Allo" in n]
    ordered_auto_results=[n for n in result_names if "Auto" in n]

    bin_size = 0.01
    #f, a = plt.subplots(4, 2)

    # making subplots
    fig, ax = plt.subplots(4, 2,figsize=(10, 10))
    fig.suptitle(sample_name + " for " + replicate +", "+ alg + " algorithm")
    ax[0, 0].set_title("Allopolyploid\n",fontsize=20)
    ax[0, 1].set_title("Autopolyploid\n",fontsize=20)

    num_sims=len(ordered_allo_results)
    for sim_idx in range(0, 4):

        allo_result_name=ordered_allo_results[sim_idx]
        auto_result_name=ordered_auto_results[sim_idx]
        ks_for_allo_result= csvfiles_by_polyploid_by_rep_by_algorthim[allo_result_name][replicate][alg]
        ks_for_auto_result= csvfiles_by_polyploid_by_rep_by_algorthim[auto_result_name][replicate][alg]

        #plot allo result
        params=params_by_polyploid[allo_result_name]
        this_ax = ax[sim_idx, 0]
        this_ax, ymax_suggestion = make_subplot(this_ax, ks_for_allo_result, bin_size,params.WGD_time_MYA, params.SPC_time_MYA,
            max_Ks_for_x_axis, False,"blue")

        this_ax.set(ylabel="simulation #" + str(sim_idx))
        if (sim_idx==3):
            this_ax.set(xlabel="Ks")

        text_string="SPC: {0} MYA\nWGD: {1} MYA".format(params.SPC_time_MYA,params.WGD_time_MYA)
        #this_ax.text(0.8, 0.8, text_string,
        #             horizontalalignment='center', verticalalignment='center',
        #             transform=this_ax.transAxes)

        #plot auto result
        params=params_by_polyploid[auto_result_name]
        this_ax = ax[sim_idx, 1]
        this_ax, ymax_suggestion = make_subplot(this_ax, ks_for_auto_result, bin_size, params.WGD_time_MYA, params.SPC_time_MYA,
            max_Ks_for_x_axis, ymax_suggestion,"blue")

        text_string="SPC: {0} MYA\nWGD: {1} MYA".format(params.SPC_time_MYA,params.WGD_time_MYA)
        #this_ax.text(0.8, 0.8, text_string,
        #             horizontalalignment='center', verticalalignment='center',
        #             transform=this_ax.transAxes)
        if (sim_idx==3):
            this_ax.set(xlabel="Ks")

    out_file_name=os.path.join(run_folder,"histogram" + "_plot_" + replicate +"_"+ alg+".png")
    if max_Ks_for_x_axis:
        out_file_name = out_file_name.replace(".png","_maxKs" +str(max_Ks_for_x_axis)+ ".png")

    plt.savefig(out_file_name)

    plt.close()


def make_subplot(this_ax, Ks_results, bin_size,WGD_time_MYA, SPC_time_MYA, max_Ks, maxY, plot_color):

    WGD_as_Ks = WGD_time_MYA * 0.01
    SPEC_as_Ks =SPC_time_MYA * 0.01

    x = Ks_results
    if max_Ks:
        bins = np.arange(0, max_Ks + 0.1, bin_size)
        n, bins, patches = this_ax.hist(x, bins=bins, facecolor=plot_color, alpha=0.25, label='ks hist. data')
    else:
        bins = np.arange(0, SPEC_as_Ks + 1, bin_size)
        n, bins, patches = this_ax.hist(x, bins=bins, facecolor=plot_color, alpha=0.25, label='ks hist. data')

    hist_maximum=max(n)
    ymax_suggestion=hist_maximum*1.5
    fit_curve_ys1, xs_for_wgd, mode,cm = curve_fitting.fit_curve_to_hist(bins, n)

    this_ax.axvline(x=WGD_as_Ks, color='b', linestyle='-', label="WGD time, "+ str(WGD_time_MYA)+ " MYA")
    this_ax.axvline(x=SPEC_as_Ks, color='r', linestyle='--', label="SPEC time, "+ str(SPC_time_MYA)+ " MYA")

    this_ax.plot(xs_for_wgd,fit_curve_ys1,
                 color='green', linestyle=':', label="lognorm fit")

    if maxY:
        yaxis_limit=maxY
        ymax_suggestion=False
    else:
        yaxis_limit= ymax_suggestion

    this_ax.scatter(cm,0.05*yaxis_limit,
                 color='darkgreen', marker='o', label="cm", s=100)

    this_ax.scatter(mode,0.05*yaxis_limit,
                 color='cyan', marker='^', label="mode",s=80)



    this_ax.legend()

    this_ax.set(ylim=[0, yaxis_limit])
    if max_Ks:
        this_ax.set(xlim=[0, max_Ks * 1.1])
    else:
        this_ax.set(xlim=[0, SPEC_as_Ks + 1])

    return this_ax,ymax_suggestion

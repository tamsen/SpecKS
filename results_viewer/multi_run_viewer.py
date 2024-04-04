import os
import unittest

from scipy import stats
from scipy.optimize import curve_fit
import numpy as np
from matplotlib import pyplot as plt

import config
import process_wrapper
from results_viewer import curve_fitting
from results_viewer.parse_aggregate_results import metric_result, get_metric_result_data_headers


#https://stackoverflow.com/questions/14770735/how-do-i-change-the-figure-size-with-subplots
class MulitRunViewerTests(unittest.TestCase):

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

                    #max_Ks = False
                    #plot_histograms_for_the_sim_runs(output_folder, plot_title,
                    #                     csvfiles_by_polyploid_by_species_rep_by_algorithm,spec,
                    #                     replicate, alg, params_by_polyploid, max_Ks, False, True)

                    #max_Ks = 0.1
                    #plot_histograms_for_the_sim_runs(output_folder, plot_title,
                    #                     csvfiles_by_polyploid_by_species_rep_by_algorithm,spec,
                    #                     replicate, alg, params_by_polyploid,max_Ks, 0.001, False )

                    #plot_histograms_for_the_sim_runs(output_folder, plot_title,
                    #                     csvfiles_by_polyploid_by_species_rep_by_algorithm,spec,
                    #                     replicate, alg, params_by_polyploid,0.5, 0.001 , False)

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

        #allo_xs_spc_times = []  # times
        #auto_xs_spc_times = []  # times
        #allo_ys_diffs = []  # metrics
        #auto_ys_diffs = []  # metrics

        #for sim_name, metric in metrics_by_result_names.items():

        #        spec_time = params_by_polyploid[sim_name].SPC_time_MYA
        #        mode = metric.mode
        #        cm = metric.cm
        #        num_paralogs=metric.num_paralogs
        #        diff = abs(mode - cm)
        #        if "Allo" in sim_name:
        #            allo_xs_spc_times.append(spec_time)
        #            allo_ys_diffs.append(diff)
        #        else:
        #            auto_xs_spc_times.append(spec_time)
        #            auto_ys_diffs.append(diff)

        #plot_allo_vs_auto_metrics(output_folder, allo_xs_spc_times, allo_ys_diffs,
        #                          auto_xs_spc_times, auto_ys_diffs, "Allo vs Auto Discrimination plot")

    def save_metrics_to_csv(self, allo_results, auto_results, output_folder, params_by_polyploid):
        out_csv = "metrics.csv"
        out_metrics_path = os.path.join(output_folder, out_csv)
        with open(out_metrics_path, 'w') as f:
            metric_result_data_headers= get_metric_result_data_headers()
            f.writelines(",".join(metric_result_data_headers) +"\n")

            for sim_name, metric in allo_results.items():
                f.writelines(metric.to_csv_string() + "\n")

            for sim_name, metric in auto_results.items():
                f.writelines(metric.to_csv_string() + "\n")

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

                    if "Kn" in csv_file:
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
        print("spec:\t" + spec)
        print("replicate:\t" + replicate)
        print("alg:\t" + alg)

        if spec in csvs_for_allo_result:
            ks_for_allo_result= csvs_for_allo_result[spec][replicate][alg]

            #plot allo result
            params=params_by_polyploid[allo_result_name]
            this_ax = ax[sim_idx, 0]
            this_ax, ymax_suggestion, lognorm_fit_metrics,gaussian_fit_metrics = make_subplot(this_ax, spec, ks_for_allo_result, bin_size,params.WGD_time_MYA, params.SPC_time_MYA,
                max_Ks_for_x_axis, False,"blue", do_curve_fit)

            if lognorm_fit_metrics:
                fit_data = fit_data_to_keep("allo",allo_result_name,
                                            params.SPC_time_MYA,params.WGD_time_MYA,
                                            lognorm_fit_metrics,gaussian_fit_metrics)
                metrics_by_result_names[allo_result_name]=fit_data

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
            this_ax, ymax_suggestion, lognorm_fit_metrics, gaussian_fit_metrics= make_subplot(this_ax,spec,
                ks_for_auto_result, bin_size, params.WGD_time_MYA, params.SPC_time_MYA,
                max_Ks_for_x_axis, ymax_suggestion,"blue", do_curve_fit)

            if lognorm_fit_metrics:
                fit_data = fit_data_to_keep("auto",auto_result_name,
                                            params.SPC_time_MYA,params.WGD_time_MYA,
                                            lognorm_fit_metrics, gaussian_fit_metrics)
                metrics_by_result_names[auto_result_name]=fit_data

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
        #bins = np.arange(0, default_xaxis_limit, bin_size)
        nBins=100
        n, bins, patches = this_ax.hist(x, bins=nBins, facecolor=plot_color, alpha=0.25,
            label='ks hist ({0} pairs)'.format(num_gene_pairs_str))

    hist_maximum=max(n)
    ymax_suggestion=hist_maximum*1.6

    #do_curve_fitting=False
    #do_curve_fitting = (spec != "outgroup" ) and (max_Ks and (max_Ks > 0.1))
    curve_fit_done=False

    lognorm_goodness_of_fit =False
    gaussian_goodness_of_fit =False
    if do_curve_fit and (spec != "outgroup" ):

        this_ax.axvline(x=WGD_as_Ks, color='b', linestyle='-', label="WGD time, "+ str(WGD_time_MYA)+ " MYA")
        this_ax.axvline(x=SPEC_as_Ks, color='r', linestyle='--', label="SPEC time, "+ str(SPC_time_MYA)+ " MYA")

        gaussian_fit_curve_ys1, xs_for_wgd, gaussian_goodness_of_fit = \
            curve_fitting.fit_curve_to_hist(bins, n,curve_fitting.wgd_gaussian)

        lognorm_fit_curve_ys1, xs_for_wgd,lognorm_goodness_of_fit =  \
            curve_fitting.fit_curve_to_hist(bins, n, curve_fitting.wgd_lognorm )


        if lognorm_fit_curve_ys1:
            curve_fit_done=True

            #https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.kstest.html
            #ks_norm_result= stats.kstest(bins,stats.norm.cdf)
            #ks_lognorm_result = stats.kstest(bins, stats.lognorm.cdf(popt))
            #https://stackoverflow.com/questions/51902996/scipy-kstest-used-on-scipy-lognormal-distrubtion
            #my_args = popt[0:3]
            #ks_lognorm_result = stats.kstest(bins, 'lognorm', args=my_args)

            #print("ks_norm_result " + str(ks_norm_result))
            #print("ks_lognorm_result " + str(ks_lognorm_result))


        if lognorm_fit_curve_ys1 and (hist_maximum>0):
            rmse_str= str(round(lognorm_goodness_of_fit.RMSE,4))
            this_ax.plot(xs_for_wgd,lognorm_fit_curve_ys1,
                 color='green', linestyle=':', label="lognorm fit (err={0})".format(rmse_str))

        if gaussian_fit_curve_ys1 and (hist_maximum>0):
            rmse_str= str(round(gaussian_goodness_of_fit.RMSE,4))
            this_ax.plot(xs_for_wgd,gaussian_fit_curve_ys1,
                 color='blue', linestyle=':', label="gaussian fit (err={0})".format(rmse_str))

        if SPC_time_MYA == 20:
            ymax_suggestion=400


    if maxY:
        yaxis_limit=maxY
        ymax_suggestion=False
    else:
        yaxis_limit= ymax_suggestion

    if curve_fit_done:
        if lognorm_fit_curve_ys1 and (hist_maximum>0):
            this_ax.scatter(lognorm_goodness_of_fit.cm,0.05*yaxis_limit,
                 color='darkgreen', marker='o', s=100)# label="cm",)

            this_ax.scatter(lognorm_goodness_of_fit.mode,0.05*yaxis_limit,
                 color='cyan', marker='^', s=80)# label="mode")



    this_ax.legend()

    this_ax.set(ylim=[0, yaxis_limit])

    if max_Ks:
        this_ax.set(xlim=[0, max_Ks * 1.1])
        if max_Ks < 0.51:
            this_ax.set(ylim=[0, 5])
    else:
        this_ax.set(xlim=[0, default_xaxis_limit])

    return this_ax,ymax_suggestion, lognorm_goodness_of_fit,gaussian_goodness_of_fit

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


class fit_data_to_keep():
    input_type=''
    sim_name=''
    spc_time=-1
    wgd_time=-1
    lognorm_fit_data=-1
    gaussian_fit_data=-1

    def __init__(self, input_type, sim_name, spc_time, wgd_time, lognorm_fit_data, gaussian_fit_data):
        self.input_type = input_type
        self.sim_name =sim_name
        self.spc_time =spc_time
        self.wgd_time = wgd_time
        self.lognorm_fit_data = lognorm_fit_data
        self.gaussian_fit_data = gaussian_fit_data
    def to_csv_string(self):

        final_data = self.to_data_list()
        return ",".join(final_data)

    def to_data_list(self):
        simple_data = [self.input_type, self.sim_name, self.spc_time, self.wgd_time]#,self.fit_data]
        fit_data_list=self.lognorm_fit_data.to_data_list() + [str(self.gaussian_fit_data.RMSE)]
        final_data = [str(d) for d in simple_data] + [str(d) for d in fit_data_list]
        return final_data

def get_metric_result_data_headers():
         return ["sim_type","sim_name","spc_time","wgd_time","mode","cm","num_paralogs",
                 "popt","lognorm_RMSE","gaussian_RMSE"]
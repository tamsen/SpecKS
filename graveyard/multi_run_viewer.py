import os
import unittest

from scipy import stats
from scipy.optimize import curve_fit
import numpy as np
from matplotlib import pyplot as plt

import config
import process_wrapper
from results_viewer import curve_fitting, run_metrics


#https://stackoverflow.com/questions/14770735/how-do-i-change-the-figure-size-with-subplots
class MulitRunViewerTests(unittest.TestCase):

    def test_multi_run_viewer(self):

        plot_title='Simulation with custom GBD model, \nwith Ne-driven allopolyploid ortholog divergence'
        #suppose you have lots of results (cvs files) with all the KS results from many specks runs,
        #and you want to see them all together on one plot.

        output_folder="/home/tamsen/Data/Specks_outout_from_mesx/sim34_log"
        batch_run_name=os.path.basename(output_folder)
        #params_by_polyploid = self.get_truth_for_1MY_sim() #self.get_truth_for_5MY_sim()
        params_by_polyploid = get_truth_for_sim31()
        print("Reading csv files..")
        csvfiles_by_polyploid_by_species_rep_by_algorithm = get_ks_data_from_folders(output_folder)
        example_sim=list(csvfiles_by_polyploid_by_species_rep_by_algorithm.keys())[0]
        species=list(csvfiles_by_polyploid_by_species_rep_by_algorithm[example_sim].keys())
        replicates=list(csvfiles_by_polyploid_by_species_rep_by_algorithm[example_sim][species[0]].keys())
        #algs=list(csvfiles_by_polyploid_by_rep_by_algorthim[example_sim][replicates[0]].keys())
        algs=["ML"]

        print("Making plots..")
        metric_by_sim_name=[]
        do_kde=False
        do_curve_fit=False
        bin_size = 0.001

        for spec in species:#['outgroup']:#species:
            print(spec)
            for replicate in replicates:
                for alg in algs:


                    max_Ks = 1.0
                    metrics_by_result_names = plot_ks_disributions_for_the_sim_runs(output_folder, plot_title,
                                   csvfiles_by_polyploid_by_species_rep_by_algorithm, spec,
                                   replicate, alg, params_by_polyploid, max_Ks, bin_size, do_curve_fit, do_kde)

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

                    out_csv = "{0}_{1}_{2}_{3}_metrics.csv".format(batch_run_name,spec,replicate,alg)
                    out_file_name = os.path.join(output_folder, out_csv)
                    run_metrics.plot_and_save_metrics(metrics_by_result_names, out_file_name)

        self.assertEqual(True,True)

def get_truth_for_Fig1_sim():
    params_by_polyploid = {}
    params_by_polyploid["Allo1"] = config.PolyploidParams(80, 75, "Allo1")
    params_by_polyploid["Allo2"] = config.PolyploidParams(60, 55, "Allo2")
    params_by_polyploid["Allo3"] = config.PolyploidParams(40, 35, "Allo3")
    params_by_polyploid["Allo4"] = config.PolyploidParams(20, 15, "Allo4")
    params_by_polyploid["Allo5"] = config.PolyploidParams(10, 5, "Allo5")
    params_by_polyploid["Allo6"] = config.PolyploidParams(5, 1, "Allo6")
    params_by_polyploid["Auto1"] = config.PolyploidParams(80, 80, "Auto1")
    params_by_polyploid["Auto2"] = config.PolyploidParams(60, 60, "Auto2")
    params_by_polyploid["Auto3"] = config.PolyploidParams(40, 40, "Auto3")
    params_by_polyploid["Auto4"] = config.PolyploidParams(20, 20, "Auto4")
    params_by_polyploid["Auto5"] = config.PolyploidParams(10, 10, "Auto5")
    params_by_polyploid["Auto6"] = config.PolyploidParams(5, 5, "Auto6")
    return params_by_polyploid

def get_truth_for_sim31():
    params_by_polyploid = {}
    params_by_polyploid["Allo1"] = config.PolyploidParams(80, 75, "Allo1")
    params_by_polyploid["Allo2"] = config.PolyploidParams(70, 65, "Allo2")
    params_by_polyploid["Allo3"] = config.PolyploidParams(60, 50, "Allo3")
    params_by_polyploid["Allo4"] = config.PolyploidParams(50, 45, "Allo4")
    params_by_polyploid["Auto1"] = config.PolyploidParams(80, 80, "Auto1")
    params_by_polyploid["Auto2"] = config.PolyploidParams(70, 70, "Auto2")
    params_by_polyploid["Auto3"] = config.PolyploidParams(60, 60, "Auto3")
    params_by_polyploid["Auto4"] = config.PolyploidParams(50, 50, "Auto4")
    return params_by_polyploid
def get_truth_for_1MY_sim():
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
def get_truth_for_5MY_sim():
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

def get_ks_data_from_folders(output_folder):
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

def plot_ks_disributions_for_the_sim_runs(run_folder, sample_name, csvfiles_by_polyploid_by_rep_by_algorthim,
                                          spec, replicate, alg, params_by_polyploid, max_Ks_for_x_axis,
                                          bin_size, do_curve_fit, do_kde):

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

            if do_kde:
                this_ax, ymax_suggestion, lognorm_fit_metrics, gaussian_fit_metrics = make_kde_subplot(this_ax,
                                                                                                             spec,
                                                                                                             ks_for_allo_result,
                                                                                                             bin_size,
                                                                                                             params.WGD_time_MYA,
                                                                                                             params.SPC_time_MYA,
                                                                                                             max_Ks_for_x_axis,
                                                                                                             False,
                                                                                                             "blue",
                                                                                                             do_curve_fit)
            else:
                this_ax, ymax_suggestion, lognorm_fit_metrics,gaussian_fit_metrics = make_histogram_subplot(this_ax, spec, ks_for_allo_result, bin_size, params.WGD_time_MYA, params.SPC_time_MYA,
                                                                                                        max_Ks_for_x_axis, False,"blue", do_curve_fit)

            if lognorm_fit_metrics:
                fit_data = run_metrics.run_metrics("allo", allo_result_name,
                                                   params.SPC_time_MYA, params.WGD_time_MYA,
                                                   lognorm_fit_metrics, gaussian_fit_metrics)
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

            if do_kde:
                this_ax, ymax_suggestion, lognorm_fit_metrics, gaussian_fit_metrics= make_kde_subplot(this_ax, spec,
                                            ks_for_auto_result, bin_size, params.WGD_time_MYA, params.SPC_time_MYA,
                                            max_Ks_for_x_axis, ymax_suggestion,"blue", do_curve_fit)
            else:
                this_ax, ymax_suggestion, lognorm_fit_metrics, gaussian_fit_metrics = make_histogram_subplot(this_ax,spec,
                                            ks_for_auto_result,
                                                                                                             bin_size,
                                                                                                             params.WGD_time_MYA,
                                                                                                             params.SPC_time_MYA,
                                                                                                             max_Ks_for_x_axis,
                                                                                                             ymax_suggestion,
                                                                                                             "blue",
                                                                                                             do_curve_fit)

            if lognorm_fit_metrics:
                fit_data = run_metrics.run_metrics("auto", auto_result_name,
                                                   params.SPC_time_MYA, params.WGD_time_MYA,
                                                   lognorm_fit_metrics, gaussian_fit_metrics)
                metrics_by_result_names[auto_result_name]=fit_data

            text_string="SPC: {0} MYA\nWGD: {1} MYA".format(params.SPC_time_MYA,params.WGD_time_MYA)
            #this_ax.text(0.8, 0.8, text_string,
            #             horizontalalignment='center', verticalalignment='center',
            #             transform=this_ax.transAxes)
            if (sim_idx==(num_sims-1)):
                this_ax.set(xlabel="Ks")

    if do_kde:
        out_file_name = os.path.join(run_folder, "kde" + "_plot_" + spec + "_" + replicate + "_" + alg + ".png")
    else:
        out_file_name=os.path.join(run_folder,"histogram" + "_plot_" + spec +"_" + replicate +"_"+ alg+".png")
    if max_Ks_for_x_axis:
        out_file_name = out_file_name.replace(".png","_maxKs" +str(max_Ks_for_x_axis)+ ".png")

    plt.savefig(out_file_name)

    plt.close()
    return metrics_by_result_names


def make_histogram_subplot(this_ax, spec, Ks_results, bin_size, WGD_time_MYA, SPC_time_MYA,
                           max_Ks, maxY, plot_color, do_curve_fit):

    WGD_as_Ks = WGD_time_MYA * 0.01 #/ 1.04 ..the peak max is about 96% off from where it should be
    SPEC_as_Ks =SPC_time_MYA * 0.01 #/ 1.04
    default_xaxis_limit =SPEC_as_Ks + 0.2
    x = Ks_results
    num_gene_pairs_str=str(len(Ks_results))
    if max_Ks:

        bins = np.arange(0, max_Ks + 0.1, bin_size)
        #bins = [i*0.01 for i in range(0,110)]
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
    this_ax.axvline(x=WGD_as_Ks, color='b', linestyle='-', label="WGD time, "+ str(WGD_time_MYA)+ " MYA")
    this_ax.axvline(x=SPEC_as_Ks, color='r', linestyle='--', label="SPEC time, "+ str(SPC_time_MYA)+ " MYA")

    if do_curve_fit and (spec != "outgroup" ):

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

    #this_ax.set(ylim=[0, yaxis_limit])
    this_ax.set(ylim=[0, 100])

    if max_Ks:
        this_ax.set(xlim=[0, max_Ks * 1.1])
        if max_Ks < 0.51:
            this_ax.set(ylim=[0, 5])
    else:
        this_ax.set(xlim=[0, default_xaxis_limit])

    return this_ax,ymax_suggestion, lognorm_goodness_of_fit,gaussian_goodness_of_fit


def make_kde_subplot(this_ax, spec, Ks_results, bin_size, WGD_time_MYA, SPC_time_MYA,
                           max_Ks, maxY, plot_color, do_curve_fit):

    WGD_as_Ks = WGD_time_MYA * 0.01 #/ 1.04 ..the peak max is about 96% off from where it should be
    SPEC_as_Ks =SPC_time_MYA * 0.01 #/ 1.04
    default_xaxis_limit =SPEC_as_Ks + 0.2
    x = Ks_results
    num_gene_pairs_str=str(len(Ks_results))
    xs_for_wgd = np.arange(0, max_Ks + 0.1, bin_size)
    kde = stats.gaussian_kde(x)
    kde_fit_curve_ys1 = kde(xs_for_wgd)
    kde_maximum=max(kde_fit_curve_ys1)
    ymax_suggestion=kde_maximum*1.6

    curve_fit_done=False
    lognorm_goodness_of_fit =False
    gaussian_goodness_of_fit =False
    if do_curve_fit and (spec != "outgroup" ):

        this_ax.axvline(x=WGD_as_Ks, color='b', linestyle='-', label="WGD time, "+ str(WGD_time_MYA)+ " MYA")
        this_ax.axvline(x=SPEC_as_Ks, color='r', linestyle='--', label="SPEC time, "+ str(SPC_time_MYA)+ " MYA")

        gaussian_fit_curve_ys1, xs_for_wgd, gaussian_goodness_of_fit = \
            curve_fitting.fit_curve_to_xs_and_ys(xs_for_wgd, kde_fit_curve_ys1,curve_fitting.wgd_gaussian)

        lognorm_fit_curve_ys1, xs_for_wgd,lognorm_goodness_of_fit =  \
            curve_fitting.fit_curve_to_xs_and_ys(xs_for_wgd, kde_fit_curve_ys1, curve_fitting.wgd_lognorm )


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


        if lognorm_fit_curve_ys1 and (kde_maximum>0):
            rmse_str= str(round(lognorm_goodness_of_fit.RMSE,4))
            this_ax.plot(xs_for_wgd,lognorm_fit_curve_ys1,
                 color='green', linestyle=':', label="lognorm fit (err={0})".format(rmse_str))

        if gaussian_fit_curve_ys1 and (kde_maximum>0):
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
        if lognorm_fit_curve_ys1 and (kde_maximum>0):
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


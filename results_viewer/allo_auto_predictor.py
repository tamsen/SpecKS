import math
import os
import random
import statistics
import unittest

import numpy as np
from scipy.optimize import curve_fit
from sklearn import linear_model
import sklearn.metrics as metrics
from matplotlib import pyplot as plt
from matplotlib.patches import Rectangle
from results_viewer import curve_fitting
from results_viewer.batch_histogrammer import get_truth_from_name_list


class AlloAutoPredictor(unittest.TestCase):
    def test_spec_time_predictor(self):

        out_folder = "/home/tamsen/Data/Specks_outout_from_mesx"
        #file1="mode_vs_spec_time.csv"
        file1="mode_vs_spec_time _without_N100.csv"
        full_path=os.path.join(out_folder,file1)
        spec_xs,mode_ys, sims = read_xs_ys_csv(full_path)

        fig, ax = plt.subplots(1, 1, figsize=(5, 5))
        plt.scatter(mode_ys, spec_xs, c='b',label="empirical data")

        fit_spec_times, xs, goodness_of_fit=curve_fitting.fit_curve_to_xs_and_ys(
            mode_ys, spec_xs, curve_fitting.linear)

        linear_fit_popt=goodness_of_fit.popt
        ax.set(xlabel="mode (Ks space)")
        ax.set(ylabel="spec time as a fxn of mode (MYA)")
        plt.plot(xs,fit_spec_times, c='k',label="line fit (model))\n"
                                                              +str(linear_fit_popt))
        example_modes=np.arange(0.2, 1.0, 0.2)
        predicted_spec_times= [curve_fitting.linear(x, *linear_fit_popt) for x in example_modes]
        plt.bar(example_modes,predicted_spec_times, width=0.002, color='pink',label="model input (mode)")
        plt.scatter(example_modes,predicted_spec_times, c='r',label="model output (est. spec time)")
        print("linear_fit_popt:\t" + linear_fit_popt)

        plot_file=full_path.replace(".csv",".new.png")
        plt.legend()
        plt.savefig(plot_file)
        plt.close()

    def test_wgd_time_predictor(self):

        out_folder = "/home/tamsen/Data/Specks_outout_from_mesx"
        file1="genes_remaining_vs_wgd_time_noN100.csv"
        full_path=os.path.join(out_folder,file1)

        wgd_xs,genes_remaining,sims = read_xs_ys_csv(full_path)
        popt, pcov = curve_fit(curve_fitting.logfit, genes_remaining, wgd_xs)

        fig, ax = plt.subplots(1, 1, figsize=(5, 5))
        plt.scatter(genes_remaining, wgd_xs, c='b',label="empirical data")
        genes_remaining.sort()
        fit_wgd_times = [curve_fitting.logfit(x, *popt) for x in genes_remaining]
        plt.plot(genes_remaining, fit_wgd_times, c='c', label="log fit (model)\n"
                                                              +str(popt))
        print("logfit_fit_popt:\t" + str(popt))

        example_genes_remaining=np.arange(500, 3000, 500)
        predicted_wgd_times= [curve_fitting.logfit(x, *popt) for x in example_genes_remaining]
        plt.bar(example_genes_remaining,predicted_wgd_times, width=10, color='pink',label="model input (# genes)")
        plt.scatter(example_genes_remaining,predicted_wgd_times, c='r',label="model output (est. wgd time)")

        ax.set(xlabel="# genes remaining")
        ax.set(ylabel="WGD time")
        plt.legend()
        plot_file=full_path.replace(".csv",".new.png").replace("shed","remaining")

        plt.savefig(plot_file)
        plt.close()

    def test_auto_vs_allo_predictor(self):

        out_folder = "/home/tamsen/Data/Specks_outout_from_mesx"
        file1="genes_remaining_vs_wgd_time_noN100.csv"
        file2="mode_vs_spec_time.csv"
        full_path1 = os.path.join(out_folder, file1)
        full_path2 = os.path.join(out_folder, file2)

        wgd_xs,genes_remaining,wgd_sims = read_xs_ys_csv(full_path1)
        genes_remaining_dict = dict(zip(wgd_sims,genes_remaining))
        wgd_popt=[ 43.89486241,574.88694572,143.86691067]
        print("wgd_popt:\t" + str(wgd_popt))

        spec_xs,mode,spec_sims = read_xs_ys_csv(full_path2)
        mode_dict = dict(zip(spec_sims, mode))
        spec_popt=[97.05019317,-0.18999441]
        print("spec_popt:\t" + str(spec_popt))

        truth_by_sim_name=get_truth_from_name_list(wgd_sims)
        allo_vs_auto_truth_by_sim={}
        allo_vs_auto_prediction_by_sim={}
        print("strating truth assemby")
        sims_names=truth_by_sim_name.keys()
        predicted_indexes_for_true_autos=[]
        predicted_indexes_for_true_allos=[]
        for sim_name in sims_names:

            print(sim_name)
            #if "uto" in sim_name:
            #    print("FOUND ONE!!!!!!!!!!!!!")

            polyploid_params=truth_by_sim_name[sim_name]
            how_auto_t=polyploid_params.SPC_time_MYA - polyploid_params.WGD_time_MYA
            truth=[polyploid_params.SPC_time_MYA,polyploid_params.WGD_time_MYA,how_auto_t]
            allo_vs_auto_truth_by_sim[sim_name]=truth
            print("true spec {0}\ttrue wgd {1}".format(polyploid_params.SPC_time_MYA,polyploid_params.WGD_time_MYA))

            genes_remaing= genes_remaining_dict[sim_name]
            wgd_prediction=curve_fitting.logfit(genes_remaing, *wgd_popt)
            mode= mode_dict[sim_name]
            spec_prediction=curve_fitting.linear(mode, *spec_popt)
            how_auto_predicition=spec_prediction - wgd_prediction
            prediction = [spec_prediction, wgd_prediction,how_auto_predicition]
            allo_vs_auto_prediction_by_sim[sim_name]=prediction
            #print("predicted index:\t" + str(how_auto_predicition))
            print("predicted spec {0}\tpredicted wgd {1}".format(spec_prediction, wgd_prediction))
            print("\n")

            if (polyploid_params.WGD_time_MYA==polyploid_params.SPC_time_MYA):
                predicted_indexes_for_true_autos.append(how_auto_predicition)
            else:
                predicted_indexes_for_true_allos.append(how_auto_predicition)

        highest_predicted_index_for_true_autos=max(predicted_indexes_for_true_autos)
        lowest_predicted_index_for_true_allos=min(predicted_indexes_for_true_allos)
        print("highest_predicted_index_for_true_autos:\t" + str(highest_predicted_index_for_true_autos))
        print("lowest_predicted_index_for_true_allos:\t" + str(lowest_predicted_index_for_true_allos))

        discrim_criteria_midpoint= 0.5 *(highest_predicted_index_for_true_autos + lowest_predicted_index_for_true_allos)
        print("midpoint:\n" + str(discrim_criteria_midpoint))

        tests=["spec","wgd","allopolyploid_index"]
        sims_names_list=list(sims_names)
        num_data_points=len(sims_names)
        for i in range(0,len(tests)):
            ava_truth=[allo_vs_auto_truth_by_sim[s][i] for s in sims_names_list]
            ava_predictions=[allo_vs_auto_prediction_by_sim[s][i] for s in sims_names_list]
            fig, ax = plt.subplots(1, 1, figsize=(5, 5))
            plt.scatter(ava_truth,ava_predictions, alpha=0.25,
                        c='k',label=tests[i] +", truth vs prediction\nn={0}".format(num_data_points))
            if i == 2:
                ax.set(xlabel="<--auto    truth (MY)   allo-->")
                ax.set(ylabel="<--auto  prediction (MY) allo-->")
                ax.set(title="Predictions of polyploid index (SPEC - WGD time, in MY)")

                ax.axhline(y=discrim_criteria_midpoint, color='r', linestyle='--', label="discimination criteria"
                       + " (y={0})".format(round(discrim_criteria_midpoint,2)))
            else:
                ax.set(xlabel="truth")
                ax.set(ylabel="prediction")
            plt.legend()
            plot_file = os.path.join(out_folder,tests[i] +"_truth_vs_predictions.png")
            plt.savefig(plot_file)
            plt.close()

            plot_data = [sims_names_list,ava_truth,ava_predictions]
            data_file = plot_file.replace("png", "csv")
            save_metrics_to_csv(plot_data,data_file)

    #https://stackoverflow.com/questions/25009284/how-to-plot-roc-curve-in-python


    def test_highN_vs_lowN_predictor(self):

        out_folder = "/home/tamsen/Data/Specks_outout_from_mesx"
        file1="metric3.csv"
        full_path1 = os.path.join(out_folder, file1)

        spc_xs,metric,wgd_sims = read_xs_ys_csv(full_path1)
        low_N_sims= ["Auto","sim37_N0p1", "sim37_N1"]
        med_N_sims=["sim37_N5", "sim35_log"]
        high_N_sims=["sim36_N10","sim37_N20"]

        truth_by_sim={}
        metrics_by_category = {"Low":[],"Medium":[],"High":[]}
        spec_by_category = {"Low":[],"Medium":[],"High":[]}
        cutoff=60
        for i in range(0,len(wgd_sims)):
            sim = wgd_sims[i]
            spc_time=spc_xs[i]
            if spc_time >= cutoff:
                continue
            if metric[i] <= 0:
                continue
            true_category= categorize_sim(low_N_sims, med_N_sims, high_N_sims, sim)
            truth_by_sim[sim]=true_category
            #print(metric[i])
            #plus_1=metric[i]+1
            truth_by_sim[sim] = true_category
            metrics_by_category[true_category].append(metric[i])
            spec_by_category[true_category].append(spc_xs[i])


        max_low_N_result=max(metrics_by_category["Low"])
        min_med_N_result=min(metrics_by_category["Medium"])
        max_med_N_result=max(metrics_by_category["Medium"])
        min_high_N_result=min(metrics_by_category["High"])

        print(statistics.stdev(metrics_by_category["Low"]))
        low_n_mean=(statistics.mean(metrics_by_category["Low"]))
        med_n_mean=(statistics.mean(metrics_by_category["Medium"]))
        high_n_mean=(statistics.mean(metrics_by_category["High"]))
        low_vs_medium_discrimination = 0.5 * sum([max_low_N_result, min_med_N_result ]) #0.01
        medium_vs_high_discrimination = 0.5 * sum([max_med_N_result, min_high_N_result]) #0.05

        fig, ax = plt.subplots(1, 1, figsize=(5, 5))
        colors_by_category= {"Low":"gray","Medium":"cyan","High":"blue"}

        accuracy_by_sim={}
        for i in range(0,len(wgd_sims)):
            sim = wgd_sims[i]
            true_category= categorize_sim(low_N_sims, med_N_sims, high_N_sims, sim)
            plt.scatter(spc_xs[i],metric[i], alpha=0.25,
                        c=colors_by_category[true_category])

            predicted_category=predict_category(
                low_vs_medium_discrimination, medium_vs_high_discrimination,metric[i])
            accuracy_by_sim[sim]=[true_category,predicted_category]

        ax.axhline(y=low_vs_medium_discrimination, color='cyan', linestyle='--', label="lvm disc. criteria"
                + " ({0})".format(round(low_vs_medium_discrimination, 2)))

        ax.axhline(y=medium_vs_high_discrimination, color='purple', linestyle='--', label="mvh disc. criteria"
                + " ({0})".format(round(medium_vs_high_discrimination, 2)))

        ax.set(yscale='log')
        plt.legend()
        plot_file = os.path.join(out_folder, "highN_vs_lowN_predictor.png")
        plt.savefig(plot_file)
        plt.close()

        means_by_category= {"Low":low_n_mean,"Medium":med_n_mean,"High":high_n_mean}
        stutter=0.5
        fig, ax = plt.subplots(1, 1, figsize=(5, 5))

        ax.add_patch(Rectangle((0, 0), low_vs_medium_discrimination, 1, color='gray', alpha=0.15))
        ax.add_patch(Rectangle((0, 0), 1, low_vs_medium_discrimination, color='gray', alpha=0.15))

        ax.add_patch(Rectangle((low_vs_medium_discrimination,0),
                               medium_vs_high_discrimination-low_vs_medium_discrimination, 1,
                               color='cyan', alpha=0.15))
        ax.add_patch(Rectangle((0, low_vs_medium_discrimination),
                               1, medium_vs_high_discrimination-low_vs_medium_discrimination,
                               color='cyan', alpha=0.15))

        ax.add_patch(Rectangle((medium_vs_high_discrimination, 0),
                               1- medium_vs_high_discrimination, 1,
                               color='blue', alpha=0.15))
        ax.add_patch(Rectangle((0, medium_vs_high_discrimination),
                               1, 1- medium_vs_high_discrimination,
                               color='blue', alpha=0.15))


        for sim,results in accuracy_by_sim.items():

            [true_category,predicted_category] = accuracy_by_sim[sim]
            x_value=means_by_category[true_category]*(1+stutter*random.random())
            y_value=means_by_category[predicted_category]*(1+stutter*random.random())

            plt.scatter(x_value, y_value, alpha=0.25,
                        c=colors_by_category[true_category])

        ax.axvline(x=low_vs_medium_discrimination, color='cyan', linestyle='--', label="lvm disc. criteria"
                                                                                       + " ({0})".format(
            round(low_vs_medium_discrimination, 2)))

        #ax.add_patch(Rectangle((0, 0), low_vs_medium_discrimination, 1), color='cyan')

        ax.axvline(x=medium_vs_high_discrimination, color='purple', linestyle='--', label="mvh disc. criteria"
                                                                                          + " ({0})".format(
            round(medium_vs_high_discrimination, 2)))
        ax.set(xlabel="<-- truth -->")
        ax.set(ylabel="<-- prediction -->")
        ax.set(xscale='log')
        ax.set(yscale='log')
        plt.legend()
        plot_file = os.path.join(out_folder, "highN_vs_lowN_accuracy.png")
        plt.savefig(plot_file)
        plt.close()

        make_box_plots(colors_by_category, metrics_by_category, out_folder)

def predict_category(low_vs_medium_discrimination, medium_vs_high_discrimination, metric):
    if metric < low_vs_medium_discrimination:
        return "Low"
    elif metric < medium_vs_high_discrimination:
        return "Medium"
    return "High"

def make_box_plots(colors_by_category, metrics_by_category, out_folder):
    fig, ax = plt.subplots(1, 1, figsize=(5, 5))
    # https://matplotlib.org/stable/gallery/statistics/boxplot_color.html#sphx-glr-gallery-statistics-boxplot-color-py
    boxplot_data = [metrics_by_category["Low"], metrics_by_category["Medium"], metrics_by_category["High"]]
    labels = ["Low", "Medium", "High"]
    bplot = ax.boxplot(boxplot_data, labels=labels, patch_artist=True)
    for patch, label in zip(bplot['boxes'], labels):
        patch.set_facecolor(colors_by_category[label])
    plot_file = os.path.join(out_folder, "highN_vs_lowN_boxplots.png")
    plt.savefig(plot_file)
    plt.close()


def categorize_sim(low_N_sims, med_N_sims, high_N_sims, sim):
    for n in low_N_sims:
        if n in sim:
            return "Low"
    for n in med_N_sims:
        if n in sim:
            return "Medium"
    for n in high_N_sims:
        if n in sim:
            return "High"
    else:
        return False


if __name__ == '__main__':
    unittest.main()

def save_metrics_to_csv(plot_data_for_sims, out_file_name):

    with open(out_file_name, 'w') as f:
        data_headers= ['sims','truth','prediction']
        f.writelines(",".join(data_headers) +"\n")
        [sims_names, ava_truth, prediction] = plot_data_for_sims
        for i in range(0,len(sims_names)):
            data_list=[sims_names[i], ava_truth[i], prediction[i]]
            data_list_string=[str(d) for d in data_list]
            f.writelines(",".join(data_list_string) +"\n")


def read_xs_ys_csv(csv_file):

    batches=[]
    sim_names=[]
    xs=[]
    ys=[]
    with open(csv_file, "r") as f:

        while True:
            line = f.readline()
            if "batch" in line:
                continue
            if len(line)==0:
                break
            if "NA" in line:
                continue
            data = line.strip().split(",")
            #batches.append(data[0])
            sim_names.append(data[1]+ "_" + data[0])
            xs.append(float(data[2]))
            ys.append(float(data[3]))

    return xs,ys,sim_names
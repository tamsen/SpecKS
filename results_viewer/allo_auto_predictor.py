import math
import os
import unittest

import numpy as np
from scipy.optimize import curve_fit
from sklearn import linear_model
from matplotlib import pyplot as plt

from results_viewer import curve_fitting
from results_viewer.batch_histogrammer import get_truth_from_name_list


class AlloAutoPredictor(unittest.TestCase):
    def test_spec_time_predictor(self):

        out_folder = "/home/tamsen/Data/Specks_outout_from_mesx"
        file1="mode_vs_spec_time.csv"

        full_path=os.path.join(out_folder,file1)
        spec_xs,mode_ys, sims = read_xs_ys_csv(full_path)

        fig, ax = plt.subplots(1, 1, figsize=(5, 5))
        plt.scatter(mode_ys, spec_xs, c='b',label="empirical data")

        fit_spec_times, xs, goodness_of_fit=curve_fitting.fit_curve_to_xs_and_ys(
            mode_ys, spec_xs, curve_fitting.linear)

        linear_fit_popt=goodness_of_fit.popt
        ax.set(xlabel="mode (Ks space)")
        ax.set(ylabel="spec time as a fxn of mode (MYA)")
        plt.plot(xs,fit_spec_times, c='k',label="line fit (model)")

        example_modes=np.arange(0.2, 1.0, 0.2)
        predicted_spec_times= [curve_fitting.linear(x, *linear_fit_popt) for x in example_modes]
        plt.bar(example_modes,predicted_spec_times, width=0.002, color='pink',label="model input (mode)")
        plt.scatter(example_modes,predicted_spec_times, c='r',label="model output (est. spec time)")


        plot_file=full_path.replace(".csv",".new.png")
        plt.legend()
        plt.savefig(plot_file)
        plt.close()

    def test_wgd_time_predictor(self):

        out_folder = "/home/tamsen/Data/Specks_outout_from_mesx"
        file1="genes_shed_vs_wgd_time.csv"
        full_path=os.path.join(out_folder,file1)
        wgd_xs,genes_shed,sims = read_xs_ys_csv(full_path)
        genes_remaining = [3000-g for g in genes_shed]

        fig, ax = plt.subplots(1, 1, figsize=(5, 5))
        plt.scatter(genes_remaining, wgd_xs, c='b',label="empirical data")


        popt, pcov = curve_fit(curve_fitting.logfit, genes_remaining, wgd_xs)
        genes_remaining.sort()
        fit_wgd_times = [curve_fitting.logfit(x, *popt) for x in genes_remaining]
        plt.plot(genes_remaining, fit_wgd_times, c='c', label="line fit (model)")


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
        file1="genes_shed_vs_wgd_time.csv"
        file2="mode_vs_spec_time.csv"
        full_path1 = os.path.join(out_folder, file1)
        full_path2 = os.path.join(out_folder, file2)

        wgd_xs,genes_shed,wgd_sims = read_xs_ys_csv(full_path1)
        genes_remaining = [3000-g for g in genes_shed]
        genes_remaining_dict = dict(zip(wgd_sims,genes_remaining))
        wgd_popt, pcov = curve_fit(curve_fitting.logfit, genes_remaining, wgd_xs)

        spec_xs,mode,spec_sims = read_xs_ys_csv(full_path2)
        mode_dict = dict(zip(spec_sims, mode))
        spec_popt, pcov = curve_fit(curve_fitting.linear, mode, wgd_xs)

        truth_by_sim_name=get_truth_from_name_list(wgd_sims)
        allo_vs_auto_truth_by_sim={}
        allo_vs_auto_prediction_by_sim={}
        print("strating truth assemby")
        sims_names=truth_by_sim_name.keys()
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

        tests=["spec","wgd","allopolyploid_index"]

        for i in range(0,len(tests)):
            ava_truth=[allo_vs_auto_truth_by_sim[s][i] for s in sims_names]
            ava_predictions=[allo_vs_auto_prediction_by_sim[s][i] for s in sims_names]
            fig, ax = plt.subplots(1, 1, figsize=(5, 5))
            plt.scatter(ava_truth,ava_predictions, c='r',label=tests[i] +", truth vs prediction")
            ax.set(xlabel="truth")
            ax.set(ylabel="prediction")
            plt.legend()
            plot_file = os.path.join(out_folder,tests[i] +"_truth_vs_predictions.png")
            plt.savefig(plot_file)
            plt.close()
if __name__ == '__main__':
    unittest.main()

def read_xs_ys_csv(csv_file):

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
            sim_names.append(data[1])
            xs.append(float(data[2]))
            ys.append(float(data[3]))

    return xs,ys,sim_names
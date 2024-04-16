import math
import os
import unittest

import numpy as np
from scipy.optimize import curve_fit
from sklearn import linear_model
from matplotlib import pyplot as plt

from results_viewer import curve_fitting


class AlloAutoPredictor(unittest.TestCase):
    def test_spec_time_predictor(self):

        out_folder = "/home/tamsen/Data/Specks_outout_from_mesx"
        file1="mode_vs_spec_time.csv"
        #genes_shed_vs_wgd_time.csv

        full_path=os.path.join(out_folder,file1)
        spec_xs,mode_ys = read_xs_ys_csv(full_path)

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
        plt.scatter(example_modes,predicted_spec_times, c='r',label="model output (inf. spec time)")


        plot_file=full_path.replace(".csv",".new.png")
        plt.legend()
        plt.savefig(plot_file)
        plt.close()

    def test_wgd_time_predictor(self):

        out_folder = "/home/tamsen/Data/Specks_outout_from_mesx"
        file1="genes_shed_vs_wgd_time.csv"

        full_path=os.path.join(out_folder,file1)
        wgd_xs,genes_shed = read_xs_ys_csv(full_path)
        genes_remaining = [3000-g for g in genes_shed]

        fig, ax = plt.subplots(1, 1, figsize=(5, 5))
        plt.scatter(genes_remaining, wgd_xs, c='b',label="empirical data")
        t=type(math.log(10))
        print(t)
        #    return m * (x**c) + b
        #param_bounds = ([0.02, -4.4, 0.9], [0.03,-4.0, 1.1])
        xs=[1000,2000,3000]
        ys=[50,20,10]
        popt, pcov = curve_fit(curve_fitting.logfit, xs, ys)
        #popt, pcov = curve_fit(curve_fitting.linear, mode_ys, spec_xs)
        #popt=[0.02914589, -4.04914737, 1.05]
        print(popt)
        #fit_spec_times, xs, goodness_of_fit=curve_fitting.fit_curve_to_xs_and_ys(
        #        mode_ys, spec_xs, curve_fitting.exponential)

        #linear_fit_popt=goodness_of_fit.popt
        genes_remaining.sort()
        fit_wgd_times = [600 - 38 * math.log(2000 * x) for x in genes_remaining]
        plt.plot(genes_remaining, fit_wgd_times, c='k', label="line fit (model)")

        fit_wgd_times = [curve_fitting.logfit(x, *popt) for x in genes_remaining]
        plt.plot(genes_remaining, fit_wgd_times, c='c', label="line fit (model)")


        #example_modes=np.arange(0, 2500, 500)
        #predicted_spec_times= [curve_fitting.exponential(x, *linear_fit_popt) for x in example_modes]
        #plt.bar(example_modes,predicted_spec_times, width=0.002, color='pink',label="model input (mode)")
        #plt.scatter(example_modes,predicted_spec_times, c='r',label="model output (inf. spec time)")

        ax.set(xlabel="# genes shed")
        ax.set(ylabel="WGD time")
        plt.legend()
        plot_file=full_path.replace(".csv",".new.png")

        plt.savefig(plot_file)
        plt.close()



if __name__ == '__main__':
    unittest.main()

def read_xs_ys_csv(csv_file):

    data_by_sim_name={}
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
            xs.append(float(data[2]))
            ys.append(float(data[3]))

    return xs,ys
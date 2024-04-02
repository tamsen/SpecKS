import os
import unittest

from scipy import stats
from scipy.optimize import curve_fit
import numpy as np
from matplotlib import pyplot as plt

import config
import process_wrapper
from results_viewer import curve_fitting
from results_viewer.multi_run_viewer import read_Ks_csv, make_subplot
from results_viewer.parse_aggregate_results import metric_result, get_metric_result_data_headers


#https://stackoverflow.com/questions/14770735/how-do-i-change-the-figure-size-with-subplots
class KsFittingTests(unittest.TestCase):

    def test_Ks_fit(self):

        csv_folder = "/home/tamsen/Data/Specks_outout_from_mesx/sim19_N1/Auto5"
        csv_file_base0 = "Auto5_S010W010_ML_rep0_Ks_by_GeneTree.csv"
        plot_title=csv_file_base0.replace(".csv",".png")
        full_csv_path=os.path.join(csv_folder,csv_file_base0)
        WGD_time_MYA = 5
        SPC_time_MYA = 10
        WGD_time_Ks = WGD_time_MYA * 0.01
        SPC_time_Ks = SPC_time_MYA * 0.01
        max_Ks_for_x_axis = 0.7

        Ks_results = read_Ks_csv(full_csv_path)

        fig, ax = plt.subplots(1, 1, figsize=(10, 10))
        fig.suptitle(plot_title)
        this_ax = ax

        x = Ks_results
        num_gene_pairs_str = str(len(Ks_results))
        #if max_Ks:
        #    bins = np.arange(0, max_Ks + 0.1, bin_size)
        #    n, bins, patches = this_ax.hist(x, bins=bins, facecolor=plot_color, alpha=0.25,
        #                                    label='ks hist ({0} pairs)'.format(num_gene_pairs_str))

        nBins = 100
        n, bins, patches = this_ax.hist(x, bins=nBins, facecolor='b', alpha=0.25,
                                            label='ks hist ({0} pairs)'.format(num_gene_pairs_str))

        this_ax.axvline(x=WGD_time_Ks, color='b', linestyle='-', label="WGD time, " + str(WGD_time_MYA) + " MYA")
        this_ax.axvline(x=SPC_time_Ks, color='r', linestyle='--',
                        label="SPEC time, " + str(SPC_time_MYA) + " MYA")

        fit_curve_ys1, xs_for_wgd, mode, cm, popt = curve_fitting.fit_curve_to_hist(bins, n)
        print("mode " + str(mode))
        print("cm " + str(cm))

        # https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.kstest.html
        ks_norm_result = stats.kstest(bins, stats.norm.cdf)
        # ks_lognorm_result = stats.kstest(bins, stats.lognorm.cdf(popt))
        # https://stackoverflow.com/questions/51902996/scipy-kstest-used-on-scipy-lognormal-distrubtion
        my_args = popt[0:3]
        ks_lognorm_result = stats.kstest(bins, 'lognorm', args=my_args)

        print("ks_norm_result " + str(ks_norm_result))
        print("ks_lognorm_result " + str(ks_lognorm_result))

        this_ax.plot(xs_for_wgd, fit_curve_ys1,color='green', linestyle=':', label="lognorm fit")


        this_ax.scatter(cm, 0.05,color='darkgreen', marker='o', s=100)  # label="cm",)
        this_ax.scatter(mode, 0.05,color='cyan', marker='^', s=80)  # label="mode")
        this_ax.legend()
        this_ax.set(xlabel="Ks")
        out_file_name = full_csv_path.replace(".csv", ".hist.png")

        out_file_name = out_file_name.replace(".png", "_Ks_fit" + str(max_Ks_for_x_axis) + ".png")
        plt.savefig(out_file_name)
        plt.close()
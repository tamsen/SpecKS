import os
import unittest

import numpy as np
from matplotlib import pyplot as plt

from results_viewer import curve_fitting


class BatchAnalyser(unittest.TestCase):

    def test_compute_metrics_for_batch(self):

        output_folder="/home/tamsen/Data/Specks_outout_from_mesx/sim35_log"
        files = os.listdir(output_folder)
        max_Ks=1
        maxY= False
        hist_data_by_file={}
        for file in files:
            if not ".hist.csv" in file:
                continue
            #print(file)
            full_path=os.path.join(output_folder,file)
            hist_data=read_hist_csv(full_path)
            hist_data_by_file[full_path]=hist_data

        for file_path,hist_data in hist_data_by_file.items():
            print(file_path)
            out_file_name = file_path.replace(".csv", "_Ks_hist_fit" + str(max_Ks) + ".png")
            WGD_time_MYA =70
            SPC_time_MYA=60
            plot_histogram(*hist_data, WGD_time_MYA, SPC_time_MYA,
                           max_Ks, maxY, 'blue', out_file_name)

            break

def read_hist_csv(csv_file):

    num_per_bin = []
    bins=[]

    with open(csv_file, "r") as f:

        while True:

            line = f.readline()
            if not line:
                break

            if "bin" in line:
                continue

            data=line.split(',')
            bins.append(float(data[0]))
            num_per_bin.append(float(data[1]))

    return [bins, num_per_bin ]



def plot_histogram(bins, n, WGD_time_MYA, SPC_time_MYA,
                           max_Ks, maxY, plot_color, out_file_name):

    fig, ax = plt.subplots(1, 1, figsize=(10, 10))
    fig.suptitle('plot_title')
    WGD_as_Ks = WGD_time_MYA * 0.01 #/ 1.04 ..the peak max is about 96% off from where it should be
    SPEC_as_Ks =SPC_time_MYA * 0.01 #/ 1.04
    plt.bar(bins, n, width=0.001)
    #plt.plot(bins,n)
    #plt.savefig(out_file_name)
    #plt.close()
    #print("saving " + out_file_name)
    #return
    num_gene_pairs_str=str(sum(n))
    hist_maximum=max(n)
    ymax_suggestion=hist_maximum*1.6

    #do_curve_fitting=False
    #do_curve_fitting = (spec != "outgroup" ) and (max_Ks and (max_Ks > 0.1))
    curve_fit_done=False

    lognorm_goodness_of_fit =False
    gaussian_goodness_of_fit =False
    ax.axvline(x=WGD_as_Ks, color='b', linestyle='-', label="WGD time, " + str(WGD_time_MYA) + " MYA")
    ax.axvline(x=SPEC_as_Ks, color='r', linestyle='--', label="SPEC time, " + str(SPC_time_MYA) + " MYA")

    gaussian_fit_curve_ys1, xs_for_wgd, gaussian_goodness_of_fit = \
            curve_fitting.fit_curve_to_xs_and_ys(bins, n,curve_fitting.wgd_gaussian)

    lognorm_fit_curve_ys1, xs_for_wgd,lognorm_goodness_of_fit =  \
            curve_fitting.fit_curve_to_xs_and_ys(bins, n, curve_fitting.wgd_lognorm )


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
            plt.plot(xs_for_wgd,lognorm_fit_curve_ys1,
                 color='green', linestyle=':', label="lognorm fit (err={0})".format(rmse_str))

    if gaussian_fit_curve_ys1 and (hist_maximum>0):
            rmse_str= str(round(gaussian_goodness_of_fit.RMSE,4))
            plt.plot(xs_for_wgd,gaussian_fit_curve_ys1,
                 color='blue', linestyle=':', label="gaussian fit (err={0})".format(rmse_str))

        #if SPC_time_MYA == 20:
        #    ymax_suggestion=400


    if maxY:
        yaxis_limit=maxY
        ymax_suggestion=False
    else:
        yaxis_limit= ymax_suggestion

    if curve_fit_done:
        if lognorm_fit_curve_ys1 and (hist_maximum>0):
            plt.scatter(lognorm_goodness_of_fit.cm,0.05*yaxis_limit,
                 color='darkgreen', marker='o', s=100)# label="cm",)

            plt.scatter(lognorm_goodness_of_fit.mode,0.05*yaxis_limit,
                 color='cyan', marker='^', s=80)# label="mode")

    plt.savefig(out_file_name)
    plt.close()

    return plt,ymax_suggestion, lognorm_goodness_of_fit,gaussian_goodness_of_fit


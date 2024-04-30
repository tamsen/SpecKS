import os.path
import unittest

import math
import os
import unittest
import pandas as pd
from scipy.optimize import curve_fit
from scipy.signal import find_peaks
from scipy.stats import poisson, bernoulli,gamma,lognorm,norm
import numpy as np
from matplotlib import pyplot as plt
from sklearn.metrics import mean_squared_error

import config
from results_viewer import batch_analyzer, batch_histogrammer, curve_fitting, batch_aggregator

class TestAgainstKnownAlloAndAutos(unittest.TestCase):

    def test_parse_coffee(self):

        output_folder = "/home/tamsen/Data/Test"
        input_folder="/home/tamsen/Data/SpecKS_input/ks_data"
        ks_file = os.path.join(input_folder,"coffea.ks.tsv")
        ks_data = parse_ks(ks_file)
        counts, bins, bars = plt.hist(ks_data)

        total_counts=sum(counts)
        print(total_counts)
    def test_parse_ks2_old_way(self):

        output_folder = "/home/tamsen/Data/Test"
        input_folder="/home/tamsen/Data/SpecKS_input/ks_data"
        ks_files = ["actinidia.ks.tsv","coffea.ks.tsv",
                    "poplar.ks.tsv","triticum.ks.tsv","mays.ks.tsv","final_ks_values_TORX.fa"]
        ks_files = ["coffea.ks.tsv"]

        full_paths=[os.path.join(input_folder,f) for f in ks_files]
        x_max = 1.0

        for ks_file in full_paths:

            if ".fa" in ks_file: #1KP file
                olea_ks_df = pd.read_csv(ks_file, sep='\t', header=0)
                olea_ks_array = olea_ks_df.loc[:, "Node Ks"]
                ks_data = [k for k in olea_ks_array.tolist() if k <= 2]
            else:    #KS rates input file
                ks_data = parse_ks(ks_file)
            species_name = ks_file.split(".")[0]

            #parameters for all except coffee:1000, 40,x_max = 1.0
            #for coffee
            x_max = 0.5
            fine_points = 1000
            course_points=40
            out_file = os.path.join(output_folder,species_name + "_xmax" + str(x_max) + ".png")
            xs = fit_fxns_to_Ks(ks_data, species_name, fine_points, course_points,
                            x_max, out_file)

    def test_parse_ks1(self):

        output_folder = "/home/tamsen/Data/Test"
        #input_folder='/home/tamsen/Data/ks_rates_output'
        #file="Actinidia_chinensis/paralog_distributions/wgd_actinidia/actinidia.ks.tsv"
        #file = "triticum.ks.tsv"
        output_folder = "/home/tamsen/Data/Test"
        input_folder="/home/tamsen/Data/SpecKS_input/ks_data"
        #ks_files = ["actinidia.ks.tsv","coffea.ks.tsv",
        #            "poplar.ks.tsv","triticum.ks.tsv","mays.ks.tsv","final_ks_values_TORX.fa"]
        ks_files = ["coffea.ks.tsv"]

        ks_files_to_reanalyze=[os.path.join(input_folder,f) for f in ks_files]

        bin_size=0.001
        right_most_ssd_peak = 0.08
        bin_sizex1000 = bin_size * 1000.0
        kernel_size = int(50 * (1.0 / bin_sizex1000))
        max_Ks = 1.0

        if not os.path.exists(output_folder):
            os.makedirs(output_folder)

        results_by_file={}
        for i in range(0,len(ks_files_to_reanalyze)):

            Ks_file = ks_files_to_reanalyze[i]
            ks_values = parse_ks(Ks_file)
            ln_ks_values = [math.log(k) for k in ks_values if k>0 ]
            sample_name = os.path.basename(Ks_file).replace(".ks.tsv","")

            #get hist data
            out_file_name = os.path.join(output_folder, "histogram" + "_plot_" +  sample_name+
                                         "_" + str(max_Ks) + ".png")
            params = config.PolyploidParams(-1, -1, sample_name)

            hist_data = plot_hist(bin_size, ks_values, max_Ks, out_file_name, params, sample_name)

            dummy_file_name = os.path.join(output_folder, "ln_histogram" + "_plot_" +  sample_name+
                                      "_" + str(max_Ks) + ".png")

            fig, ax = plt.subplots(1, 1, figsize=(10, 10))
            fig.suptitle(sample_name)
            n, bins, patches = ax.hist(ln_ks_values , bins=100,  alpha=0.25)
            plt.savefig(dummy_file_name)
            plt.close()

            hist_data2 = [hist_data[0][1:], hist_data[1][1:]]
            csv_file_name = os.path.join(output_folder, sample_name + ".hist.csv")
            batch_histogrammer.save_hist_to_csv(hist_data2, csv_file_name)

            #analyze hist data
            hist_data3 =  batch_analyzer.read_hist_csv(csv_file_name)
            out_fit_png = os.path.join(output_folder, "fit" + "_plot_" +  sample_name+
                                         "_"+ str(max_Ks) + ".png")

            fit_results = batch_analyzer.analyze_histogram(*hist_data3, params.WGD_time_MYA, params.SPC_time_MYA,
                                            kernel_size, False, right_most_ssd_peak, out_fit_png)

            results_by_file[Ks_file]=fit_results
            print(Ks_file + " analyzed")

        out_csv = "1KP_metrics.csv"
        csv_file_name = os.path.join(output_folder, out_csv)
        #if reanalyze:
        batch_analyzer.plot_and_save_metrics(results_by_file, csv_file_name)

        #well_behaved_wgd_histograms=analyze_metric_results(csv_file_name,output_folder)

def plot_hist(bin_size, ks_values, max_Ks, out_file_name, params, sample_name):
    fig, ax = plt.subplots(1, 1, figsize=(10, 10))
    fig.suptitle(sample_name)
    hist_data = batch_histogrammer.make_histogram_subplot(ax, sample_name, ks_values,
                                                          bin_size, params,
                                                          max_Ks, False, "blue")
    plt.savefig(out_file_name)
    plt.close()
    return hist_data


if __name__ == '__main__':
    unittest.main()

def parse_ks(input_file):
    ks_df=pd.read_csv(input_file, sep='\t', header=0)
    ks_array = ks_df.loc[:,"Ks"]
    ks_list = [k for k in ks_array.tolist() if k <= 2]
    #print(ks_list[1:10])
    return ks_list


def plot_ks_histogram(ks_data, plot_title, num_bins, max_value):
    # force all plots to have the exact same bins!
    my_bins = np.arange(0, max_value + (max_value / num_bins), max_value / num_bins)

    # counts, bins, bars  = plt.hist(ks_data, num_bins,
    counts, bins, bars = plt.hist(ks_data, bins=my_bins,
                                  density=1,
                                  color='green',
                                  alpha=0.3)

    plt.title(plot_title)
    plt.xlabel('Ks')
    plt.ylabel('density')
    #plt.show()
    return counts, bins


def bins_to_xs(bins):
    xs = [0.5 * (bins[i] + bins[i + 1]) for i in range(0, len(bins) - 1)]
    return xs


def limit_to_range(xs, ys, xrange):
    x_start = xrange[0]
    x_end = xrange[1]
    xs_in_range = []
    ys_in_range = []

    for i in range(0, len(xs)):
        if (xs[i] >= x_start) and (xs[i] <= x_end):
            xs_in_range.append(xs[i])
            ys_in_range.append(ys[i])

    return xs_in_range, ys_in_range


def smooth_data_with_kernel(ys, kernel_size):
    kernel = np.ones(kernel_size) / kernel_size
    smoothed_trace = np.convolve(ys, kernel, mode='same')
    return (smoothed_trace)


def find_local_minima(xs, ys):
    # threshold,height=threshold,
    min_distance_between_peaks = 1
    min_peak_width = 0.1

    # https://plotly.com/python/peak-finding/
    upside_down_ys = [-1.0 * y for y in ys]
    indices = find_peaks(upside_down_ys,
                         distance=min_distance_between_peaks,
                         width=min_peak_width)[0]

    minima = [(xs[i], ys[i], i) for i in indices]
    # sort, left-most first.
    minima.sort(key=lambda x: x[0], reverse=False)
    return minima


def find_local_maxima(xs, ys):
    # threshold,height=threshold,
    min_distance_between_peaks = 1
    min_peak_width = 0.1

    # https://plotly.com/python/peak-finding/
    indices = find_peaks(ys,
                         distance=min_distance_between_peaks,
                         width=min_peak_width)[0]

    minima = [(xs[i], ys[i], i) for i in indices]
    # sort, left-most first.
    minima.sort(key=lambda x: x[0], reverse=False)
    return minima


def ssd_func(x, C0, C1):
    return C0 * np.exp(- x * C1)


def wgd_gaussian(x, amp, mu, sig):
    return amp * norm.pdf(x, mu, sig)


def wgd_lognorm(x, amp, scale, x_shift, skew):
    return amp * lognorm.pdf(scale * x + x_shift, skew)


def fit_fxn_to_data(xs, ys, xrange, fit_fxn):
    xs_in_range, ys_in_range = limit_to_range(xs, ys, xrange)
    popt, pcov = curve_fit(fit_fxn, xs_in_range, ys_in_range)
    fit_curve_ys = [fit_fxn(x, *popt) for x in xs_in_range]
    return xs_in_range, fit_curve_ys, popt


def simple_smooth(x, side):
    # print(len(x))
    front_flank = [x[0] for i in range(0, side)]
    back_flank = [x[len(x) - 1] for i in range(0, side)]
    dummy = front_flank + list(x) + back_flank
    # print(dummy)
    window_size = side * 2 + 1
    result = []

    for i in range(0, len(x)):
        subseq = dummy[i:i + window_size]
        new_value = sum(subseq) / float(window_size)
        result.append(new_value)

    return result


def fit_fxns_to_Ks(ks_data, sci_name, num_fine_points,
                   num_course_points, x_max, out_file):
    max_ks_value = 2

    # find extremum
    ys_smooth, bins_smooth = plot_ks_histogram(ks_data, sci_name, num_course_points, max_ks_value)
    xs_smooth = bins_to_xs(bins_smooth)
    minima = find_local_minima(xs_smooth, ys_smooth)
    minima_xs = [m[0] for m in minima]
    minima_ys = [m[1] for m in minima]
    maxima = find_local_maxima(xs_smooth, ys_smooth)
    maxima_xs = [m[0] for m in maxima]
    maxima_ys = [m[1] for m in maxima]

    plt.scatter(minima_xs, minima_ys, color="blue", marker="X", s=5)
    plt.scatter(maxima_xs, maxima_ys, color="red", marker="X", s=5)

    ys_400, bins_400 = plot_ks_histogram(ks_data, sci_name, num_fine_points, max_ks_value)
    xs_400 = bins_to_xs(bins_400)
    bin_size = xs_400[1] - xs_400[0]
    len_ks_data=len(ks_data)
    AUC=sum(ys_400)
    print("len_ks_data\t" + str(len_ks_data))
    print("AUC\t" + str(AUC))
    print("scaled total AUC\t" + str(AUC*bin_size))

    plt.hist(ks_data, num_fine_points, density=1,
             color='green',
             alpha=0.3)

    plt.scatter(minima_xs, minima_ys, color="blue", marker="X", s=5)
    plt.scatter(maxima_xs, maxima_ys, color="red", marker="X", s=5)

    # check for very recent WGD hiding in the front..
    ys_through_window = simple_smooth(ys_400, 1)
    closest_extremum = min(minima_xs[0], maxima_xs[0])
    minima_through_window = find_local_minima(xs_400, ys_through_window)
    maxima_through_window = find_local_maxima(xs_400, ys_through_window)
    early_minima = [m for m in minima_through_window[0:1] if m[0] < closest_extremum]
    early_maxima = [m for m in maxima_through_window[0:1] if m[0] < closest_extremum]
    plt.scatter([m[0] for m in early_minima],
                [m[1] for m in early_minima], color="blue", s=20)
    plt.scatter([m[0] for m in early_maxima],
                [m[1] for m in early_maxima], color="red", s=20)

    # ssd trace:
    if minima_xs[0] < maxima_xs[0]:
        ssd_range = [0, minima[0][0]]
        print("SSD range:" + str(ssd_range) + "\n")
        xs_ssd, ys_ssd, popt_ssd = fit_fxn_to_data(xs_400, ys_400, ssd_range, ssd_func)
        full_ssd_ys = [ssd_func(x, *popt_ssd) for x in xs_400]
        ssd_AUC=sum(full_ssd_ys)*bin_size
        print("scaled ssd AUC\t" + str(ssd_AUC))
        # xs_ssd,ys_ssd,popt_ssd=fit_fxn_to_data(xs_400,ys_400,ssd_range,wgd_lognorm)
        # full_ssd_ys = [wgd_lognorm(x, *popt_ssd) for x in xs_400]
    else:
        full_ssd_ys = [0 for x in xs_400]
        minima = [[xs_smooth[0], ys_smooth[0]]] + minima
        print("SSD range:" + "unclear...\n")

    # TODO - check if there is a lognorm component hidden in the SSD range.

    plt.plot(xs_400, full_ssd_ys, color="gold", label="SSD")
    mix_model_so_far = full_ssd_ys

    # add in those early, sketch minima, to see what we find
    ssd_ys_after_early_min = [full_ssd_ys[i] for i in range(0, len(xs_400)) if xs_400[i] >= early_maxima[0][0]]
    sketchy_WGD = False

    if (early_maxima[0][1] > max(ssd_ys_after_early_min)):
        minima = early_minima + minima
        maxima = early_maxima + maxima
        minima_xs = [m[0] for m in minima]
        minima_ys = [m[1] for m in minima]
        maxima_xs = [m[0] for m in maxima]
        maxima_ys = [m[1] for m in maxima]
        sketchy_WGD = True

    # WGD trace:
    num_minima = min(len(minima), 4)
    # num_minima=len(minima)
    wgd_fits = []
    wgd_parameters = []
    colors = ["blue", "orange", "red"]
    if sketchy_WGD:
        colors = ["gray", "blue", "orange", "red"]
        print("first WGD is a rescue attempt")

    for i in range(0, num_minima - 1):
        # print("i"+ str(i))

        wgd_range = [minima[i][0] + bin_size, minima[i + 1][0] - bin_size]

        print("WGD#" + str(i + 1) + " range:" + str(wgd_range))
        remainder = [ys_400[i] - mix_model_so_far[i] for i in range(0, len(xs_400))]
        xss, yss= limit_to_range(xs_400,remainder , wgd_range)

        #plt.plot(xss, yss,
        #         color=colors[i], linewidth=3, label="fit WGD#" + str(i+1), )
        try:
            xs_wgd, ys_wgd, popt_wgd = fit_fxn_to_data(xs_400, remainder, wgd_range, wgd_lognorm)
            AUC_wgd = sum(ys_wgd)*bin_size
            print("scaled wgd AUC\t" + str(AUC_wgd))
            plt.plot(xs_wgd,[6 for x in xs_wgd], color=colors[i], linewidth=3, label="WGD range" + str(i+1), )
            full_wgd_ys = [wgd_lognorm(x, *popt_wgd) for x in xs_400]
            center_of_mass, x_value_of_ymax = curve_fitting.get_mode_and_cm(full_wgd_ys, xs_400)
            metric3=center_of_mass-x_value_of_ymax
            metric9 = curve_fitting.get_asymmetry_ratio(full_wgd_ys, xs_400, x_value_of_ymax)
            if metric3 <=0:
                ln_metric3=-10
            else:
                ln_metric3=math.log(metric3)
            xss, fit_in_ranged = limit_to_range(xs_400, full_wgd_ys, wgd_range)
            lg_rms = mean_squared_error(yss, fit_in_ranged , squared=False)
            #label = ("LG WGD#" + str(i + 1) + " RMSE=" + str(round(lg_rms, 3)))
            label = ("LG WGD#" + str(i + 1) +
                     " \nAUC=" + str(round(AUC_wgd, 3)))# +
                     #" \nRMSE=" + str(round(lg_rms, 3)) +
                     #"\nM9=" + str(round(metric9, 3)) +
                     #"\nlnM3=" + str(round(ln_metric3, 3))
                     #)

            plt.plot(xs_400, full_wgd_ys, label=label, color=colors[i])

            #if (i + 1) ==2:
            #    plt.plot(xs_400, remainder, label="remainder"+ str(i + 1) , color='k')

            ga_xs_wgd, ga_ys_wgd, ga_popt_wgd = fit_fxn_to_data(xs_400, remainder, wgd_range,
                                                       curve_fitting.wgd_gaussian)
            if ga_xs_wgd:
                ga_full_wgd_ys = [curve_fitting.wgd_gaussian(x, *ga_popt_wgd) for x in xs_400]
                xss, ga_fit_in_ranged = limit_to_range(xs_400, ga_full_wgd_ys, wgd_range)
                ga_rms = mean_squared_error(yss,ga_fit_in_ranged, squared=False)
                label = ("GA WGD#" + str(i + 1) + " RMSE=" + str(round(ga_rms, 3)) +
                     "\nM9="+str(round(metric9, 3)) +
                     "\nlnM3="+str(round(ln_metric3, 3))
                     )
                plt.plot(xs_400, ga_full_wgd_ys, label=label, color=colors[i])
            wgd_fits.append(full_wgd_ys)
            wgd_parameters.append(popt_wgd)
            mix_model_so_far = [mix_model_so_far[i] + full_wgd_ys[i] for i in range(0, len(xs_400))]
            print("WGD#" + str(i + 1) + " parameters:" + str(popt_wgd))
        except Exception as e:

            if (i == 0 and sketchy_WGD):
                try:
                    smoothed_remainer = simple_smooth(ys_400, 1)
                    xss, yss = limit_to_range(xs_400, smoothed_remainer, wgd_range)
                    plt.plot(xss, yss, color=colors[i])
                    xs_wgd, ys_wgd, popt_wgd = fit_fxn_to_data(xs_400, smoothed_remainer,
                                                               wgd_range, wgd_lognorm)
                    full_wgd_ys = [wgd_lognorm(x, *popt_wgd) for x in xs_400]
                    center_of_mass, x_value_of_ymax = curve_fitting.get_mode_and_cm(full_wgd_ys, xs_400)
                    metric3 = center_of_mass - x_value_of_ymax
                    metric9 = curve_fitting.get_asymmetry_ratio(full_wgd_ys, xs_400, x_value_of_ymax)
                    if metric3 <= 0:
                        ln_metric3 = -10
                    else:
                        ln_metric3 = math.log(metric3)
                    label = ("GA WGD#" + str(i + 1) + "? RMSE=" + str(round(ga_rms, 3)) +
                             "\nM9=" + str(round(metric9, 3)) +
                             "\nlnM3=" + str(round(ln_metric3, 3))
                             )
                    plt.plot(xs_400, full_wgd_ys, label=label, color=colors[i])
                    wgd_fits.append(full_wgd_ys)
                    wgd_parameters.append(popt_wgd)
                    mix_model_so_far = [mix_model_so_far[i] for i in range(0, len(xs_400))]
                    print("WGD#" + str(i + 1) + " parameters:" + str(popt_wgd))

                except Exception as e:
                    print(str(e))
                    print("no clear WGD #" + str(i + 1) + " on range " + str(wgd_range))

                    plt.plot(xss, yss, color=colors[i], label="fit WGD#" + str(i + 1) + " failed fit", )
            else:
                print(str(e))
                print("no clear WGD #" + str(i + 1) + " on range " + str(wgd_range))

                plt.plot(xss, yss, color=colors[i], label="fit WGD#" + str(i + 1) + " failed fit", )

        finally:
            print("testing for WGD#" + str(i + 1) + " complete.\n")

    # sum of traces (mixture model result)
    plt.plot(xs_400, mix_model_so_far, color='k', linestyle=":", label="mixture model")

    # for i in range(0,len(popt_wgd)):
    #    #plt.text(0.5,8-i*0.25,str(popt_wgd[i])
    #    plt.text(1.5,6-i*0.25,str(popt_wgd[i]))

    # make it pretty
    plt.title(sci_name + " Ks Fit")
    plt.xlabel('Ks')
    plt.ylabel('density')
    plt.legend()
    plt.xlim(0, x_max)
    plt.ylim(0, 10)
    plt.savefig(out_file)
    #plt.show()
    # plt.clf()
    plt.close()
    return xs_400, wgd_parameters

#get_low_to_medium_threshold= -4.53
#get_medium_to_high_threshold= -3.12
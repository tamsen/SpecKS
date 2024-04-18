import os
import unittest

import numpy as np
from matplotlib import pyplot as plt
from scipy.signal import find_peaks
from results_viewer import curve_fitting
from results_viewer.batch_histogrammer import get_truth_from_name_list
from results_viewer.run_metrics import run_metrics, plot_and_save_metrics


class BatchAnalyser(unittest.TestCase):

    def test_compute_metrics_for_batch(self):

        batch_name= "sim37_N100"
        base_output_folder = "/home/tamsen/Data/Specks_outout_from_mesx/"
        compute_metrics_for_batch(batch_name,base_output_folder )

    def test_compute_metrics_for_csv(self):
        csv_file = "Allo2_S070W060.hist.csv"  # "Allo6_S030W010.hist.csv"
        max_Ks = 1.0
        output_folder = "/home/tamsen/Data/Specks_outout_from_mesx/sim37_N1/analysis"
        # /home/tamsen/Data/Specks_outout_from_mesx/sim36_N10/analysis/Allo1_S080W075.hist.csv
        full_path = os.path.join(output_folder, csv_file)
        hist_data = read_hist_csv(full_path)
        polyploid_params_by_name = get_truth_from_name_list([csv_file])
        params = polyploid_params_by_name[csv_file]
        out_file_name = full_path.replace(".csv", "_Ks_hist_fit" + str(max_Ks) + ".png")
        fit_results = analyze_histogram(*hist_data, params.WGD_time_MYA, params.SPC_time_MYA,
                                        max_Ks, False, 'tan', out_file_name)

        results_by_file = {"test": fit_results}
        out_csv = "{0}_metrics.csv".format("test")
        out_file_name = os.path.join(output_folder, out_csv)
        plot_and_save_metrics(results_by_file, out_file_name)


def compute_metrics_for_batch(batch_name,base_output_folder):

    analysis_folder="analysis"
    output_folder=os.path.join(base_output_folder,batch_name,analysis_folder)

    files = os.listdir(output_folder)
    max_Ks=2.0
    maxY= False
    hist_data_by_file={}
    polyploid_names=[]
    results_by_file={}

    #gather data
    for file in files:
        if not ".hist.csv" in file:
            continue
        #print(file)
        full_path=os.path.join(output_folder,file)
        hist_data=read_hist_csv(full_path)
        hist_data_by_file[file]=hist_data
        polyploid_names.append(file)#file.replace("hist.cvs",""))

    #run analysis
    polyploid_params_by_name = get_truth_from_name_list(polyploid_names)
    for file,hist_data in hist_data_by_file.items():
        file_path = os.path.join(output_folder, file)
        params =polyploid_params_by_name[file]
        print("starting " + file_path)
        out_file_name = file_path.replace(".csv", "_Ks_hist_fit" + str(max_Ks) + ".png")
        fit_results = analyze_histogram(*hist_data, params.WGD_time_MYA, params.SPC_time_MYA,
                                        max_Ks, maxY, 'blue', out_file_name)
        results_by_file[file]=fit_results
        print(file_path + " analyzed")

    out_csv = "{0}_metrics.csv".format(batch_name)
    out_file_name = os.path.join(output_folder, out_csv)
    plot_and_save_metrics(results_by_file, out_file_name)

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



def analyze_histogram(bins, n, WGD_time_MYA, SPC_time_MYA,
                      max_Ks, maxY, hist_plot_color, out_file_name):

    polyploid_name=os.path.basename(out_file_name).replace(".png","")
    fig, ax = plt.subplots(1, 1, figsize=(10, 10))
    fig.suptitle('Histogram analysis for ' + polyploid_name)
    WGD_as_Ks = WGD_time_MYA * 0.01 #/ 1.04 ..the peak max is about 96% off from where it should be
    SPEC_as_Ks =SPC_time_MYA * 0.01 #/ 1.04
    #plt.bar(bins, n, width=0.001, color=hist_plot_color, label="hist",alpha=0.5)


    num_gene_pairs_str=str(sum(n))
    hist_maximum=max(n)

    ax.axvline(x=WGD_as_Ks, color='b', linestyle='-', label="WGD time, " + str(WGD_time_MYA) + " MYA")
    ax.axvline(x=SPEC_as_Ks, color='r', linestyle='--', label="SPEC time, " + str(SPC_time_MYA) + " MYA")

    kernel_size=50
    smoothed_color='gray'
    smoothed_ys = smooth_data(kernel_size, n)
    plt.plot(bins, smoothed_ys,color=smoothed_color, label="smoothed data")


    smoothed_minima = find_local_minima(bins, smoothed_ys)
    smoothed_maxima = find_global_maxima(bins, smoothed_ys, 0.075)
    plt.scatter([m[0] for m in smoothed_minima], [m[1] for m in smoothed_minima], color="red", label="minima", marker='*')

    ssd_end, next_min =smallest_min(smoothed_minima,smoothed_maxima)
    ssds_xs_to_subtract, ssds_ys_to_subtract = estimate_overlap(bins, ssd_end, next_min)

    plt.scatter(ssd_end[0],3*ssd_end[1],color="black", label="wgd_start",marker='*',s=50)
    #super_smoothed_ys=smooth_data(200, n)
    ssd_xs, ssd_ys, wgd_xs, wgd_ys = sort_ssds_and_wgds(bins, n, ssd_end, ssds_ys_to_subtract)

    wgd_maxima = find_global_maxima(wgd_xs, wgd_ys, 0.075)
    plt.scatter(wgd_maxima[0], -0.05 * maxY,
                color='blue', marker='^', s=80, label="wgd max")

    wgd_max_d=get_loc_of_max_derivative(kernel_size, wgd_ys, wgd_xs)
    raw_cm, raw_x_value_of_ymax = curve_fitting.get_mode_and_cm(wgd_ys, wgd_xs)
    ax.set(xlabel="Ks")
    ax.set(ylabel="# paralogs")
    plt.bar(ssd_xs, ssd_ys, width=0.001, color="lightgray", label="ssd")
    plt.bar(wgd_xs, wgd_ys, width=0.001, color="gray", label="wgd")
    #plt.bar(ssds_xs_to_subtract, ssds_ys_to_subtract, width=0.001, color="lightgray", label="ssd under wgd")

    gaussian_fit_curve_ys1, xs_for_wgd, gaussian_goodness_of_fit = \
            curve_fitting.fit_curve_to_xs_and_ys(wgd_xs, wgd_ys, curve_fitting.wgd_gaussian)

    lognorm_fit_curve_ys1, xs_for_wgd,lognorm_goodness_of_fit =  \
            curve_fitting.fit_curve_to_xs_and_ys(wgd_xs, wgd_ys, curve_fitting.wgd_lognorm )


    if lognorm_fit_curve_ys1 and (hist_maximum>0):
            rmse_str= str(round(lognorm_goodness_of_fit.RMSE,4))
            pser_str= str(round(lognorm_goodness_of_fit.pearsons_corr_result[2],4))
            plt.plot(xs_for_wgd,lognorm_fit_curve_ys1,
                 color='green', linestyle=':',
                     label="lognorm fit \n(RMSE={0})\n(PRSE={1})".format(rmse_str,pser_str))

    if gaussian_fit_curve_ys1 and (hist_maximum>0):
            rmse_str= str(round(gaussian_goodness_of_fit.RMSE,4))
            pser_str= str(round(lognorm_goodness_of_fit.pearsons_corr_result[2],4))
            plt.plot(xs_for_wgd,gaussian_fit_curve_ys1,
                 color='blue', linestyle=':',
                     label="gaussian fit \n(RMSE={0})\n(PRSE={1})".format(rmse_str,pser_str))


    if lognorm_fit_curve_ys1 and (hist_maximum>0):
            plt.scatter(lognorm_goodness_of_fit.cm,-0.05*maxY,
                 color='darkgreen', marker='o', s=100, label="cm",)

            plt.scatter(lognorm_goodness_of_fit.mode,-0.05*maxY,
                 color='cyan', marker='^', s=80, label="mode")

    plt.scatter(wgd_maxima[0], -0.05 * maxY,
                color='green', marker='^', s=80, label="wgd abs max")

    plt.scatter(wgd_max_d[0], -0.05 * maxY,
                color='k', marker='^', s=80, label="wgd max derivative")

    #plt.bar(ssd_xs, ssd_ys, width=0.001, color="lightgray", label="ssd")
    #plt.bar(wgd_xs, wgd_ys, width=0.001, color="gray", label="wgd")
    #plt.bar(ssds_xs_to_subtract, ssds_ys_to_subtract, width=0.001, color="lightgray", label="ssd under wgd")

    plt.legend()
    plt.savefig(out_file_name)
    plt.close()

    fit_data = run_metrics(polyploid_name[0:4], polyploid_name,
                           SPC_time_MYA, WGD_time_MYA,
                           raw_cm, raw_x_value_of_ymax,
                           wgd_maxima[0],wgd_max_d[0],
                           lognorm_goodness_of_fit, gaussian_goodness_of_fit)

    return fit_data


def estimate_overlap(bins, ssd_end, next_min):

    ssds_ys_to_subtract = [0 for x in bins]
    ssds_xs_to_subtract = [x for x in bins]

    if next_min:
        ssd_injection = ssd_end[1] - next_min[1]
        distance_in_ind = next_min[2] - ssd_end[2]
        slope = ssd_injection / distance_in_ind
        ssds_ys_to_subtract = [0 for x in bins]
        ssds_xs_to_subtract = [x for x in bins]
        if ssd_injection > 0:
            for i in range(ssd_end[2], next_min[2]):
                #ssds_xs_to_subtract[i](bins[i])
                ssds_ys_to_subtract[i]=(0.5 * (ssd_end[1] - slope * i))

    return ssds_xs_to_subtract, ssds_ys_to_subtract


def sort_ssds_and_wgds(bins, n, ssd_end, overlap_ssds):

    ssd_xs = [];
    ssd_ys = [];
    wgd_xs = [];
    wgd_ys = [];
    for i in range(0, len(bins)):
        if i < ssd_end[2]:
            ssd_xs.append(bins[i])
            ssd_ys.append(n[i])
            #wgd_xs.append(0)
            #wgd_ys.append(0)
        else:
            wgd_xs.append(bins[i])
            wgd_ys.append(n[i]-overlap_ssds[i])
            #ssd_xs.append(0)
            #ssd_ys.append(0)

    #bin_size=bins[1]-bins[0]
    #print(bin_size)
    #for i in range(0, 100):
    #    wgd_xs= [-1*bin_size*i]+wgd_xs
    #    wgd_ys= [0]+wgd_ys

    return ssd_xs, ssd_ys, wgd_xs, wgd_ys
def smallest_min(minima,peak_max):

    if peak_max:
        smallest_min_y=peak_max[1]
    else:
        smallest_min_y = 1000

    best_min=False
    next_min=False
    for i in range(0, len(minima) - 1):
        m = minima[i]

        # presuming the peak max is the WGD and not the SSD
        # dont look left of the peak maximum
        if peak_max and (peak_max[0] > 0.1):
            if ( m[0] > peak_max[0]):
                break

        if smallest_min_y > m[1]:
            smallest_min_y=m[1]
            best_min=m
            next_min=minima[i+1]

    if next_min and (next_min[1] >= best_min[1]):
        next_min = False

    if not best_min:
        return smallest_min(minima,False)

    return best_min, next_min

def smooth_data(kernel_size, ys):

        kernel = np.ones(kernel_size) / kernel_size
        smoothed_ys = np.convolve(ys, kernel, mode='same')
        return smoothed_ys

def get_loc_of_max_derivative(kernel_size, wgd_ys, wgd_xs):

    smooth_wgd=smooth_data(kernel_size,wgd_ys)
    d=find_derivative(smooth_wgd)
    max_d=find_global_maxima(wgd_xs[0:len(d)], d, 0)
    return max_d

def find_derivative(ys):

    d=[]
    for i in range(0,len(ys)-1):
        d.append(ys[i+1]-ys[i])
    return d
def find_global_maxima(xs,ys,keep_right_of_ssds):

    maxima=0
    maxima_idx=-1
    for i in range(0,len(xs)):

        xi=xs[i]
        yi=ys[i]

        if xi < keep_right_of_ssds:
            continue

        if yi>maxima:
            maxima=yi
            maxima_idx=i
    return [xs[maxima_idx],ys[maxima_idx],maxima_idx]
def find_local_minima(xs,ys):

    #threshold,height=threshold,
    min_distance_between_peaks = 1
    min_peak_width = 1

    # https://plotly.com/python/peak-finding/
    upside_down_ys=[-1.0*y for y in ys]
    indices = find_peaks(upside_down_ys,
                         distance=min_distance_between_peaks,
                         width=min_peak_width)[0]

    minima = [(xs[i], ys[i], i) for i in indices]
    # sort, left-most first.
    minima.sort(key=lambda x: x[0], reverse=False)
    if (len(minima) < 2) and (len(minima) >0):
        minima.append((2.0,minima[-1][1],minima[-1][2]+1))

    #print("minima: " + str(minima))
    return minima

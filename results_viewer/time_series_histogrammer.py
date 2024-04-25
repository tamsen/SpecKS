import os
import unittest
import numpy as np
from matplotlib import pyplot as plt
import batch_histogrammer
import config
from results_viewer.batch_analyzer import smooth_data


class TimeSeriesHistogrammer(unittest.TestCase):

    def test_make_time_series_histograms_for_batch(self):

        batch_name="sim38_test" ##"sim37_N20" #sim37_N0p1,sim37_N5
        plot_title=("Time series for sim ({0})".format(batch_name))

        input_folder=os.path.join( "/home/tamsen/Data/Specks_outout_from_mesx",batch_name)
        output_folder=os.path.join(input_folder,"analysis")
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)

        print("Reading csv files..")
        csvfiles_by_polyploid_by_species_rep_by_algorithm = batch_histogrammer.get_ks_data_from_folders(input_folder)
        params_by_polyploid = batch_histogrammer.get_truth_from_name_dict(csvfiles_by_polyploid_by_species_rep_by_algorithm)
        alg="ML"

        print("Making plots..")
        bin_size = 0.001
        max_Ks = 0.75
        max_Y= 250

        for spec in ['polyploid']:#species:
            get_time_series_histograms_for_runs_in_batch(output_folder, plot_title,
                                                     csvfiles_by_polyploid_by_species_rep_by_algorithm, spec,
                                                     "rep0", alg, params_by_polyploid,
                                                      max_Y, max_Ks, bin_size)


def get_time_series_histograms_for_runs_in_batch(out_folder, sample_name, csvfiles_by_polyploid_by_rep_by_algorthim,
                                     spec, replicate, alg, params_by_polyploid,max_Y, max_Ks_for_x_axis,
                                     bin_size):

    result_names=(list(csvfiles_by_polyploid_by_rep_by_algorthim.keys()))
    result_names.sort()
    ordered_allo1_results=[n for n in result_names if "Allo1" in n]
    ordered_auto_results=[n for n in result_names if "Auto" in n]
    metrics_by_result_names= {}

    # making subplots
    spec_times= [params_by_polyploid[auto_name].SPC_time_MYA for auto_name in ordered_auto_results]
    wgd_offsets= [params_by_polyploid[allo_name].SPC_time_MYA - params_by_polyploid[allo_name].WGD_time_MYA
                  for allo_name in ordered_allo1_results] + [0]
    num_spec_times=len(spec_times)

    num_autos=len(ordered_allo1_results)
    num_wgd_times = num_autos + 1

    fig, ax = plt.subplots(1, num_wgd_times,figsize=(10,5))
    fig.suptitle(sample_name + " for " + replicate +", "+ alg + " algorithm")
    ax[0].set_title("Allopolyploid\n",fontsize=15)
    ax[1].set_title("Autopolyploid\n",fontsize=15)

    for sim_idx in range(0, num_spec_times):

        allo_group="Allo"+str(sim_idx+1)
        ordered_allo_results=[n for n in result_names if allo_group in n]
        ordered_allo_results.sort()
        results_for_this_sim_index=ordered_allo_results + [ordered_auto_results[sim_idx]]
        num_omitted_sims= num_wgd_times - len(results_for_this_sim_index)
        for i in range(0, num_omitted_sims):
            results_for_this_sim_index = [False] + results_for_this_sim_index

        for result_idx in range(0,num_wgd_times):

            allo_result_name=results_for_this_sim_index[result_idx]
            if not allo_result_name:
                continue
            this_ax = ax[result_idx]
            csvs_for_allo_result= csvfiles_by_polyploid_by_rep_by_algorthim[allo_result_name]
            ks_for_allo_result= csvs_for_allo_result[spec][replicate][alg]
            params=params_by_polyploid[allo_result_name]
            hist_data=make_time_series_histogram_subplot(this_ax, spec, ks_for_allo_result,
                                                         bin_size, params,
                                                         max_Ks_for_x_axis, max_Y,"blue")

            out_file_name = os.path.join(out_folder, allo_result_name + ".hist.csv")
            save_hist_to_csv(hist_data, out_file_name)


            this_ax.set(xlabel="<-- ks as MYA -->")

    ax[0].set(ylabel="<- # paralog pairs -> ")
    plt.tight_layout()
    out_file_name=os.path.join(out_folder, "histogram" + "_plot_" + spec +
                               "_" + replicate + "_" + str(max_Ks_for_x_axis) + ".png")
    plt.savefig(out_file_name)
    plt.close()
    return metrics_by_result_names


def make_time_series_histogram_subplot(this_ax, spec, Ks_results, bin_size, params,
                                       max_Ks, maxY, plot_color):

    x = Ks_results
    #kernel_size=1
    bins = np.arange(0, max_Ks + 0.1, bin_size)
    n, bins, patches = this_ax.hist(x, bins=bins, facecolor='none', alpha=0.25)
    trimmed_bins=bins[0:len(n)]
    #smoothed_ys = smooth_data(kernel_size, n)
    this_ax.plot(trimmed_bins,n,
                 label='SPC/WGD time: {0}/{1} MYA'.format(params.SPC_time_MYA, params.WGD_time_MYA)                 )
    hist_result=[n, bins]
    hist_maximum=max(n)
    ymax_suggestion=hist_maximum*1.6

    if maxY:
        yaxis_limit=maxY
        ymax_suggestion=False
    else:
        yaxis_limit= ymax_suggestion

    xticks=this_ax.get_xticks()
    as_my=[x/config.SpecKS_config.Ks_per_Myr for x in xticks]
    tick_labels=[str(round(my,3)) for my in as_my]
    #this_ax.set_xticks(ticks)
    this_ax.set_xticklabels(tick_labels)

    this_ax.legend()
    this_ax.set(ylim=[0, yaxis_limit])
    #this_ax.set(ylim=[0, 150])

    this_ax.set(xlim=[0, max_Ks])


    return hist_result


def save_hist_to_csv(hist_result, out_file_name):

    [n, bins]=hist_result
    with open(out_file_name, 'w') as f:
        data_headers= ['bin','# in bin']
        f.writelines(",".join(data_headers) +"\n")

        for i in range(0,len(n)):
            data_list=[str(bins[i]),str(n[i])]
            f.writelines(",".join(data_list) +"\n")

import os
import unittest

from matplotlib import pyplot as plt
import config
import kp_reader
from results_viewer import batch_analyzer, batch_histogrammer, curve_fitting


class Test1KP(unittest.TestCase):
    def test_1KP_analysis(self):
        data_directory = "/home/tamsen/Data/1KP_classfier"
        Ks_file_1 = "final_ks_values_CSUV.fa"
        kp_directory = os.path.join(data_directory, "1KP_final_ks_files")
        out_data_directory = os.path.join(data_directory, "1KP_classifier_out_Apr21")

        lookup_file = "1KP-Sample-List.csv"
        sample_lookup = kp_reader.get_sample_code_lookup(lookup_file)

        if not os.path.exists(out_data_directory):
            os.makedirs(out_data_directory)

        original_bin_size=2.0/50
        print(original_bin_size)
        bin_size = 0.02
        max_Ks = 2.0
        KS_data_files = [Ks_file_1 ] + [f for f in os.listdir(kp_directory) if ".fa" in f]
        max_to_process=5
        i=0
        curated_samples= ['OBUY']#known_auto_examples+known_allo_examples+required

        for Ks_file in KS_data_files:

            if i > max_to_process:
                break
            Ks_file_path = os.path.join(kp_directory, Ks_file)
            (species_code, species) = kp_reader.get_species_code(Ks_file, sample_lookup)
            print("classifying " + str(i) + " of " + str(max_to_process))
            ks_values = kp_reader.get_Ks_from_file(Ks_file_path)

            #get hist data
            sample_name=species
            fig, ax = plt.subplots(1, 1, figsize=(10, 10))
            fig.suptitle(sample_name)
            fig.suptitle(sample_name)
            ax.set_title("Allopolyploid\n", fontsize=20)
            ax.set_title("Autopolyploid\n", fontsize=20)

            params=config.PolyploidParams(-1, -1, sample_name)
            hist_data=batch_histogrammer.make_histogram_subplot(ax,  sample_name,ks_values,
                                        bin_size, params,
                                        max_Ks, False,"blue")
            #hist_data2=[hist_data[0],hist_data[1][0:len(hist_data[1])-1]]
            out_file_name = os.path.join(out_data_directory, "histogram" + "_plot_" +  sample_name+
                                         "_" + species_code + "_"+ str(max_Ks) + ".png")
            plt.savefig(out_file_name)
            plt.close()

            out_file_name = os.path.join(out_data_directory, sample_name + ".hist.csv")
            batch_histogrammer.save_hist_to_csv(hist_data, out_file_name)

            #analyze hist data
            hist_data2 =  batch_analyzer.read_hist_csv(out_file_name)
            ys=hist_data2[0][0:len(hist_data2[0])-1]
            ns = hist_data2[1][0:len(hist_data2[1]) - 1]
            hist_data3=[ys,ns]
            out_fit_png = out_file_name.replace(".csv", "_Ks_hist_fit" + str(max_Ks) + ".png")
            fit_results = analyze_histogram2(*hist_data3, params.WGD_time_MYA, params.SPC_time_MYA,
                                            max_Ks, False, 'tan', out_fit_png)

            #results_by_file = {"test": fit_results}
            #out_csv = "{0}_metrics.csv".format("test")
            #out_file_name = os.path.join(out_data_directory, out_csv)
            #batch_analyzer.plot_and_save_metrics(results_by_file, out_file_name)


def analyze_histogram2(bins, n, WGD_time_MYA, SPC_time_MYA,
                      max_Ks, maxY, hist_plot_color, out_file_name):

    polyploid_name=os.path.basename(out_file_name).replace(".png","")
    fig, ax = plt.subplots(1, 1, figsize=(10, 10))
    fig.suptitle('Histogram analysis for ' + polyploid_name)
    WGD_as_Ks = WGD_time_MYA * 0.01 #/ 1.04 ..the peak max is about 96% off from where it should be
    SPEC_as_Ks =SPC_time_MYA * 0.01 #/ 1.04
    #plt.bar(bins, n, width=0.001, color=hist_plot_color, label="hist",alpha=0.5)
    width=bins[2]-bins[1]

    num_gene_pairs_str=str(sum(n))
    hist_maximum=max(n)

    ax.axvline(x=WGD_as_Ks, color='b', linestyle='-', label="WGD time, " + str(WGD_time_MYA) + " MYA")
    ax.axvline(x=SPEC_as_Ks, color='r', linestyle='--', label="SPEC time, " + str(SPC_time_MYA) + " MYA")

    kernel_size=2
    smoothed_color='gray'
    smoothed_ys = batch_analyzer.smooth_data(kernel_size, n)
    plt.plot(bins, smoothed_ys,color=smoothed_color, label="smoothed data")


    smoothed_minima = batch_analyzer.find_local_minima(bins, smoothed_ys)
    smoothed_maxima = batch_analyzer.find_global_maxima(bins, smoothed_ys, 0.075)
    plt.scatter([m[0] for m in smoothed_minima], [m[1] for m in smoothed_minima], color="red", label="minima", marker='*')

    ssd_end, next_min =batch_analyzer.smallest_min(smoothed_minima,smoothed_maxima)
    ssds_xs_to_subtract, ssds_ys_to_subtract = batch_analyzer.estimate_overlap(bins, ssd_end, next_min)

    plt.scatter(ssd_end[0],3*ssd_end[1],color="black", label="wgd_start",marker='*',s=50)
    #super_smoothed_ys=smooth_data(200, n)
    ssd_xs, ssd_ys, wgd_xs, wgd_ys = batch_analyzer.sort_ssds_and_wgds(bins, n, ssd_end, ssds_ys_to_subtract)

    wgd_maxima = batch_analyzer.find_global_maxima(wgd_xs, wgd_ys, 0.075)
    plt.scatter(wgd_maxima[0], -0.05 * maxY,
                color='blue', marker='^', s=80, label="wgd max")

    wgd_max_d=batch_analyzer.get_loc_of_max_derivative(kernel_size, wgd_ys, wgd_xs)
    raw_cm, raw_x_value_of_ymax = curve_fitting.get_mode_and_cm(wgd_ys, wgd_xs)
    ax.set(xlabel="Ks")
    ax.set(ylabel="# paralogs")
    plt.bar(ssd_xs, ssd_ys, width=width, color="lightgray", label="ssd")
    plt.bar(wgd_xs, wgd_ys, width=width, color="gray", label="wgd")
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

    fit_data = batch_analyzer.run_metrics(polyploid_name[0:4], polyploid_name,
                           SPC_time_MYA, WGD_time_MYA,
                           raw_cm, raw_x_value_of_ymax,
                           wgd_maxima[0],wgd_max_d[0],
                           lognorm_goodness_of_fit, gaussian_goodness_of_fit)

    return fit_data

if __name__ == '__main__':
    unittest.main()

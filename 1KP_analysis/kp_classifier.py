import math
import os
import statistics

import seaborn
import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit

import curve_equations
import fit_curve_to_histogram
import kp_data_to_process
from fit_curve_to_histogram import find_local_minima


class KP_classifier:

    def __init__(self, max_Ks, num_bins, kernel_size,known_auto_examples,known_allo_examples):
        self.max_Ks = max_Ks
        self.num_bins = num_bins
        self.bin_width = max_Ks / num_bins
        self.kernel_size = kernel_size
        self.allo_vs_auto_threshhold=0.5
        self.ssd_vs_wgd_threshhold=0.05
        self.processed_results_by_species_code = {}
        self.known_auto_examples=known_auto_examples
        self.known_allo_examples=known_allo_examples

    def classify(self,species_data_to_classify,out_data_directory):

        species_data_to_classify.ks_to_histogram_data(self.max_Ks,
                                                      self.num_bins, self.bin_width)

        kernel = np.ones(self.kernel_size) / self.kernel_size
        smoothed_ys = np.convolve(species_data_to_classify.ys, kernel, mode='same')
        minima = find_local_minima(species_data_to_classify.xs, smoothed_ys)

        try:
            species_data_to_classify.ssd_fit_results = fit_curve_to_histogram.do_kp_SSD_fit(species_data_to_classify.xs,
                                                                                species_data_to_classify.ys,
                                                                                self.bin_width)


        except Exception as inst:
            print(type(inst))    # the exception type
            print(inst.args)     # arguments stored in .args
            print(inst)          # __str__ allows args to be printed directly
            species_data_to_classify.ssd_fit_results = [0]

        species_data_to_classify.smoothed_ys = smoothed_ys

        if len(minima)==0:
            print(species_data_to_classify.species + " has no local mimima.")
        else:
            WGD_detected = False
            WGD_detected, WGD_minima = self.find_minima(WGD_detected, minima, species_data_to_classify.ys)
            print("WGD detected? " + str(WGD_detected))

            fits_to_try=[curve_equations.wgd_lognorm,curve_equations.wgd_gaussian]

            fit_range=self.determine_fit_range(WGD_minima, self.bin_width,
                                           species_data_to_classify.xs, species_data_to_classify.ys)



            for fit_to_try in fits_to_try:

                WGD_detected, fit_result = fit_curve_to_histogram.do_kp_WGD_fit(
                    fit_to_try, fit_range, WGD_detected,
                    species_data_to_classify.xs, species_data_to_classify.ys,
                    species_data_to_classify.ssd_fit_results,
                    self.ssd_vs_wgd_threshhold, )

                species_data_to_classify.add_fit_results(WGD_detected, fit_to_try, fit_result)

            self.score_data(species_data_to_classify)

        species_data_to_classify.plot_species_results_histogram(out_data_directory)
        #species_data_to_classify.plot_species_results_SSD_histogram(out_data_directory)
        #species_data_to_classify.plot_species_results_WGD_histogram(out_data_directory)
        self.processed_results_by_species_code[species_data_to_classify.species_code] = \
            species_data_to_classify

    def score_data(self,species_data_to_classify):

        #species_data_to_classify.best_wgd_fit_result = species_data_to_classify.all_wgd_fit_results[0]
        if species_data_to_classify.gaussian_error and species_data_to_classify.lognorm_error:

            ratio = species_data_to_classify.gaussian_error / species_data_to_classify.lognorm_error
            diff = species_data_to_classify.gaussian_error-species_data_to_classify.lognorm_error
            species_data_to_classify.classification_score = \
                    round(max(min(ratio - 1,1.0),0.0),2)

            if species_data_to_classify.classification_score < self.allo_vs_auto_threshhold:
                species_data_to_classify.classification = "AUTO"
            else:
                species_data_to_classify.classification = "ALLO"

            #if #diff < self.allo_vs_auto_threshhold*species_data_to_classify.lognorm_error:
            #if (species_data_to_classify.gaussian_error < self.allo_vs_auto_threshholdspecies_data_to_classify.lognorm_error):
            #    species_data_to_classify.classification = "AUTO"

        else:
            if species_data_to_classify.gaussian_error:
                species_data_to_classify.classification_score =0.0
                species_data_to_classify.classification = "AUTO"

            if species_data_to_classify.gaussian_error:
                species_data_to_classify.classification_score =1.0
                species_data_to_classify.classification = "ALLO"

        if species_data_to_classify.classification == "ALLO" :
            species_data_to_classify.best_wgd_fit_result = species_data_to_classify.all_wgd_fit_results[curve_equations.wgd_lognorm]
        elif  species_data_to_classify.classification == "AUTO":
            species_data_to_classify.best_wgd_fit_result = species_data_to_classify.all_wgd_fit_results[curve_equations.wgd_gaussian]


        return
    def print_summary_file(self, out_data_directory):

        results_file = os.path.join(out_data_directory,
                                    "results.txt")
        lines_for_file = [kp_data_to_process.results_to_string_headers()]

        for species_code in self.processed_results_by_species_code:

            data_to_print =self.processed_results_by_species_code[species_code].results_to_string()
            lines_for_file = lines_for_file + data_to_print

        with open(results_file, "w") as f:
            f.writelines(lines_for_file)

    def plot_summary_files(self, out_data_directory):
        for i in range(0,4):
            self.plot_summary_file(i, out_data_directory)

        self.plot_fit_errors(out_data_directory)

    def plot_summary_file(self, metric_idx, out_data_directory):

        metrics=["AUC","nAUC","nAUC2","AUH_WGD_to_SSD"]
        yLabel = metrics[metric_idx]
        AUC_plot_file = os.path.join(out_data_directory,"Ks_vs_" +yLabel +"_plot.png")

        plt.close()
        plt.clf()
        for species_code in self.processed_results_by_species_code:
            species_data_results = self.processed_results_by_species_code[species_code ]
            if (species_data_results.classification != "NO WGD") :

                best_fit = species_data_results.best_wgd_fit_result
                Ks = best_fit.xForKsPeak

                if metric_idx==0:
                    data_to_report = best_fit.AUC
                if metric_idx==1:
                    data_to_report = best_fit.nAUC
                if metric_idx==2:
                    data_to_report = best_fit.nAUC
                if metric_idx==3:
                    data_to_report = best_fit.AUH_WGD_to_SSD

                color = (species_data_results.classification_score, 0.5, 0.5, 1.0)
                size = species_data_results.classification_score*50+5

                if species_code in self.known_allo_examples:
                    plt.scatter(Ks, data_to_report, c="red", s=size+40)

                if species_code in self.known_auto_examples:
                    plt.scatter(Ks, data_to_report, c="blue", s=size+40)

                plt.scatter(Ks, data_to_report, c=color, s=size)
                label=(species_code+ "\n" +  species_data_results.classification +", " +
                       str(species_data_results.classification_score))



                #plt.text(Ks, data_to_report, label, c=(0.8, 0.8, 0.8, 0.8))

        plt.ylabel(yLabel)
        plt.xlabel("Ks")
        plt.title("Ks vs " + yLabel)
        plt.savefig(AUC_plot_file)
        plt.close()
        plt.clf()

    def plot_fit_errors(self, out_data_directory):

        yLabel = "error (rms Yi-Fi)"

        plot_file_1 = os.path.join(out_data_directory, "Error_box_plot.png")
        plot_file_2 = os.path.join(out_data_directory, "Error_violin.png")
        plot_file_3 = os.path.join(out_data_directory, "Error_scatter.png")
        plot_file_4 = os.path.join(out_data_directory, "Error_hist.png")

        plt.close()
        plt.clf()
        errors_by_fit_fxn={}
        for species_code in self.processed_results_by_species_code:
            species_data_results = self.processed_results_by_species_code[species_code]
            if (species_data_results.classification != "NO WGD"):

                fit_fxns = list(species_data_results.all_wgd_fit_results.keys())

                for i in range(0, len(fit_fxns)):
                    fit_fxn=fit_fxns[i]
                    if fit_fxn not in errors_by_fit_fxn:
                        errors_by_fit_fxn[fit_fxn] = []
                    fit_result = species_data_results.all_wgd_fit_results[fit_fxn]
                    errors_by_fit_fxn[fit_fxn].append(fit_result.fit_error)

                    #error_to_print = str(round(fit_result.fit_error, 2))
                    #label = str(fit_result.fit_fxn_name) + ", error: " + error_to_print

        groups=errors_by_fit_fxn.keys()

        nice_labels=[str(g).split()[1] for g in groups]
        data_to_plot = [errors_by_fit_fxn[fit_fxn] for fit_fxn in groups]
        num_datapoints_per_group = len(data_to_plot[0])
        diffs = [data_to_plot[1][r]-data_to_plot[0][r] for r in range(0,num_datapoints_per_group)]
        data_to_plot.append(diffs)
        nice_labels.append("diffs")

        means = [statistics.mean(data) for data in data_to_plot]
        print("Error means:" + str(means))
        ticks=[1, 2, 3]
        plt.boxplot(data_to_plot)
        plt.xticks(ticks, nice_labels)

        for i in range(0,len(ticks)):
            plt.text(ticks[i], -1.5, str(round(means[i],2)), c="k")

        plt.ylabel(yLabel)
        plt.xlabel("fit functions")
        plt.title("Errors by Fit Function")
        plt.savefig(plot_file_1)
        plt.close()
        plt.clf()

        #stats.kstest(stats.uniform.rvs(size=100, random_state=rng),
        #             stats.norm.cdf)

        seaborn.violinplot(data_to_plot)
        plt.xticks([0, 1, 2], nice_labels)
        plt.ylabel(yLabel)
        plt.xlabel("fit functions")
        plt.title("Errors by Fit Function")
        plt.savefig(plot_file_2)
        plt.close()
        plt.clf()

        noise_level=0.1
        plt.scatter(np.random.normal(0,noise_level,len(data_to_plot[0])),data_to_plot[0])
        plt.scatter(np.random.normal(1,noise_level,len(data_to_plot[1])),data_to_plot[1])
        plt.scatter(np.random.normal(2,noise_level,len(data_to_plot[2])),data_to_plot[2])
        plt.xticks([0, 1, 2], nice_labels)
        plt.ylabel(yLabel)
        plt.xlabel("fit functions")
        plt.title("Errors by Fit Function")
        plt.savefig(plot_file_3)
        plt.close()
        plt.clf()

        plt.hist(diffs, bins=50,color="gray", ec="black")
        plt.ylabel("counts")
        plt.xlabel("fit functions")
        plt.title("Distribution of improvements using Lognorm")
        plt.savefig(plot_file_4)
        plt.close()
        plt.clf()

    def find_minima(self, WGD_detected, minima, ys):
        biggest_y = 0
        best_interval = 0
        local_minima = []
        if len(minima) > 1:

            for i in range(0, len(minima) - 1):

                start_idx = minima[i][2]
                stop_idx = minima[i + 1][2]
                ys_in_range = [ys[j] for j in range(start_idx, stop_idx)]
                y_max = max(ys_in_range)
                if y_max > biggest_y:
                    biggest_y = y_max
                    best_interval = i

            local_minima = [minima[best_interval][0], minima[best_interval + 1][0]]
            WGD_detected = True
        return WGD_detected, local_minima

    def determine_fit_range(self, local_minima,interval, xs, ys):

        wgd_range = local_minima  # [0.05, 0.55]#
        xs_for_fit = []
        ys_for_fit = []
        is_for_fit = []
        extended_range = [wgd_range[0] - interval, wgd_range[1] + interval]
        for i in range(0, len(xs)):
            x = xs[i]
            y = ys[i]
            if (x <= (extended_range[1]) and (x > extended_range[0])):
                xs_for_fit = xs_for_fit + [x]
                ys_for_fit = ys_for_fit + [y]
                is_for_fit = is_for_fit + [i]
        return xs_for_fit, ys_for_fit, is_for_fit

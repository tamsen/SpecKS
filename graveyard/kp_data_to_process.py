import math
import os

import numpy as np
from matplotlib import pyplot as plt

import curve_equations


class data_to_process:

    def __init__(self, species, species_code, ks_values, max_Ks):

        self.max_Ks = max_Ks
        self.species= species
        self.species_code = species_code
        self.ks_values = ks_values
        self.bin_width = -1
        self.num_bins = -1
        self.xs = []#xs
        self.ys = []#ys
        self.smoothed_ys = []#ys
        self.all_wgd_fit_results = {}
        self.best_wgd_fit_result = False
        self.ssd_fit_results = []
        self.gaussian_error=False
        self.lognorm_error=False
        self.classification_score=False
        self.classification="NO WGD"

    def add_fit_results(self, WGD_detected, fit_fxn, fit_result):

        self.all_wgd_fit_results[fit_fxn]=fit_result

        if WGD_detected:
            if (fit_fxn == curve_equations.wgd_gaussian):
                self.gaussian_error = fit_result.fit_error
            if (fit_fxn == curve_equations.wgd_lognorm):
                self.lognorm_error = fit_result.fit_error

    def ks_to_histogram_data(self, max_Ks, num_bins,step_size):

        bin_sides = np.arange(0,max_Ks + step_size, step_size)
        bins = [0 for i in range(0, num_bins)]
        for ks in self.ks_values:

            if ks >= max_Ks:
                continue
            ith_bin = math.floor(ks / step_size)
            bins[ith_bin] = bins[ith_bin] + 1.0

        self.xs = bin_sides[0:num_bins]
        self.ys = bins
        self.bin_width = step_size
        self.num_bins = num_bins

    def plot_species_results_histogram(self, out_data_directory):

        plot_file = os.path.join(out_data_directory,
                                 self.species_code + "_" + self.species + ".png")

        plt.bar(self.xs, self.ys,
                width=self.bin_width, color='gray', alpha=0.5,label='Ks histogram')

        colors=["pink","red"]
        fit_fxns=list(self.all_wgd_fit_results.keys())

        for i in range(0, len(fit_fxns)):

            fit_result =self.all_wgd_fit_results[fit_fxns[i]]
            error_to_print = str(round(fit_result.fit_error, 2))
            label = str(fit_result.fit_fxn_name) + ", error: " + error_to_print
            plt.plot(self.xs, fit_result.fit_ys, color=colors[i],linestyle = '-',
                     label=label) #marker=":")

        #plt.plot(self.xs, self.ssd_fit_results, color="yellow",linestyle = "-",
        #         label="ssd fit") #marker=":")

        #plt.plot(self.xs[2:-1], self.smoothed_ys[2:-1], color="gray",linestyle = "-",
        #         label="smoothed fxn", linewidth=0.5)

        score_to_print= str(round(100.0*self.classification_score,0))
        plt.legend()

        #plt.ylim(0, 550)
        plt.xlim(-.1, 4.0)
        if self.best_wgd_fit_result:
            max_y_at_ks_peak=max(self.best_wgd_fit_result.fit_ys)
            plt.ylim(0,max_y_at_ks_peak*2.0)

        plt.xlabel("Ks")
        plt.ylabel("Count in Bin")
        plt.title("Ks histogram for " + self.species +
                  ", last ~" + str(self.max_Ks * 100) + " MY\n" +
                  self.classification + ",  "
                  "Allopolyploid %" + score_to_print)
        plt.savefig(plot_file)
        plt.close()
        plt.clf()

    def plot_species_results_WGD_histogram(self, out_data_directory):

        plot_file = os.path.join(out_data_directory,
                                 self.species_code + "_" + self.species + "_WGD.png")

        colors = ["red", "black"]
        fit_fxns = list(self.all_wgd_fit_results.keys())
        #for i in range(0, len(fit_fxns)):
        for i in [1]:
            fit_result = self.all_wgd_fit_results[fit_fxns[i]]
            error_to_print = str(round(fit_result.fit_error, 2))
            label = str(fit_result.fit_fxn_name) + ", error: " + error_to_print
            plt.plot(self.xs, fit_result.fit_ys, color=colors[i], linestyle='-',
                     label=label)  # marker=":")\

        bar_under_ssd = [max(0,(self.ys[i]- self.ssd_fit_results[i])) for i in range(0, len(self.ys))]
        plt.bar(self.xs, bar_under_ssd,
                width=self.bin_width, color='gray', alpha=0.5, label='Ks histogram')
        score_to_print = str(round(100.0 * self.classification_score, 0))
        # plt.legend()

        if self.best_wgd_fit_result:
            max_y_at_ks_peak = max(self.best_wgd_fit_result.fit_ys)
            plt.ylim(0, max_y_at_ks_peak * 2.0)

        plt.ylim(0, 550)
        plt.xlim(-.1, 4.0)
        plt.xlabel("Ks")
        plt.ylabel("Count in Bin")
        plt.title("Ks histogram for " + self.species +
                  ", last ~" + str(self.max_Ks * 100) + " MY\n" +
                  self.classification + ",  "
                                        "Allopolyploid %" + score_to_print)
        plt.savefig(plot_file)
        plt.close()
        plt.clf()

    def plot_species_results_SSD_histogram(self, out_data_directory):

        plot_file = os.path.join(out_data_directory,
                                 self.species_code + "_" + self.species + "_SSD.png")

        #plt.bar(self.xs, self.ys,
        #        width=self.bin_width, color='gray', alpha=0.5,label='Ks histogram')

        plt.plot(self.xs, self.ssd_fit_results, color="yellow",linestyle = "-",
                 label="ssd fit") #marker=":")

        bar_under_ssd=[min(self.ys[i], self.ssd_fit_results[i]) for i in range(0,len(self.ys)) ]
        plt.bar(self.xs, bar_under_ssd,
                width=self.bin_width, color='gray', alpha=0.5,label='Ks histogram')
        score_to_print= str(round(100.0*self.classification_score,0))
        #plt.legend()

        if self.best_wgd_fit_result:
            max_y_at_ks_peak=max(self.best_wgd_fit_result.fit_ys)
            plt.ylim(0,max_y_at_ks_peak*2.0)

        plt.ylim(0, 550)
        plt.xlim(-.1, 4.0)

        plt.xlabel("Ks")
        plt.ylabel("Count in Bin")
        plt.title("Ks histogram for " + self.species +
                  ", last ~" + str(self.max_Ks * 100) + " MY\n" +
                  self.classification + ",  "
                  "Allopolyploid %" + score_to_print)
        plt.savefig(plot_file)
        plt.close()
        plt.clf()


    def results_to_string(self):

        results = []

        for fit_result in self.all_wgd_fit_results.values():

            data_list = [self.species_code,self.species,
                         fit_result.fit_fxn_name,
                         fit_result.xForKsPeak,
                         fit_result.fit_error,
                         fit_result.AUC,
                         fit_result.nAUC,
                         self.classification_score,
                         self.classification ]
            data_string="\t".join([str(d) for d in data_list]) + "\n"
            results.append(data_string)

        return results

def results_to_string_headers():
    #return "code\tspecies\tfit_fxn\tKs\terr\tAUC\tnAUC\tsAUC\tsnAUC\tclassification\tWGD?\n"
    return "code\tspecies\tfit_fxn\tKs\terr\tAUC\tnAUC\tclassification\tWGD?\n"
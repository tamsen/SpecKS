import os
import unittest

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

from fitting_to_known_allos_and_autos.fit_to_ks_rates_data import parse_ks
from results_viewer import batch_analyzer, batch_histogrammer, curve_fitting, batch_aggregator

class MyTestCase(unittest.TestCase):
    def test_Single_histogram(self):

        hist_comparison_out_folder = "/home/tamsen/Data/SpecKS_output/hist_comparison"
        specks_out_folder="/home/tamsen/Data/SpecKS_output/" + \
                    "SpecKS_m04d25y2024_h17m47s14/Allo_Maize/8_final_results"
        ksrates_out_folder = "/home/tamsen/Data/SpecKS_input/ks_data"
        specks_csv_file="Allo_Maize_ML_rep0_LCA_to_Ortholog_Ks_by_GeneTree.csv"
        ksrates_csv_file="mays.ks.tsv"

        specks_full_path=os.path.join(specks_out_folder,specks_csv_file)
        specks_ks_results = batch_histogrammer.read_Ks_csv(specks_full_path)

        real_full_path=os.path.join(ksrates_out_folder,ksrates_csv_file)
        real_ks_results = parse_external_ksfile(real_full_path)

        bin_size=0.002
        max_Ks=0.5
        color='blue'

        if not os.path.exists(hist_comparison_out_folder):
            os.makedirs(hist_comparison_out_folder)

        out_png = os.path.join(hist_comparison_out_folder, "specks_out.png")
        make_simple_histogram(specks_ks_results, bin_size, color, max_Ks, out_png)

        out_png = os.path.join(hist_comparison_out_folder, "real_out.png")
        make_simple_histogram(real_ks_results, bin_size, color, max_Ks, out_png)

def make_simple_histogram(Ks_results, bin_size, color, max_Ks, out_png):


    fig = plt.figure(figsize=(10, 10), dpi=100)
    x = Ks_results
    # print(PAML_hist_out_file)
    if max_Ks:
        bins = np.arange(bin_size, max_Ks + 0.1, bin_size)
        n, bins, patches = plt.hist(x, bins=bins, facecolor=color, alpha=0.25, label='histogram data')
        plt.xlim([0, max_Ks * (1.1)])
    # plt.ylim([0, max_y])
    # plt.axvline(x=WGD_as_Ks, color='b', linestyle='-', label="WGD time as Ks")
    # plt.axvline(x=SPEC_as_Ks, color='r', linestyle='--', label="SPEC time as Ks")
    plt.legend()
    plt.xlabel("Ks")
    plt.ylabel("Count in Bin")
    # plt.title("Ks histogram for " + species_name + ", last ~" + str(max_Ks * 100) + " MY\n" +
    #          "algorithm: PAML " + alg_name)
    plt.savefig(out_png)
    plt.clf()
    plt.close()


def parse_external_ksfile(ks_file):

    if ".fa" in ks_file: #1KP file
        olea_ks_df = pd.read_csv(ks_file, sep='\t', header=0)
        olea_ks_array = olea_ks_df.loc[:, "Node Ks"]
        ks_data = [k for k in olea_ks_array.tolist() if k <= 2]
    else:    #KS rates input file
        ks_data = parse_ks(ks_file)
    return ks_data

if __name__ == '__main__':
    unittest.main()

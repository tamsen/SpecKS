import os
import unittest

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from fitting_to_known_allos_and_autos.fit_to_ks_rates_data import parse_ks, fit_fxns_to_Ks
from results_viewer import batch_analyzer, batch_histogrammer, curve_fitting, batch_aggregator

class MyTestCase(unittest.TestCase):

    def test_coffee_histogram(self):

        hist_comparison_out_folder = "/home/tamsen/Data/SpecKS_output/hist_comparison"
        ksrates_out_folder = "/home/tamsen/Data/SpecKS_input/ks_data"

        specks_out_folder="/home/tamsen/Data/Specks_outout_from_mesx/sim41_coffee"
        specks_csv_file = "Allo_Coffea5_ML_rep0_LCA_to_Ortholog_Ks_by_GeneTree.csv"
        ksrates_csv_file="coffea.ks.tsv"

        splat=specks_csv_file.split("_")
        species_run_name=splat[0]+splat[1]
        specks_full_path=os.path.join(specks_out_folder,specks_csv_file)
        real_full_path=os.path.join(ksrates_out_folder,ksrates_csv_file)

        bin_size=0.002
        max_Ks=0.2
        color='blue'
        wgd_ks=0.002
        density = 1

        make_both_histograms(bin_size, color, hist_comparison_out_folder,
                             wgd_ks, max_Ks, density, real_full_path,
                         species_run_name, specks_full_path)

    def test_poplar_histogram(self):
        hist_comparison_out_folder = "/home/tamsen/Data/SpecKS_output/hist_comparison"
        ksrates_out_folder = "/home/tamsen/Data/SpecKS_input/ks_data"

        specks_out_folder = "/home/tamsen/Data/Specks_outout_from_mesx/sim41_poplar"
        specks_csv_file = "Allo_Poplar4_ML_rep0_LCA_to_Ortholog_Ks_by_GeneTree.csv"
        ksrates_csv_file = "poplar.ks.tsv"

        splat = specks_csv_file.split("_")
        species_run_name = splat[0] + splat[1]
        specks_full_path = os.path.join(specks_out_folder, specks_csv_file)
        real_full_path = os.path.join(ksrates_out_folder, ksrates_csv_file)

        bin_size = 0.002
        max_Ks = 0.5
        color = 'blue'
        wgd_ks=0.18
        density = 1

        make_both_histograms(bin_size, color, hist_comparison_out_folder, wgd_ks,
                             max_Ks, density, real_full_path,
                             species_run_name, specks_full_path)
    def test_maize_histogram(self):

        hist_comparison_out_folder = "/home/tamsen/Data/SpecKS_output/hist_comparison"
        ksrates_out_folder = "/home/tamsen/Data/SpecKS_input/ks_data"

        specks_out_folder="/home/tamsen/Data/Specks_outout_from_mesx/sim41_maize"
        #specks_csv_file = "Allo2_Maize_ML_rep0_LCA_to_Ortholog_Ks_by_GeneTree.csv"
        specks_csv_file = "Allo5_Maize_ML_rep0_LCA_to_Ortholog_Ks_by_GeneTree.csv"
        ksrates_csv_file="mays.ks.tsv"

        splat=specks_csv_file.split("_")
        species_run_name=splat[0]+splat[1]
        specks_full_path=os.path.join(specks_out_folder,specks_csv_file)
        real_full_path=os.path.join(ksrates_out_folder,ksrates_csv_file)

        bin_size=0.002
        max_Ks=0.3
        color='blue'
        wgd_ks=0.125
        density = 1

        make_both_histograms(bin_size, color, hist_comparison_out_folder, wgd_ks,
                         max_Ks, density, real_full_path,
                         species_run_name, specks_full_path)
    def test_olive_histogram(self):

        hist_comparison_out_folder = "/home/tamsen/Data/SpecKS_output/hist_comparison"
        ksrates_out_folder = "/home/tamsen/Data/SpecKS_input/ks_data"

        specks_out_folder="/home/tamsen/Data/SpecKS_output/" + \
                    "SpecKS_m04d26y2024_h13m45s40/Auto_Olive/8_final_results"

        specks_csv_file = "Auto_Olive_ML_rep0_LCA_to_Ortholog_Ks_by_GeneTree.csv"
        ksrates_csv_file="final_ks_values_TORX.fa"

        splat=specks_csv_file.split("_")
        species_run_name=splat[0]+splat[1]
        specks_full_path=os.path.join(specks_out_folder,specks_csv_file)
        real_full_path=os.path.join(ksrates_out_folder,ksrates_csv_file)

        bin_size=0.002
        max_Ks=0.5
        color='blue'

        make_both_histograms(bin_size, color, hist_comparison_out_folder, max_Ks, real_full_path,
                         species_run_name, specks_full_path)
    def test_Single_histogram(self):

        hist_comparison_out_folder = "/home/tamsen/Data/SpecKS_output/hist_comparison"
        specks_out_folder="/home/tamsen/Data/SpecKS_output/" + \
                    "SpecKS_m04d25y2024_h17m47s14/Allo_Maize/8_final_results"

        specks_out_folder="/home/tamsen/Data/SpecKS_output/" + \
                    "SpecKS_m04d26y2024_h12m30s41/Auto_Poplar/8_final_results"
        #specks_out_folder="/home/tamsen/Data/Specks_outout_from_mesx/sim41_maize"

        specks_out_folder="/home/tamsen/Data/SpecKS_output/" + \
                    "SpecKS_m04d26y2024_h13m45s40/Auto_Olive/8_final_results"

        ksrates_out_folder = "/home/tamsen/Data/SpecKS_input/ks_data"
        #specks_csv_file="Allo_Maize_ML_rep0_LCA_to_Ortholog_Ks_by_GeneTree.csv"
        #specks_csv_file = "Auto_Poplar_ML_rep0_LCA_to_Ortholog_Ks_by_GeneTree.csv"
        specks_csv_file = "Auto_Poplar_ML_rep0_LCA_to_Ortholog_Ks_by_GeneTree.csv"
        ksrates_csv_file="mays.ks.tsv"
        ksrates_csv_file="poplar.ks.tsv"

        splat=specks_csv_file.split("_")
        species_run_name=splat[0]+splat[1]
        specks_full_path=os.path.join(specks_out_folder,specks_csv_file)

        real_full_path=os.path.join(ksrates_out_folder,ksrates_csv_file)

        specks_ks_results = batch_histogrammer.read_Ks_csv(specks_full_path)
        real_ks_results = parse_external_ksfile(real_full_path)

        bin_size=0.002
        max_Ks=0.3
        WGD_ks=0.001
        color='blue'
        density = None

        make_both_histograms(bin_size, color, hist_comparison_out_folder, WGD_ks,max_Ks, density,real_ks_results,
                                  species_run_name, specks_ks_results)

def make_both_histograms(bin_size, color, hist_comparison_out_folder, WGD_ks, max_Ks, density, real_full_path,
                         species_run_name, specks_full_path):

    specks_ks_results = batch_histogrammer.read_Ks_csv(specks_full_path)
    real_ks_results = parse_external_ksfile(real_full_path)

    if not os.path.exists(hist_comparison_out_folder):
        os.makedirs(hist_comparison_out_folder)
    out_png1 = os.path.join(hist_comparison_out_folder, "specks_" + species_run_name + "_out.png")
    out_png2 = os.path.join(hist_comparison_out_folder, "specks_" + species_run_name + "_fit.png")
    specks_hist_data =make_simple_histogram(specks_ks_results, species_run_name, bin_size, color, WGD_ks,
                          max_Ks, density, out_png1)
    #fit_fxns_to_Ks(specks_ks_results, species_run_name,5000,400,
    #               max_Ks, out_png2)

    out_png1 = os.path.join(hist_comparison_out_folder, "real_" + species_run_name + "_out.png")
    out_png2 = os.path.join(hist_comparison_out_folder, "real_" + species_run_name + "_fit.png")
    real_hist_data = make_simple_histogram(real_ks_results, species_run_name, bin_size, color, WGD_ks,
                          max_Ks, density, out_png1)

    out_png3 = os.path.join(hist_comparison_out_folder, "overlay_" + species_run_name + "_fit.png")
    overlay_histogram(species_run_name, bin_size, color,
                      [specks_hist_data,real_hist_data], WGD_ks, out_png3)

def make_simple_histogram(Ks_results, species_name, bin_size, color,WGD_ks, max_Ks, density, out_png):

    fig = plt.figure(figsize=(10, 10), dpi=100)
    x = Ks_results
    # print(PAML_hist_out_file)
    label="hist for " + os.path.basename(out_png).replace("_out.png","")
    if max_Ks:
        bins = np.arange(bin_size, max_Ks + 0.1, bin_size)
        n, bins, patches = plt.hist(x, bins=bins, facecolor=color, alpha=0.25,
                                    label=label, density=density)
        plt.xlim([0, max_Ks * (1.1)])


    plt.axvline(x=WGD_ks, color='b', linestyle='-', label="WGD paralog start")
    num_pairs=sum(n)
    num_after_wgd=sum([n[i] for i in range(0,len(n)) if bins[i] > WGD_ks])
    plt.legend()
    plt.xlabel("Ks")
    plt.ylabel("Count in Bin")
    plt.title("Ks histogram for {0}.\n{1} pairs of genes. ~{2} retained from WGD.".format(
        species_name,num_pairs,num_after_wgd))
    plt.savefig(out_png)
    plt.clf()
    plt.close()

    return [n,bins]


def overlay_histogram(species_name, bin_size, color, list_of_hist_data, WGD_ks,out_png):
    fig = plt.figure(figsize=(10, 10), dpi=100)
    colors = ['blue','green']
    labels = ['specks','truth']

    for i in range(0,len(list_of_hist_data)):
        hist_data =list_of_hist_data[i]
        [n, bins]=hist_data
        width=(bins[2]-bins[1])/2
        xs=[b + i*width for b in bins[0:len(bins)-1]]
        plt.bar(xs,n,width=width,
                color=colors[i], alpha=0.5, label=labels[i])
        #plt.bar(bins[0:len(bins)-1],n,width=width,
        #        color=colors[i], alpha=0.5, label=labels[i])

    plt.axvline(x=WGD_ks, color='b', linestyle='-', label="WGD paralog start")
    num_pairs = sum(n)
    num_after_wgd = sum([n[i] for i in range(0, len(n)) if bins[i] > WGD_ks])
    plt.legend()
    plt.xlabel("Ks")
    plt.ylabel("Count in Bin")
    plt.title("Ks histogram for {0}.\n{1} pairs of genes. ~{2} retained from WGD.".format(
        species_name, num_pairs, num_after_wgd))
    plt.savefig(out_png)
    plt.clf()
    plt.close()

    return n, bins

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

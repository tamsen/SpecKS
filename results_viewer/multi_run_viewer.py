import os
import unittest

import numpy as np
from matplotlib import pyplot as plt

#https://stackoverflow.com/questions/14770735/how-do-i-change-the-figure-size-with-subplots
class MulitRunViewerTests(unittest.TestCase):
    def test_multi_run_viewer(self):

        #suppose you have lots of results (cvs files) with all the KS results from many specks runs,
        #and you want to see them all together on one plot.

        output_folder="/home/tamsen/Data/SpecKS_mesx_data/mesx_sim1_no_genebirth_or_death"
        csv_file=os.path.join(output_folder,"ML_rep0_Ks_by_GeneTree.csv")
        ks_results=read_Ks_csv(csv_file)
        ks_results_by_scenario={}
        species="Allo1_S050W025"
        ks_results_by_scenario["Sim1"]=ks_results
        ks_results_by_scenario["Sim2"]=ks_results
        ks_results_by_scenario["Sim3"]=ks_results

        plot_traces_for_the_sample(output_folder, 'sample_name',
                                    ks_results_by_scenario)

        self.assertEqual(True,True)


def read_Ks_csv(csv_file):

    ks_results = []
    with open(csv_file, "r") as f:

        reading_header=True
        while True:
            line = f.readline()
            if not line:
                break
            if len(line)==0:
                break
            if reading_header:
                reading_header=False
                continue
            data = line.split(",")
            #print(data)
            ks_value=float(data[2])
            ks_results.append(ks_value)

    return ks_results

def plot_traces_for_the_sample(run_folder, sample_name, ks_results_by_scenario):

    #fig = plt.figure(figsize=(10, 10))
    Ks_results=ks_results_by_scenario["Sim1"]
    sim_labels=["Sim0","Sim1","Sim2","Sim3"]
    max_Ks = 2
    bin_size = 0.01
    #f, a = plt.subplots(4, 2)

    # making subplots
    fig, ax = plt.subplots(4, 2,figsize=(10, 10))
    ax[0, 0].set_title("Allopolyploid\n\n",fontsize=20)
    ax[0, 1].set_title("Autopolyploid\n\n",fontsize=20)

    num_sims=len(sim_labels)
    for sim_idx in range(0, num_sims):

        this_ax = ax[sim_idx, 0]
        make_subplot(this_ax , Ks_results, bin_size,
            "some text",max_Ks ,  "blue")

        this_ax.set(ylabel=sim_labels[sim_idx])

        text_string="SPC: {0} MYA\nWGD: {1} MYA".format(50,25)
        this_ax.text(0.8, 0.8, text_string,
                     horizontalalignment='center', verticalalignment='center',
                     transform=this_ax.transAxes)

    plot_title = sample_name
    plt.savefig(run_folder + "/" + sample_name + "_histogram" + "_plot" + ".png")

    plt.close()


def make_subplot(this_ax, Ks_results, bin_size, text_to_write, max_Ks, plot_color):

    bins = np.arange(0, max_Ks + 0.1, bin_size)
    #nBins=50
    x = Ks_results
    n, bins, patches = this_ax.hist(x, bins=bins, facecolor=plot_color, alpha=0.25, label='histogram data')

    #this_ax= plt.gca()
    #this_ax.text(0.05, 0.95, text_to_write, horizontalalignment='left',
    #              verticalalignment='top', transform=ax.transAxes)

    plt.xlim([0, max_Ks * (1.1)])

    return this_ax
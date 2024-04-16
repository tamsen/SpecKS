import os
import unittest
from sklearn import linear_model
from matplotlib import pyplot as plt

from results_viewer import curve_fitting


class AlloAutoPredictor(unittest.TestCase):
    def test_allo_auto_predictor(self):

        out_folder = "/home/tamsen/Data/Specks_outout_from_mesx"
        file1="mode_vs_spec_time.csv"
        #genes_shed_vs_wgd_time.csv

        full_path=os.path.join(out_folder,file1)
        xs,ys = read_xs_ys_csv(full_path)

        fig, ax = plt.subplots(1, 1, figsize=(5, 5))
        plt.scatter(xs, ys, c='b')
        fit_ys, xs, goodness_of_fit=curve_fitting.fit_curve_to_xs_and_ys(
            xs, ys, curve_fitting.linear)

        plt.plot(xs,fit_ys, c='k')

        plot_file=full_path.replace(".csv",".new.png")
        plt.savefig(plot_file)
        plt.close()



if __name__ == '__main__':
    unittest.main()

def read_xs_ys_csv(csv_file):

    data_by_sim_name={}
    xs=[]
    ys=[]
    with open(csv_file, "r") as f:

        while True:
            line = f.readline()
            if "batch" in line:
                continue
            if len(line)==0:
                break
            if "NA" in line:
                continue
            data = line.strip().split(",")
            xs.append(float(data[2]))
            ys.append(float(data[3]))

    return xs,ys
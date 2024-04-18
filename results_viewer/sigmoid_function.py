import os
import unittest

import numpy as np
from matplotlib import pyplot as plt

from results_viewer import curve_fitting


class TestSigmoid(unittest.TestCase):

    def test_sigmoid(self):

        out_folder = "/home/tamsen/Data/Specks_outout_from_mesx"
        xs = np.arange(-8, 8, 0.001)

        slopes =[1,2,5]
        fig, ax = plt.subplots(1, 1, figsize=(5, 5))
        colors=['r','b', 'c']

        for i in range(0,len(slopes)):
            c=colors[i]
            a=slopes[i]
            ys = [curve_fitting.sigmoid(x,a) for x in xs]
            dx=xs[1]-xs[0]
            dsdx=[ (1/dx)*(curve_fitting.sigmoid(xs[i+1],a)- curve_fitting.sigmoid(xs[i],a))
                   for i in range(0,len(xs)-1)]
            max_dsdx=round(max(dsdx),3)

            plt.plot(xs, ys, c=colors[i], label="sigmoid {0}".format(a),linestyle="-")
            plt.plot(xs[0:len(xs)-1], dsdx, alpha=0.25,c=colors[i],linestyle=":",
                     label="ds({0})/dx (max: {1})".format(a,max_dsdx))

        ax.set(yscale='log')
        ax.set(xlabel="x")
        ax.set(ylabel="y")
        plt.legend()
        plot_file = os.path.join(out_folder, "sigmoid.png")
        plt.savefig(plot_file)
        plt.close()

    def test_coalescent_t(self):

        out_folder = "/home/tamsen/Data/Specks_outout_from_mesx"
        xs = np.arange(-8, 8, 0.001)

        Ns =[1,2,5]
        fig, ax = plt.subplots(1, 1, figsize=(5, 5))
        colors=['r','b', 'c']

        for i in range(0,len(Ns)):
            c=colors[i]
            N=Ns[i]

            ys = [curve_fitting.coal(x,N) for x in xs]
            dx=xs[1]-xs[0]
            dsdx=[ (1/dx)*(curve_fitting.coal(xs[i+1],N)- curve_fitting.coal(xs[i],N))
                   for i in range(0,len(xs)-1)]
            max_dsdx=round(max(dsdx),3)

            plt.plot(xs, ys, c=colors[i], label="coal {0}".format(N),linestyle="-")
            plt.plot(xs[0:len(xs)-1], dsdx, alpha=0.25,c=colors[i],linestyle=":",
                     label="dc({0})/dx (max: {1})".format(N,max_dsdx))

        ax.set(yscale='log')
        ax.set(xlabel="x")
        ax.set(ylabel="y")
        plt.legend()
        plot_file = os.path.join(out_folder, "coal.png")
        plt.savefig(plot_file)
        plt.close()

    def test_coaloid(self):

        out_folder = "/home/tamsen/Data/Specks_outout_from_mesx"
        xs = np.arange(-8, 16, 0.001)
        #xs = np.arange(0, 8, 0.001)
        do_log_plot=False
        Ns =[10]
        slopes =[1,2,3]
        #gray_zone_length=[16,4,2]
        fig, ax = plt.subplots(1, 1, figsize=(5, 5))
        colors=['r','b', 'c']

        for i in range(0,len(Ns)):

            N=Ns[i]
            for i in range(0, len(slopes)):
                a = slopes[i]
                c = 5
                gray_area_time = 2 * c / a
                #ys1 = [curve_fitting.sigmoid(x, a, c)*curve_fitting.coal(x,N)  for x in xs]
                ys2 = [curve_fitting.sigmoid(x, a, 0) for x in xs]
                ys2b = [curve_fitting.sigmoid(x, a, c) for x in xs]
                #ys3 = [curve_fitting.coal(x,N)  for x in xs]


                #plt.plot(xs, ys1, c=colors[0], label="coaloid N={0},a={1}".format(N,a),linestyle="-")
                #plt.plot(xs, ys2, c=colors[i], label="c0 sigmoid {0}".format(a),linestyle=":")
                plt.plot(xs, ys2b, c=colors[i], label="c1 sigmoid {0}".format(gray_area_time),linestyle="-")
                plt.plot([i*0.001 for i in range(0,1000)],
                    [-1 for i in range(0,1000)], c=colors[i],
                label="gray area",linestyle="-")
                #plt.plot(xs, ys3, c=colors[2], label="coales {0}".format(N),linestyle="-")
                #plt.plot(xs[0:len(xs)-1], dsdx, alpha=0.25,c=colors[i],linestyle=":",
                #         label="dc({0})/dx (max: {1})".format(N,max_dsdx))

        plot_name="coaloid.png"
        if do_log_plot:
            ax.set(yscale='log')
            plot_name=plot_name.replace(".png","log.png")
        ax.set(xlabel="x")
        ax.set(ylabel="y")
        plt.legend()
        plot_file = os.path.join(out_folder,plot_name)
        plt.savefig(plot_file)
        plt.close()
if __name__ == '__main__':
    unittest.main()

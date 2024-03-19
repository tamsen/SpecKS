import math
import os
import statistics
import unittest
from scipy.stats import norm,lognorm
from matplotlib import pyplot as plt
from scipy import stats

from results_viewer import curve_fitting
import numpy as np


class LognormTestCase(unittest.TestCase):

    #https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.lognorm.html
    def test_a_lognorm(self):
        test_out = "test_out"
        if not os.path.exists(test_out):
            os.makedirs(test_out)

        max_X = 50.0
        bin_size = 0.01
        xs = np.arange(0, max_X + 0.01, bin_size)
        fig, ax = plt.subplots(1, 1, figsize=(10, 10))

        title='lognorm'
        fig.suptitle(title)

        #for ~0.5 MY spread
        #90% < 0.55
        shape_parameter=0.8
        xscale = 0.2

        #for 10 MY spread
        #90% < 10 MY
        shape_parameter=0.5
        xscale = 5.27
        #cdfs = [lognorm.cdf(xscale*x,shape_parameter) for x in xs]
        ys = [lognorm.pdf(x,shape_parameter, scale=xscale ) for x in xs]
        cdfs = [lognorm.cdf(x,shape_parameter, scale=xscale ) for x in xs]
        ax.set(xlabel="MYA")

        center_of_mass, x_value_of_ymax, xs_lt_90, ys_lt_90 = self.get_mode_and_cm(cdfs, xs, ys)

        print("mode:" + str(x_value_of_ymax ))
        print("cm" + str(center_of_mass))
        print("90% within x>=" + str(max(xs_lt_90)))

        ax.scatter(x_value_of_ymax, 0.01,
                            color='darkgreen', marker='^', label="mode", s=100)
        #ax.scatter(x_value_of_cdf_gt_90, 0.01,
        #                    color='pink', marker='^', label="cdf > .90 ", s=100)
        plt.plot(xs,ys,label='pdf')
        plt.plot(xs,cdfs,color='gold',label='cdf')
        #plt.bar(xs,ys,width=bin_size/2,color='gray',alpha=0.5)
        plt.bar(xs_lt_90,ys_lt_90, width=bin_size / 2, color='gray', alpha=0.5, label="90% of the distribution\n"+
                                                                                      "within x <= " + str(max(xs_lt_90)))


        ax.set(xlabel="MYA")

        #plt.text(0.8, 0.9, 'mode at x=' + str(x_value_of_ymax),
        #         horizontalalignment='center',
        #         verticalalignment='center',
        #         transform=ax.transAxes)

        ax.legend()
        out_file_name = os.path.join(test_out, title+"_s" + str(shape_parameter) + ".png")
        plt.savefig(out_file_name)
        plt.close()

        xaxis_limit=max_X + 0.1
        r = lognorm.rvs(shape_parameter, size=1000, scale=xscale)
        fig, ax = plt.subplots(1, 1, figsize=(10, 10))
        bins = np.arange(0, xaxis_limit, bin_size)
        ax.hist(r, density=True, bins=bins, alpha=0.2, label='simulated draws')
        plt.plot(xs,ys,label='original pdf')
        #ax.set_xlim([x[0], x[-1]])
        #ax.legend(loc='best', frameon=False)
        out_file_name = os.path.join(test_out,"Testing"+"_s" + str(shape_parameter) + ".png")
        ax.set(xlim=[0, xaxis_limit])
        ax.legend()
        plt.savefig(out_file_name)
        plt.close()

        self.assertEqual(True, False)  # add assertion here

    def get_mode_and_cm(self, cdfs, xs, ys):
        # get mode & center of mass
        ymax = max(ys)
        xs_of_ymax = []
        xs_of_cdf_gt_90 = []
        weighted_mass = 0
        weights = 0
        xs_lt_90 = []
        ys_lt_90 = []
        for i in range(0, len(xs)):
            x = xs[i]
            y = ys[i]
            cdf = cdfs[i]
            weighted_mass = weighted_mass + (x * y)
            weights = weights + y
            if y == ymax:
                xs_of_ymax.append(x)
            if cdf > .90:
                xs_of_cdf_gt_90.append(x)
            else:
                xs_lt_90.append(x)
                ys_lt_90.append(y)
        x_value_of_ymax = sum(xs_of_ymax) / len(xs_of_ymax)
        #x_value_of_cdf_gt_90 = min(xs_of_cdf_gt_90)
        center_of_mass = weighted_mass / weights
        return center_of_mass, x_value_of_ymax, xs_lt_90, ys_lt_90


if __name__ == '__main__':
    unittest.main()

import os
import unittest

import numpy as np
from matplotlib import pyplot as plt
from scipy import stats


#https://stackoverflow.com/questions/53177243/fitting-a-lognormal-distribution-to-the-data-and-performing-kolmogorov-smirnov-t


class BatchAnalyser(unittest.TestCase):

    def test_ks_test_example(self):
        n = 10000
        # reject null hypothesis if p value is less than significance.
        # https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.kstest.html
        random_draws_from_distribution = stats.norm.rvs(loc=0, size=n, scale=1)
        ks_norm_result = stats.kstest(random_draws_from_distribution, stats.norm.cdf)
        print("accept:\t" + str(ks_norm_result))

        random_draws_from_distribution = stats.norm.rvs(loc=10, size=n, scale=1)
        ks_norm_result = stats.kstest(random_draws_from_distribution, stats.norm.cdf)
        print("reject:\t" + str(ks_norm_result))

        random_draws_from_distribution = stats.norm.rvs(loc=10, size=n, scale=1)
        ks_norm_result = stats.kstest([r - 10 for r in random_draws_from_distribution], stats.norm.cdf)
        print(ks_norm_result)
        s = sum(random_draws_from_distribution)
        m = s / len(random_draws_from_distribution)

        print("mean:" + str(m))
        print("sum:" + str(m))
        random_draws_from_distribution = stats.norm.rvs(loc=10, size=n, scale=1)
        normalized = (np.array(random_draws_from_distribution) - m) / s
        ks_norm_result = stats.kstest([r - 10 for r in random_draws_from_distribution], stats.norm.cdf)
        print(ks_norm_result)
        print("normalized:" + str(normalized))
        ks_norm_result = stats.kstest(normalized, stats.norm.cdf)
        print(ks_norm_result)
        random_draws_from_distribution = stats.norm.rvs(loc=10, size=n, scale=1)
        ks_norm_result = stats.kstest([r - 10 for r in random_draws_from_distribution], stats.norm.cdf)
        print(ks_norm_result)
        # ks_lognorm_result = stats.kstest(Ks_results, lambda x: wgd_lognorm_cdf(x, popt, 0,10))

        # https://stackoverflow.com/questions/51902996/scipy-kstest-used-on-scipy-lognormal-distrubtion
        # my_args = popt[0:3]
        # ks_lognorm_result = stats.kstest(bins, 'lognorm', args=my_args)
        # print("ks_norm_result " + str(ks_norm_result))
        # print("ks_lognorm_result " + str(ks_lognorm_result))

        # https://seaborn.pydata.org/generated/seaborn.kdeplot.html

        # ff1 = normalize([f1], norm="l1")
        # ff2 = normalize([f2], norm="l1")
        # print(ff1[0])
        # chi2=  chisquare(ff1[0], f_exp=ff2[0], ddof=2, axis=0)
        # chi2 = chisquare(f1, f_exp=f2, ddof=(len(f1)-1), axis=0)
        # chi2 = chisquare([1,2,3], f_exp=[1,.1,2.2,2.9], ddof=2, axis=0)
        # print(chi2)
        # f_exp = np.array([44, 24, 29, 3]) / 100 * 189
        # f_obs = np.array([43, 52, 54, 40])
        # chi2 = chisquare(f_obs=f_obs, f_exp=f_exp)

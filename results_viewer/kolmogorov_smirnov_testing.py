import os
import unittest

import numpy as np
from scipy.stats import lognorm, kstest, pearsonr
from scipy import stats


#https://stackoverflow.com/questions/53177243/fitting-a-lognormal-distribution-to-the-data-and-performing-kolmogorov-smirnov-t
#https://www.rdocumentation.org/packages/dgof/versions/1.4/topics/ks.test

class BatchAnalyser(unittest.TestCase):

    #https://forecastegy.com/posts/correlation-between-two-time-series-python/#:~:text=NumPy%20is%20the%20most%20popular,corrcoef%20function.&text=This%20function%20calculates%20the%20Pearson%20correlation%20coefficient.
    # https://pieriantraining.com/pearson-correlation-coefficient-with-scipy-pearsonr/#:~:text=To%20find%20the%20Pearson%20Correlation,value%20is%20the%20p%2Dvalue.

    def test_ks_test_from_web(self):
        x = [1, 2, 3, 4, 5]
        y = [5, 4, 3, 2, 1]

        # Finding Pearson Correlation Coefficient
        corr_coef, p_value = pearsonr(x, y)
        print("Correlation Coefficient:", corr_coef)
    def test_ks_test_from_web(self):

        # reject null hypothesis if p value is less than significance.
        # Ie, high p-value, means it is from the distribution we think it is.

        # https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.lognorm.html
        #Specifically, lognorm.pdf(x, s, loc, scale) is identically equivalent to
        # lognorm.pdf(y, s) / scale
        # with y = (x - loc) / scale

        n=1000
        x = stats.lognorm.rvs(loc=0, size=n, scale=1, s=1)
        print("x:\t" +str(x))
        #x_mean=sum(x)/n
        #print(x_mean)
        #x = [341, 291, 283, 155, 271, 270, 250, 272, 209, 236,
        #     295, 214, 443, 632, 310, 334, 376, 305, 216, 339]
        #x_shifted= [0 for xi in x ]
        sigma, loc, scale = lognorm.fit(x, floc=0)
        #sigma, loc, scale = lognorm.fit(x, loc=0, scale=1)

        mu = np.log(scale)

        print("mu    = %9.5f" % mu)
        print("sigma = %9.5f" % sigma)

        stat, p = kstest(x, 'lognorm', args=(sigma, 0, scale), alternative='two-sided')
        print("KS Test:")
        print("stat    = %9.5f" % stat)
        print("p-value = %9.5f" % p)
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

import math
import os
import unittest

from scipy import stats, integrate
from scipy.optimize import curve_fit
import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
from scipy.stats import lognorm, norm
import config
import process_wrapper
#from results_viewer import curve_fitting
from results_viewer.multi_run_viewer import read_Ks_csv, make_subplot
from results_viewer.collect_aggregate_metrics import metric_result, get_metric_result_data_headers


#https://stackoverflow.com/questions/14770735/how-do-i-change-the-figure-size-with-subplots
class KsFittingTests(unittest.TestCase):

    def test_Ks_fit(self):

        #csv_folder = "/home/tamsen/Data/Specks_outout_from_mesx/sim20_log/Allo4"
        #csv_file_base0 = "Allo4_S020W015_ML_rep0_Ks_by_GeneTree.csv"

        csv_folder = "/home/tamsen/Data/Specks_outout_from_mesx/sim20_log/Auto4"
        csv_file_base0 = "Auto4_S020W020_ML_rep0_Ks_by_GeneTree.csv"

        plot_title=csv_file_base0.replace(".csv",".png")
        full_csv_path=os.path.join(csv_folder,csv_file_base0)
        WGD_time_MYA = 15
        SPC_time_MYA = 20
        WGD_time_Ks = WGD_time_MYA * 0.01
        SPC_time_Ks = SPC_time_MYA * 0.01
        max_Ks_for_x_axis = 0.7

        Ks_results = read_Ks_csv(full_csv_path)

        fig, ax = plt.subplots(1, 1, figsize=(10, 10))
        fig.suptitle(plot_title)
        this_ax = ax

        x = Ks_results
        num_gene_pairs_str = str(len(Ks_results))
        #if max_Ks:
        #    bins = np.arange(0, max_Ks + 0.1, bin_size)
        #    n, bins, patches = this_ax.hist(x, bins=bins, facecolor=plot_color, alpha=0.25,
        #                                    label='ks hist ({0} pairs)'.format(num_gene_pairs_str))

        nBins = 100
        n, bins, patches = this_ax.hist(x, bins=nBins, facecolor='b', alpha=0.25,
                                            label='ks hist ({0} pairs)'.format(num_gene_pairs_str))

        this_ax.axvline(x=WGD_time_Ks, color='b', linestyle='-', label="WGD time, " + str(WGD_time_MYA) + " MYA")
        this_ax.axvline(x=SPC_time_Ks, color='r', linestyle='--',
                        label="SPEC time, " + str(SPC_time_MYA) + " MYA")

        #using std error / covariance
        #https://education.molssi.org/python-data-analysis/03-data-fitting/index.html
        fit_curve_ys1, xs_for_wgd, mode, cm, popt, standard_error = fit_curve_to_hist(bins, n)
        print("mode " + str(mode))
        print("cm " + str(cm))
        print("standard_error " + str(standard_error))

        #reject null hypothesis if p value is less than significance.

        # https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.kstest.html
        # ks_norm_result = stats.kstest(Ks_results, stats.norm.cdf)
        # ks_lognorm_result = stats.kstest(Ks_results, lambda x: wgd_lognorm_cdf(x, popt, 0,10))

        # https://stackoverflow.com/questions/51902996/scipy-kstest-used-on-scipy-lognormal-distrubtion
        #my_args = popt[0:3]
        #ks_lognorm_result = stats.kstest(bins, 'lognorm', args=my_args)
        # print("ks_norm_result " + str(ks_norm_result))
        # print("ks_lognorm_result " + str(ks_lognorm_result))

        label="lognorm fit"# + str(popt)
        this_ax.plot(xs_for_wgd, fit_curve_ys1,color='green', linestyle=':', label=label)


        this_ax.scatter(cm, 0.05,color='darkgreen', marker='o', s=100)  # label="cm",)
        this_ax.scatter(mode, 0.05,color='cyan', marker='^', s=80)  # label="mode")
        this_ax.legend()
        this_ax.set(xlabel="Ks")
        this_ax.set(xlim=[0, 0.5])
        this_ax.set(ylim=[0, 250])
        out_file_name = full_csv_path.replace(".csv", ".hist.png")

        out_file_name = out_file_name.replace(".png", "_Ks_fit" + str(max_Ks_for_x_axis) + ".png")
        plt.savefig(out_file_name)
        plt.close()

def wgd_lognorm(x, amp, scale, x_shift,skew):
    return amp * lognorm.pdf(scale * x + x_shift, skew)

def wgd_lognorm_cdf(x,popt,range_start,range_end):
    values = []
    for value in x:
        integral = integrate.quad(lambda k: wgd_lognorm(k,*popt),range_start,value)[0]
        normalized = integral/integrate.quad(lambda k: wgd_lognorm(k,*popt),range_start,range_end)[0]
        values.append(normalized)
    return np.array(values)

#pdf(x, s, loc=0, scale=1)
#def wgd_lognorm(x, scale, x_shift,skew):
#    return lognorm.pdf(scale * x + x_shift, skew)
#def wgd_lognorm(x, s, loc, scale):
#    return lognorm.pdf(x,s,loc=loc, scale=scale)
def get_xs_from_histogram(Bins,N):

    Num=len(list(N))
    interval=Bins[1]-Bins[0]
    #print(interval)

    xs=[(Bins[i]+Bins[i+1])*0.5 for i in range(0,Num)]
    ys=list(N)
    last_x=0
    while last_x <= 2:
        last_x=xs[-1]+interval
        xs=xs+[last_x]
        ys=ys+[0]

    #plt.plot(xs,ys)
    return [xs,ys]

def fit_curve_to_hist(bins, n):
    fit_fxn = wgd_lognorm  # wgd_skewnorm
    [xs_for_wgd, ys_for_wgd] = get_xs_from_histogram(bins, n)


    try:
        popt, pcov = curve_fit(fit_fxn, xs_for_wgd, ys_for_wgd)
    except Exception as inst:
        print(type(inst))  # the exception type
        print(inst.args)  # arguments stored in .args
        print(inst)  # __str__ allows args to be printed directly,
        return False, xs_for_wgd, 0,0, False

    fit_curve_ys = [fit_fxn(x, *popt) for x in xs_for_wgd]
    SE = np.sqrt(np.diag(pcov))

    RMSE_to_sum = [(fit_curve_ys[i] - ys_for_wgd[i])* (fit_curve_ys[i] - ys_for_wgd[i]) for i in range(0,len(xs_for_wgd))]
    RMSE = math.sqrt( sum(RMSE_to_sum) / len(RMSE_to_sum))

    #get mode & center of mass
    ymax=max(fit_curve_ys)
    xs_of_ymax=[]
    weighted_mass=0
    weights=0
    for i in range(0,len(xs_for_wgd)):
        x=xs_for_wgd[i]
        y=fit_curve_ys[i]
        weighted_mass = weighted_mass+(x*y)
        weights = weights + y
        if y==ymax:
            xs_of_ymax.append(x)


    x_value_of_ymax=sum(xs_of_ymax)/len(xs_of_ymax)
    center_of_mass=weighted_mass/weights

    return fit_curve_ys, xs_for_wgd, x_value_of_ymax,center_of_mass,popt, SE
import math

import numpy as np
from sklearn.metrics import mean_squared_error
from scipy.optimize import curve_fit
from scipy.stats import lognorm, norm, chisquare
from sklearn.preprocessing import normalize
def wgd_gaussian(x, amp, mu, sig):
    return amp * norm.pdf(x, mu, sig)

def wgd_lognorm(x, amp, scale, x_shift,skew):
    return amp * lognorm.pdf(scale * x + x_shift, skew)

#TODO - I cant rememver why this was so complicated...
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

def fit_curve_to_hist(bins, n, fit_fxn ):

    #fit_fxn = wgd_lognorm  # wgd_skewnorm

    minimum_x=0.1 #to keep us away from SSDs
    [xs_for_wgd, ys_for_wgd] = get_xs_from_histogram(bins, n)
    return fit_curve_to_xs_and_ys(xs_for_wgd, ys_for_wgd, fit_fxn)

def fit_curve_to_xs_and_ys(xs_for_wgd, ys_for_wgd, fit_fxn ):

    try:
        popt, pcov = curve_fit(fit_fxn, xs_for_wgd, ys_for_wgd)
    except Exception as inst:
        print(type(inst))  # the exception type
        print(inst.args)  # arguments stored in .args
        print(inst)  # __str__ allows args to be printed directly,
        return False, xs_for_wgd, False

    fit_curve_ys = [fit_fxn(x, *popt) for x in xs_for_wgd]
    RMSE_to_sum = [(fit_curve_ys[i] - ys_for_wgd[i])* (fit_curve_ys[i] - ys_for_wgd[i]) for i in range(0,len(xs_for_wgd))]
    RMSE = math.sqrt( sum(RMSE_to_sum) / len(RMSE_to_sum))
    rms2 = mean_squared_error(ys_for_wgd, fit_curve_ys, squared=False)
    print(rms2)

    chi2 = do_chi2(fit_curve_ys, ys_for_wgd)

    ff1=normalize([f1], norm="l1")
    ff2=normalize([f2], norm="l1")
    #print(ff1[0])
    #chi2=  chisquare(ff1[0], f_exp=ff2[0], ddof=2, axis=0)
    #chi2 = chisquare(f1, f_exp=f2, ddof=(len(f1)-1), axis=0)
    #chi2 = chisquare([1,2,3], f_exp=[1,.1,2.2,2.9], ddof=2, axis=0)
    #print(chi2)
    #f_exp = np.array([44, 24, 29, 3]) / 100 * 189
    #f_obs = np.array([43, 52, 54, 40])
    #chi2 = chisquare(f_obs=f_obs, f_exp=f_exp)
    print(chi2)

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
    num_paralogs = sum(ys_for_wgd)
    goodness_of_fit = curve_fit_metrics(x_value_of_ymax,center_of_mass,num_paralogs,popt, RMSE)
    return fit_curve_ys, xs_for_wgd, goodness_of_fit


def do_chi2(fit_curve_ys, ys_for_wgd):
    f_obs = np.array(ys_for_wgd)
    f_exp = np.array(fit_curve_ys)
    scale = sum(f_obs) / sum(f_exp)
    f_exp_normalized = f_exp * scale
    chi2 = chisquare(f_obs=f_obs, f_exp=f_exp_normalized)
    return chi2


class curve_fit_metrics():

    mode=-1
    cm=-1
    popt=[]
    RMSE = 0
    def __init__(self,mode,cm,num_paralogs,popt, RMSE):
        self.mode = float(mode)
        self.cm = float(cm)
        self.num_paralogs = int(num_paralogs)
        self.popt = popt
        self.RMSE = RMSE

    def to_data_list(self):
        return [self.mode,self.cm,self.num_paralogs,self.popt,self.RMSE]
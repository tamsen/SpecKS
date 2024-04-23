import math

import numpy as np
from scipy.optimize import curve_fit
from scipy.stats import lognorm, norm, chisquare, pearsonr


def coal(x, N):

    k = 1 / N
    return k * math.exp(k*(1 - x))
def sigmoid(x, a, c):
    return 1 / (1 + math.exp(-a*x +c))
def logfit(x, m, b, c):
    return b - m * np.log(c*x)

def logfit2(x, m, b, c, d):
    return b - m * np.log(c*x) -d*x
def linear(x, m, b):
    return m * x + b
def wgd_gaussian(x, amp, mu, sig):
    return amp * norm.pdf(x, mu, sig)

def wgd_lognorm(x, amp, scale, x_shift,skew):
    return amp * lognorm.pdf(scale * x + x_shift, skew)

def logfit_lognorm(x, m, b, c, amp, scale, x_shift,skew):
    return b - m * np.log(c*x) + amp * lognorm.pdf(scale * x + x_shift, skew)
def my_exp(x, Ao,k,b,c):
    return Ao * math.exp(-k*(x-c)) + b
def my_line(x, m, b):
    return m*x+b

def my_decay(x, Ao,k, r):
    return Ao*((1-r)**(k*x))

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

    #rms2 = mean_squared_error(ys_for_wgd, fit_curve_ys, squared=False)
    #chi2 = do_chi2(fit_curve_ys, ys_for_wgd)

    pearson_result=do_pearsons_corr(fit_curve_ys, ys_for_wgd)
    center_of_mass, x_value_of_ymax = get_mode_and_cm(fit_curve_ys, xs_for_wgd)

    num_paralogs = sum(ys_for_wgd)
    goodness_of_fit = curve_fit_metrics(x_value_of_ymax,center_of_mass,num_paralogs,popt, RMSE, pearson_result)
    return fit_curve_ys, xs_for_wgd, goodness_of_fit


def get_mode_and_cm(fit_curve_ys, xs_for_wgd):
    ymax = max(fit_curve_ys)
    xs_of_ymax = []
    weighted_mass = 0
    weights = 0
    for i in range(0, len(xs_for_wgd)):
        x = xs_for_wgd[i]
        y = fit_curve_ys[i]
        weighted_mass = weighted_mass + (x * y)
        weights = weights + y
        if y == ymax:
            xs_of_ymax.append(x)
    x_value_of_ymax = sum(xs_of_ymax) / len(xs_of_ymax)
    center_of_mass = weighted_mass / weights
    return center_of_mass, x_value_of_ymax


def do_chi2(fit_curve_ys, ys_for_wgd):
    f_obs = np.array(ys_for_wgd)
    f_exp = np.array(fit_curve_ys)
    scale = sum(f_obs) / sum(f_exp)
    f_exp_normalized = f_exp * scale
    chi2 = chisquare(f_obs=f_obs, f_exp=f_exp_normalized)
    return chi2

def do_pearsons_corr(fit_curve_ys, ys_for_wgd):

    #Pearson correlation
    #The P-value is the probability that you would have found the current result
    # if the correlation coefficient were in fact zero (null hypothesis).
    # If this probability is lower than the conventional 5% (P<0.05)
    # the correlation coefficient is called statistically significant.

    # (ie lower p-value is better)
    corr_coef, p_value = pearsonr(fit_curve_ys, ys_for_wgd)
    se=standard_error_from_pearsons_corr_r(corr_coef, len(ys_for_wgd))

    return np.array([corr_coef, p_value, se])
def standard_error_from_pearsons_corr_r(r,n):

    r2= r*r
    numerator=1.0-r2
    denominator=n-2.0
    SEr = math.sqrt(numerator/denominator)
    return SEr

class curve_fit_metrics():

    mode=-1
    cm=-1
    popt=[]
    RMSE = 0
    pearsons_corr_result = [] # (corr_coef, p_value, se)
    def __init__(self,mode,cm,num_paralogs,popt, RMSE, pearson_result):
        self.mode = float(mode)
        self.cm = float(cm)
        self.num_paralogs = int(num_paralogs)
        self.popt = popt
        self.RMSE = RMSE

        res = isinstance(pearson_result, str)
        if res:
            splat=pearson_result.replace("[","").replace("]","").split(" ")
            data = [float(s) for s in splat if not s == '']
            self.pearsons_corr_result = data
        else:
            self.pearsons_corr_result = pearson_result

    def to_data_list(self):
        return [self.mode,self.cm,self.num_paralogs,self.popt,self.RMSE,self.pearsons_corr_result]

    def get_pearsons_corr_coef(self):
        return self.pearsons_corr_result[0]

    def get_pearsons_p_value(self):
        return self.pearsons_corr_result[1]

    def get_pearsons_std_error(self):
        return self.pearsons_corr_result[2]
def col_headers_for_curve_fit_metrics():
        return ["mode","cm","num_paralogs","popt","RMSE","pearsons_corr_result (ccoef pv prse)"]
from scipy.optimize import curve_fit
from scipy.stats import lognorm, norm

def wgd_gaussian(x, amp, mu, sig):
    return amp * norm.pdf(x, mu, sig)

def wgd_lognorm(x, amp, scale, x_shift,skew):
    return amp * lognorm.pdf(scale * x + x_shift, skew)


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
    fit_fxn = wgd_gaussian  # wgd_skewnorm
    [xs_for_wgd, ys_for_wgd] = get_xs_from_histogram(bins, n)



    popt, pcov = curve_fit(fit_fxn, xs_for_wgd, ys_for_wgd)
    fit_curve_ys = [fit_fxn(x, *popt) for x in xs_for_wgd]

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

    return fit_curve_ys, xs_for_wgd, x_value_of_ymax,center_of_mass
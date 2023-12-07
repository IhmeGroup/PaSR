import numpy as np
import pandas as pd
from scipy import integrate
from scipy import interpolate
from scipy import stats
from scipy.special import gamma
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
import matplotlib.font_manager
plt.style.use('default')
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": ["Computer Modern Roman"],
})

XSMALL_SIZE = 12
SMALL_SIZE = 14
MEDIUM_SIZE = 16
BIGGER_SIZE = 18

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=XSMALL_SIZE)   # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

mu_arr = 10**(np.linspace(-4, -2, 20))
sk_arr = np.linspace(0.01, 3.0, 20)

data = pd.DataFrame()

for mu in mu_arr:
    for sk in sk_arr:
        lam = mu / (sk / 3)**2
        scale = 1 / lam
        mean, var, skew, kurt = stats.invgauss.stats(mu/lam, scale=lam, moments='mvsk') 

        age = np.linspace(1.0e-5,
                          5 * mu, 1000)
        
        pdf = stats.invgauss.pdf(age, mu/lam, scale=lam)

        bin_edges = np.concatenate((age, [age[-1] + (age[-1] - age[-2])]))
        hist = np.concatenate((pdf, [np.nan]))
        data['bin_edges'] = bin_edges
        data['hist'] = hist
        data.to_csv("hists/hist_mu_{0:.4e}_var_{1:.4e}_skew_{2:.4e}.csv".format(mean, var, skew), header=False, index=False)

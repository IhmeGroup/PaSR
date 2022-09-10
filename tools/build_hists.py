import numpy as np
import pandas as pd
from scipy import integrate
from scipy import interpolate
from scipy import stats
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

mu_arr = 10**(np.linspace(-4, -1, 20))
sk_arr = np.linspace(0.1, 1.1, 10)

data = pd.DataFrame()

for mu in mu_arr:
    for sk in sk_arr:
        alpha = (2.0 / sk)**2
        beta = alpha / mu
        scale = 1 / beta

        mean, var, skew, kurt = stats.gamma.stats(alpha, 0, scale, moments='mvsk')
        # print("Mean: {0:.4e}".format(mean))
        # print("Var: {0:.4e}".format(var))
        # print("Skew: {0:.4e}".format(skew))
        # print("Kurt: {0:.4e}".format(kurt))

        age = np.linspace(1.0e-5,
                        5 * mu, 1000)
        pdf = stats.gamma.pdf(age, alpha, 0, scale)

        # pdf[pdf < 1.0e-300] = 0.0 # smallest double for c++

        bin_edges = np.concatenate((age, [age[-1] + (age[-1] - age[-2])]))
        hist = np.concatenate((pdf, [np.nan]))
        data['bin_edges'] = bin_edges
        data['hist'] = hist
        data.to_csv("hists/hist_mu_{0:.4e}_var_{1:.4e}_skew_{2:.4e}.csv".format(mean, var, skew), header=False, index=False)

    #     plt.plot(age, pdf, linewidth=2, label=r"$\sigma^2 = {0:.2e}$".format(var))
    # plt.xlabel(r"$\tau_{res}$")
    # plt.ylabel(r"$PDF(\tau_{res})$")
    # plt.legend()
    # plt.show()
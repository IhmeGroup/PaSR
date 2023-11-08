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

# mu_arr = [1.8e-3]
# # var = 1.0e-7
# var_arr = 10**np.linspace(-5, -9, 5)
# var_arr = np.linspace(10**-8, 10**-6, 10)
# sk = 10.0
# sk_arr = np.linspace(0.1, 10.0, 5)

data = pd.DataFrame()
# fig, ax = plt.subplots(figsize=[4, 3.2])

for mu in mu_arr:
    for sk in sk_arr:
    # for var in var_arr:

        # Gamma distribution
        # alpha = (2.0 / sk)**2
        # beta = alpha / mu
        # scale = 1 / beta
        # mean, var_test, skew, kurt = stats.gamma.stats(alpha, 0, scale, moments='mvsk')

        # Lognormal distribution
        # var = np.log((sk**2 + np.sqrt(sk**4 + 4*sk**2) + 2)**(1/3) /
        #              2**(1/3) + 2**(1/3)/(sk**2 + np.sqrt(sk**4 + 4*sk**2) + 2)**(1/3) - 1)
        # s = np.sqrt(var)
        # mu_log = np.log(mu) - var / 2
        # scale = np.exp(mu_log)
        # mean, var_test, skew, kurt = stats.lognorm.stats(s, 0, scale, moments='mvsk')

        # Frechet distribution
        # f = lambda a : ((gamma(1 - 3/a) - 3*gamma(1-2/a)*gamma(1-1/a) + 2*(gamma(1-1/a))**3) /
        #                 np.sqrt((gamma(1-2/a) - (gamma(1-1/a))**2)**3)) - sk
        # alpha = fsolve(f, 5)[0]
        # s = np.sqrt(var / (gamma(1-2/alpha) - (gamma(1-1/alpha))**2))
        # m = mu - s * gamma(1-1/alpha)
        # mean, var_test, skew, kurt = stats.invweibull.stats(alpha, m, s, moments='mvsk')

        # Inverse Gauss distribution
        lam = mu / (sk / 3)**2
        scale = 1 / lam
        mean, var, skew, kurt = stats.invgauss.stats(mu/lam, scale=lam, moments='mvsk') 

        # print("Mean: {0:.4e}".format(mean))
        # print("Var: {0:.4e}".format(var))
        # print("Skew: {0:.4e}".format(skew))
        # print("Kurt: {0:.4e}".format(kurt))

        age = np.linspace(1.0e-5,
                          5 * mu, 1000)
        
        # pdf = stats.gamma.pdf(age, alpha, 0, scale)
        # pdf = stats.lognorm.pdf(age, s, 0, scale)
        # pdf = stats.invweibull.pdf(age, alpha, m, s)
        pdf = stats.invgauss.pdf(age, mu/lam, scale=lam)

        # pdf[pdf < 1.0e-300] = 0.0 # smallest double for c++

        bin_edges = np.concatenate((age, [age[-1] + (age[-1] - age[-2])]))
        hist = np.concatenate((pdf, [np.nan]))
        data['bin_edges'] = bin_edges
        data['hist'] = hist
        data.to_csv("hists/hist_mu_{0:.4e}_var_{1:.4e}_skew_{2:.4e}.csv".format(mean, var, skew), header=False, index=False)

    #     # ax.plot(age * 1.0e3, pdf, linewidth=2, label=r"$\sigma^2 = {0:.2e}$".format(var))
    #     ax.plot(age * 1.0e3, pdf / np.amax(pdf), linewidth=2, label=r"$\sigma^2 = {0:.2e}$".format(var))
    # ax.set_xlabel(r"$\tau_{res}$ [ms]")
    # # ax.set_ylabel(r"$\mathcal{P}(\tau_{res})$")
    # ax.set_ylabel(r"$\mathcal{P}(\tau_{res}) / \hat{\mathcal{P}}(\tau_{res})$")
    # plt.legend()
    # plt.tight_layout()
    # plt.savefig("hists.png", bbox_inches='tight', dpi=300)
    # plt.show()

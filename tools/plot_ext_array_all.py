import os
import glob
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

parent_dir_pattern = "./Pe_*"
sim_dir_pattern = "./sim_*"
data_file = "particle_data.csv"
stats_dir = "stats"
Pe_plot = [0.1, 0.15, 0.2]
scale = 1.0
N_levels = 55
T_extinct = 1000.0
T_bounds = [800, 1500]
figsize = [10, 3.2]

mu_les = 1.727e-3
skew_les = 1.823

x_ticks = [1.0e-4, 1.0e-3, 1.0e-2]

levels = np.linspace(T_bounds[0], T_bounds[1], N_levels)
levels_label = np.arange(T_bounds[0], T_bounds[1]+100, 100)

levels_01 = np.linspace(0, 1, 101)
levels_01_label = np.linspace(0, 1, 5)

def parseValue(filename, key):
    if not key in filename:
        raise ValueError("Key " + key + " not found in filename")
    start_index = filename.find(key + "_") + len(key) + 1
    end_index = filename.find("_", start_index)
    if end_index == -1:
        end_index = len(filename)
    return float(filename[start_index:end_index])

parent_dirs = np.array(glob.glob(parent_dir_pattern))

Pe_arr = np.array([parseValue(parent_dir, "Pe") for parent_dir in parent_dirs])

# Sort
i_sorted = np.argsort(Pe_arr)
Pe_arr = Pe_arr[i_sorted]
parent_dirs = parent_dirs[i_sorted]

# Filter
i_filter = np.isin(Pe_arr, Pe_plot)
Pe_arr = Pe_arr[i_filter]
parent_dirs = parent_dirs[i_filter]

n_parents = len(parent_dirs)


fig_T, axs_T = plt.subplots(1, n_parents + 1,
                            figsize=figsize,
                            gridspec_kw={'width_ratios': [1, 1, 1, 0.05]})
fig_p, axs_p = plt.subplots(1, n_parents + 1,
                            figsize=figsize,
                            gridspec_kw={'width_ratios': [1, 1, 1, 0.05]})

for ip, parent_dir in enumerate(parent_dirs):

    os.chdir(parent_dir)

    sim_dirs = glob.glob(sim_dir_pattern)
    n_cases = len(sim_dirs)

    mu = np.zeros(n_cases)
    var = np.zeros(n_cases)
    skew = np.zeros(n_cases)
    # T_fmean = np.zeros(n_cases)
    T_fmean = [[] for i in range(n_cases)]
    T_fmean_mean = np.zeros(n_cases)
    p_lit = np.zeros(n_cases)

    for i, sim in enumerate(sim_dirs):
        print("Reading case " + sim + "...")

        mu[i] = parseValue(sim, "mu")
        var[i] = parseValue(sim, "var")
        skew[i] = parseValue(sim, "skew")

        cases = os.listdir(sim)
        for case in cases:
            stats_file = os.path.join(sim, case, stats_dir, "stats_fmean.csv")

            if not os.path.exists(stats_file):
                T_fmean[i].append(0.0)
                continue

            stats_fmean = pd.read_csv(stats_file)
            if len(stats_fmean) == 0:
                continue
            T_fmean[i].append(stats_fmean['T'].iloc[-1])
        
        T_arr = np.array(T_fmean[i])
        p_lit[i] = np.sum(T_arr >= T_extinct) / len(T_arr)
        T_fmean_mean[i] = np.mean(T_fmean[i])
    
    T_fmean_mean[np.isnan(T_fmean_mean)] = 0.0
    p_lit[np.isnan(p_lit)] = 0.0
    p_lit[p_lit > 0.99] = 0.99
    p_lit[p_lit < 0.01] = 0.01

    os.chdir("../")

    im = axs_T[ip].tricontourf(mu*scale, skew, T_fmean_mean, levels, extend='both')
    axs_T[ip].tricontour(mu*scale, skew, T_fmean_mean, levels=[T_extinct], colors=['w'], linestyles='dashed', linewidths=2)
    axs_T[ip].scatter(mu*scale, skew, c='k', s=10, marker='X')
    axs_T[ip].scatter(mu_les, skew_les, c='r', s=150, marker='X')
    axs_T[ip].set_xscale('log')
    axs_T[ip].set_title(r"$Pe = {0:.2f}$".format(Pe_arr[ip]))
    axs_T[ip].set_xlabel(r"$\overline{\tau}_{res}$ [s]")
    axs_T[ip].set_xticks(x_ticks)
    axs_T[ip].get_xticklabels()[0].set_horizontalalignment('left')
    axs_T[ip].get_xticklabels()[-1].set_horizontalalignment('right')
    if ip == 0:
        axs_T[ip].set_ylabel(r"$\gamma_1$ [-]")
    else:
        axs_T[ip].set_yticklabels([])
    if ip == n_parents - 1:
        cbar = plt.colorbar(im, cax=axs_T[-1], label=r"$\overline{T}$ [-]", ticks=levels_label)
        cbar.ax.set_yticklabels(levels_label)

    im = axs_p[ip].tricontourf(mu*scale, skew, p_lit, levels_01)
    axs_p[ip].scatter(mu*scale, skew, c='k', s=10, marker='X')
    axs_p[ip].scatter(mu_les, skew_les, c='r', s=150, marker='X')
    axs_p[ip].set_xscale('log')
    axs_p[ip].set_title(r"$Pe = {0:.2f}$".format(Pe_arr[ip]))
    axs_p[ip].set_xlabel(r"$\overline{\tau}_{res}$ [s]")
    if ip == 0:
        axs_p[ip].set_ylabel(r"$\gamma_1$ [-]")
    else:
        axs_p[ip].set_yticklabels([])
    if ip == n_parents - 1:
        cbar = plt.colorbar(im, cax=axs_p[-1], label=r"$P(ig)$", ticks=levels_01_label)
        cbar.ax.set_yticklabels(levels_01_label)

fig_T.subplots_adjust(wspace=0.1)
# fig_T.tight_layout()
# fig_T.savefig("T_fmean_array_skew.png", bbox_inches='tight', dpi=300)
fig_T.savefig("T_fmean_array_skew.pdf", bbox_inches='tight')

fig_p.subplots_adjust(wspace=0.1)
# fig_p.tight_layout()
# fig_p.savefig("P_ig_array_skew.png", bbox_inches='tight', dpi=300)
fig_p.savefig("P_ig_array_skew.pdf", bbox_inches='tight')

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

hist_dir = "../hists"
parent_dir_pattern = "./Da_*"
sim_dir_pattern = "./sim_*"
data_file = "particle_data.csv"
stats_dir = "stats"
scale = 1.0
N_levels = 55
T_extinct = 800.0
figsize = [5, 4]

def parseValue(filename, key):
    if not key in filename:
        raise ValueError("Key " + key + " not found in filename")
    start_index = filename.find(key + "_") + len(key) + 1
    end_index = filename.find("_", start_index)
    if end_index == -1:
        end_index = len(filename)
    return float(filename[start_index:end_index])

parent_dirs = glob.glob(parent_dir_pattern)
n_parents = len(parent_dirs)

Da_arr = np.zeros(n_parents)

for ip, parent_dir in enumerate(parent_dirs):

    os.chdir(parent_dir)

    Da_arr[ip] = parseValue(parent_dir, "Da")

    sim_dirs = glob.glob(sim_dir_pattern)
    n_cases = len(sim_dirs)

    mu = np.zeros(n_cases)
    var = np.zeros(n_cases)
    skew = np.zeros(n_cases)
    mode = np.zeros(n_cases)
    # T_fmean = np.zeros(n_cases)
    T_fmean = [[] for i in range(n_cases)]
    T_fmean_mean = np.zeros(n_cases)
    p_lit = np.zeros(n_cases)

    for i, sim in enumerate(sim_dirs):
        print("Reading case " + sim + "...")

        mu[i] = parseValue(sim, "mu")
        var[i] = parseValue(sim, "var")
        skew[i] = parseValue(sim, "skew")

        hist_file = "hist_" + os.path.split(sim)[-1][4:] + ".csv"
        hist = pd.read_csv(os.path.join(hist_dir, hist_file),
                           names=['edges', 'pdf'], header=None)
        centers = (hist['edges'].to_numpy()[1:] + hist['edges'].to_numpy()[:-1]) / 2
        mode[i] = centers[np.argmax(hist['pdf'])]

        cases = os.listdir(sim)
        for case in cases:
            stats_fmean = pd.read_csv(os.path.join(sim, case, stats_dir, "stats_fmean.csv"))
            T_fmean[i].append(stats_fmean['T'].iloc[-1])
        
        T_arr = np.array(T_fmean[i])
        p_lit[i] = np.sum(T_arr >= T_extinct) / len(T_arr)
        T_fmean_mean[i] = np.mean(T_fmean[i])
    
    p_lit[p_lit > 0.99] = 0.99
    p_lit[p_lit < 0.01] = 0.01

    os.chdir("../")

    fig, ax = plt.subplots(figsize=figsize)
    ax.scatter(mu * 1.0e3, T_fmean_mean)
    ax.set_xlabel(r"$\mu$ [ms]")
    ax.set_ylabel(r"$\widetilde{T}$ (K)")
    plt.tight_layout()
    plt.savefig("scatter_T_mu_Da_{0:.4e}.png".format(Da_arr[ip]), bbox_inches='tight', dpi=300)

    fig, ax = plt.subplots(figsize=figsize)
    ax.scatter(var, T_fmean_mean)
    ax.set_xscale('log')
    ax.set_xlabel(r"$\sigma^2$")
    ax.set_ylabel(r"$\widetilde{T}$ (K)")
    plt.tight_layout()
    plt.savefig("scatter_T_var_Da_{0:.4e}.png".format(Da_arr[ip]), bbox_inches='tight', dpi=300)

    fig, ax = plt.subplots(figsize=figsize)
    ax.scatter(skew, T_fmean_mean)
    ax.set_xlabel(r"$\tilde{\mu}_3$")
    ax.set_ylabel(r"$\widetilde{T}$ (K)")
    plt.tight_layout()
    plt.savefig("scatter_T_skew_Da_{0:.4e}.png".format(Da_arr[ip]), bbox_inches='tight', dpi=300)

    fig, ax = plt.subplots(figsize=figsize)
    ax.scatter(mode * 1.0e3, T_fmean_mean)
    ax.set_xlabel(r"$\hat{\tau_{res}}$ [ms]")
    ax.set_ylabel(r"$\widetilde{T}$ (K)")
    plt.tight_layout()
    plt.savefig("scatter_T_mode_Da_{0:.4e}.png".format(Da_arr[ip]), bbox_inches='tight', dpi=300)
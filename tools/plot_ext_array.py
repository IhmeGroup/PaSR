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

sim_dir_pattern = "sim_*"
data_file = "particle_data.csv"
stats_dir = "stats"
scale = 1.0e3

def parseValue(filename, key):
    if not key in filename:
        raise ValueError("Key " + key + " not found in filename")
    start_index = filename.find(key + "_") + len(key) + 1
    end_index = filename.find("_", start_index)
    if end_index == -1:
        end_index = len(filename)
    return float(filename[start_index:end_index])

sim_dirs = glob.glob(sim_dir_pattern)
n_cases = len(sim_dirs)

mu = np.zeros(n_cases)
var = np.zeros(n_cases)
skew = np.zeros(n_cases)
T_fmean = np.zeros(n_cases)

for i, sim in enumerate(sim_dirs):
    print("Reading case " + sim + "...")

    mu[i] = parseValue(sim, "mu")
    var[i] = parseValue(sim, "var")
    skew[i] = parseValue(sim, "skew")

#     data = pd.read_csv(os.path.join(sim, data_file))
#     data['C'] = data['Y_CO'] + data['Y_CO2'] + data['Y_H2'] + data['Y_H2O']
#     plot_data = data.iloc[-10000:]

    stats_fmean = pd.read_csv(os.path.join(sim, stats_dir, "stats_fmean.csv"))
    T_fmean[i] = stats_fmean['T'][-1:]
#     T_fmean[i] = np.mean(stats_fmean['T'][-10:])

fig, ax = plt.subplots()
im = ax.tricontourf(mu*scale, skew, T_fmean)
# ax.tricontour(mu*scale, skew, T_fmean, levels=[0], colors=['k'], linestyles='dashed', linewidths=2)
plt.colorbar(im, label=r"$\widetilde{T}$ (K)")
ax.scatter(mu*scale, skew, c='k', s=10, marker='X')
ax.set_xscale('log')
ax.set_xlabel(r"$\widetilde{\tau_{res}}$ (ms)")
# ax.set_ylabel(r"$\widetilde{\tau_{res}^{''2}}$")
ax.set_ylabel(r"$\widetilde{Skew\left[\tau_{res}\right]}$")
ax.set_title(r"Extinction Map")

plt.tight_layout()
plt.savefig("T_fmean_array.png", bbox_inches='tight', dpi=300)

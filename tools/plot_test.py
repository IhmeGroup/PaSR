import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.font_manager
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

data_file = "./sample/particle_data.csv"

data = pd.read_csv(data_file)

data['C'] = data['Y_CO'] + data['Y_CO2'] + data['Y_H2'] + data['Y_H2O']

plot_data = data.iloc[-10000:]

plt.figure()
plt.scatter(plot_data['Z'], plot_data['C'], s=0.1)

plt.figure()
plt.hist(plot_data['tau_res'])

plt.show()
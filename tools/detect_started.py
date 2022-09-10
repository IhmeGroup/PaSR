import os
import glob

sims = glob.glob("*/sim_*")

for sim in sims:
    started = False
    files = os.listdir(sim)
    for file in files:
        if ".out" in file:
            started = True
            break
    if not started:
        print("Not started: " + sim)


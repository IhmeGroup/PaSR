import os
import shutil
import argparse
import numpy as np
import toml

ref_file_name = "./input.toml"
job_file_name = "./job.slurm"
hist_dir = "./hists"
# Pe_arr = np.array([0.125, 0.25, 0.5])
Pe_arr = np.array([0.06, 0.07, 0.08])
N_samples = 5
sim_dir_prefix = "sim_"
overwrite = True
write_only = False

def parseValue(filename, key):
    if not key in filename:
        raise ValueError("Key " + key + " not found in filename")
    start_index = filename.find(key + "_") + len(key) + 1
    end_index = filename.find("_", start_index)
    if end_index == -1:
        end_index = len(filename)-4
    return float(filename[start_index:end_index])

input_file = toml.load(ref_file_name)
hist_files = np.array(os.listdir(hist_dir))

for Pe in Pe_arr:
    print("Writing cases for Pe = {0:.4e}...".format(Pe))

    Pe_dir_name = "Pe_{0:.4e}".format(Pe)

    if overwrite and os.path.exists(Pe_dir_name):
        shutil.rmtree(Pe_dir_name)
    os.makedirs(Pe_dir_name)

    os.chdir(Pe_dir_name)

    for i, hist_file in enumerate(hist_files):
        id = os.path.splitext(hist_file)[0][5:]
        mu = parseValue(id, "mu")

        sim_dir_name = sim_dir_prefix + id

        if overwrite and os.path.exists(sim_dir_name):
            shutil.rmtree(sim_dir_name)
        os.makedirs(sim_dir_name)

        os.chdir(sim_dir_name)

        print("Writing cases " + sim_dir_name + "...")

        for i_s in range(N_samples):
            sample_name = "{0:02d}".format(i_s) 
            os.makedirs(sample_name)
            os.chdir(sample_name)

            shutil.copyfile(
                os.path.join('../../../', ref_file_name),
                os.path.join('./', ref_file_name))
            input_file = toml.load(ref_file_name)
            input_file['Conditions']['tau_res'] = r'{}'.format(
                os.path.join("../../../", hist_dir, hist_file))
            input_file['Conditions']['tau_mix'] = float(mu * Pe)
            input_file['Conditions']['pilot_t_stop'] = 10.0 * mu
            input_file['Numerics']['t_stop'] = 20.0 * mu

            with open(ref_file_name, "w") as file:
                toml.dump(input_file, file)

            shutil.copyfile(
                os.path.join('../../../', job_file_name),
                os.path.join('./', job_file_name))

            if not write_only:
                os.system("sbatch job.slurm")
            
            os.chdir("../")
        os.chdir("../")
    os.chdir("../")

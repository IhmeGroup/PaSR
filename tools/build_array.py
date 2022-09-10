import os
import shutil
import argparse
import numpy as np
import toml

ref_file_name = "./input.toml"
job_file_name = "./job.slurm"
hist_dir = "../hists"
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
n_files = len(hist_files)

for i, hist_file in enumerate(hist_files):
    id = os.path.splitext(hist_file)[0][5:]
    sim_dir_name = sim_dir_prefix + id
    mu = parseValue(id, "mu")

    if overwrite and os.path.exists(sim_dir_name):
        shutil.rmtree(sim_dir_name)
    os.makedirs(sim_dir_name)

    os.chdir(sim_dir_name)

    shutil.copyfile(
        os.path.join('../', ref_file_name),
        os.path.join('./', ref_file_name))
    input_file = toml.load(ref_file_name)
    input_file['Conditions']['tau_res'] = r'{}'.format(os.path.join("../", hist_dir, hist_file))
    input_file['Conditions']['tau_mix'] = 0.1 * mu
    input_file['Numerics']['t_stop'] = 10.0 * mu

    with open(ref_file_name, "w") as file:
        toml.dump(input_file, file)

    shutil.copyfile(
        os.path.join('../', job_file_name),
        os.path.join('./', job_file_name))
    
#     shutil.copyfile(
#         os.path.join('../', hist_dir, hist_file),
#         os.path.join('./', hist_file))
    
    if not write_only:
        os.system("sbatch job.slurm")
    
    os.chdir("../")

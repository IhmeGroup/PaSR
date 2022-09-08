import os
import shutil
import argparse
import numpy as np
import toml

ref_file_name = "./input.toml"
job_file_name = "./job.slurm"
hist_dir = "./hists"
overwrite = True
write_only = True

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

    if overwrite and os.path.exists(id):
        shutil.rmtree(id)
    os.makedirs(id)

    os.chdir(id)

    shutil.copyfile(
        os.path.join('../', ref_file_name),
        os.path.join('./', ref_file_name))
    input_file = toml.load(ref_file_name)
    input_file['Conditions']['tau_res'] = r'{}'.format(hist_file)

    with open(ref_file_name, "w") as file:
        toml.dump(input_file, file)

    shutil.copyfile(
        os.path.join('../', job_file_name),
        os.path.join('./', job_file_name))
    
    shutil.copyfile(
        os.path.join('../', hist_dir, hist_file),
        os.path.join('./', hist_file))
    
    if not write_only:
        os.system("sbatch job.slurm")
    
    os.chdir("../")


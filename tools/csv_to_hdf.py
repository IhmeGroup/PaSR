import argparse
from tqdm import tqdm
import numpy as np
import pandas as pd

CHUNK_SIZE = 5000000

parser = argparse.ArgumentParser()
parser.add_argument('--datafile', '-d', required=True,
                    type=str, default='particle_data.csv')
args = parser.parse_args()

iter_csv = pd.read_csv(
    args.datafile,
    dtype=np.float32, encoding='utf-8',
    low_memory=True)

iter_csv.to_hdf(args.datafile.replace('.csv', '.h5'),
                key='data', mode='w')


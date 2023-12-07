import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument('-d', '--datafile', type=str, required=True,
                    help='datafile to extract last n lines')
parser.add_argument('-n', '--n_particles', type=int, required=True,
                    help='number of lines to extract')
parser.add_argument('-o', '--output', type=str, required=False,
                    default='restart.csv',
                    help='output file')
args = parser.parse_args()

os.system(f'head -n 1 {args.datafile} > header.tmp')
os.system(f'tail -n {args.n_particles} {args.datafile} > data.tmp')
os.system(f'cat header.tmp data.tmp > {args.output}')
os.system(f'rm header.tmp data.tmp')
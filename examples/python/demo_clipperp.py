import numpy as np
import networkx as nx
import matplotlib.pyplot as plt

import sys, os, pathlib, time, argparse
sys.path.append(os.path.abspath('../../build/bindings/python/clipperpluspy/'))

import clipperpluspy


parser = argparse.ArgumentParser()
parser.add_argument('--folder', type=str, default='../data')
parser.add_argument('--name', type=str, default='7-Scenes')
parser.add_argument('--sequence', type=str, default='kitchen')
parser.add_argument('idx1', type=int)
parser.add_argument('idx2', type=int)
args = parser.parse_args()


FOLDER_PATH = pathlib.Path(args.folder) / args.name / args.sequence / f'scans_{args.idx1}_{args.idx2}'
print(f'loading scans {args.idx1} and {args.idx2} of in {args.name} dataset.')


adj = np.loadtxt(FOLDER_PATH / 'adj.txt')
print(f'Loaded adj matrix of size {adj.shape}')


#NOTE: matlab indices start from 1
omega_gt_idx = np.loadtxt(FOLDER_PATH / 'omega_gt_idx.txt').astype(int) - 1

print('CLIPPER+ start')
start_time = time.perf_counter()
clique_size, clique, certified = clipperpluspy.clipperplus_clique(adj)
end_time = time.perf_counter()
print(f'CLIPPER+ end took {end_time - start_time: 0.03f}s')

intersect = set(clique) & set(omega_gt_idx)

print(f'CLIPPER+ found clique of size {clique_size}');
print(f'ground truth clique size is {len(omega_gt_idx)}.');

print()
print(f'Omega ratio: {len(clique) / len(omega_gt_idx):0.02f}')
print(f'Precision: {len(intersect) / len(clique):0.02f} ')
print(f'Recall: {len(intersect) / len(omega_gt_idx):0.02f}')


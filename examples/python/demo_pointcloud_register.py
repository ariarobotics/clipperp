import numpy as np
import open3d as o3d
from scipy.spatial.transform import Rotation


import sys, os, pathlib, time, argparse
sys.path.append(os.path.abspath('../../build/bindings/python/clipperpluspy/'))

import clipperpluspy


parser = argparse.ArgumentParser()
parser.add_argument('--folder', type=str, default='../data')
parser.add_argument('--name', type=str, default='7-Scenes')
parser.add_argument('--sequence', type=str, default='kitchen')
parser.add_argument('idx1', type=int)
parser.add_argument('idx2', type=int)
parser.add_argument('--threshold', type=float, default=0.05)
args = parser.parse_args()


THREDHOLD = args.threshold
FOLDER_PATH = pathlib.Path(args.folder) / args.name / args.sequence / f'scans_{args.idx1}_{args.idx2}'
print(f'loading scans {args.idx1} and {args.idx2} of in {args.name} dataset.')

source_points = np.loadtxt(FOLDER_PATH / 'feat2.txt')
target_points = np.loadtxt(FOLDER_PATH / 'feat1.txt')

#NOTE: matlab indices start from 1
associations = np.loadtxt(FOLDER_PATH / 'assoc_putative.txt').astype(int) - 1
source_indices, target_indices = associations[:, 1], associations[:, 0]

true_associations = np.loadtxt(FOLDER_PATH / 'assoc_inliers.txt').astype(int) - 1

source_distance_matrix = np.linalg.norm(source_points[source_indices, None] - source_points[source_indices], axis=2)
target_distance_matrix = np.linalg.norm(target_points[target_indices, None] - target_points[target_indices], axis=2)
consistency_matrix = np.abs(source_distance_matrix - target_distance_matrix)

adj = (consistency_matrix <= THREDHOLD).astype(float)
for x in np.unique(source_indices):
    i, = np.where(source_indices == x)
    adj[np.ix_(i, i)] = 0


for x in np.unique(target_indices):
    j, = np.where(target_indices == x)
    adj[np.ix_(j, j)] = 0


X = np.arange(len(target_indices))
adj[X, X] = 0


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



source_indices = source_indices[clique]
target_indices = target_indices[clique]

correspondance = np.vstack([source_indices, target_indices])
correspondance = o3d.utility.Vector2iVector(correspondance.T)


source = o3d.geometry.PointCloud()
source.colors = o3d.utility.Vector3dVector([[1, 0, 0] for i in range(len(source_points))])
source.points = o3d.utility.Vector3dVector(source_points)

target = o3d.geometry.PointCloud()
target.points = o3d.utility.Vector3dVector(target_points)
target.colors = o3d.utility.Vector3dVector([[0, 0, 1] for i in range(len(target_points))])


ground_transformation = np.loadtxt(FOLDER_PATH / 'transf_gt.txt')
gt_R, gt_t = ground_transformation[:3, :3], ground_transformation[:3, 3]


T = o3d.pipelines.registration.TransformationEstimationPointToPoint(False).compute_transformation(
    source,
    target,
    correspondance
)
R, t = T[:3, :3], T[:3, 3]

print(f'Translation error = {np.linalg.norm(t - gt_t):0.02f}m')

rotation_diff = Rotation.from_matrix(R.T @ gt_R)
rotation_diff = rotation_diff.as_euler('xyz', degrees=True)
print(f'Rotation error = {np.linalg.norm(rotation_diff):0.02f} deg')


n = len(source_indices)
points = np.concatenate((source_points[source_indices] @ R.T + t, target_points[target_indices]), axis=0)
lines = []
for i in range(n):
    lines.append([i, i + n])

colors = [[0, 1, 0] for i in range(len(lines))]

line_set = o3d.geometry.LineSet(
    points=o3d.utility.Vector3dVector(points),
    lines=o3d.utility.Vector2iVector(lines),
)
line_set.colors = o3d.utility.Vector3dVector(colors)

o3d.visualization.draw_geometries([source.transform(T), target, line_set])


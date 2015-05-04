import ase
import numpy as np

tol = 1e-9

cell = np.array([
    [1., 0., 0.],
    [0.5, np.sqrt(3) / 2, 0.],
    [0., 0., 1.],
    ]) * 10

pos = np.dot(np.array([
    [0.0, 0.0, 0.0],
    [0.5, 0.5, 0.5],
    [0.2, 0.2, 0.2],
    ]), cell)

a = ase.Atoms('C3', pos, cell=cell, pbc=True)

# get_distance()
assert abs(a.get_distance(0, 1) - 10.0) < tol
assert abs(a.get_distance(0, 1, mic=True) - 5 * np.sqrt(2)) < tol

# get_distance(mic=True)
assert abs(a.get_distance(0, 2) - 4.0) < tol
assert abs(a.get_distance(0, 2, mic=True) - 4.0) < tol

# get_distance(vector=True)
assert all(abs(a.get_distance(0, 1, vector=True)
    - np.array([7.5, np.sqrt(18.75), 5.0])) < tol)
assert all(abs(a.get_distance(0, 2, vector=True)
    - np.array([3., np.sqrt(3.), 2.0])) < tol)

# get_distance(mic=True, vector=True)
assert np.all(abs(a.get_distance(0, 1, mic=True, vector=True)
    - np.array([-2.5, np.sqrt(18.75), -5.0])) < tol)
assert np.all(abs(a.get_distance(0, 2, mic=True, vector=True)
    - np.array([3., np.sqrt(3.), 2.0])) < tol)

# get_all_distances()
all_dist = a.get_all_distances()
assert abs(all_dist[0, 1] - 10.0) < tol
assert abs(all_dist[0, 2] - 4.0) < tol
assert all(abs(np.diagonal(all_dist)) < tol)

# get_all_distances(mic=True)
all_dist_mic = a.get_all_distances(mic=True)
assert abs(all_dist_mic[0, 1] - 5 * np.sqrt(2)) < tol
assert abs(all_dist_mic[0, 2] - 4.0) < tol
assert all(abs(np.diagonal(all_dist)) < tol)

# get_distances()
assert all(abs(a.get_distances(0, [0, 1, 2]) - all_dist[0]) < tol)
assert all(abs(a.get_distances(1, [0, 1, 2]) - all_dist[1]) < tol)
assert all(abs(a.get_distances(2, [0, 1, 2]) - all_dist[2]) < tol)

# get_distances(mic=True)
assert all(abs(a.get_distances(0, [0, 1, 2], mic=True)
    - all_dist_mic[0]) < tol)
assert all(abs(a.get_distances(1, [0, 1, 2], mic=True)
    - all_dist_mic[1]) < tol)
assert all(abs(a.get_distances(2, [0, 1, 2], mic=True)
    - all_dist_mic[2]) < tol)

# get_distances(vec=True)
assert np.all(abs(a.get_distances(0, [0, 1, 2], vector=True)
    - np.array([a.get_distance(0, i, vector=True) for i in [0, 1, 2]])) < tol)
assert np.all(abs(a.get_distances(1, [0, 1, 2], vector=True)
    - np.array([a.get_distance(1, i, vector=True) for i in [0, 1, 2]])) < tol)
assert np.all(abs(a.get_distances(2, [0, 1, 2], vector=True)
    - np.array([a.get_distance(2, i, vector=True) for i in [0, 1, 2]])) < tol)

# get_distances(mic=True, vec=True)
assert np.all(abs(a.get_distances(0, [0, 1, 2], mic=True, vector=True)
    - np.array([a.get_distance(0, i, mic=True, vector=True)
        for i in [0, 1, 2]])) < tol)
assert np.all(abs(a.get_distances(1, [0, 1, 2], mic=True, vector=True)
    - np.array([a.get_distance(1, i, mic=True, vector=True)
        for i in [0, 1, 2]])) < tol)
assert np.all(abs(a.get_distances(2, [0, 1, 2], mic=True, vector=True)
    - np.array([a.get_distance(2, i, mic=True, vector=True)
        for i in [0, 1, 2]])) < tol)

# set_distance
a.set_distance(0, 1, 11.)
assert abs(a.get_distance(0, 1) - 11.) < tol
assert abs(a.get_distance(0, 1, mic=True) - np.sqrt(46)) < tol

# set_distance(mic=True)
a.set_distance(0, 1, 3.)
assert abs(a.get_distance(0, 1, mic=True) - 3.) < tol

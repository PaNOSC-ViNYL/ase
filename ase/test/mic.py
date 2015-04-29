import ase
import numpy as np

tol = 1e-9

cell = np.array([
    [1., 0., 0.],
    [0.5, np.sqrt(3)/2, 0.],
    [0., 0., 1.],
    ]) * 10

pos = np.dot(np.array([[0., 0., 0.], [0.5, 0.5, 0.5]]), cell)

a = ase.Atoms('CC', pos, cell=cell, pbc=True)

assert abs(a.get_distance(0, 1) - 10.0) < tol
assert abs(a.get_distance(0, 1, mic=True) - 5 * np.sqrt(2)) < tol

a.set_distance(0, 1, 5, mic=True)

#assert abs(a.get_distance(0, 1) - 5) < tol
assert abs(a.get_distance(0, 1, mic=True) - 5) < tol

a.set_distance(0, 1, 1.5)

assert abs(a.get_distance(0, 1) - 1.5) < tol
assert abs(a.get_distance(0, 1, mic=True) - 1.5) < tol

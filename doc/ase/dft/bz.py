# creates: cubic.svg, fcc.svg, bcc.svg, tetragonal.svg, orthorhombic.svg
# creates: hexagonal.svg, monoclinic.svg
import numpy as np
import matplotlib.pyplot as plt

from ase.dft.kpoints import (get_special_points, special_paths,
                             parse_path_string)
from ase.dft.bz import bz3d_plot


for X, cell in [
    ('cubic', np.eye(3)),
    ('fcc', [[0, 1, 1], [1, 0, 1], [1, 1, 0]]),
    ('bcc', [[-1, 1, 1], [1, -1, 1], [1, 1, -1]]),
    ('tetragonal', [[1, 0, 0], [0, 1, 0], [0, 0, 1.3]]),
    ('orthorhombic', [[1, 0, 0], [0, 1.2, 0], [0, 0, 1.4]]),
    ('hexagonal', [[1, 0, 0], [-0.5, 3**0.5 / 2, 0], [0, 0, 1]]),
    ('monoclinic', [[1, 0, 0], [0, 1, 0], [0, 0.2, 1]])]:

    icell = np.linalg.inv(cell)
    print(cell, X)
    special_points = get_special_points(cell, X)
    paths = []
    for names in parse_path_string(special_paths[X]):
        points = []
        for name in names:
            points.append(np.dot(icell, special_points[name]))
        paths.append((names, points))

    if X == 'bcc':
        scale = 0.6
        elev = 0.24
        # pi / 13
    else:
        scale = 1
        elev = None

    bz3d_plot(cell=cell, paths=paths, elev=elev, scale=scale)
    plt.savefig(X + '.svg')

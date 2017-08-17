# creates: cubic.svg, fcc.svg, bcc.svg, tetragonal.svg, orthorhombic.svg
# creates: hexagonal.svg, monoclinic.svg
from math import pi, sin, cos

import numpy as np
import matplotlib.pyplot as plt

from ase.dft.kpoints import (get_special_points, special_paths,
                             parse_path_string)


def bz_vertices(cell):
    from scipy.spatial import Voronoi
    icell = np.linalg.inv(cell)
    I = np.indices((3, 3, 3)).reshape((3, 27)) - 1
    G = np.dot(icell, I).T
    vor = Voronoi(G)
    bz1 = []
    for vertices, points in zip(vor.ridge_vertices, vor.ridge_points):
        if -1 not in vertices and 13 in points:
            normal = G[points].sum(0)
            normal /= (normal**2).sum()**0.5
            bz1.append((vor.vertices[vertices], normal))
    return bz1


def plot(cell, paths, elev=None, scale=1):
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    Axes3D  # silence pyflakes

    fig = plt.figure(figsize=(5, 5))
    ax = fig.gca(projection='3d')

    azim = pi / 5
    elev = elev or pi / 6
    x = sin(azim)
    y = cos(azim)
    view = [x * cos(elev), y * cos(elev), sin(elev)]

    bz1 = bz_vertices(cell)

    for points, normal in bz1:
        if np.dot(normal, view) < 0:
            ls = ':'
        else:
            ls = '-'
        x, y, z = np.concatenate([points, points[:1]]).T
        ax.plot(x, y, z, c='k', ls=ls)

    txt = ''
    for names, points in paths:
        x, y, z = np.array(points).T
        ax.plot(x, y, z, c='b', ls='-')

        for name, point in zip(names, points):
            x, y, z = point
            if name == 'G':
                name = '\\Gamma'
            elif len(name) > 1:
                name = name[0] + '_' + name[1]
            ax.text(x, y, z, '$' + name + '$',
                    ha='center', va='bottom', color='r')
            txt += '`' + name + '`-'

        txt = txt[:-1] + '|'

    print(txt[:-1])

    ax.set_axis_off()
    ax.autoscale_view(tight=True)
    s = np.array(paths[0][1]).max() / 0.5 * 0.45 * scale
    ax.set_xlim(-s, s)
    ax.set_ylim(-s, s)
    ax.set_zlim(-s, s)
    ax.set_aspect('equal')

    ax.view_init(azim=azim / pi * 180, elev=elev / pi * 180)


for X, cell in [
    ('cubic', np.eye(3)),
    ('fcc', [[0, 1, 1], [1, 0, 1], [1, 1, 0]]),
    ('bcc', [[-1, 1, 1], [1, -1, 1], [1, 1, -1]]),
    ('tetragonal', [[1, 0, 0], [0, 1, 0], [0, 0, 1.3]]),
    ('orthorhombic', [[1, 0, 0], [0, 1.2, 0], [0, 0, 1.4]]),
    ('hexagonal', [[1, 0, 0], [-0.5, 3**0.5 / 2, 0], [0, 0, 1]]),
    ('monoclinic', [[1, 0, 0], [0, 1, 0], [0, 0.2, 1]])]:

    icell = np.linalg.inv(cell)
    special_points = get_special_points(cell, X)
    paths = []
    for names in parse_path_string(special_paths[X]):
        points = []
        for name in names:
            points.append(np.dot(icell, special_points[name]))
        paths.append((names, points))

    if X == 'bcc':
        scale = 0.6
        elev = pi / 13
    else:
        scale = 1
        elev = None

    plot(cell, paths, elev, scale)
    plt.savefig(X + '.svg')

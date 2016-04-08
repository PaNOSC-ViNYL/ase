# creates: cubic.svg, fcc.svg
from math import pi, sin, cos

import numpy as np

from ase.dft.kpoints import high_symm_path, ibz_points


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
    
            
def plot(cell, points, names):
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    Axes3D  # silence pyflakes
    
    s = np.array(points)[:, 0].max() / 0.5 * 0.45

    bz1 = bz_vertices(cell)

    fig = plt.figure()
    ax = fig.gca(projection='3d')
            
    x, y, z = np.array(points).T
    ax.plot(x, y, z, c='b', ls='-')

    for name, point in zip(names, points):
        x, y, z = point
        if name == 'Gamma':
            name = '\\Gamma'
        ax.text(x, y, z, '$' + name + '$',
                ha='center', va='bottom', color='r')
        
    azim = pi / 6
    elev = pi / 6
    x = sin(azim)
    y = cos(azim)
    view = [x * cos(elev), y * cos(elev), sin(elev)]
    for points, normal in bz1:
        if np.dot(normal, view) < 0:
            ls = ':'
        else:
            ls = '-'
        x, y, z = np.concatenate([points, points[:1]]).T
        ax.plot(x, y, z, c='k', ls=ls)
    
    ax.set_axis_off()
    ax.autoscale_view(tight=True)
    ax.set_xlim(-s, s)
    ax.set_ylim(-s, s)
    ax.set_zlim(-s, s)
    ax.set_aspect('equal')
    
    ax.view_init(azim=azim / pi * 180, elev=elev / pi * 180)

    
import matplotlib.pyplot as plt

for X, cell in [('cubic', np.eye(3)),
                ('fcc', [[0, 1, 1], [1, 0, 1], [1, 1, 0]])]:
    icell = np.linalg.inv(cell)
    points = []
    names = high_symm_path[X]
    for name in names:
        points.append(np.dot(icell, ibz_points[X][name]))
    plot(cell, points, names)
    plt.savefig(X + '.svg')

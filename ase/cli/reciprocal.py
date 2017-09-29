from math import pi, sin, cos

import numpy as np
import matplotlib.pyplot as plt

from ase.io import read
from ase.geometry import crystal_structure_from_cell as csfc
from ase.dft.kpoints import (get_special_points, special_paths,
                             parse_path_string)


class CLICommand:
    short_description = 'Show the reciprocal space'

    @staticmethod
    def add_arguments(parser):
        add = parser.add_argument
        add('name', metavar='input-file')
        add('output', nargs='?')
        add('-v', '--verbose', action='store_true')
        add('--band-path', action='store_true')

    @staticmethod
    def run(args, parser):
        # read from file
        atoms = read(args.name)
        
        # cell
        cell = atoms.get_cell()
        icell = np.linalg.inv(cell)
        cryst = csfc(cell)

        # show info
        if args.verbose:
            print(icell, cryst)

        # band path
        if args.band_path:
            paths = []
            special_points = get_special_points(cell, cryst)
            for names in parse_path_string(special_paths[cryst]):
                points = []
                for name in names:
                    points.append(np.dot(icell, special_points[name]))
                paths.append((names, points))
        else:
            paths = None

        plot(cell, paths)

        if args.output:
            plt.savefig(args.output)
        else:
            plt.show()



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

    if paths is not None:
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

    ax.set_axis_off()
    ax.autoscale_view(tight=True)
    s = points[:][0].max() / 0.5 * 0.45 * scale
    ax.set_xlim(-s, s)
    ax.set_ylim(-s, s)
    ax.set_zlim(-s, s)
    ax.set_aspect('equal')

    ax.view_init(azim=azim / pi * 180, elev=elev / pi * 180)

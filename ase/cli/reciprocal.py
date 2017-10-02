from math import pi, sin, cos

import numpy as np

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
        add('--vectors', action='store_true',
            help="Add reciprocal vectors")
        add('--band-path', action='store_true',
            help="Add the band path")
        kp = parser.add_mutually_exclusive_group(required=False)
        kp.add_argument('--k-points', action='store_true',
                        help="Add k-points of the calculator")
        kp.add_argument('--i-k-points', action='store_true',
                        help="Add irreducible k-points of the calculator")

    @staticmethod
    def run(args, parser):
        # read from file
        atoms = read(args.name)
        
        # cell
        cell = atoms.get_cell()
        icell = atoms.get_reciprocal_cell()
        cryst = csfc(cell)

        # show info
        if args.verbose:
            print('Crystal: ' + cryst)
            print(icell)

        # band path
        if args.band_path:
            paths = []
            special_points = get_special_points(cell, cryst)
            for names in parse_path_string(special_paths[cryst]):
                points = []
                for name in names:
                    points.append(np.dot(icell.T, special_points[name]))
                paths.append((names, points))
        else:
            paths = None

        # k points
        points = None
        if atoms.calc is not None:
            if args.k_points:
                points = atoms.calc.get_bz_k_points()
            elif args.i_k_points:
                points = atoms.calc.get_ibz_k_points()
            if points is not None:
                for i in range(len(points)):
                    points[i] = np.dot(icell.T, points[i])

        # get the correct backend
        if not args.output:
            import matplotlib
            matplotlib.use('Qt4Agg')
        import matplotlib.pyplot as plt

        bz3d_plot(plt, cell, vectors=args.vectors, paths=paths, points=points)

        if args.output:
            plt.savefig(args.output)
        else:
            plt.show()


def bz_vertices(icell):
    from scipy.spatial import Voronoi
    I = (np.indices((3, 3, 3)) - 1).reshape((3, 27))
    G = np.dot(icell.T, I).T
    vor = Voronoi(G)
    bz1 = []
    for vertices, points in zip(vor.ridge_vertices, vor.ridge_points):
        if -1 not in vertices and 13 in points:
            normal = G[points].sum(0)
            normal /= (normal**2).sum()**0.5
            bz1.append((vor.vertices[vertices], normal))
    return bz1


def bz3d_plot(plt, cell, vectors=False, paths=None, points=None,
              elev=None, scale=1):
    from mpl_toolkits.mplot3d import Axes3D
    from mpl_toolkits.mplot3d import proj3d
    from matplotlib.patches import FancyArrowPatch
    Axes3D  # silence pyflakes

    class Arrow3D(FancyArrowPatch):
        def __init__(self, xs, ys, zs, *args, **kwargs):
            FancyArrowPatch.__init__(self, (0, 0), (0, 0), *args, **kwargs)
            self._verts3d = xs, ys, zs

        def draw(self, renderer):
            xs3d, ys3d, zs3d = self._verts3d
            xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, renderer.M)
            self.set_positions((xs[0], ys[0]), (xs[1], ys[1]))
            FancyArrowPatch.draw(self, renderer)

    icell = np.linalg.inv(cell).T
    kpoints = points
    fig = plt.figure(figsize=(5, 5))
    ax = fig.gca(projection='3d')

    azim = pi / 5
    elev = elev or pi / 6
    x = sin(azim)
    y = cos(azim)
    view = [x * cos(elev), y * cos(elev), sin(elev)]

    bz1 = bz_vertices(icell)

    maxp = 0.0
    for points, normal in bz1:
        if np.dot(normal, view) < 0:
            ls = ':'
        else:
            ls = '-'
        x, y, z = np.concatenate([points, points[:1]]).T
        ax.plot(x, y, z, c='k', ls=ls)
        maxp = max(maxp, points.max())

    if vectors:
        ax.add_artist(Arrow3D([0, icell[0, 0]],
                              [0, icell[0, 1]],
                              [0, icell[0, 2]],
                              mutation_scale=20, lw=1,
                              arrowstyle="-|>", color="k"))
        ax.add_artist(Arrow3D([0, icell[1, 0]],
                              [0, icell[1, 1]],
                              [0, icell[1, 2]],
                              mutation_scale=20, lw=1,
                              arrowstyle="-|>", color="k"))
        ax.add_artist(Arrow3D([0, icell[2, 0]],
                              [0, icell[2, 1]],
                              [0, icell[2, 2]],
                              mutation_scale=20, lw=1,
                              arrowstyle="-|>", color="k"))

    if paths is not None:
        txt = ''
        for names, points in paths:
            x, y, z = np.array(points).T
            ax.plot(x, y, z, c='r', ls='-')

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

    if kpoints is not None:
        for p in kpoints:
            # as long it is transposed
            ax.scatter(p[0], p[1], p[2], c='b')

    ax.set_axis_off()
    ax.autoscale_view(tight=True)
    s = maxp / 0.5 * 0.45 * scale
    ax.set_xlim(-s, s)
    ax.set_ylim(-s, s)
    ax.set_zlim(-s, s)
    ax.set_aspect('equal')

    ax.view_init(azim=azim / pi * 180, elev=elev / pi * 180)


def bz2d_plot(plt, icell, vectors=False, paths=None, points=None):
    pass


def bz1d_plot(plt, icell, vectors=False, paths=None, points=None):
    pass

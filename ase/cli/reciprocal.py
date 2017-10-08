import numpy as np

from ase.io import read
from ase.geometry import crystal_structure_from_cell as csfc
from ase.dft.kpoints import (get_special_points, special_paths,
                             parse_path_string)
from ase.bz import bz1d_plot, bz2d_plot, bz3d_plot


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

        bz3d_plot(plt, cell, vectors=args.vectors, paths=paths, points=points,
                  interactive=True)

        if args.output:
            plt.savefig(args.output)
        else:
            plt.show()

from __future__ import print_function
import numpy as np

from ase.io import read
from ase.geometry import crystal_structure_from_cell
from ase.dft.kpoints import (get_special_points, special_paths,
                             parse_path_string, labels_from_kpts,
                             get_monkhorst_pack_size_and_offset)
from ase.dft.bz import bz1d_plot, bz2d_plot, bz3d_plot


class CLICommand:
    short_description = 'Show the reciprocal space'

    @staticmethod
    def add_arguments(parser):
        add = parser.add_argument
        add('name', metavar='input-file')
        add('output', nargs='?')
        add('-v', '--verbose', action='store_true')
        add('-p', '--path', nargs='?', type=str, const='default',
            help='Add a band path.  Example: "GXL".')
        add('-d', '--dimension', type=int, default=3,
            help='Dimension of the cell.')
        add('--no-vectors', action='store_true',
            help="Don't show reciprocal vectors.")
        kp = parser.add_mutually_exclusive_group(required=False)
        kp.add_argument('-k', '--k-points', action='store_true',
                        help='Add k-points of the calculator.')
        kp.add_argument('-i', '--ibz-k-points', action='store_true',
                        help='Add irreducible k-points of the calculator.')

    @staticmethod
    def run(args, parser):
        atoms = read(args.name)

        cell = atoms.get_cell()
        icell = atoms.get_reciprocal_cell()

        try:
            cs = crystal_structure_from_cell(cell)
        except ValueError:
            cs = None

        if args.verbose:
            if cs:
                print('Crystal:', cs)
                print('Special points:', special_paths[cs])
            print('Lattice vectors:')
            for i, v in enumerate(cell):
                print('{}: ({:16.9f},{:16.9f},{:16.9f})'.format(i + 1, *v))
            print('Reciprocal vectors:')
            for i, v in enumerate(icell):
                print('{}: ({:16.9f},{:16.9f},{:16.9f})'.format(i + 1, *v))

        # band path
        if args.path:
            if args.path == 'default':
                args.path = special_paths[cs]
            paths = []
            special_points = get_special_points(cell)
            for names in parse_path_string(args.path):
                points = []
                for name in names:
                    points.append(np.dot(icell.T, special_points[name]))
                paths.append((names, points))
        else:
            paths = None

        # k points
        points = None
        if atoms.calc is not None and hasattr(atoms.calc, 'get_bz_k_points'):
            bzk = atoms.calc.get_bz_k_points()
            if args.path is None:
                try:
                    size, offset = get_monkhorst_pack_size_and_offset(bzk)
                except ValueError:
                    # This was not a MP-grid.  Must be a path in the BZ:
                    args.path = ''.join(labels_from_kpts(bzk, cell)[2])

            if args.k_points:
                points = bzk
            elif args.ibz_k_points:
                points = atoms.calc.get_ibz_k_points()
            if points is not None:
                for i in range(len(points)):
                    points[i] = np.dot(icell.T, points[i])

        # get the correct backend
        if not args.output:
            import matplotlib
            matplotlib.use('Qt4Agg')
        import matplotlib.pyplot as plt

        kwargs = {'cell': cell,
                  'vectors': not args.no_vectors,
                  'paths': paths,
                  'points': points}
        if args.dimension == 1:
            bz1d_plot(**kwargs)
        elif args.dimension == 2:
            bz2d_plot(**kwargs)
        else:
            bz3d_plot(interactive=True, **kwargs)

        if args.output:
            plt.savefig(args.output)
        else:
            plt.show()

from __future__ import print_function

import numpy as np

from ase.io import read
from ase.dft.kpoints import (get_monkhorst_pack_size_and_offset, interpolate,
                             bandpath)
from ase.dft.band_structure import BandStructure


class CLICommand:
    short_description = 'Calculate band-structure'

    @staticmethod
    def add_arguments(parser):
        parser.add_argument('calculation')
        parser.add_argument('-k', '--path')
        parser.add_argument('-p', '--plot', action='store_true')

    @staticmethod
    def run(args):
        main(args)


def main(args):
    atoms = read(args.calculation)
    calc = atoms.calc
    bzkpts = calc.get_bz_k_points()
    ibzkpts = calc.get_ibz_k_points()
    nibz = len(ibzkpts)
    print(nibz)
    eps = np.array([[calc.get_eigenvalues(kpt=k, spin=s)
                     for k in range(nibz)]
                    for s in range(1)])
    try:
        size, offset = get_monkhorst_pack_size_and_offset(bzkpts)
    except:
        path = ibzkpts
    else:
        print(size, offset)
        bz2ibz = calc.get_bz_to_ibz_map()
        path = bandpath(args.path, atoms.cell, 100)[0]
        eps = interpolate(path, eps, 42, bz2ibz, size, offset)
    bs = BandStructure(atoms.cell, path, eps)
    bs.plot()

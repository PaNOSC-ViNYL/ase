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
    kpts = calc.get_bz_k_points()
    size, offset = get_monkhorst_pack_size_and_offset(kpts)
    print(size, offset)
    bz2ibz = calc.get_bz_to_ibz_map()
    nkpts = bz2ibz.max() + 1
    print(nkpts)
    eps = np.array([calc.get_eigenvalues(kpt=k) for k in range(nkpts)])
    path = bandpath(args.path, atoms.cell, 100)
    eps = interpolate(path, eps, 42, bz2ibz, size, offset)
    bs = BandStructure(atoms.cell, path, eps)
    bs.plot()

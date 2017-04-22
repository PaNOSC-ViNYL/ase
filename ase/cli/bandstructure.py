from __future__ import print_function

from ase.io import read
from ase.dft.kpoints import mp


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
    kpts = calc.bz()
    size, offset = mp(kpts)
    np.fft(eps)
    
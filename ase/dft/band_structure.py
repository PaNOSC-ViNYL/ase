import pickle

import numpy as np

from ase.dft.kpoints import xaxis_from_kpts
from ase.parallel import paropen


class BandStructure:
    def __init__(self, atoms=None, calc=None, filename=None, labels=None):
        if filename:
            self.read(filename)
        else:
            atoms = atoms or calc.atoms
            calc = calc or atoms.calc
            self.cell = atoms.cell
            self.kpts = calc.get_ibz_k_points()
            energies = []
            for s in range(calc.get_number_of_spins()):
                energies.append([calc.get_eigenvalues(kpt=k, spin=s)
                                 for k in range(len(self.kpts))])
            self.energies = np.array(energies)
            self.fermilevel = calc.get_fermi_level()
            self.xcoords, self.label_xcoords, self.labels = xaxis_from_kpts(
                self.kpts, self.cell)

        if labels:
            self.labels = labels

    def write(self, filename):
        data = {key: getattr(self, key) for key in
                ['cell', 'kpts', 'energies', 'fermilevel', 'labels',
                 'xcoord', 'label_xcoords']}
        with paropen(filename, 'wb') as f:
            pickle.dump(data, f, protocol=2)  # Python 2+3 compatible?

    def read(self, filename):
        with paropen(filename, 'rb') as f:
            data = pickle.load(f)
        self.__dict__.update(data)

    def plot(self, spin=None, emax=10.0, filename=None, ax=None, show=True):
        if ax is None:
            import matplotlib.pyplot as plt
            ax = plt.gca()

        def pretty(kpt):
            if kpt == 'G':
                kpt = r'\Gamma'
            elif len(kpt) == 2:
                kpt = kpt[0] + '_' + kpt[1]
            return '$' + kpt + '$'

        if spin is None:
            e_skn = self.energies
        else:
            e_skn = self.energies[spin][None]

        for spin, e_kn in enumerate(e_skn):
            color = 'br'[spin]
            for e_k in e_kn.T:
                ax.plot(self.xcoords, e_k, color=color)

        for x in self.label_xcoords[1:-1]:
            ax.axvline(x, color='0.5')

        ax.set_xticks(self.label_xcoords)
        ax.set_xticklabels([pretty(name) for name in self.labels])
        ax.axis(xmin=0, xmax=self.xcoords[-1], ymax=self.fermilevel + emax)
        ax.set_ylabel('eigenvalues [eV]')
        ax.grid(axis='x')
        ax.axhline(self.fermilevel, color='k')

        if filename:
            plt.savefig(filename)

        if show:
            plt.show()

        return ax

import pickle
import sys

import numpy as np

from ase.dft.kpoints import labels_from_kpts
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
                 'xcoords', 'label_xcoords']}
        with paropen(filename, 'wb') as f:
            pickle.dump(data, f, protocol=2)  # Python 2+3 compatible

    def read(self, filename):
        with paropen(filename, 'rb') as f:
            if sys.version_info[0] == 2:
                data = pickle.load(f)
            else:
                data = pickle.load(f, encoding='latin1')
        self.__dict__.update(data)

    def plot(self, spin=None, emax=None, filename=None, ax=None, show=True):
        import matplotlib.pyplot as plt
        if ax is None:
            ax = plt.gca()

        def pretty(kpt):
            if kpt == 'G':
                kpt = r'\Gamma'
            elif len(kpt) == 2:
                kpt = kpt[0] + '_' + kpt[1]
            return '$' + kpt + '$'

        if emax is not None:
            emax = emax + self.fermilevel

        if spin is None:
            e_skn = self.energies
        else:
            e_skn = self.energies[spin][None]

        labels = [pretty(name) for name in self.labels]
        i = 1
        while i < len(labels):
            if self.label_xcoords[i - 1] == self.label_xcoords[i]:
                labels[i - 1] = labels[i - 1][:-1] + ',' + labels[i][1:]
                labels[i] = ''
            i += 1

        for spin, e_kn in enumerate(e_skn):
            color = 'br'[spin]
            for e_k in e_kn.T:
                ax.plot(self.xcoords, e_k, color=color)

        for x in self.label_xcoords[1:-1]:
            ax.axvline(x, color='0.5')

        ax.set_xticks(self.label_xcoords)
        ax.set_xticklabels(labels)
        ax.axis(xmin=0, xmax=self.xcoords[-1], ymax=emax)
        ax.set_ylabel('eigenvalues [eV]')
        ax.axhline(self.fermilevel, color='k')
        plt.tight_layout()

        if filename:
            plt.savefig(filename)

        if show:
            plt.show()

        return ax

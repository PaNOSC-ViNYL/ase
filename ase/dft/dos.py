from math import pi, sqrt

import numpy as npy


class DOS:
    def __init__(self, calc, width=0.1, window=None, npts=201):
        """Electronic Density Of States object"""

        self.npts = npts
        
        self.w_k = calc.get_k_point_weights()

        self.nspins = calc.get_number_of_spins()

        self.e_skn = npy.array([[calc.get_eigenvalues(kpt=k, spin=s)
                                 for k in range(len(self.w_k))]
                                for s in range(self.nspins)])

        self.e_skn -= calc.get_fermi_level()

        if window is None:
            emin = self.e_skn.min() - 5 * self.width
            emax = self.e_skn.max() + 5 * self.width
        else:
            emin, emax = window

        self.energies = npy.linspace(emin, emax, npts)

    def get_energies(self):
        return self.energies

    def delta(self, energy):
        """Return a delta-function centered at 'energy'."""
        x = -((self.energies - energy) / self.width)**2
        return npy.exp(x) / (sqrt(pi) * self.width)

    def get_dos(self, spin=None):
        """ """
        
        if spin is None:
            if self.nspins == 2:
                # Spin-polarized calculation, but no spin specified -
                # return the total DOS:
                return self.get_dos(spin=0) + self.get_dos(spin=1)
            else:
                spin = 0
        
        dos = npy.zeros(self.npts)
        for w, e_n in zip(self.w_k, self.e_skn[spin]):
            for e in e_n:
                dos += w * self.delta(e)
        return dos

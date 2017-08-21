import functools
from math import pi, sqrt

import numpy as np


class DOS:
    def __init__(self, calc, width=0.1, window=None, npts=201):
        """Electronic Density Of States object.

        calc: calculator object
            Any ASE compliant calculator object.
        width: float
            Width of guassian smearing.
        window: tuple of two float
            Use ``window=(emin, emax)``.  If not specified, a window
            big enough to hold all the eigenvalues will be used.
        npts: int
            Number of points.

        """

        self.npts = npts
        self.width = width
        self.w_k = calc.get_k_point_weights()
        self.nspins = calc.get_number_of_spins()
        self.e_skn = np.array([[calc.get_eigenvalues(kpt=k, spin=s)
                                for k in range(len(self.w_k))]
                               for s in range(self.nspins)])
        self.e_skn -= calc.get_fermi_level()

        if window is None:
            emin = self.e_skn.min() - 5 * self.width
            emax = self.e_skn.max() + 5 * self.width
        else:
            emin, emax = window

        self.energies = np.linspace(emin, emax, npts)

    def get_energies(self):
        """Return the array of energies used to sample the DOS.

        The energies are reported relative to the Fermi level.
        """
        return self.energies

    def delta(self, energy):
        """Return a delta-function centered at 'energy'."""
        x = -((self.energies - energy) / self.width)**2
        return np.exp(x) / (sqrt(pi) * self.width)

    def get_dos(self, spin=None):
        """Get array of DOS values.

        The *spin* argument can be 0 or 1 (spin up or down) - if not
        specified, the total DOS is returned.
        """

        if spin is None:
            if self.nspins == 2:
                # Spin-polarized calculation, but no spin specified -
                # return the total DOS:
                return self.get_dos(spin=0) + self.get_dos(spin=1)
            else:
                spin = 0

        dos = np.zeros(self.npts)
        for w, e_n in zip(self.w_k, self.e_skn[spin]):
            for e in e_n:
                dos += w * self.delta(e)
        return dos


def tint(energies, cell, dos, kpts, E):
    zero = energies[0]
    de = energies[1] - zero
    for e in E.T:
        i = e.argsort()
        k = kpts[i, :, np.newaxis]
        e0, e1, e2, e3 = ee = e[i]
        dedk = (np.dot(cell.T, e[1:] - e[0])**2).sum()**0.5
        for j in range(3):
            m = int((ee[j] - zero) / de) + 1
            n = int((ee[j + 1] - zero) / de)
            if n > m:
                v = energies[m:n]
                if j == 0:
                    k1 = (k[0] * (e1 - v) + k[1] * (v - e0)) / (e1 - e0)
                    k2 = (k[0] * (e2 - v) + k[2] * (v - e0)) / (e2 - e0) - k1
                    k3 = (k[0] * (e3 - v) + k[3] * (v - e0)) / (e3 - e0) - k1
                    dos[m:n] += (np.cross(k2, k3, 0, 0)**2).sum(1)**0.5 / dedk
                elif j == 1:
                    k1 = (k[1] * (e2 - v) + k[2] * (v - e1)) / (e2 - e1)
                    k2 = (k[0] * (e2 - v) + k[2] * (v - e0)) / (e2 - e0) - k1
                    k3 = (k[0] * (e3 - v) + k[3] * (v - e0)) / (e3 - e0) - k1
                    k4 = (k[1] * (e3 - v) + k[3] * (v - e1)) / (e3 - e1) - k1
                    dos[m:n] += (np.cross(k2, k3, 0, 0)**2).sum(1)**0.5 / dedk
                    dos[m:n] += (np.cross(k4, k3, 0, 0)**2).sum(1)**0.5 / dedk
                else:
                    k0 = (k[0] * (e3 - v) + k[3] * (v - e0)) / (e3 - e0)
                    k1 = (k[1] * (e3 - v) + k[3] * (v - e1)) / (e3 - e1) - k0
                    k2 = (k[2] * (e3 - v) + k[3] * (v - e2)) / (e3 - e2) - k0
                    dos[m:n] += (np.cross(k1, k2, 0, 0)**2).sum(1)**0.5 / dedk


def tetrahedron(cell, eigs, energies):
    from scipy.spatial import Delaunay

    if 1:
        energies = np.linspace(-2, 3, 2000)
        #eigs = np.array([[[[0, 2], [0, 2]], [[0, 2], [1.0, 1]]]])[:,:,:,:]
        #cell = np.eye(3)

    print(eigs.shape)
    B = np.linalg.inv(cell).T
    I, J, K = eigs.shape[:3]

    indices = np.array([[i, j, k]
                        for i in [0, 1] for j in [0, 1] for k in [0, 1]])
    dt = Delaunay(np.dot(indices, B))

    dos = np.zeros_like(energies)
    integrate = functools.partial(tint, energies, cell, dos)

    for s in dt.simplices:
        kpts = dt.points[s]
        for i in range(I):
            for j in range(J):
                for k in range(K):
                    E = np.array([eigs[(i + a) % I, (j + b) % J, (k + c) % K]
                                  for a, b, c in indices[s]])
                    integrate(kpts, E)
    import matplotlib.pyplot as plt
    plt.plot(energies, dos)
    plt.show()
    asdgf
    return dos

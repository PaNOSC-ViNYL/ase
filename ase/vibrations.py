import pickle
from math import sin, pi, sqrt
from os import remove
from os.path import isfile

import numpy as npy

import ase.units as units
from ase.data import atomic_masses
from ase.io.trajectory import PickleTrajectory


class Vibrations:
    def __init__(self, atoms=None, indices=None, name='vib', delta=0.01):
	self.atoms = atoms
        if indices is None:
            indices = range(len(atoms))
        self.indices = npy.asarray(indices)
        self.name = name
        self.delta = delta
        self.H = None

    def run(self):
        p = self.atoms.positions.copy()
        for a in self.indices:
            for i in range(3):
                for sign in [-1, 1]:
                    filename = '%s.%d%s%s.pckl' % (self.name, a,
                                                   'xyz'[i], ' +-'[sign])
                    if isfile(filename):
                        continue
                    fd = open(filename, 'w')
                    self.atoms.positions[a, i] = p[a, i] + sign * self.delta
                    forces = self.atoms.get_forces()
                    pickle.dump(forces.ravel(), fd)
                    self.atoms.positions[a, i] = p[a, i]
        self.atoms.set_positions(p)

    def clean(self):
        for a in self.indices:
            for i in 'xyz':
                for sign in '-+':
                    name = '%s.%d%s%s.pckl' % (self.name, a, i, sign)
                    if isfile(name):
                        remove(name)
        
    def read(self):
        n = 3 * len(self.indices)
        H = npy.empty((n, n))
        r = 0
        for a in self.indices:
            for i in 'xyz':
                name = '%s.%d%s' % (self.name, a, i)
                H[r] = (pickle.load(open(name + '-.pckl')) -
                        pickle.load(open(name + '+.pckl'))) / (4 * self.delta)
                r += 1
        H += H.copy().T
        self.H = H
        m = self.atoms.get_masses()
        if m is None:
            m = atomic_masses[self.atoms.get_atomic_numbers()]
        self.im = npy.repeat(m[self.indices]**-0.5, 3)
        Q = npy.diag(self.im)
        omega2, modes = npy.linalg.eigh(npy.dot(Q, npy.dot(H, Q)))
        self.modes = modes.T.copy()
        s = units._hbar * 1e10 / sqrt(units._e * units._amu)
        self.hnu = s * omega2.astype(complex)**0.5

    def get_energies(self):
        if self.H is None:
            self.read()
        return self.hnu

    def get_frequencies(self):
        s = 0.01 * units._e / units._c / units._hplanck
        return s * self.get_energies()

    def summary(self):
        hnu = self.get_energies()
        s = 0.01 * units._e / units._c / units._hplanck
        print '---------------------'
        print '  #    meV     cm^-1'
        print '---------------------'
        for n, e in enumerate(hnu):
            if e.imag != 0:
                c = 'i'
                e = e.imag
            else:
                c = ' '
            print '%3d %6.1f%s  %7.1f%s' % (n, 1000 * e, c, s * e, c)
        print '---------------------'
        print 'Zero-point energy: %.3f eV' % self.get_zero_point_energy()
        print

    def get_zero_point_energy(self):
        return 0.5 * self.hnu.real.sum()

    def get_mode(self, n):
        mode = npy.zeros((len(self.atoms), 3))
        mode[self.indices] = (self.modes[n] * self.im).reshape((-1, 3))
        return mode

    def write_mode(self, n, kT=units.kB * 300, nimages=30):
        mode = self.get_mode(n) * sqrt(kT / self.hnu[n])
        p = self.atoms.positions.copy()
        print p
        n %= 3 * len(self.indices)
        traj = PickleTrajectory('%s.%d.traj' % (self.name, n), 'w')
        for x in npy.linspace(0, 2 * pi, nimages, endpoint=False):
            self.atoms.set_positions(p + sin(x) * mode)
            traj.write(self.atoms)
            ## -calc XXXXX
        self.atoms.set_positions(p)
        traj.close()

# -*- coding: utf-8 -*-

"""Resonant Raman intensities"""
from __future__ import print_function
import pickle
import os
from math import sin, pi, sqrt, exp, log
from sys import stdout

import numpy as np

import ase.units as units
from ase.io.trajectory import PickleTrajectory
from ase.parallel import rank, barrier, parprint, paropen
from ase.vibrations import Vibrations
from ase.utils.timing import Timer

class ResonantRaman(Vibrations):
    """Class for calculating vibrational modes and 
    resonant Raman intensities using finite difference.

    atoms: 
      Atoms object
    Excitations:
      Class to calculate the excitations. The class object is
      initialized as Excitations(atoms.get_calculator()) or
      by reading form a file as Excitations('filename', **exkwargs). 
      The file is written by calling the method
      Excitations.write('filename').
    """
    def __init__(self, atoms, Excitations,
                 indices=None, 
                 gsname='rraman', # name for ground state calculations
                 exname=None,     # name for excited state calculations
                 delta=0.01, 
                 nfree=2, 
                 directions=None,
                 exkwargs={},   # kwargs to be passed to Excitations
             ):
        assert(nfree == 2)
        Vibrations.__init__(self, atoms, indices, gsname, delta, nfree)
        self.name = gsname + '-d%.3f' % delta
        if exname is None:
            exname = gsname
        self.exname = exname + '-d%.3f' % delta

        if directions is None:
            self.directions = np.asarray([0, 1, 2])
        else:
            self.directions = np.asarray(directions)

        self.exobj = Excitations
        self.exkwargs = exkwargs

        self.timer = Timer()
        if rank > 0:
            self.txt = open('/dev/null', 'w')
        else:
            self.txt = stdout 

    def calculate(self, filename, fd):
        self.timer.start('Ground state')
        forces = self.atoms.get_forces()
        if rank == 0:
            pickle.dump(forces, fd)
            fd.close()
        self.timer.stop('Ground state')
        self.timer.start('Excitations')
        basename, _ = os.path.splitext(filename)
        excitations = self.exobj(self.atoms.get_calculator())
        excitations.write(basename + '.excitations')
        self.timer.stop('Excitations')

    def get_intensity_tensor(self, omega, gamma=0.1):
        if not hasattr(self, 'modes'):
            self.read()

        if not hasattr(self, 'kss0'):
            self.kss0 = self.exobj(
                self.exname + '.eq.excitations', **self.exkwargs)
            self.kssminus = []
            self.kssplus = []

        ndof = 3 * len(self.indices)
        H = np.empty((ndof, ndof))
        amplitudes = np.zeros((ndof, 3, 3), dtype=complex)
        pre = 1. / (2 * self.delta)
        
        def kappa(exl_f, exl_e, omega, gamma, form='v', eu=units.Hartree):
            """Kappa tensor after Profeta and Mauri
            PRB 63 (2001) 245415"""
            result = np.zeros((3,3), dtype=complex)
            for ex_f, ex_e in zip(exl_f, exl_e):
                # XXX can we savely assume the same ordering ?
                me = ex_f.get_dipole_me(form=form)
                result += (np.outer(me, me.conj()) / 
                           (ex_e.energy * eu - omega - 1j * gamma))
                result += (np.outer(me, me.conj()) / 
                           (ex_e.energy * eu + omega + 1j * gamma))
            return result
            
        r = 0
        for a in self.indices:
            for i in 'xyz':
                try:
                    kssminus = self.kssminus[r]
                    kssplus = self.kssplus[r]
                except IndexError:
                    name = '%s.%d%s' % (self.exname, a, i)
                    self.kssminus.append(
                        self.exobj(name + '-.excitations', 
                                   **self.exkwargs))
                    self.kssplus.append(
                        self.exobj(name + '+.excitations', 
                                   **self.exkwargs))
                    kssminus = self.kssminus[r]
                    kssplus = self.kssplus[r]
                amplitudes[r] = pre * (
                    kappa(kssplus, self.kss0, omega, gamma) - 
                    kappa(kssminus, self.kss0, omega, gamma))
                r += 1
        
        # map to modes
        am = np.dot(amplitudes.T, self.modes.T).T
        return omega**4 * (am * am.conj()).real

    def get_intensities(self, omega, gamma=0.1):
        return self.get_intensity_tensor(omega, gamma).sum(axis=1).sum(axis=1)

    def get_spectrum(self, omega, gamma=0.1, 
                     start=200, end=4000, npts=None, width=4, 
                     type='Gaussian', method='standard', direction='central', 
                     intensity_unit='????', normalize=False):
        """Get resonant Raman spectrum.

        The method returns wavenumbers in cm^-1 with corresponding 
        absolute infrared intensity.
        Start and end point, and width of the Gaussian/Lorentzian should 
        be given in cm^-1.
        normalize=True ensures the integral over the peaks to give the 
        intensity.
        """

        self.type = type.lower()
        assert self.type in ['gaussian', 'lorentzian']

        if not npts: 
            npts = (end - start) / width * 10 + 1
        frequencies = self.get_frequencies(method, direction).real
        intensities = self.get_intensities(omega, gamma)
        prefactor = 1 
        if type == 'lorentzian':
            intensities = intensities * width * pi / 2.
            if normalize:
                prefactor = 2. / width / pi
        else:
            sigma = width / 2. / sqrt(2. * log(2.))
            if normalize:
                prefactor = 1. / sigma / sqrt(2 * pi)
        #Make array with spectrum data
        spectrum = np.empty(npts,np.float)
        energies = np.empty(npts,np.float)
        ediff = (end-start)/float(npts-1)
        energies = np.arange(start, end+ediff/2, ediff)
        for i, energy in enumerate(energies):
            energies[i] = energy
            if type == 'lorentzian':
                spectrum[i] = (intensities * 0.5 * width / pi / (
                        (frequencies - energy)**2 + 0.25 * width**2)).sum()
            else:
                spectrum[i] = (intensities * 
                               np.exp(-(frequencies - energy)**2 / 
                                       2. / sigma**2)).sum()
        return [energies, prefactor * spectrum]

    def write_spectra(self, omega, gamma,
                      out='resonant-raman-spectra.dat', 
                      start=200, end=4000, 
                      npts=None, width=10, 
                      type='Gaussian', method='standard', 
                      direction='central'):
        """Write out spectrum to file.

        First column is the wavenumber in cm^-1, the second column the 
        absolute infrared intensities, and
        the third column the absorbance scaled so that data runs 
        from 1 to 0. Start and end 
        point, and width of the Gaussian/Lorentzian should be given 
        in cm^-1."""
        energies, spectrum = self.get_spectrum(omega, gamma,
                                               start, end, npts, width, 
                                               type, method, direction)

        #Write out spectrum in file. First column is absolute intensities. 
        #Second column is absorbance scaled so that data runs from 1 to 0
        spectrum2 = 1. - spectrum / spectrum.max()
        outdata = np.empty([len(energies), 3])
        outdata.T[0] = energies
        outdata.T[1] = spectrum
        fd = open(out, 'w')
        fd.write('# Resonat Raman spectrum\n')
        fd.write('# omega={0:g} eV, gamma={1:g} eV\n'.format(omega, gamma))
        fd.write('# %s folded, width=%g cm^-1\n' % (type.title(), width))
        fd.write('# [cm^-1]  [a.u.]\n')

        for row in outdata:
            fd.write('%.3f  %15.5g\n' % 
                     (row[0], row[1]))
        fd.close()

    def summary(self, omega, gamma=0.1,
                method='standard', direction='central', 
                intensity_unit='(D/A)2/amu', log=stdout):
        """Print summary for given omega [eV]"""
        hnu = self.get_energies(method, direction)
        s = 0.01 * units._e / units._c / units._hplanck
        intensities = self.get_intensities(omega, gamma)

        if isinstance(log, str):
            log = paropen(log, 'a')

        parprint('-------------------------------------', file=log)
        parprint(' excitation at ' + str(omega) + ' eV', file=log)
        parprint(' gamma ' + str(gamma) + ' eV\n', file=log)
        parprint(' Mode    Frequency        Intensity', file=log)
        parprint('  #    meV     cm^-1      [a.u.]', file=log)
        parprint('-------------------------------------', file=log)
        for n, e in enumerate(hnu):
            if e.imag != 0:
                c = 'i'
                e = e.imag
            else:
                c = ' '
                e = e.real
            parprint(('%3d %6.1f%s  %7.1f%s  %9.3g') % 
                     (n, 1000 * e, c, s * e, c, intensities[n]),
                     file=log)
        parprint('-------------------------------------', file=log)
        parprint('Zero-point energy: %.3f eV' % self.get_zero_point_energy(), 
                 file=log)

    def __del__(self):
        self.timer.write(self.txt)

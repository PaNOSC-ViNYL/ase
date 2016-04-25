# -*- coding: utf-8 -*-

"""Resonant Raman intensities"""

from __future__ import print_function, division
import pickle
import os
import sys

import numpy as np

import ase.units as units
from ase.parallel import rank, parprint, paropen
from ase.vibrations import Vibrations
from ase.vibrations.franck_condon import FranckCondonOverlap
from ase.utils.timing import Timer

# XXX remove gpaw dependence
from gpaw.output import get_txt


class ResonantRaman(Vibrations):
    """Class for calculating vibrational modes and
    resonant Raman intensities using finite difference.

    atoms:
        Atoms object
    Excitations:
        Class to calculate the excitations. The class object is
        initialized as::
            
            Excitations(atoms.get_calculator())
            
        or by reading form a file as::
            
            Excitations('filename', **exkwargs)
            
        The file is written by calling the method
        Excitations.write('filename').
    """
    def __init__(self, atoms, Excitations,
                 indices=None,
                 gsname='rraman',  # name for ground state calculations
                 exname=None,      # name for excited state calculations
                 delta=0.01,
                 nfree=2,
                 directions=None,
                 approximation='Profeta',
                 exkwargs={},      # kwargs to be passed to Excitations
                 exext='.ex.gz',   # extension for Excitation names
                 txt='-',
                 verbose=False,):
        assert(nfree == 2)
        Vibrations.__init__(self, atoms, indices, gsname, delta, nfree)
        self.name = gsname + '-d%.3f' % delta
        if exname is None:
            exname = gsname
        self.exname = exname + '-d%.3f' % delta
        self.exext = exext

        if directions is None:
            self.directions = np.array([0, 1, 2])
        else:
            self.directions = np.array(directions)

        self.approximation = approximation
        self.exobj = Excitations
        self.exkwargs = exkwargs

        self.timer = Timer()
        self.txt = get_txt(txt, rank)

        self.verbose = verbose

    @staticmethod
    def m2(z):
        return (z * z.conj()).real

    def log(self, message, pre='# ', end='\n'):
        if self.verbose:
            self.txt.write(pre + message + end)
            self.txt.flush()

    def calculate(self, filename, fd):
        """Call ground and excited state calculation"""
        self.timer.start('Ground state')
        forces = self.atoms.get_forces()
        if rank == 0:
            pickle.dump(forces, fd)
            fd.close()
        self.timer.stop('Ground state')

        self.timer.start('Excitations')
        basename, _ = os.path.splitext(filename)
        excitations = self.exobj(
            self.atoms.get_calculator(), **self.exkwargs)
        excitations.write(basename + self.exext)
        self.timer.stop('Excitations')

    def read_excitations(self):
        self.timer.start('read excitations')
        self.timer.start('really read')
        self.log('reading ' + self.exname + '.eq' + self.exext)
        ex0_object = self.exobj(self.exname + '.eq' + self.exext,
                                **self.exkwargs)
        self.timer.stop('really read')
        self.timer.start('index')
        matching = frozenset(ex0_object)
        self.timer.stop('index')

        def append(lst, exname, matching):
            self.timer.start('really read')
            self.log('reading ' + exname, end=' ')
            exo = self.exobj(exname, **self.exkwargs)
            lst.append(exo)
            self.timer.stop('really read')
            self.timer.start('index')
            matching = matching.intersection(exo)
            self.log('len={0}, matching={1}'.format(len(exo),
                                                    len(matching)), pre='')
            self.timer.stop('index')
            return matching

        exm_object_list = []
        exp_object_list = []
        for a in self.indices:
            for i in 'xyz':
                name = '%s.%d%s' % (self.exname, a, i)
                matching = append(exm_object_list,
                                  name + '-' + self.exext, matching)
                matching = append(exp_object_list,
                                  name + '+' + self.exext, matching)
        self.ndof = 3 * len(self.indices)
        self.nex = len(matching)
        self.timer.stop('read excitations')

        self.timer.start('select')

        def select(exl, matching):
            mlst = [ex for ex in exl if ex in matching]
            assert(len(mlst) == len(matching))
            return mlst
        ex0 = select(ex0_object, matching)
        exm = []
        exp = []
        r = 0
        for a in self.indices:
            for i in 'xyz':
                exm.append(select(exm_object_list[r], matching))
                exp.append(select(exp_object_list[r], matching))
                r += 1
        self.timer.stop('select')

        self.timer.start('me and energy')

        eu = units.Hartree
        self.ex0E_p = np.array([ex.energy * eu for ex in ex0])
        self.ex0m_pc = np.array(
            [ex.get_dipole_me(form='v') for ex in ex0])
        self.exF_rp = []
        self.exmm_rpc = []
        self.expm_rpc = []
        r = 0
        for a in self.indices:
            for i in 'xyz':
                self.exF_rp.append(
                    [(ep.energy - em.energy)
                     for ep, em in zip(exp[r], exm[r])])
                self.exmm_rpc.append(
                    [ex.get_dipole_me(form='v') for ex in exm[r]])
                self.expm_rpc.append(
                    [ex.get_dipole_me(form='v') for ex in exp[r]])
                r += 1
        self.exF_rp = np.array(self.exF_rp) * eu / 2 / self.delta
        self.exmm_rpc = np.array(self.exmm_rpc)
        self.expm_rpc = np.array(self.expm_rpc)

        self.timer.stop('me and energy')

    def read(self, method='standard', direction='central'):
        """Read data from a pre-performed calculation."""
        if not hasattr(self, 'modes'):
            self.timer.start('read vibrations')
            Vibrations.read(self, method, direction)
            # we now have:
            # self.H     : Hessian matrix
            # self.im    : 1./sqrt(masses)
            # self.modes : Eigenmodes of the mass weighted H
            self.om_r = self.hnu.real    # energies in eV
            self.timer.stop('read vibrations')
        if not hasattr(self, 'ex0E_p'):
            self.read_excitations()

    def get_Huang_Rhys_factors(self, forces_r):
        """Evaluate Huang-Rhys factors derived from forces."""
        self.timer.start('Huang-Rhys')
        assert(len(forces_r.flat) == self.ndof)

        # mass weighted quantities
        Hm_rr = self.im[:, None] * self.H * self.im
        Fm_r = forces_r.flat * self.im
        
        # solve the matrix equation for the equilibrium displacements
        # XXX why are the forces mass weighted ???
        X_r = np.linalg.solve(self.im[:, None] * self.H * self.im,
                              forces_r.flat * self.im)
        d_r = np.dot(self.modes, X_r)

        # Huang-Rhys factors S
        s = 1.e-20 / units.kg / units.C / units._hbar**2  # SI units
        self.timer.stop('Huang-Rhys')
        return s * d_r**2 * self.om_r / 2.

    def get_matrix_element_AlbrechtA(self, omega, gamma=0.1, ml=range(10)):
        """Evaluate Albrecht A term."""
        self.read()
        
        self.timer.start('AlbrechtA')
        
        self.fco = FranckCondonOverlap()

        # excited state forces
        F_pr = self.exF_rp.T
        
        m_rcc = np.zeros((self.ndof, 3, 3), dtype=complex)
        for p, energy in enumerate(self.ex0E_p):
            S_r = self.get_Huang_Rhys_factors(F_pr[p])
            me_cc = np.outer(self.ex0m_pc[p], self.ex0m_pc[p].conj())

            for m in ml:
                self.timer.start('0mm1')
                fco_r = self.fco.direct0mm1(m, S_r)
                self.timer.stop('0mm1')
                self.timer.start('einsum')
                m_rcc += np.einsum('a,bc->abc',
                    fco_r / (energy + m * self.om_r - omega - 1j * gamma),
                                   me_cc)
                m_rcc += np.einsum('a,bc->abc',
                    fco_r / (energy + (m - 1) * self.om_r + omega + 1j * gamma),
                                   me_cc)
                self.timer.stop('einsum')

        self.timer.stop('AlbrechtA')
        return m_rcc

    def get_matrix_element_AlbrechtBC(self, omega, gamma=0.1, ml=range(10),
                                      term='BC'):
        """Evaluate Albrecht B and/or C term(s)."""
        self.read()
        
        self.timer.start('AlbrechtBC')
        
        self.fco = FranckCondonOverlap()

        # excited state forces
        F_pr = self.exF_rp.T

        m_rcc = np.zeros((self.ndof, 3, 3), dtype=complex)
        pre = 1. / (2 * self.delta)
        for p, energy in enumerate(self.ex0E_p):
            S_r = self.get_Huang_Rhys_factors(F_pr[p])
            
            for m in ml:
                self.timer.start('Franck-Condon overlaps')
                fc1mm1_r = self.fco.direct(1, m, S_r)
                fc0mm02_r = self.fco.direct(0, m, S_r)
                fc0mm02_r += np.sqrt(2) * self.fco.direct0mm2(m, S_r)
                self.timer.stop('Franck-Condon overlaps')

                self.timer.start('me dervivatives')
                dm_rc = []
                r = 0
                for a in self.indices:
                    for i in 'xyz':
                        dm_rc.append(
                            (self.expm_rpc[r][p] - self.exmm_rpc[r][p])
                            * self.im[r])
                        r += 1
                dm_rc = np.array(dm_rc) / (2 * self.delta)
                self.timer.stop('me dervivatives')

                self.timer.start('map to modes')
                dm_rc = np.dot(dm_rc.T, self.modes.T).T
                self.timer.stop('map to modes')

                self.timer.start('multiply')
                me_cc = np.outer(self.ex0m_pc[p], self.ex0m_pc[p].conj())
                for r in range(self.ndof):
                    if 'B' in term:
                        denom = (1. /
                    (energy + m * self.om_r[r] - omega - 1j * gamma))
                        m_rcc[r] += (np.outer(dm_rc[r], self.ex0m_pc[p].conj())
                                     * fc1mm1_r[r] * denom)
                        m_rcc[r] += (np.outer(self.ex0m_pc[p], dm_rc[r].conj())
                                     * fc0mm02_r[r] * denom)
                    if 'C' in term:
                        denom = (1. /
                    (energy + (m - 1) * self.om_r[r] + omega + 1j * gamma))
                        m_rcc[r] += (np.outer(self.ex0m_pc[p], dm_rc[r].conj())
                                     * fc1mm1_r[r] * denom)
                        m_rcc[r] += (np.outer(dm_rc[r], self.ex0m_pc[p].conj())
                                     * fc0mm02_r[r] * denom)
                self.timer.stop('multiply')

        self.timer.start('pre_r')
        pre_r = np.where(self.om_r > 0,
                         np.sqrt(units.hbar**2 / 2. / self.om_r), 0)
        for r, p in enumerate(pre_r):
            m_rcc[r] *= p
        self.timer.stop('pre_r')
        self.timer.stop('AlbrechtBC')
        return m_rcc

    def get_matrix_element_Profeta(self, omega, gamma=0.1):
        """Evaluate Albrecht B+C term in Profeta and Mauri approximation"""
        self.read()

        self.timer.start('amplitudes')

        self.timer.start('init')
        V_rcc = np.zeros((self.ndof, 3, 3), dtype=complex)
        pre = 1. / (2 * self.delta)
        self.timer.stop('init')
        
        def kappa(me_pc, e_p, omega, gamma, form='v'):
            """Kappa tensor after Profeta and Mauri
            PRB 63 (2001) 245415"""
            me_ccp = np.empty((3, 3, len(e_p)), dtype=complex)
            for p, me_c in enumerate(me_pc):
                me_ccp[:, :, p] = np.outer(me_pc[p], me_pc[p].conj())
            kappa_ccp = (me_ccp / (e_p - omega - 1j * gamma) +
                         me_ccp.conj() / (e_p + omega + 1j * gamma))
            return kappa_ccp.sum(2)

        self.timer.start('kappa')
        r = 0
        for a in self.indices:
            for i in 'xyz':
                V_rcc[r] = pre * self.im[r] * (
                    kappa(self.expm_rpc[r], self.ex0E_p, omega, gamma) -
                    kappa(self.exmm_rpc[r], self.ex0E_p, omega, gamma))
                r += 1
        self.timer.stop('kappa')

        self.timer.stop('amplitudes')
        
        # map to modes
        self.timer.start('pre_r')
        with np.errstate(divide='ignore'):
            pre_r = np.where(self.om_r > 0,
                             np.sqrt(units.hbar**2 / 2. / self.om_r), 0)
        V_rcc = np.dot(V_rcc.T, self.modes.T).T
        for r, p in enumerate(pre_r):
            V_rcc[r] *= p
        self.timer.stop('pre_r')
        return V_rcc

    def get_matrix_element(self, omega, gamma):
        self.read()
        V_rcc = np.zeros((self.ndof, 3, 3), dtype=complex)
        if self.approximation.lower() == 'profeta':
            V_rcc += self.get_matrix_element_Profeta(omega, gamma)
        elif self.approximation.lower() == 'albrecht a':
            V_rcc += self.get_matrix_element_AlbrechtA(omega, gamma)
        elif self.approximation.lower() == 'albrecht b':
            V_rcc += self.get_matrix_element_AlbrechtBC(omega, gamma, term='B')
        elif self.approximation.lower() == 'albrecht c':
            V_rcc += self.get_matrix_element_AlbrechtBC(omega, gamma, term='C')
        elif self.approximation.lower() == 'albrecht bc':
            V_rcc += self.get_matrix_element_AlbrechtBC(omega, gamma)
        elif self.approximation.lower() == 'albrecht':
            raise NotImplementedError('not yet')
            V_rcc += self.get_matrix_element_AlbrechtA(omega, gamma)
            V_rcc += self.get_matrix_element_AlbrechtBC(omega, gamma)
        else:
            raise NotImplementedError(
                'Approximation {0} not implemented. '.format(
                    self.approximation) +
                'Please use "Profeta", "Albrecht A/B/C/BC", ' +
                'or "Albrecht".')

        return V_rcc

    def get_intensity_tensor(self, omega, gamma):
        V_rcc = self.get_matrix_element(omega, gamma)
        m2 = ResonantRaman.m2
        return omega**4 * m2(V_rcc)

    def get_intentsites_Porezag(self, omega, gamma=0.1):
        """Intensity after Porezag PRB 54 (1996)7830"""
        alpha_rcc = self.get_intensity_tensor(omega, gamma)
        alpha_r = (alpha_rcc[:, 0, 0] + alpha_rcc[:, 1, 1] + 
                   alpha_rcc[:, 2, 2]) / 3.
        m2 = ResonantRaman.m2
        beta2_r = (m2(alpha_rcc[:, 0, 0] - alpha_rcc[:, 1, 1]) +
                   m2(alpha_rcc[:, 0, 0] - alpha_rcc[:, 2, 2]) +
                   m2(alpha_rcc[:, 1, 1] - alpha_rcc[:, 2, 2]) +
                   6. * (m2(alpha_rcc[:, 0, 1]) + m2(alpha_rcc[:, 0, 2]) +
                         m2(alpha_rcc[:, 1, 2]))) / 2.
        return 45 * m2(alpha_r) + 7 * beta2_r

    def get_intensities(self, omega, gamma=0.1, 
#                        method='sum',
                        method='Porezag',
    ):
        if method.lower() == 'sum':
            return self.get_intensity_tensor(omega, gamma).sum(
                axis=1).sum(axis=1)
        elif method.lower() == 'porezag':
            """Intensity after Porezag PRB 54 (1996)7830"""
            alpha_rcc = self.get_intensity_tensor(omega, gamma)
            alpha_r = (alpha_rcc[:, 0, 0] + alpha_rcc[:, 1, 1] + 
                       alpha_rcc[:, 2, 2]) / 3.
            m2 = ResonantRaman.m2
            beta2_r = (m2(alpha_rcc[:, 0, 0] - alpha_rcc[:, 1, 1]) +
                       m2(alpha_rcc[:, 0, 0] - alpha_rcc[:, 2, 2]) +
                       m2(alpha_rcc[:, 1, 1] - alpha_rcc[:, 2, 2]) +
                       6. * (m2(alpha_rcc[:, 0, 1]) + m2(alpha_rcc[:, 0, 2]) +
                             m2(alpha_rcc[:, 1, 2]))) / 2.
            return 45 * m2(alpha_r) + 7 * beta2_r
        else:
            raise NotImplementedError('Method ' + method +
                                      'not implemented')

    def get_spectrum(self, omega, gamma=0.1,
                     start=200.0, end=4000.0, npts=None, width=4.0,
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
            npts = int((end - start) / width * 10 + 1)
        frequencies = self.get_frequencies(method, direction).real
        intensities = self.get_intensities(omega, gamma)
        prefactor = 1
        if type == 'lorentzian':
            intensities = intensities * width * np.pi / 2.
            if normalize:
                prefactor = 2. / width / np.pi
        else:
            sigma = width / 2. / np.sqrt(2. * np.log(2.))
            if normalize:
                prefactor = 1. / sigma / np.sqrt(2 * np.pi)
        # Make array with spectrum data
        spectrum = np.empty(npts)
        energies = np.linspace(start, end, npts)
        for i, energy in enumerate(energies):
            energies[i] = energy
            if type == 'lorentzian':
                spectrum[i] = (intensities * 0.5 * width / np.pi /
                               ((frequencies - energy)**2 +
                                0.25 * width**2)).sum()
            else:
                spectrum[i] = (intensities *
                               np.exp(-(frequencies - energy)**2 /
                                      2. / sigma**2)).sum()
        return [energies, prefactor * spectrum]

    def write_spectrum(self, omega, gamma,
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

        # Write out spectrum in file. First column is absolute intensities.
        outdata = np.empty([len(energies), 3])
        outdata.T[0] = energies
        outdata.T[1] = spectrum
        fd = open(out, 'w')
        fd.write('# Resonant Raman spectrum\n')
        fd.write('# omega={0:g} eV, gamma={1:g} eV\n'.format(omega, gamma))
        fd.write('# %s folded, width=%g cm^-1\n' % (type.title(), width))
        fd.write('# [cm^-1]  [a.u.]\n')

        for row in outdata:
            fd.write('%.3f  %15.5g\n' %
                     (row[0], row[1]))
        fd.close()

    def summary(self, omega, gamma=0.1,
                method='standard', direction='central',
                log=sys.stdout):
        """Print summary for given omega [eV]"""
        hnu = self.get_energies(method, direction)
        s = 0.01 * units._e / units._c / units._hplanck
        intensities = self.get_intensities(omega, gamma)

        if isinstance(log, str):
            log = paropen(log, 'a')

        parprint('-------------------------------------', file=log)
        parprint(' excitation at ' + str(omega) + ' eV', file=log)
        parprint(' gamma ' + str(gamma) + ' eV', file=log)
        parprint(' approximation:', self.approximation, '\n', file=log)
        parprint(' Mode    Frequency        Intensity', file=log)
        parprint('  #    meV     cm^-1      [A^4/amu]', file=log)
        parprint('-------------------------------------', file=log)
        for n, e in enumerate(hnu):
            if e.imag != 0:
                c = 'i'
                e = e.imag
            else:
                c = ' '
                e = e.real
            parprint('%3d %6.1f%s  %7.1f%s  %9.3g' %
                     (n, 1000 * e, c, s * e, c, intensities[n]),
                     file=log)
        parprint('-------------------------------------', file=log)
        parprint('Zero-point energy: %.3f eV' % self.get_zero_point_energy(),
                 file=log)

    def __del__(self):
        self.timer.write(self.txt)

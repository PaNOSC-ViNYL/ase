# -*- coding: utf-8 -*-

"""Resonant Raman intensities"""

from __future__ import print_function, division
import pickle
import os
import sys

import numpy as np

import ase.units as u
from ase.parallel import rank, parprint, paropen
from ase.vibrations import Vibrations
from ase.vibrations.franck_condon import FranckCondonOverlap
from ase.utils.timing import Timer
from ase.utils import convert_string_to_fd

# XXX remove gpaw dependence
from gpaw.lrtddft.spectrum import polarizability


def overlap(calc1, calc2):
    """GPAW overlaps"""
    n1 = calc1.get_number_of_bands()
#    spin1 = calc1.get_number_of_spins() == 2
    n2 = calc2.get_number_of_bands()
#    spin2 = calc2.get_number_of_spins() == 2
    overlap_nn = np.zeros((n1, n2), dtype=float)
    for i1 in range(n1):
        psi1 = calc1.wfs.kpt_u[0].psit_nG[i1]
        for i2 in range(n2):
            psi2 = calc2.wfs.kpt_u[0].psit_nG[i2]
            overlap_nn[i1, i2] = calc1.wfs.gd.integrate(psi1 * psi2)
            parprint(i1, i2, '{0:7.2f}'.format(overlap_nn[i1, i2]))
    return overlap_nn
            

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

        Excitations should work like a list of ex obejects, where:
            ex.get_dipole_me(form='v'):
                gives the velocity form dipole matrix element in
                units |e| * Angstrom
            ex.energy:
                is the transition energy in Hartrees
    """
    def __init__(self, atoms, Excitations,
                 indices=None,
                 gsname='rraman',  # name for ground state calculations
                 exname=None,      # name for excited state calculations
                 delta=0.01,
                 nfree=2,
                 directions=None,
                 approximation='Profeta',
                 observation={'geometry': '-Z(XX)Z'},
                 form='v',         # form of the dipole operator
                 exkwargs={},      # kwargs to be passed to Excitations
                 exext='.ex.gz',   # extension for Excitation names
                 txt='-',
                 verbose=False,
                 overlap=False,
    ):
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
        self.observation = observation
        self.exobj = Excitations
        self.exkwargs = exkwargs
        self.dipole_form = form

        self.timer = Timer()
        self.txt = convert_string_to_fd(txt)

        self.verbose = verbose
        self.overlap = overlap
        
    @staticmethod
    def m2(z):
        return (z * z.conj()).real

    def log(self, message, pre='# ', end='\n'):
        if self.verbose:
            self.txt.write(pre + message + end)
            self.txt.flush()

    def run(self):
        if self.overlap:
            # XXXX stupid way to make a copy
            self.atoms.get_potential_energy()
            self.eq_calculator = self.atoms.get_calculator()
            fname = 'eq.gpw'
            self.eq_calculator.write(fname, 'all')
            self.eq_calculator = self.eq_calculator.__class__(fname)
        Vibrations.run(self)

    def calculate(self, filename, fd):
        """Call ground and excited state calculation"""
        self.timer.start('Ground state')
        forces = self.atoms.get_potential_energy()
        if rank == 0:
            pickle.dump(forces, fd, protocol=2)
            fd.close()
        if self.overlap:
            ov_nn = overlap(self.atoms.get_calculator(), self.eq_calculator)
            np.save(filename + '.ov', ov_nn)
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

        eu = u.Hartree
        self.ex0E_p = np.array([ex.energy * eu for ex in ex0])
        self.ex0m_pc = (np.array(
            [ex.get_dipole_me(form=self.dipole_form) for ex in ex0])
            * u.Bohr)
        exmE_rp = []
        expE_rp = []
        exF_rp = []
        exmm_rpc = []
        expm_rpc = []
        r = 0
        for a in self.indices:
            for i in 'xyz':
                exmE_rp.append([em.energy for em in exm[r]])
                expE_rp.append([ep.energy for ep in exp[r]])
                exF_rp.append(
                    [(em.energy - ep.energy)
                     for ep, em in zip(exp[r], exm[r])])
                exmm_rpc.append(
                    [ex.get_dipole_me(form=self.dipole_form)
                     for ex in exm[r]])
                expm_rpc.append(
                    [ex.get_dipole_me(form=self.dipole_form)
                     for ex in exp[r]])
                r += 1
        # indicees: r=coordinate, p=excitation
        # energies in eV
        self.exmE_rp = np.array(exmE_rp) * eu
        self.expE_rp = np.array(expE_rp) * eu
        # forces in eV / Angstrom
        self.exF_rp = np.array(exF_rp) * eu / 2 / self.delta
        # matrix elements in e * Angstrom
        self.exmm_rpc = np.array(exmm_rpc) * u.Bohr
        self.expm_rpc = np.array(expm_rpc) * u.Bohr

        self.timer.stop('me and energy')

    def read(self, method='standard', direction='central'):
        """Read data from a pre-performed calculation."""
        if not hasattr(self, 'modes'):
            self.timer.start('read vibrations')
            Vibrations.read(self, method, direction)
            # we now have:
            # self.H     : Hessian matrix
            # self.im    : 1./sqrt(masses)
            # self.modes : Eigenmodes of the mass weighted Hessian
            self.om_Q = self.hnu.real    # energies in eV
            # pre-factors for one vibrational excitation
            with np.errstate(divide='ignore'):
                self.vib01_Q = np.where(self.om_Q > 0,
                                        1. / np.sqrt(self.om_Q), 0)
            # -> sqrt(amu) * Angstrom
            self.vib01_Q *= np.sqrt(u._me * u.kg * u.Ha) * u.Bohr
            self.timer.stop('read vibrations')
        if not hasattr(self, 'ex0E_p'):
            self.read_excitations()

    def get_Huang_Rhys_factors(self, forces_r):
        """Evaluate Huang-Rhys factors derived from forces."""
        self.timer.start('Huang-Rhys')
        assert(len(forces_r.flat) == self.ndof)

        # solve the matrix equation for the equilibrium displacements
        X_q = np.linalg.solve(self.im[:, None] * self.H * self.im,
                              forces_r.flat * self.im)
        d_Q = np.dot(self.modes, X_q)

        # Huang-Rhys factors S
        s = 1.e-20 / u.kg / u.C / u._hbar**2  # SI units
        self.timer.stop('Huang-Rhys')
        return s * d_Q**2 * self.om_Q / 2.

    def me_AlbrechtA(self, omega, gamma=0.1, ml=range(16)):
        """Evaluate Albrecht A term.

        Unit: 1/eV
        """
        self.read()

        self.timer.start('AlbrechtA')

        if not hasattr(self, 'fco'):
            self.fco = FranckCondonOverlap()

        # excited state forces
        F_pr = self.exF_rp.T

        m_Qcc = np.zeros((self.ndof, 3, 3), dtype=complex)
        for p, energy in enumerate(self.ex0E_p):
            S_Q = self.get_Huang_Rhys_factors(F_pr[p])
            me_cc = np.outer(self.ex0m_pc[p], self.ex0m_pc[p].conj())

            for m in ml:
                self.timer.start('0mm1')
                fco_Q = self.fco.direct0mm1(m, S_Q)
                self.timer.stop('0mm1')
                self.timer.start('einsum')
                m_Qcc += np.einsum('a,bc->abc',
                                   fco_Q / (energy + m * self.om_Q - omega -
                                            1j * gamma),
                                   me_cc)
                m_Qcc += np.einsum('a,bc->abc',
                                   fco_Q / (energy + (m - 1) * self.om_Q +
                                            omega + 1j * gamma),
                                   me_cc)
                self.timer.stop('einsum')

        self.timer.stop('AlbrechtA')
        return m_Qcc

    def get_matrix_element_AlbrechtBC(self, omega, gamma=0.1, ml=[1],
                                      term='BC'):
        """Evaluate Albrecht B and/or C term(s)."""
        self.read()

        self.timer.start('AlbrechtBC')

        if not hasattr(self, 'fco'):
            self.fco = FranckCondonOverlap()

        # excited state forces
        F_pr = self.exF_rp.T

        m_rcc = np.zeros((self.ndof, 3, 3), dtype=complex)
        for p, energy in enumerate(self.ex0E_p):
            S_r = self.get_Huang_Rhys_factors(F_pr[p])

            for m in ml:
                self.timer.start('Franck-Condon overlaps')
                fc1mm1_r = self.fco.direct(1, m, S_r)
                fc0mm02_r = self.fco.direct(0, m, S_r)
                fc0mm02_r += np.sqrt(2) * self.fco.direct0mm2(m, S_r)
                # XXXXX
                fc1mm1_r[-1] = 1
                fc0mm02_r[-1] = 1
                print(m, fc1mm1_r[-1], fc0mm02_r[-1])
                self.timer.stop('Franck-Condon overlaps')

                self.timer.start('me dervivatives')
                dm_rc = []
                r = 0
                for a in self.indices:
                    for i in 'xyz':
                        dm_rc.append(
                            (self.expm_rpc[r, p] - self.exmm_rpc[r, p]) *
                            self.im[r])
                        print('pm=', self.expm_rpc[r, p], self.exmm_rpc[r, p])
                        r += 1
                dm_rc = np.array(dm_rc) / (2 * self.delta)
                self.timer.stop('me dervivatives')

                self.timer.start('map to modes')
                # print('dm_rc[2], dm_rc[5]', dm_rc[2], dm_rc[5])
                print('dm_rc=', dm_rc)
                dm_rc = np.dot(dm_rc.T, self.modes.T).T
                print('dm_rc[-1][2]', dm_rc[-1][2])
                self.timer.stop('map to modes')

                self.timer.start('multiply')
                # me_cc = np.outer(self.ex0m_pc[p], self.ex0m_pc[p].conj())
                for r in range(self.ndof):
                    if 'B' in term:
                        # XXXX
                        denom = (1. /
                                 (energy + m * 0 * self.om_Q[r] -
                                  omega - 1j * gamma))
                        # ok print('denom=', denom)
                        m_rcc[r] += (np.outer(dm_rc[r],
                                              self.ex0m_pc[p].conj()) *
                                     fc1mm1_r[r] * denom)
                        if r == 5:
                            print('m_rcc[r]=', m_rcc[r][2, 2])
                        m_rcc[r] += (np.outer(self.ex0m_pc[p],
                                              dm_rc[r].conj()) *
                                     fc0mm02_r[r] * denom)
                    if 'C' in term:
                        denom = (1. /
                                 (energy + (m - 1) * self.om_Q[r] +
                                  omega + 1j * gamma))
                        m_rcc[r] += (np.outer(self.ex0m_pc[p],
                                              dm_rc[r].conj()) *
                                     fc1mm1_r[r] * denom)
                        m_rcc[r] += (np.outer(dm_rc[r],
                                              self.ex0m_pc[p].conj()) *
                                     fc0mm02_r[r] * denom)
                self.timer.stop('multiply')
        print('m_rcc[-1]=', m_rcc[-1][2, 2])

        self.timer.start('pre_r')
        with np.errstate(divide='ignore'):
            pre_r = np.where(self.om_Q > 0,
                             np.sqrt(u._hbar**2 / 2. / self.om_Q), 0)
            # print('BC: pre_r=', pre_r)
        for r, p in enumerate(pre_r):
            m_rcc[r] *= p
        self.timer.stop('pre_r')
        self.timer.stop('AlbrechtBC')
        return m_rcc

    def electronic_me_profeta_rcc(self, omega, gamma=0.1,
                                  energy_derivative=False):
        """Evaluate Albrecht B+C term in Profeta and Mauri approximation"""
        self.read()

        self.timer.start('amplitudes')

        self.timer.start('init')
        V_rcc = np.zeros((self.ndof, 3, 3), dtype=complex)
        pre = 1. / (2 * self.delta)
        pre *= u.Hartree * u.Bohr  # e^2Angstrom^2/Ha -> Angstrom^3
        self.timer.stop('init')

        def kappa_cc(me_pc, e_p, omega, gamma, form='v'):
            """Kappa tensor after Profeta and Mauri
            PRB 63 (2001) 245415"""
            k_cc = np.zeros((3, 3), dtype=complex)
            for p, me_c in enumerate(me_pc):
                me_cc = np.outer(me_c, me_c.conj())
                k_cc += (me_cc / (e_p[p] - omega - 1j * gamma) +
                         me_cc.conj() / (e_p[p] + omega + 1j * gamma))
            return k_cc
        
        self.timer.start('kappa')
        r = 0
        for a in self.indices:
            for i in 'xyz':
                if not energy_derivative < 0:
                    V_rcc[r] += pre * (
                        kappa_cc(self.expm_rpc[r], self.ex0E_p,
                                 omega, gamma, self.dipole_form) -
                        kappa_cc(self.exmm_rpc[r], self.ex0E_p,
                                 omega, gamma, self.dipole_form))
                if energy_derivative:
                    V_rcc[r] += pre * (
                        kappa_cc(self.ex0m_pc, self.expE_rp[r],
                                 omega, gamma, self.dipole_form) -
                        kappa_cc(self.ex0m_pc, self.exmE_rp[r],
                                 omega, gamma, self.dipole_form))
                r += 1
        self.timer.stop('kappa')
        self.timer.stop('amplitudes')

        return V_rcc

    def electronic_me_Qcc(self, omega, gamma):
        self.read()
        Vel_rcc = np.zeros((self.ndof, 3, 3), dtype=complex)
        Vel_Qcc = None
        if self.approximation.lower() == 'profeta':
            Vel_rcc += self.electronic_me_profeta_rcc(omega, gamma)
        elif self.approximation.lower() == 'placzek':
            Vel_rcc += self.electronic_me_profeta_rcc(omega, gamma, True)
        elif self.approximation.lower() == 'p-p':
            Vel_rcc += self.electronic_me_profeta_rcc(omega, gamma, -1)
        elif self.approximation.lower() == 'albrecht a':
            Vel_Qcc = self.me_AlbrechtA(omega, gamma)
            # divide through pre-factor
            with np.errstate(divide='ignore'):
                Vel_Qcc *= np.where(self.vib01_Q > 0,
                                    1. / self.vib01_Q, 0)[:, None, None]
        else:
            raise NotImplementedError(
                'Approximation {0} not implemented. '.format(
                    self.approximation) +
                'Please use "Profeta" or "Placzek".')
        
        # map to modes if needed
        if Vel_Qcc is None:
            self.timer.start('map R2Q')
            V_qcc = (Vel_rcc.T * self.im).T  # units Angstrom^2 / sqrt(amu)
            Vel_Qcc = np.dot(V_qcc.T, self.modes.T).T
            self.timer.stop('map R2Q')

        return Vel_Qcc

    def matrix_element(self, omega, gamma):
        self.read()
        V_Qcc = np.zeros((self.ndof, 3, 3), dtype=complex)
        if self.approximation.lower() in ['profeta', 'placzek',
                                          'p-p', 'placzekalpha']:
            me_Qcc = self.electronic_me_Qcc(omega, gamma)
            for Q, vib01 in enumerate(self.vib01_Q):
                V_Qcc[Q] = me_Qcc[Q] * vib01
        elif self.approximation.lower() == 'albrecht a':
            V_Qcc += self.get_matrix_element_AlbrechtA(omega, gamma)
        elif self.approximation.lower() == 'albrecht b':
            raise NotImplementedError('not working')
            V_Qcc += self.get_matrix_element_AlbrechtBC(omega, gamma, term='B')
        elif self.approximation.lower() == 'albrecht c':
            raise NotImplementedError('not working')
            V_Qcc += self.get_matrix_element_AlbrechtBC(omega, gamma, term='C')
        elif self.approximation.lower() == 'albrecht bc':
            raise NotImplementedError('not working')
            V_Qcc += self.get_matrix_element_AlbrechtBC(omega, gamma)
        elif self.approximation.lower() == 'albrecht':
            raise NotImplementedError('not working')
            V_Qcc += self.get_matrix_element_AlbrechtA(omega, gamma)
            V_Qcc += self.get_matrix_element_AlbrechtBC(omega, gamma)
        elif self.approximation.lower() == 'albrecht+profeta':
            V_Qcc += self.get_matrix_element_AlbrechtA(omega, gamma)
            V_Qcc += self.get_matrix_element_Profeta(omega, gamma)
        else:
            raise NotImplementedError(
                'Approximation {0} not implemented. '.format(
                    self.approximation) +
                'Please use "Profeta", "Placzek", "Albrecht A/B/C/BC", ' +
                'or "Albrecht".')

        return V_Qcc

    def get_intensities(self, omega, gamma=0.1):
        m2 = ResonantRaman.m2
        alpha_Qcc = self.matrix_element(omega, gamma)
        if not self.observation:  # XXXX remove
            """Simple sum, maybe too simple"""
            return m2(alpha_Qcc).sum(axis=1).sum(axis=1)
        # XXX enable when appropriate
        #        if self.observation['orientation'].lower() != 'random':
        #            raise NotImplementedError('not yet')

        # random orientation of the molecular frame
        # Woodward & Long,
        # Guthmuller, J. J. Chem. Phys. 2016, 144 (6), 64106
        m2 = ResonantRaman.m2
##        print('alpha_Qcc[-1,:,:]=', alpha_Qcc[-1,:,:].diagonal().sum()/3)
        alpha2_r = m2(alpha_Qcc[:, 0, 0] + alpha_Qcc[:, 1, 1] +
                      alpha_Qcc[:, 2, 2]) / 9.
##        print('alpha2_r**1/2=', np.sqrt(alpha2_r))
        delta2_r = 3 / 4. * (
            m2(alpha_Qcc[:, 0, 1] - alpha_Qcc[:, 1, 0]) +
            m2(alpha_Qcc[:, 0, 2] - alpha_Qcc[:, 2, 0]) +
            m2(alpha_Qcc[:, 1, 2] - alpha_Qcc[:, 2, 1]))
        gamma2_r = (3 / 4. * (m2(alpha_Qcc[:, 0, 1] + alpha_Qcc[:, 1, 0]) +
                              m2(alpha_Qcc[:, 0, 2] + alpha_Qcc[:, 2, 0]) +
                              m2(alpha_Qcc[:, 1, 2] + alpha_Qcc[:, 2, 1])) +
                    (m2(alpha_Qcc[:, 0, 0] - alpha_Qcc[:, 1, 1]) +
                     m2(alpha_Qcc[:, 0, 0] - alpha_Qcc[:, 2, 2]) +
                     m2(alpha_Qcc[:, 1, 1] - alpha_Qcc[:, 2, 2])) / 2)

        if self.observation['geometry'] == '-Z(XX)Z':  # Porto's notation
            return (45 * alpha2_r + 5 * delta2_r + 4 * gamma2_r) / 45.
        elif self.observation['geometry'] == '-Z(XY)Z':  # Porto's notation
            return gamma2_r / 15.
        elif self.observation['scattered'] == 'Z':
            # scattered light in direction of incoming light
            return (45 * alpha2_r + 5 * delta2_r + 7 * gamma2_r) / 45.
        elif self.observation['scattered'] == 'parallel':
            # scattered light perendicular and
            # polarization in plane
            return 6 * gamma2_r / 45.
        elif self.observation['scattered'] == 'perpendicular':
            # scattered light perendicular and
            # polarization out of plane
            return (45 * alpha2_r + 5 * delta2_r + 7 * gamma2_r) / 45.
        else:
            raise NotImplementedError

    def absolute_intensity(self, omega, gamma=0.1):
        """Absolute Raman intensity or Raman scattering factor

        Unit: Ang**4/amu
        """
        m2 = ResonantRaman.m2
        alpha_Qcc = self.electronic_me_Qcc(omega, gamma)
        alpha2_r = m2(alpha_Qcc[:, 0, 0] + alpha_Qcc[:, 1, 1] +
                      alpha_Qcc[:, 2, 2]) / 9.
##        print('alpha2_r**1/2=', np.sqrt(alpha2_r))
        delta2_r = 3 / 4. * (
            m2(alpha_Qcc[:, 0, 1] - alpha_Qcc[:, 1, 0]) +
            m2(alpha_Qcc[:, 0, 2] - alpha_Qcc[:, 2, 0]) +
            m2(alpha_Qcc[:, 1, 2] - alpha_Qcc[:, 2, 1]))
        gamma2_r = (3 / 4. * (m2(alpha_Qcc[:, 0, 1] + alpha_Qcc[:, 1, 0]) +
                              m2(alpha_Qcc[:, 0, 2] + alpha_Qcc[:, 2, 0]) +
                              m2(alpha_Qcc[:, 1, 2] + alpha_Qcc[:, 2, 1])) +
                    (m2(alpha_Qcc[:, 0, 0] - alpha_Qcc[:, 1, 1]) +
                     m2(alpha_Qcc[:, 0, 0] - alpha_Qcc[:, 2, 2]) +
                     m2(alpha_Qcc[:, 1, 1] - alpha_Qcc[:, 2, 2])) / 2)

        return 45 * alpha2_r + 0. * delta2_r + 7 * gamma2_r

    def get_cross_sections(self, omega, gamma=0.1):
        """Returns Raman cross sections for each vibration."""
        I_r = self.get_intensities(omega, gamma)
        pre = 1. / 16 / np.pi**2 / u._eps0**2 / u._c**4
        # frequency of scattered light
        omS_r = omega - self.hnu.real
        return pre * omega * omS_r**3 * I_r

    def get_spectrum(self, omega, gamma=0.1,
                     start=200.0, end=4000.0, npts=None, width=4.0,
                     type='Gaussian', method='standard', direction='central',
                     intensity_unit='????', normalize=False):
        """Get resonant Raman spectrum.

        The method returns wavenumbers in cm^-1 with corresponding
        Raman cross section.
        Start and end point, and width of the Gaussian/Lorentzian should
        be given in cm^-1.
        """

        self.type = type.lower()
        assert self.type in ['gaussian', 'lorentzian']

        if not npts:
            npts = int((end - start) / width * 10 + 1)
        frequencies = self.get_frequencies(method, direction).real
        intensities = self.get_cross_sections(omega, gamma)
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

    def summary(self, omega=0, gamma=0,
                method='standard', direction='central',
                log=sys.stdout):
        """Print summary for given omega [eV]"""
        hnu = self.get_energies(method, direction)
        s = 0.01 * u._e / u._c / u._hplanck
        intensities = self.absolute_intensity(omega, gamma)

        if isinstance(log, str):
            log = paropen(log, 'a')

        parprint('-------------------------------------', file=log)
        parprint(' excitation at ' + str(omega) + ' eV', file=log)
        parprint(' gamma ' + str(gamma) + ' eV', file=log)
        parprint(' approximation:', self.approximation, file=log)
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
            parprint('%3d %6.1f%s  %7.1f%s  %9.1f' %
                     (n, 1000 * e, c, s * e, c, intensities[n]),
                     file=log)
        parprint('-------------------------------------', file=log)
        parprint('Zero-point energy: %.3f eV' % self.get_zero_point_energy(),
                 file=log)

    def __del__(self):
        self.timer.write(self.txt)


class Placzek(ResonantRaman):
    """Raman spectra within the Placzek approximation."""
    def __init__(*args, **kwargs):
        # XXX check for approximation
        kwargs['approximation'] = 'PlaczekAlpha'
        ResonantRaman.__init__(*args, **kwargs)

    def read_excitations(self):
        self.timer.start('read excitations')
        self.exm_r = []
        self.exp_r = []
        r = 0
        for a in self.indices:
            for i in 'xyz':
                exname = '%s.%d%s-' % (self.exname, a, i) + self.exext
                self.log('reading ' + exname)
                self.exm_r.append(self.exobj(exname, **self.exkwargs))
                exname = '%s.%d%s+' % (self.exname, a, i) + self.exext
                self.log('reading ' + exname)
                self.exp_r.append(self.exobj(exname, **self.exkwargs))
                r += 1
        self.ndof = 3 * len(self.indices)
        self.timer.stop('read excitations')

    def electronic_me_Qcc(self, omega, gamma=0):
        self.read()
        
        self.timer.start('init')
        V_rcc = np.zeros((self.ndof, 3, 3), dtype=complex)
        pre = 1. / (2 * self.delta)
        pre *= u.Hartree * u.Bohr  # e^2Angstrom^2/Ha -> Angstrom^3

        om = omega
        if gamma:
            om += 1j * gamma
        self.timer.stop('init')
        
        self.timer.start('alpha derivatives')
        r = 0
        for a in self.indices:
            for i in 'xyz':
                V_rcc[r] = pre * (
                    polarizability(self.exp_r[r], om,
                                   form=self.dipole_form, tensor=True) -
                    polarizability(self.exm_r[r], om,
                                   form=self.dipole_form, tensor=True))
                r += 1
        self.timer.stop('alpha derivatives')
 
        # map to modes
        V_qcc = (V_rcc.T * self.im).T  # units Angstrom^2 / sqrt(amu)
        V_Qcc = np.dot(V_qcc.T, self.modes.T).T
        return V_Qcc


class LrResonantRaman(ResonantRaman):
    """Resonant Raman for linear response

    Quick and dirty approach to enable loading of LrTDDFT calculations
    """
    def read_excitations(self):
        self.timer.start('read excitations')
        self.timer.start('really read')
        self.log('reading ' + self.exname + '.eq' + self.exext)
        ex0_object = self.exobj(self.exname + '.eq' + self.exext,
                                **self.exkwargs)
        self.timer.stop('really read')
        self.timer.start('index')
        matching = frozenset(ex0_object.kss)
        self.timer.stop('index')

        def append(lst, exname, matching):
            self.timer.start('really read')
            self.log('reading ' + exname, end=' ')
            exo = self.exobj(exname, **self.exkwargs)
            lst.append(exo)
            self.timer.stop('really read')
            self.timer.start('index')
            matching = matching.intersection(exo.kss)
            self.log('len={0}, matching={1}'.format(len(exo.kss),
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
        self.timer.stop('read excitations')

        self.timer.start('select')

        def select(exl, matching):
            exl.diagonalize(**self.exkwargs)
            mlst = [ex for ex in exl]
#            mlst = [ex for ex in exl if ex in matching]
#            assert(len(mlst) == len(matching))
            return mlst
        ex0 = select(ex0_object, matching)
        self.nex = len(ex0)
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

        eu = u.Hartree
        self.ex0E_p = np.array([ex.energy * eu for ex in ex0])
#        self.exmE_p = np.array([ex.energy * eu for ex in exm])
#        self.expE_p = np.array([ex.energy * eu for ex in exp])
        self.ex0m_pc = (np.array(
            [ex.get_dipole_me(form=self.dipole_form) for ex in ex0])
            * u.Bohr)
        self.exF_rp = []
        exmE_rp = []
        expE_rp = []
        exmm_rpc = []
        expm_rpc = []
        r = 0
        for a in self.indices:
            for i in 'xyz':
                exmE_rp.append([em.energy for em in exm[r]])
                expE_rp.append([ep.energy for ep in exp[r]])
                self.exF_rp.append(
                    [(em.energy - ep.energy)
                     for ep, em in zip(exp[r], exm[r])])
                exmm_rpc.append(
                    [ex.get_dipole_me(form=self.dipole_form) for ex in exm[r]])
                expm_rpc.append(
                    [ex.get_dipole_me(form=self.dipole_form) for ex in exp[r]])
                r += 1
        self.exmE_rp = np.array(exmE_rp) * eu
        self.expE_rp = np.array(expE_rp) * eu
        self.exF_rp = np.array(self.exF_rp) * eu / 2 / self.delta
        self.exmm_rpc = np.array(exmm_rpc) * u.Bohr
        self.expm_rpc = np.array(expm_rpc) * u.Bohr

        self.timer.stop('me and energy')

# -*- coding: utf-8 -*-

from __future__ import print_function, division
import sys
import numpy as np

import ase.units as u
from ase.parallel import parprint, paropen
from ase.vibrations.resonant_raman import ResonantRaman
from ase.vibrations.franck_condon import FranckCondonOverlap


class Albrecht(ResonantRaman):
    def __init__(self, *args, **kwargs):
        """
        Parameters
        ----------
        all from ResonantRaman.__init__
        combinations: int
            Combinations to consider for multiple excitations.
            Default is 1, possible 2
        skip: int
            Number of first transitions to exclude. Default 0,
            recommended: 5 for linear molecules, 6 for other molecules
        """
        self.combinations = kwargs.pop('combinations', 1)
        self.skip = kwargs.pop('skip', 0)
        ResonantRaman.__init__(self, *args, **kwargs)
        
    def read(self, method='standard', direction='central'):
        ResonantRaman.read(self, method, direction)

        # single transitions and their occupation
        om_Q = self.om_Q[self.skip:]
        om_v = om_Q
        ndof = len(om_Q)
        n_vQ = np.eye(ndof, dtype=int)
        if self.combinations > 1:
            # double transitions and their occupation
            om_v = list(om_v)
            n_vQ = list(n_vQ)
            for i in range(ndof):
                n_Q = np.zeros(ndof, dtype=int)
                n_Q[i] = 1
                for j in range(i, ndof):
                    n_vQ.append(n_Q.copy())
                    n_vQ[-1][j] += 1
                    om_v.append(np.dot(n_vQ[-1], om_Q))
            om_v = np.array(om_v)
            n_vQ = np.array(n_vQ)
            
        self.om_v = om_v
        self.n_vQ = n_vQ

    def Huang_Rhys_factors(self, forces_r):
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

    def omegaLS(self,  omega, gamma):
        omL = omega + 1j * gamma
        omS_Q = omL - self.om_Q
        return omL, omS_Q

    def meA(self, omega, gamma=0.1, ml=range(16)):
        """Evaluate Albrecht A term.

        Returns
        -------
        Full Albrecht A matrix element. Unit: e^2 Angstrom^2 / eV
        """
        self.read()

        self.timer.start('AlbrechtA')

        if not hasattr(self, 'fco'):
            self.fco = FranckCondonOverlap()

        om = omega + 1j * gamma
        # excited state forces
        F_pr = self.exF_rp.T

        m_Qcc = np.zeros((self.ndof, 3, 3), dtype=complex)
        for p, energy in enumerate(self.ex0E_p):
            S_Q = self.Huang_Rhys_factors(F_pr[p])
            energy_Q = energy - self.om_Q * S_Q
            me_cc = np.outer(self.ex0m_pc[p], self.ex0m_pc[p].conj())

            wm_Q = np.zeros((self.ndof), dtype=complex)
            wp_Q = np.zeros((self.ndof), dtype=complex)
            for m in ml:
                self.timer.start('0mm1')
                fco_Q = self.fco.direct0mm1(m, S_Q)
                self.timer.stop('0mm1')
                
                self.timer.start('weight_Q')
                wm_Q += fco_Q / (energy_Q + m * self.om_Q - om)
                wp_Q += fco_Q / (energy_Q + (m - 1) * self.om_Q + om)
                self.timer.stop('weight_Q')
            self.timer.start('einsum')
            m_Qcc += np.einsum('a,bc->abc', wm_Q, me_cc)
            m_Qcc += np.einsum('a,bc->abc', wp_Q, me_cc.T.conj())
            self.timer.stop('einsum')
                
        self.timer.stop('AlbrechtA')
        return m_Qcc

    def meBC(self, omega, gamma=0.1, ml=range(16),
             term='BC'):
        """Evaluate Albrecht BC term.

        Returns
        -------
        Full Albrecht BC matrix element.
        Unit: e^2 Angstrom / eV / sqrt(amu)
        """
        self.read()

        self.timer.start('AlbrechtBC')

        if not hasattr(self, 'fco'):
            self.fco = FranckCondonOverlap()

        omL = omega + 1j * gamma
        omS_Q = omL - self.om_Q

        # excited state forces
        F_pr = self.exF_rp.T
        # derivatives after normal coordinates
        dmdq_qpc = (self.exdmdr_rpc.T * self.im).T  # unit e / sqrt(amu)
        dmdQ_Qpc = np.dot(dmdq_qpc.T, self.modes.T).T # unit e / sqrt(amu)
##        print('dmdQ_Qpc', dmdQ_Qpc[-1])

        me_Qcc = np.zeros((self.ndof, 3, 3), dtype=complex)
        for p, energy in enumerate(self.ex0E_p):
            S_Q = self.Huang_Rhys_factors(F_pr[p])
            # relaxed excited state energy
##            n_vQ = np.where(self.n_vQ > 0, 1, 0)
##            energy_v = energy - n_vQ.dot(self.om_Q * S_Q)
            energy_Q = energy - self.om_Q * S_Q
            
            ##me_cc = np.outer(self.ex0m_pc[p], self.ex0m_pc[p].conj())
            m_c = self.ex0m_pc[p]
            dmdQ_Qc = dmdQ_Qpc[:, p]

            wBLS_Q = np.zeros((self.ndof), dtype=complex)
            wBSL_Q = np.zeros((self.ndof), dtype=complex)
            wCLS_Q = np.zeros((self.ndof), dtype=complex)
            wCSL_Q = np.zeros((self.ndof), dtype=complex)
            for m in ml:
                self.timer.start('0mm1/2')
                f0mmQ1_Q = (self.fco.directT0(m, S_Q) +
                            np.sqrt(2) * self.fco.direct0mm2(m, S_Q))
                f0Qmm1_Q = self.fco.direct(1, m, S_Q)
##                if (self.n_vQ > 1).any():
##                    fco2_Q = self.fco.direct0mm2(m, S_Q)
                self.timer.stop('0mm1/2')
                
                self.timer.start('weight_Q')
                em_Q = energy_Q + m * self.om_Q
                wBLS_Q += f0mmQ1_Q / (em_Q - omL)
                wBSL_Q += f0Qmm1_Q / (em_Q - omL)
                wCLS_Q += f0mmQ1_Q / (em_Q + omS_Q)
                wCSL_Q += f0Qmm1_Q / (em_Q + omS_Q)
                self.timer.stop('weight_Q')
            self.timer.start('einsum')
            # unit e^2 Angstrom / sqrt(amu)
            mdmdQ_Qcc = np.einsum('a,bc->bac', m_c, dmdQ_Qc.conj())
            dmdQm_Qcc = np.einsum('ab,c->abc', dmdQ_Qc, m_c.conj())
            if 'B' in term:
                me_Qcc += np.multiply(wBLS_Q, mdmdQ_Qcc.T).T
                me_Qcc += np.multiply(wBSL_Q, dmdQm_Qcc.T).T
            if 'C' in term:
                me_Qcc += np.multiply(wCLS_Q, mdmdQ_Qcc.T).T
                me_Qcc += np.multiply(wCSL_Q, dmdQm_Qcc.T).T
            self.timer.stop('einsum')
                
        self.timer.stop('AlbrechtBC')
        return me_Qcc  # unit e^2 Angstrom / eV / sqrt(amu)

    def electronic_me_Qcc(self, omega, gamma):
        """Evaluate an electronic matric element."""
        if self.approximation.lower() == 'albrecht a':
            Vel_Qcc = self.meA(omega, gamma)  # e^2 Angstrom^2 / eV 
##            print('A Vel_Qcc=', Vel_Qcc[-1])
            # divide through pre-factor
            with np.errstate(divide='ignore'):
                Vel_Qcc *= np.where(self.vib01_Q > 0,
                                    1. / np.sqrt(self.vib01_Q), 0)[:, None, None]
            # -> e^2 Angstrom / eV / sqrt(amu)
        elif self.approximation.lower() == 'albrecht bc':
            Vel_Qcc = self.meBC(omega, gamma)  # e^2 Angstrom / eV / sqrt(amu)
##            print('BC Vel_Qcc=', Vel_Qcc[-1])
        elif self.approximation.lower() == 'albrecht b':
            Vel_Qcc = self.meBC(omega, gamma, term='B')
        elif self.approximation.lower() == 'albrecht c':
            Vel_Qcc = self.meBC(omega, gamma, term='C')
        else:
            raise NotImplementedError(
                'Approximation {0} not implemented. '.format(
                    self.approximation) +
                'Please use "Albrecht A/B/C".')

        Vel_Qcc *= u.Hartree * u.Bohr  # e^2 Angstrom^2 / eV -> Angstrom^3

        return Vel_Qcc  # Angstrom^2 / sqrt(amu)

    def matrix_element(self, omega, gamma):
        self.read()
        V_Qcc = np.zeros((self.ndof, 3, 3), dtype=complex)
        if self.approximation.lower() in ['profeta', 'placzek',
                                          'p-p', 'placzekalpha']:
            me_Qcc = self.electronic_me_Qcc(omega, gamma)
            for Q, vib01 in enumerate(self.vib01_Q):
                V_Qcc[Q] = me_Qcc[Q] * vib01
        elif self.approximation.lower() == 'albrecht a':
            V_Qcc += self.me_AlbrechtA(omega, gamma)
        elif self.approximation.lower() == 'albrecht b':
            V_Qcc += self.me_AlbrechtBC(omega, gamma, term='B')
        elif self.approximation.lower() == 'albrecht c':
            V_Qcc += self.me_AlbrechtBC(omega, gamma, term='C')
        elif self.approximation.lower() == 'albrecht bc':
            V_Qcc += self.me_AlbrechtBC(omega, gamma)
        elif self.approximation.lower() == 'albrecht':
            V_Qcc += self.me_AlbrechtA(omega, gamma)
            V_Qcc += self.me_AlbrechtBC(omega, gamma)
        elif self.approximation.lower() == 'albrecht+profeta':
            raise NotImplementedError('not useful')
            V_Qcc += self.get_matrix_element_AlbrechtA(omega, gamma)
            V_Qcc += self.get_matrix_element_Profeta(omega, gamma)
        else:
            raise NotImplementedError(
                'Approximation {0} not implemented. '.format(
                    self.approximation) +
                'Please use "Albrecht A/B/C/BC" ' +
                'or "Albrecht".')

        return V_Qcc
    
    def summary(self, omega=0, gamma=0,
                method='standard', direction='central',
                log=sys.stdout):
        """Print summary for given omega [eV]"""
        intensities = self.absolute_intensity(omega, gamma)[self.skip:]

        if isinstance(log, str):
            log = paropen(log, 'a')

        parprint('-------------------------------------', file=log)
        parprint(' excitation at ' + str(omega) + ' eV', file=log)
        parprint(' gamma ' + str(gamma) + ' eV', file=log)
        parprint(' approximation:', self.approximation, file=log)
        parprint(' Mode    Frequency        Intensity', file=log)
        parprint('  #    meV     cm^-1      [A^4/amu]', file=log)
        parprint('-------------------------------------', file=log)
        for n, e in enumerate(self.om_v):
            if e.imag != 0:
                c = 'i'
                e = e.imag
            else:
                c = ' '
                e = e.real
            parprint('%3d %6.1f   %7.1f%s  %9.1f' %
                     (n, 1000 * e, e / u.invcm, c, intensities[n]),
                     file=log)
        parprint('-------------------------------------', file=log)
        parprint('Zero-point energy: %.3f eV' % self.get_zero_point_energy(),
                 file=log)

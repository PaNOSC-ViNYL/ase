# -*- coding: utf-8 -*-

from __future__ import print_function, division
import numpy as np

import ase.units as u
from ase.vibrations.resonant_raman import ResonantRaman

# XXX remove gpaw dependence
from gpaw.lrtddft.spectrum import polarizability


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
        pre *= u.Hartree * u.Bohr  # e^2Angstrom^2 / eV -> Angstrom^3

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

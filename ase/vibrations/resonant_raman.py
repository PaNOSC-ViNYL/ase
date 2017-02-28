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
from ase.utils.timing import Timer
from ase.utils import convert_string_to_fd


def overlap(calc1, calc2):
    """GPAW overlaps"""
    n1 = calc1.get_number_of_bands()
#    spin1 = calc1.get_number_of_spins() == 2
    n2 = calc2.get_number_of_bands()
#    spin2 = calc2.get_number_of_spins() == 2
    overlap_nn = np.zeros((n1, n2), dtype=float)
    for i1 in range(n1):
        psi1 = calc1.wfs.kpt_u[0].psit_nG[i1]
        norm1 = calc1.wfs.gd.integrate(psi1.conj() * psi1)
        for i2 in range(n2):
            psi2 = calc2.wfs.kpt_u[0].psit_nG[i2]
            norm2 = calc2.wfs.gd.integrate(psi2.conj() * psi2)
            overlap_nn[i1, i2] = (calc1.wfs.gd.integrate(psi1.conj() * psi2) /
                                  np.sqrt(norm1 * norm2))
            if 0 and abs(overlap_nn[i1, i2]) > 0.5:
                parprint(i1, i2,
                         '{0:7.2f} ({1:4.2f} {2:4.2f})'.format(
                             overlap_nn[i1, i2], norm1, norm2))
    return overlap_nn
            

class ResonantRaman(Vibrations):
    """Resonant Raman intensities using finite differences."""

    def __init__(self, atoms, Excitations,
                 indices=None,
                 gsname='rraman',  # name for ground state calculations
                 exname=None,      # name for excited state calculations
                 delta=0.01,
                 nfree=2,
                 directions=None,
                 observation={'geometry': '-Z(XX)Z'},
                 form='v',         # form of the dipole operator
                 exkwargs={},      # kwargs to be passed to Excitations
                 exext='.ex.gz',   # extension for Excitation names
                 txt='-',
                 verbose=False,
                 overlap=False,
                 minoverlap=0.1):
        """
        Parameters
        ----------
        atoms: ase Atoms object
        Excitations: class
            Type of the excitation list object. The class object is
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
        indices: list
        gsname: string
            name for ground state calculations
        exname: string
            name for excited state calculations
        delta: float
            Finite difference displacement in Angstrom.
        nfree: float
        directions:
        approximation: string
            Level of approximation used.
        observation: dict
            Polarization settings
        form: string
            Form of the dipole operator, 'v' for velocity form (default)
            and 'r' for length form.
        exkwargs: dict
            Arguments given to the Excitations objects in reading.
        exext: string
            Extension for filenames of Excitation lists.
        txt:
            Output stream
        verbose:
            Verbosity level of output
        overlap: bool/float
            Use wavefunction overlaps.
        minoverlap: float
            Minimal absolute overlap to consider. Defaults to 0.1 to avoid
            numerical garbage.
        """
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

        self.observation = observation
        self.exobj = Excitations
        self.exkwargs = exkwargs
        self.dipole_form = form

        self.timer = Timer()
        self.txt = convert_string_to_fd(txt)

        self.verbose = verbose
        self.overlap = overlap
        self.minoverlap = minoverlap
        
    @property
    def approximation(self):
        return self._approx

    @approximation.setter
    def approximation(self, value):
        self.check_approximation(value)

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
            fname = self.exname + '.eq.gpw'
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
            self.timer.start('Overlap')
            ov_nn = overlap(self.atoms.get_calculator(), self.eq_calculator)
            if rank == 0:
                np.save(filename + '.ov', ov_nn)
            self.timer.stop('Overlap')
        self.timer.stop('Ground state')

        self.timer.start('Excitations')
        basename, _ = os.path.splitext(filename)
        excitations = self.exobj(
            self.atoms.get_calculator(), **self.exkwargs)
        excitations.write(basename + self.exext)
        self.timer.stop('Excitations')

    def read_excitations(self):
        """Read all finite difference excitations and select matching."""
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
            [ex.get_dipole_me(form=self.dipole_form) for ex in ex0]) *
            u.Bohr)
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

    def read_excitations_overlap(self):
        """Read all finite difference excitations and wf overlaps."""
        self.timer.start('read excitations')
        self.timer.start('really read')
        self.log('reading ' + self.exname + '.eq' + self.exext)
        ex0 = self.exobj(self.exname + '.eq' + self.exext,
                         **self.exkwargs)
        self.timer.stop('really read')

        def load(name, pm):
            self.log('reading ' + name + pm + self.exext)
            ex_p = self.exobj(name + pm + self.exext, **self.exkwargs)
            self.log('reading ' + name + pm + '.pckl.ov.npy')
            ov_nn = np.load(name + pm + '.pckl.ov.npy')
            # remove numerical garbage
            ov_nn = np.where(np.abs(ov_nn) > self.minoverlap, ov_nn, 0)
            ov_pp = ex_p.overlap(ov_nn)
#            print(ov_pp)
            return ex_p, ov_pp
            
        exm = []
        ovm = []
        exp = []
        ovp = []
        for a in self.indices:
            for i in 'xyz':
                name = '%s.%d%s' % (self.exname, a, i)
                ex, ov = load(name, '-')
                exm.append(ex)
                ovm.append(ov)
                ex, ov = load(name, '+')
                exp.append(ex)
                ovp.append(ov)
        self.ndof = 3 * len(self.indices)
        self.timer.stop('read excitations')

        self.timer.start('me and energy')

        eu = u.Hartree
        self.ex0E_p = np.array([ex.energy * eu for ex in ex0])
        self.ex0m_pc = (np.array(
            [ex.get_dipole_me(form=self.dipole_form) for ex in ex0]) *
            u.Bohr)

        def rotate(ex_p, ov_pp):
            em_p = np.array(
                [ex.get_dipole_me(form=self.dipole_form) for ex in ex_p])
            return ov_pp.dot(em_p)

        def z(arr):
            return np.where(abs(arr) > 0.05, arr, 0)

        exmE_rp = []
        expE_rp = []
        exF_rp = []
        exmm_rpc = []
        expm_rpc = []
        exdmdr_rpc = []
        r = 0
        for a in self.indices:
            for i in 'xyz':
                exmE_rp.append([em.energy for em in exm[r]])
                expE_rp.append([ep.energy for ep in exp[r]])
                exF_rp.append(
                    [(em.energy - ep.energy)
                     for ep, em in zip(exp[r], exm[r])])
                exmm_rpc.append(rotate(exm[r], ovm[r]))
                expm_rpc.append(rotate(exp[r], ovp[r]))
#                print('m,p=', z(expm_rpc[-1]), z(exmm_rpc[-1]))
                exdmdr_rpc.append(expm_rpc[-1] - exmm_rpc[-1])
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
        # matrix element derivatives in e
        self.exdmdr_rpc = np.array(exdmdr_rpc) * u.Bohr / 2 / self.delta

        self.timer.stop('me and energy')

    def read(self, method='standard', direction='central'):
        """Read data from a pre-performed calculation."""

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
                                    1. / np.sqrt(2 * self.om_Q), 0)
        # -> sqrt(amu) * Angstrom
        self.vib01_Q *= np.sqrt(u.Ha * u._me / u._amu) * u.Bohr
        self.timer.stop('read vibrations')

        if not hasattr(self, 'ex0E_p'):
            if self.overlap:
                self.read_excitations_overlap()
            else:
                self.read_excitations()

    def me_Qcc(self, omega, gamma):
        """Full matrix element

        Returns
        -------
        Matrix element in e^2 Angstrom^2 / eV
        """
        # Angstrom^2 / sqrt(amu)
        elme_Qcc = self.electronic_me_Qcc(omega, gamma)
        # Angstrom^3 -> e^2 Angstrom^2 / eV
        elme_Qcc /= u.Hartree * u.Bohr  # e^2 Angstrom / eV / sqrt(amu)
        return elme_Qcc * self.vib01_Q[:, None, None]

    def intensity(self, omega, gamma=0.1):
        """Raman intensity

        Returns
        -------
        unit e^4 Angstrom^4 / eV^2"""
        m2 = ResonantRaman.m2
        alpha_Qcc = self.me_Qcc(omega, gamma)
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
        alpha2_r = m2(alpha_Qcc[:, 0, 0] + alpha_Qcc[:, 1, 1] +
                      alpha_Qcc[:, 2, 2]) / 9.
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
        I_r = self.intensity(omega, gamma)
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
                     (n, 1000 * e, c, e / u.invcm, c, intensities[n]),
                     file=log)
        parprint('-------------------------------------', file=log)
        parprint('Zero-point energy: %.3f eV' % self.get_zero_point_energy(),
                 file=log)

    def __del__(self):
        self.timer.write(self.txt)


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
            [ex.get_dipole_me(form=self.dipole_form) for ex in ex0]) *
            u.Bohr)
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

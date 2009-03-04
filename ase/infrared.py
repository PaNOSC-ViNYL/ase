
# -*- coding: utf-8 -*-

"""Infrared spectroscopy"""

import pickle
from math import sin, pi, sqrt, exp, log
from os import remove
from os.path import isfile

import numpy as np

import ase.units as units
from ase.io.trajectory import PickleTrajectory
from ase.parallel import rank, barrier
from ase.vibrations import Vibrations

# Atomic mass unit to rest mass of electron
amutome = 1822.888485
# Conversion factors from atomic units to (D/Angstrom)^2/amu
# and electron x Angstrom to Debye
conv = 42055.44331
qeangtod = 4.803204272

class InfraRed(Vibrations):
    """Class for calculating vibrational modes and infrared intensities
    using finite difference.

    The vibrational modes are calculated from a finite difference
    approximation of the Hessian matrix and the IR intensities from
    a finite difference approximation of the gradient of the dipole
    momen. The method is described in:

      D. Porezag, M. R. Peterson:
      "Infrared intensities and Raman-scattering activities within
      density-functional theory",
      Phys. Rev. B 54, 7830 (1996)

    The calculator object (calc) linked to the Atoms object (atoms) must have the 
    attribute: calc.get_dipole_moment(atoms)

    The *summary*, *get_energies()* and *get_frequencies()*
    methods all take an optional *method* keyword.  Use
    method='Frederiksen' to use the method described in:

      T. Frederiksen, M. Paulsson, M. Brandbyge, A. P. Jauho:
      "Inelastic transport theory from first-principles: methodology
      and applications for nanoscale devices", 
      Phys. Rev. B 75, 205413 (2007) 

    atoms: Atoms object
        The atoms to work on.
    indices: list of int
        List of indices of atoms to vibrate.  Default behavior is
        to vibrate all atoms.
    name: str
        Name to use for files.
    delta: float
        Magnitude of displacements.
    nfree: int
        Number of displacements per degree of freedom, 2 or 4 are
        supported. Default is 2 which will displace each atom +delta
        and -delta for each cartesian coordinate.
    directions: list of int
        Cartesian coordinates to calculate the gradient of the dipole moment in. 
        For example directions = 2 only dipole moment in the z-direction will
        be considered, whereas for directions = [0, 1] only the dipole
        moment in the xy-plane will be considered. Default behavior is to
        consider the dipole moment in all directions.

    """
    def __init__(self, atoms, indices=None, name='ir', delta=0.01, nfree=2, directions=None):
        assert nfree in [2, 4]
        self.atoms = atoms
        if atoms.constraints:
            print "WARNING! \n Your Atoms object is constrained. Some forces may be unintended set to zero. \n"
        self.calc = atoms.get_calculator()
        if indices is None:
            indices = range(len(atoms))
        self.indices = np.asarray(indices)
        self.nfree = nfree
        self.name = name+'-d%.3f' % delta
        self.delta = delta
        self.H = None
        if directions is None:
            self.directions = np.asarray([0, 1, 2])
        else:
            self.directions = np.asarray(directions)
        self.ir = True

    def read(self, method='standard', direction='central'):
        self.method = method.lower()
        self.direction = direction.lower()
        assert self.method in ['standard', 'frederiksen']
        if direction != 'central':
            raise NotImplementedError('Only central difference is implemented at the moment.')

        # Get "static" dipole moment and forces
        name = '%s.eq.pckl' % self.name
        [forces_zero, dipole_zero] = pickle.load(open(name))
        self.dipole_zero = (sum(dipole_zero**2)**0.5)*qeangtod
        self.force_zero = max([sum((forces_zero[j])**2)**0.5 for j in self.indices])

        ndof = 3 * len(self.indices)
        H = np.empty((ndof, ndof))
        dpdx = np.empty((ndof, 3))
        r = 0
        for a in self.indices:
            for i in 'xyz':
                name = '%s.%d%s' % (self.name, a, i)
                [fminus, dminus] = pickle.load(open(name + '-.pckl'))
                [fplus, dplus] = pickle.load(open(name + '+.pckl'))
                if self.nfree == 4:
                    [fminusminus, dminusminus] = pickle.load(open(name + '--.pckl'))
                    [fplusplus, dplusplus] = pickle.load(open(name + '++.pckl'))
                if self.method == 'frederiksen':
                    fminus[a] += -fminus.sum(0)
                    fplus[a] += -fplus.sum(0)
                    if self.nfree == 4:
                        fminusminus[a] += -fminus.sum(0)
                        fplusplus[a] += -fplus.sum(0)
                if self.nfree == 2:
                    H[r] = (fminus - fplus)[self.indices].ravel() / 2.0
                    dpdx[r] = (dminus - dplus)
                if self.nfree == 4:
                    H[r] = (-fminusminus+8*fminus-8*fplus+fplusplus)[self.indices].ravel() / 12.0
                    dpdx[r] = (-dplusplus + 8*dplus - 8*dminus +dminusminus) / 6.0
                H[r] /= 2 * self.delta
                dpdx[r] /= 2 * self.delta
                for n in range(3):
                    if n not in self.directions:
                        dpdx[r][n] = 0
                        dpdx[r][n] = 0
                r += 1
        # Calculate eigenfrequencies and eigenvectors
        m = self.atoms.get_masses()
        H += H.copy().T
        self.H = H
        m = self.atoms.get_masses()
        self.im = np.repeat(m[self.indices]**-0.5, 3)
        omega2, modes = np.linalg.eigh(self.im[:, None] * H * self.im)
        self.modes = modes.T.copy()

        # Calculate intensities
        dpdq = np.array([dpdx[j]/sqrt(m[self.indices[j/3]]*amutome) for j in range(ndof)])
        dpdQ = np.dot(dpdq.T, modes)
        dpdQ = dpdQ.T
        intensities = np.array([sum(dpdQ[j]**2) for j in range(ndof)])
        # Conversion factor:
        s = units._hbar * 1e10 / sqrt(units._e * units._amu)
        self.hnu = s * omega2.astype(complex)**0.5
        self.intensities = intensities*conv

    def summary(self, method='standard', direction='central'):
        hnu = self.get_energies(method, direction)
        s = 0.01 * units._e / units._c / units._hplanck
        print '-------------------------------------'
        print ' Mode    Frequency        Intensity'
        print '  #    meV     cm^-1   (D/Å)^2 amu^-1'
        print '-------------------------------------'
        for n, e in enumerate(hnu):
            if e.imag != 0:
                c = 'i'
                e = e.imag
            else:
                c = ' '
            print '%3d %6.1f%s  %7.1f%s  %9.4f' % (n, 1000 * e, c, s * e, c, self.intensities[n])
        print '-------------------------------------'
        print 'Zero-point energy: %.3f eV' % self.get_zero_point_energy()
        print 'Static dipole moment: %.3f D' % self.dipole_zero
        print 'Maximum force on atom in eqiulibrium: %.4f eV/Å' % self.force_zero
        print

    def write_spectra(self, out='ir-spectra.dat', start=800, end=4000, npts=None, width=4, type='Gaussian', method='standard', direction='central'):
        """Write out infrared spectrum to file.

        Start and end point, and width of the Gaussian/Lorentzian should be given in cm^-1."""
        self.type = type.lower()
        assert self.type in ['gaussian', 'lorentzian']
        if not npts:
            npts = (end-start)/width*5+1
        frequencies = self.get_frequencies(method, direction).real
        intensities=self.intensities
        if type == 'lorentzian':
            lineshape = 1
            intensities = intensities*width*pi/2.
        else:
            lineshape = 0
            sigma = width/2./sqrt(2.*log(2.))
        #Make array with spectrum data
        spectrum=np.zeros(npts,np.float)
        ediff = (end-start)/float(npts-1)
        for i in range(1,npts):
            energy = end - float(i)*ediff
            for j in range(len(frequencies)):
                if lineshape:
                    spectrum[i] = spectrum[i]+intensities[j]*0.5*width/pi/((energy-frequencies[j])**2+0.25*width**2)
                else:
                    spectrum[i] = spectrum[i]+intensities[j]*exp(-(energy-frequencies[j])**2/2./sigma**2)
        #Write out spectrum in file. First column is just intensities. 
        #Second column is absorbance scaled so that data runs from 1 to 0
        spectrumfile = open(out,"w")
        for i in range(1,npts):
            energy = end - float(i)*ediff
            spectrumfile.write("%f %15.5e %15.5e\n" % (energy,spectrum[i],1.-spectrum[i]/spectrum.max()))
        spectrumfile.close()

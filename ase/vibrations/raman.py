class Raman(Vibrations):
    """Base class for calculating vibrational modes and
    Raman intensities using finite difference.

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
                gives the dipole matrix element in |e| * Angstrom
            ex.energy:
                is the transition energy in Hartrees
    """
    def __init__(self, atoms,
                 indices=None,
                 name='raman',
                 delta=0.01,
                 nfree=2,
                 directions=None,
                 observation={'geometry': '-Z(XX)Z'},
                 txt='-',
                 verbose=False,):
        assert(nfree == 2)
        Vibrations.__init__(self, atoms, indices, name, delta, nfree)
        
        if directions is None:
            self.directions = np.array([0, 1, 2])
        else:
            self.directions = np.array(directions)

        self.observation = observation

        self.timer = Timer()
        self.txt = convert_string_to_fd(txt)

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

    def get_intensities(self, omega, gamma=0.1):
        """Calculate intensities depending on polarization settings""" 
        m2 = self.m2
        alpha_rcc = self.get_matrix_element(omega, gamma)
        if not self.observation:  # XXXX remove
            """Simple sum, maybe too simple"""
            return m2(alpha_rcc).sum(axis=1).sum(axis=1)
        # XXX enable when appropriate
        #        if self.observation['orientation'].lower() != 'random':
        #            raise NotImplementedError('not yet')

        # pre-factor from vibrational overlaps
        pre_r = 1. / self.hnu

        # random orientation of the molecular frame
        # Woodward & Long,
        # Guthmuller, J. J. Chem. Phys. 2016, 144 (6), 64106
        m2 = ResonantRaman.m2
        alpha2_r = m2(alpha_rcc[:, 0, 0] + alpha_rcc[:, 1, 1] +
                      alpha_rcc[:, 2, 2]) / 9.
        delta2_r = 3 / 4. * (
            m2(alpha_rcc[:, 0, 1] - alpha_rcc[:, 1, 0]) +
            m2(alpha_rcc[:, 0, 2] - alpha_rcc[:, 2, 0]) +
            m2(alpha_rcc[:, 1, 2] - alpha_rcc[:, 2, 1]))
        gamma2_r = (3 / 4. * (m2(alpha_rcc[:, 0, 1] + alpha_rcc[:, 1, 0]) +
                              m2(alpha_rcc[:, 0, 2] + alpha_rcc[:, 2, 0]) +
                              m2(alpha_rcc[:, 1, 2] + alpha_rcc[:, 2, 1])) +
                    (m2(alpha_rcc[:, 0, 0] - alpha_rcc[:, 1, 1]) +
                     m2(alpha_rcc[:, 0, 0] - alpha_rcc[:, 2, 2]) +
                     m2(alpha_rcc[:, 1, 1] - alpha_rcc[:, 2, 2])) / 2)

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

    def get_cross_sections(self, omega, gamma=0.1):
        I_r = self.get_intensities(omega, gamma)
        pre = 1. / 16 / np.pi**2 / units.eps0**2 / units.c**4
        # frequency of scattered light
        omS_r = omega - self.hnu
        return pre * omega * omS_r**3 * I_r

    def get_spectrum(self, omega, gamma=0.1,
                     start=200.0, end=4000.0, npts=None, width=4.0,
                     type='Gaussian', method='standard', direction='central',
                     intensity_unit='????', normalize=False):
        """Get Raman spectrum.

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
        parprint(' approximation:', self.approximation, file=log)
        parprint(' observation:', self.observation, '\n', file=log)
        parprint(' Mode    Frequency        Intensity', file=log)
        parprint('  #    meV     cm^-1      [e^4A^4/eV^2]', file=log)
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


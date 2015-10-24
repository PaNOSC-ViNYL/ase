from __future__ import print_function
"""This module defines an ASE interface to SIESTA.

http://www.uam.es/departamentos/ciencias/fismateriac/siesta
"""
import os
import sys
from os.path import join, isfile, islink, getmtime
from cmath import exp
import array
import string
import numpy as np

from ase.data import chemical_symbols, atomic_numbers
from ase.units import Rydberg, fs, Bohr
from ase.io.siesta import read_rho, read_fdf, read_struct_out
from ase.io.cube import read_cube_data
from ase.calculators.siesta.basis_set import BasisSet, DZP
from ase.calculators.calculator import FileIOCalculator, all_changes
from ase.calculators.calculator import LockedParameters

class SpinType(object):pass
class UNPOLARIZED(SpinType):
    @classmethod
    def write_fdf(cls, f):
        f.write(format_fdf('SpinPolarized', False))

class COLLINEAR(SpinType):
    @classmethod
    def write_fdf(cls, f):
        f.write(format_fdf('SpinPolarized', True))

class FULL(SpinType):
    @classmethod
    def write_fdf(cls, f):
        f.write(format_fdf('SpinPolarized', True))
        f.write(format_fdf('NonCollinearSpin', True))

class SiestaXX:
    """Class for doing SIESTA calculations.

    The default parameters are the one of the SIESTA
    Fortran code use.

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Use the set_siesta_param method from the inp.siesta
    class to set the FDF parameters:
      calc.MS.set_siesta_param('PAO.EnergyShift', 0.01, unit='Ry')
    """
    def __init__(self,
                 label='siesta',
                 ghosts=[],
                 mpirun=False,
                 np=1,
                 write_fdf=True,
                 ):
        """Construct SIESTA-calculator object.

        Parameters
        ==========
        label: str
            Prefix to use for filenames (label.fdf, label.txt, ...).
            Default is 'siesta'.
        write_fdf: bool
            Use write_fdf=False to use your own fdf-file.
        ghost: list
        mpirun: bool
          run siesta in parallel
        np: integer
          number of processeur to use if mpirun =True

        Examples
        ========
        Use default values:

        >>>from ase.calculators.siesta import Siesta
        >>>from ase import Atoms
        >>>import numpy as np

        >>>bud = Atoms('CH4', np.array([[0.000000,  0.000000,  0.000000],
                                        [0.682793,  0.682793,  0.682793],
                                        [-0.682793, -0.682793,  0.68279],
                                        [-0.682793,  0.682793, -0.682793],
                                        [0.682793, -0.682793, -0.682793]]))
        >>>calc = Siesta(label='ch4-d')
        >>>calc.MS.set_siesta_param('SystemName', 'METHANE')
        >>>calc.MS.set_siesta_param('PAO.BasisSize', 'SZ')
        >>>calc.MS.set_siesta_param('LatticeConstant', 10, unit='Ang')
        >>>calc.MS.set_siesta_param('AtomicCoordinatesFormat', 'Ang')
        >>>calc.MS.set_siesta_param('AtomCoorFormatOut', 'Ang')
        >>>calc.MS.set_siesta_param('AtomicCoordinatesOrigin',
                                      np.array([0.127, 0.745, -0.33]))
        >>>calc.MS.set_siesta_param('XC.Functional', 'LDA')
        >>>calc.MS.set_siesta_param('XC.Authors', 'CA')
        >>>calc.MS.set_siesta_param('MeshCutoff', 30, unit='Ry')
        >>>calc.MS.set_siesta_param('MaxSCFIterations', 50)
        >>>calc.MS.set_siesta_param('DM.MixingWeight', 0.15)
        >>>calc.MS.set_siesta_param('DM.NumberPulay', 3)
        >>>calc.MS.set_siesta_param('DM.Tolerance', 1E-4)
        >>>calc.MS.set_siesta_param('SolutionMethod', 'diagon')
        >>>calc.MS.set_siesta_param('ElectronicTemperature', 25, unit = 'meV')
        >>>calc.MS.set_siesta_param('LongOutput', False)
        >>>bud.set_calculator(calc)
        >>>bud.get_potential_energy()
        """
        self.name = 'Siesta'
        self.label = label# ################### != out
        self.ghosts = ghosts
        # self.write_tddft_file = write_tddft
        # self.tddft_label = tddft_label#################### != out
        # self.run_calc = run_calc
        self.MS = SiestaParameters()
        # self.MS.fname = xyz_file
        self.write_fdf_file = write_fdf
        # self.run_raman = run_raman
        # self.tddft_inp = inp.tddft_inp()
        # self.Raman = Raman(self.MS)
        self.atoms = None
        self.mpirun = mpirun
        self.np = np

        # self.data_path = './'
        # self.pol_name_re = 'dipol_inter_iter_krylov_re.txt'
        # self.pol_name_im = 'dipol_inter_iter_krylov_im.txt'
        # self.fname_pol_re = self.data_path + self.pol_name_re
        # self.fname_pol_im = self.data_path + self.pol_name_im

        self.converged = False
        self.fdf = {}
        self.e_fermi = None
        self.results = {}  # calculated properties (energy, forces, ...)

        # self.nb_step = 0
        # self.get_old_calc = get_old_calc
        # self.siesta_run = False
        # self.tddft_run = False

#    def prep_inp(self, atoms, raman=False):
#        import pyiter.Read_data as RD
#
#        folder = 'calc_{0}'.format(self.nb_step)
#        RD.run('mkdir ' + folder)
#        self.write_fdf(atoms)
#        print(self.MS.list_param)
#        if raman:
#          self.tddft_inp.write_tddft_inp()
#
#        RD.run('cp * ' + folder)
#        RD.run('rm ' + self.MS.param['SystemLabel'] + '.*')
#
#        self.nb_step = self.nb_step + 1
#
#
    def update(self, atoms):
        self.atoms = atoms
        if (not self.converged or
            len(self.numbers) != len(atoms) or
            (self.numbers != atoms.get_atomic_numbers()).any()):
            self.initialize(atoms)
            self.calculate(atoms)

        elif ((self.positions != atoms.get_positions()).any() or
              (self.pbc != atoms.get_pbc()).any() or
              (self.cell != atoms.get_cell()).any()):
            self.calculate(atoms)

    def initialize(self, atoms):
        self.numbers = atoms.get_atomic_numbers().copy()
        self.species = []
        for a, Z in enumerate(self.numbers):
            if a in self.ghosts:
                Z = -Z
            if Z not in self.species:
                self.species.append(Z)

        if 'SIESTA_PP_PATH' in os.environ:
            pppaths = os.environ['SIESTA_PP_PATH'].split(':')
        else:
            pppaths = []

        for Z in self.species:
            symbol = chemical_symbols[abs(Z)]
            name = symbol + '.vps'
            name1 = symbol + '.psf'
            found = False
            for path in pppaths:
                filename = join(path, name)
                filename1 = join(path, name1)
                if isfile(filename) or islink(filename):
                    found = True
                    if path != '.':
                        if islink(name) or isfile(name):
                            os.remove(name)
                        os.symlink(filename, name)

                elif isfile(filename1) or islink(filename1):
                    found = True
                    if path != '.':
                        if islink(name1) or isfile(name1):
                            os.remove(name1)
                        os.symlink(filename1, name1)
            if not found:
                raise RuntimeError('No pseudopotential for %s!' % symbol)

        self.converged = False

    def get_potential_energy(self, atoms, force_consistent=False):
        self.update(atoms)

        if force_consistent:
            return self.efree
        else:
            # Energy extrapolated to zero Kelvin:
            return  (self.etotal + self.efree) / 2

    def get_forces(self, atoms):
        self.update(atoms)
        return self.forces.copy()

    def get_stress(self, atoms):
        self.update(atoms)
        return self.stress.copy()

    def get_dipole_moment(self, atoms):
        """Returns total dipole moment of the system."""
        self.update(atoms)
        return self.dipole

    def get_polarizability(self, atoms):
        """Returns total dipole moment of the system."""
        self.update(atoms)
        return self.polarizability

    def get_frequencies(self, atoms):
        """Returns frequency use for the TDDFT calculations."""
        self.update(atoms)
        return self.frequencies

    def read_dipole(self):
        dipolemoment = np.zeros([1, 3])
        for line in open(self.label + '.txt', 'r'):
            if line.rfind('Electric dipole (Debye)') > -1:
                dipolemoment = np.array([float(f) for f in line.split()[5:8]])
        # debye to e*Ang (the units of VASP)
        dipolemoment = dipolemoment * 0.2081943482534
        return dipolemoment

#    def read_polarizability(self):
#      P = np.loadtxt(self.fname_pol_re)[:, 2:11] + complex(0.0, 1.0)*
# np.loadtxt(self.fname_pol_im)[:, 2:11]*(Bohr**3) # convert to Ang**3
#      return P
#
#    def read_frequencies(self):
#      freq = np.loadtxt(self.fname_pol_re)[:, 0]
#      return freq

    def get_pseudo_density(self, spin=None, pad=True):
        """Return pseudo-density array.

        If *spin* is not given, then the total density is returned.
        Otherwise, the spin up or down density is returned (spin=0 or 1).
        """
        filename = self.label + '.RHO'
        if not isfile(filename):
            raise RuntimeError('Could not find rho-file (make sure '
            'to add fdf-option "SaveRho=True" to your calculation)')

        rho = read_rho(filename)

        if spin is None:
            return rho.sum(axis=3)
        elif rho.shape[3] != 2:
            raise RuntimeError('Explicit spin-value requested. '
                               'Only total density is available.')
        elif spin == 0 or spin == 1:
            return rho[:,:,:, spin]
        else:
            raise RuntimeError('Invalid spin-value requested. '
                               'Expected 0 or 1, got %s' % spin)

    def get_pseudo_wave_function(self, band=0, kpt=0, spin=None):
        """Return pseudo-wave-function array.

        The method is limited to the gamma point, and is implemented
        as a wrapper to denchar (a tool shipped with siesta);
        denchar must be available in the command path.

        When retrieving a p_w_f from a non-spin-polarized calculation,
        spin must be None (default), and for spin-polarized
        calculations, spin must be set to either 0 (up) or 1 (down).

        As long as the necessary files are present and named
        correctly, old p_w_fs can be read as long as the
        calculator label is set. E.g.

        >>> c = Siesta(label='name_of_old_calculation')
        >>> pwf = c.get_pseudo_wave_function()

        The broadcast and pad options are not implemented.
        """

        # Not implemented: kpt=0, broadcast=True, pad=True
        # kpoint must be Gamma
        assert kpt == 0, \
            "siesta.get_pseudo_wave_function is unfortunately limited " \
            "to the gamma point only. kpt must be 0."

        # In denchar, band numbering starts from 1
        assert isinstance(band, int) and band >= 0
        band = band+1

        if spin is None:
            spin_name = ""
        elif spin == 0:
            spin_name = ".UP"
        elif spin == 1:
            spin_name = ".DOWN"

        label = self.label
        # If <label>.WF<band>.cube already exist and is newer than <label>.fdf,
        # just return it
        fn_wf = label+('.WF%i%s.cube'%(band, spin_name))
        fn_fdf = label+'.fdf'
        if isfile(fn_wf) and isfile(fn_fdf) and \
            (getmtime(fn_wf) > getmtime(fn_fdf)):
            x, _ = read_cube_data(fn_wf)
            return x

        if not isfile(fn_fdf):
            raise RuntimeError('Could not find the fdf-file. '
              'It is required as part of the input for denchar.')

        fdf_mtime = getmtime(fn_fdf)
        for suf in ['.WFS', '.PLD', '.DM', '.DIM']:
            if not isfile(label+suf):
                raise RuntimeError('Could not find file "%s%s" which is '
                              'required when extracting wave functions '
                              '(make sure the fdf options "WriteDenchar" is '
                              'True, and WaveFuncKpoints is [0.0 0.0 0.0]")' %
                                   (label, suf))
            if not getmtime(label+suf) > fdf_mtime:
                # This should be handled in a better way, e.g. by implementing
                # a "calculation_required() and calculate()"
                raise RuntimeError('The calculation is not up to date.')

        # Simply read the old fdf-file and pick some meta info from there.
        # However, strictly it's not always neccesary
        fdf = read_fdf(fn_fdf)
        if 'latticeconstant' in fdf:
            const = float(fdf['latticeconstant'][0])
            unit =  fdf['latticeconstant'][1]
        else:
            const = 1.0
            unit = 'Ang'

        if 'latticevectors' in fdf:
            cell = np.array(fdf['latticevectors'], dtype='d')
        else:
            raise RuntimeError('Failed to find the lattice vectors in the fdf-file.')

        if 'spinpolarized' in fdf and \
                fdf['spinpolarized'][0].lower() in \
                ['yes', 'true', False, 'T', '']:
            if spin is None:
                raise RuntimeError('The calculation was spin polarized, pick either '
                                   'spin=0 or 1.')
        else:
            if not spin is None:
                raise RuntimeError('The calculation was not spin polarized, '
                                   'spin argument must be None.')

        denc_fdf = open(fn_fdf).readlines()
        denc_fdf.append('Denchar.TypeOfRun 3D\n')
        denc_fdf.append('Denchar.PlotWaveFunctions T\n')
        for dim, dir in zip(cell.transpose(), ['X', 'Y', 'Z']):
            # Naive square box limits to denchar
            denc_fdf.append('Denchar.Min%s %f %s\n' % (dir, const*dim.min(), unit))
            denc_fdf.append('Denchar.Max%s %f %s\n' % (dir, const*dim.max(), unit))

        # denchar rewinds stdin and fails if stdin is a pipe
        denc_fdf_file = open(label+'.denchar.fdf', 'w')
        denc_fdf_file.write(''.join(denc_fdf))
        denc_fdf_file.close()

        try:
            from subprocess import Popen, PIPE
            p = Popen('denchar', shell=True, stdin=open(label+'.denchar.fdf'),
                      stdout=PIPE, stderr=PIPE, close_fds=True)
            exitcode = p.wait()
        except ImportError:
            raise RuntimeError('get_pseudo_wave_function implemented only with subprocess.')

        if exitcode == 0:
            if not isfile(fn_wf):
                raise RuntimeError('Could not find the requested file (%s)'%fn_wf)
            x, _ = read_cube_data(fn_wf)
            return x
        elif exitcode == 127:
            raise RuntimeError('No denchar executable found. Make sure it is in the path.')
        else:
            import sys
            print(''.join(p.stderr.readlines()), file=sys.stderr)
            raise RuntimeError('Execution of denchar failed!')


    def calculate(self, atoms):

        self.positions = atoms.get_positions().copy()
        self.cell = atoms.get_cell().copy()
        self.pbc = atoms.get_pbc().copy()

        if self.write_fdf_file:
            self.write_fdf(atoms)

        # if self.run_calc and not self.siesta_run:
        if self.mpirun:
            siesta = os.environ['SIESTA_MPI_SCRIPT']
            locals = {'np': self.np, 'label': self.label}
            exec(compile(open(siesta).read(), siesta, 'exec'), {}, locals)
        else:
            siesta = os.environ['SIESTA_SCRIPT']
            locals = {'label': self.label}
            exec(compile(open(siesta).read(), siesta, 'exec'), {}, locals)
        exitcode = locals['exitcode']
        self.siesta_run = True
        if exitcode != 0:
            raise RuntimeError(('Siesta exited with exit code: %d.  ' +
                                'Check %s.out for more information.') %
                               (exitcode, self.label))

        self.dipole = self.read_dipole()
        self.read()

        atoms_structout = read_struct_out('%s.STRUCT_OUT' % self.label)
        atoms.cell = atoms_structout.cell
        atoms.positions = atoms_structout.positions


        self.converged = True

    def set_fdf(self, key, value):
        """Set FDF parameter."""
        self.fdf[key] = value

#    def write_tddft_lr(self):
#      self.tddft_inp.write_tddft_inp()


    def write_fdf(self, atoms):
        """Write input parameters to fdf-file using the inp.siesta class."""
        self.MS.param['SystemLabel'] = self.label
        self.MS.param['SystemName'] = self.label

        magmoms = atoms.get_initial_magnetic_moments()
        if magmoms.any():
            self.MS.param['SpinPolarized'] = True
            self.MS.list_param.append('SpinPolarized')

        self.MS.write_siesta_file(atoms=atoms)

        if magmoms.any():
            fh = open(self.MS.param['SystemLabel'] + '.fdf', 'a')
            fh.write('\n')
            fh.write('%block InitSpin\n')
            for n, M in enumerate(magmoms):
                if M != 0:
                    fh.write('%d %.14f\n' % (n + 1, M))
                fh.write('%endblock InitSpin\n')
            fh.close()


    def calculation_required(self, atoms, properties):
        system_changes = self.check_state(atoms)
        if system_changes:
            return True

        for name in properties:
            if name not in self.results:
                return True

        return False


    def check_state(self, atoms, tol=1e-15):
        """Check for system changes since last calculation."""
        from .calculator import equal
        if self.atoms is None:
            system_changes = all_changes
        else:
            system_changes = []
            if not equal(self.atoms.positions, atoms.positions, tol):
                system_changes.append('positions')
            if not equal(self.atoms.numbers, atoms.numbers):
                system_changes.append('numbers')
            if not equal(self.atoms.cell, atoms.cell, tol):
                system_changes.append('cell')
            if not equal(self.atoms.pbc, atoms.pbc):
                system_changes.append('pbc')
            if not equal(self.atoms.get_initial_magnetic_moments(),
                         atoms.get_initial_magnetic_moments(), tol):
                system_changes.append('initial_magmoms')
            if not equal(self.atoms.get_initial_charges(),
                         atoms.get_initial_charges(), tol):
                system_changes.append('initial_charges')

        return system_changes


    def read(self):
        """Read results from SIESTA's text-output file."""
        text = open(self.label + '.txt', 'r').read().lower()
        assert 'error' not in text
        lines = iter(text.split('\n'))

        # Get the number of grid points used:
        for line in lines:
            if line.startswith('initmesh: mesh ='):
                self.grid = [int(word) for word in line.split()[3:8:2]]
                break

        # Stress (fixed so it's compatible with a MD run from siesta):
        for line in lines:
            if line.startswith('siesta: stress tensor '):
                stress = np.empty((3, 3))
                for i in range(3):
                    tmp = lines.next().split()
                    if len(tmp) == 4:
                        stress[i] = [float(word) for word in tmp[1:]]
                    else:
                        stress[i] = [float(word) for word in tmp]
                self.stress = np.array(
                    [stress[0, 0], stress[1, 1], stress[2, 2],
                     stress[1, 2], stress[0, 2], stress[0, 1]])
                break
        else:
            raise RuntimeError

        text = open(self.label + '.txt', 'r').read().lower()
        lines = iter(text.split('\n'))

        for line in lines:
            if line.startswith('siesta: etot    ='):
                self.etotal = float(line.split()[-1])
                self.efree = float(lines.next().split()[-1])
                break
        else:
            raise RuntimeError

        # Forces (changed so forces smaller than -999eV/A can be fetched):
        lines = open(self.label + '.FA', 'r').readlines()
        assert int(lines[0]) == len(self.numbers)
        assert len(lines) == len(self.numbers) + 1
        lines = lines[1:]
        self.forces = np.zeros((len(lines), 3))
        for i in range(len(lines)):
            self.forces[i, 0] = float(lines[i][6:18].strip())
            self.forces[i, 1] = float(lines[i][18:30].strip())
            self.forces[i, 2] = float(lines[i][30:42].strip())

    def read_eig(self):
        if self.e_fermi is not None:
            return

        assert os.access(self.label + '.EIG', os.F_OK)
        assert os.access(self.label + '.KP', os.F_OK)

        # Read k point weights
        text = open(self.label + '.KP', 'r').read()
        lines = text.split('\n')
        n_kpts = int(lines[0].strip())
        self.weights = np.zeros((n_kpts,))
        for i in range(n_kpts):
            l = lines[i + 1].split()
            self.weights[i] = float(l[4])

        # Read eigenvalues and fermi-level
        text = open(self.label+'.EIG', 'r').read()
        lines = text.split('\n')
        self.e_fermi = float(lines[0].split()[0])
        tmp = lines[1].split()
        self.n_bands = int(tmp[0])
        n_spin_bands = int(tmp[1])
        self.spin_pol = n_spin_bands == 2
        lines = lines[2:-1]
        lines_per_kpt = (self.n_bands * n_spin_bands / 10 +
                         int((self.n_bands * n_spin_bands) % 10 != 0))
        self.eig = dict()
        for i in range(len(self.weights)):
            tmp = lines[i * lines_per_kpt:(i + 1) * lines_per_kpt]
            v = [float(v) for v in tmp[0].split()[1:]]
            for l in tmp[1:]:
                v.extend([float(t) for t in l.split()])
            if self.spin_pol:
                self.eig[(i, 0)] = np.array(v[0:self.n_bands])
                self.eig[(i, 1)] = np.array(v[self.n_bands:])
            else:
                self.eig[(i, 0)] = np.array(v)

    def get_k_point_weights(self):
        self.read_eig()
        return self.weights

    def get_fermi_level(self):
        self.read_eig()
        return self.e_fermi

    def get_eigenvalues(self, kpt=0, spin=0):
        self.read_eig()
        return self.eig[(kpt, spin)]

    def get_number_of_spins(self):
        self.read_eig()
        if self.spin_pol:
            return 2
        else:
            return 1

    def get_temperature(self, unit='K'):
        T = self.MS.param['ElectronicTemperature']

        if T['unit'] != unit:
        # convert
            raise ValueError('conversion not yet implemented, please take'\
                            ' the same temperature unit (kelvin in preference)')
            conv_fact = 1 #need to be implemented
            T['Valeur'] = conv_fact*T['Valeur']
            return T['Valeur']
        else:
            return T['Valeur']

    def read_hs(self, filename, is_gamma_only=False, magnus=False):
        """Read the Hamiltonian and overlap matrix from a Siesta
           calculation in sparse format.

        Parameters
        ==========
        filename: str
            The filename should be on the form jobname.HS
        is_gamma_only: {False, True), optional
            Is it a gamma point calculation?
        magnus: bool
            The fileformat was changed by Magnus in Siesta at some
            point around version 2.xxx.
            Use mangus=False, to use the old file format.

        Note
        ====
        Data read in is put in self._dat.

        Examples
        ========
            >>> calc = Siesta()
            >>> calc.read_hs('jobname.HS')
            >>> print calc._dat.fermi_level
            >>> print 'Number of orbitals: %i' % calc._dat.nuotot
        """
        assert not magnus, 'Not implemented; changes by Magnus to file io'
        assert not is_gamma_only, 'Not implemented. Only works for k-points.'
        class Dummy:
            pass
        self._dat = dat = Dummy()
        # Try to read supercell and atom data from a jobname.XV file
        filename_xv = filename[:-2] + 'XV'
        # assert isfile(filename_xv), 'Missing jobname.XV file'
        if isfile(filename_xv):
            print('Reading supercell and atom data from ' + filename_xv)
            fd = open(filename_xv, 'r')
            dat.cell = np.zeros((3, 3)) # Supercell
            for a_vec in dat.cell:
                a_vec[:] = np.array(fd.readline().split()[:3], float)
            dat.rcell = 2 * np.pi * np.linalg.inv(dat.cell.T)
            dat.natoms = int(fd.readline().split()[0])
            dat.symbols = []
            dat.pos_ac = np.zeros((dat.natoms, 3))
            for a in range(dat.natoms):
                line = fd.readline().split()
                dat.symbols.append(chemical_symbols[int(line[1])])
                dat.pos_ac[a,:] = [float(line[i]) for i in range(2, 2 + 3)]
        # Read in the jobname.HS file
        fileobj = file(filename, 'rb')
        fileobj.seek(0)
        dat.fermi_level = float(open(filename[:-3] + '.EIG', 'r').readline())
        dat.is_gammay_only = is_gamma_only
        dat.nuotot, dat.ns, dat.mnh = getrecord(fileobj, 'l')
        nuotot, ns, mnh = dat.nuotot, dat.ns, dat.mnh
        print('Number of orbitals found: %i' % nuotot)
        dat.numh = numh = np.array([getrecord(fileobj, 'l')
                                    for i in range(nuotot)], 'l')
        dat.maxval = max(numh)
        dat.listhptr = listhptr = np.zeros(nuotot, 'l')
        listhptr[0] = 0
        for oi in range(1, nuotot):
            listhptr[oi] = listhptr[oi - 1] + numh[oi - 1]
        dat.listh = listh = np.zeros(mnh, 'l')

        print('Reading sparse info')
        for oi in range(nuotot):
            for mi in range(numh[oi]):
                listh[listhptr[oi] + mi] = getrecord(fileobj, 'l')

        dat.nuotot_sc = max(listh)
        dat.h_sparse = h_sparse = np.zeros((mnh, ns), float)
        dat.s_sparse = s_sparse = np.zeros(mnh, float)
        print('Reading H')
        for si in range(ns):
            for oi in range(nuotot):
                for mi in range(numh[oi]):
                    h_sparse[listhptr[oi] + mi, si] = getrecord(fileobj, 'd')
        print('Reading S')
        for oi in range(nuotot):
            for mi in range(numh[oi]):
                s_sparse[listhptr[oi] + mi] = getrecord(fileobj, 'd')

        dat.qtot, dat.temperature = getrecord(fileobj, 'd')
        if not is_gamma_only:
            print('Reading X')
            dat.xij_sparse = xij_sparse = np.zeros([3, mnh], float)
            for oi in range(nuotot):
                for mi in range(numh[oi]):
                    xij_sparse[:, listhptr[oi] + mi] = getrecord(fileobj, 'd')
        fileobj.close()

    def get_hs(self, kpt=(0, 0, 0), spin=0, remove_pbc=None, kpt_scaled=True):
        """Hamiltonian and overlap matrices for an arbitrary k-point.

        The default values corresponds to the Gamma point for
        spin 0 and periodic boundary conditions.

        Parameters
        ==========
        kpt : {(0, 0, 0), (3,) array_like}, optional
            k-point in scaled or absolute coordinates.
            For the latter the units should be Bohr^-1.
        spin : {0, 1}, optional
            Spin index
        remove_pbc : {None, ({'x', 'y', 'z'}, basis)}, optional
            Use remove_pbc to truncate h and s along a cartesian
            axis.
        basis: {str, dict}
            The basis specification as either a string or a dictionary.
        kpt_scaled : {True, bool}, optional
            Use kpt_scaled=False if `kpt` is in absolute units (Bohr^-1).

        Note
        ====
        read_hs should be called before get_hs gets called.

        Examples
        ========
        >>> calc = Siesta()
        >>> calc.read_hs('jobname.HS')
        >>> h, s = calc.get_hs((0.0, 0.375, 0.375))
        >>> h -= s * calc._dat.fermi_level # fermi level is now at 0.0
        >>> basis = 'szp'
        >>> h, s = calc.get_hs((0.0, 0.375, 0.375), remove_pbc=('x', basis))
        >>> basis = {'Au:'sz}', 'C':'dzp', None:'szp'}
        >>> h, s = calc.get_hs((0.0, 0.375, 0.375), remove_pbc=('x', basis))

        """
        if not hasattr(self, '_dat'):# XXX Crude check if data is avail.
            print('Please read in data first by calling the method read_hs.')
            return None, None
        dot = np.dot
        dat = self._dat
        kpt_c = np.array(kpt, float)
        if kpt_scaled:
            kpt_c = dot(kpt_c, dat.rcell)

        h_MM = np.zeros((dat.nuotot, dat.nuotot), complex)
        s_MM = np.zeros((dat.nuotot, dat.nuotot), complex)
        h_sparse, s_sparse = dat.h_sparse, dat.s_sparse
        x_sparse = dat.xij_sparse
        numh, listhptr, listh = dat.numh, dat.listhptr, dat.listh
        indxuo = np.mod(np.arange(dat.nuotot_sc), dat.nuotot)

        for iuo in range(dat.nuotot):
            for j in range(numh[iuo]):
                ind =  listhptr[iuo] + j
                jo = listh[ind] - 1
                juo = indxuo[jo]
                kx = dot(kpt_c, x_sparse[:, ind])
                phasef = exp(1.0j * kx)
                h_MM[iuo, juo] += phasef * h_sparse[ind, spin]
                s_MM[iuo, juo] += phasef * s_sparse[ind]

        if remove_pbc is not None:
            direction, basis = remove_pbc
            centers_ic = get_bf_centers(dat.symbols, dat.pos_ac, basis)
            d = 'xyz'.index(direction)
            cutoff = dat.cell[d, d] * 0.5
            truncate_along_axis(h_MM, s_MM, direction, centers_ic, cutoff)

        h_MM *= complex(Rydberg)
        return h_MM, s_MM


class WithUnit:
    allowed = tuple()
    def __init__(value, unit):
        self.__value = value
        self.setUnit(unit)

    def setUnit(self, unit):
        assert unit in self.allowed
        self.__unit = unit

class Mass(WithUnit): allowed = ['Kg', 'g', 'amu']
class Length(WithUnit): allowed = ['m', 'cm', 'nm', 'Ang', 'Bohr']
class Time(WithUnit): allowed=['s', 'fs', 'ps', 'ns', 'mins', 'hours', 'days']
class Force(WithUnit): allowed=['N', 'eV/Ang', 'Ry/Bohr']
class Charge(WithUnit): allowed = ['C', 'e']
class Dipole(WithUnit): allowed = ['C*m', 'D', 'debye', 'e*Bohr', 'e*Ang']
class InertialMoment(WithUnit): allowed=['Kg*m**2', 'Ry*fs**2']
class EField(WithUnit): allowed=['V/m', 'V/nm', 'V/Ang', 'V/Bohr', 'Ry/Bohr/e', 'Har/Bohr/e']
class Angle(WithUnit): allowed=['deg', 'rad']
class Energy(WithUnit):
    allowed=['J', 'erg', 'eV', 'meV', 'Ry', 'mRy',
             'Hartree', 'mHartree', 'K', 'Kcal/mol', 'KJ/mol',
             'Hz', 'THz', 'cm-1', 'cm**-1', 'cm^-1']
class Pressure(WithUnit):
    allowed=['Pa', 'MPa', 'GPa', 'atm', 'bar', 'Kbar', 'Mbar', 'Ry/Bohr**3', 'eV/Ang**3']
class Torque(WithUnit): allowed=['eV/deg', 'eV/rad', 'Ry/deg', 'Ry/rad', 'meV/deg',
                                 'meV/rad', 'mRy/deg', 'mRy/rad']



class SiestaParameters(LockedParameters):
    def write_fdf(self, f):
        for key, value in self.iteritems():
            key = self.prefix() + '.' + key
            f.write(format_fdf(key, value))

class SolutionMethod(SiestaParameters):
    def identitier(self):
        raise NotImplementedError

    def write_fdf(self, f):
        f.write(format_fdf('SolutionMethod', self.identifier()))
        SiestaParameters.write_fdf(self, f)

class Diag(SolutionMethod):
    def prefix(self):
        return 'Diag'

    def identifier(self):
        return 'diagon'

    def __init__(self,
        DivideAndConquer=False,
        AllInOne=False,
        NoExpert=False,
        PreRotate=False,
        Use2D=False,
        Memory=1.0,
        ParallelOverK=False,
        ):
        kwargs = locals()
        kwargs.pop('self')
        SolutionMethod.__init__(self, **kwargs)

class OrderN(SolutionMethod):
    def prefix(self):
        return 'ON'

    def identifier(self):
        return 'ON'

    def __init__(self,
        functional='Kim',
        MaxNumIter=1000,
        etol=1e-8,
        eta=(0.0, 'eV'),
        eta_alpha=(0.0, 'eV'),
        eta_beta=(0.0, 'eV'),
        RcLWF=(9.5, 'Bohr'),
        ChemicalPotential=False,
        ChemicalPotentialUse=False,
        ChemicalPotentialRc=(9.5, 'Bohr'),
        ChemicalPotentialTemperature=(0.05,'Ry'),
        ChemicalPotentialOrder=100,
        LowerMemory=False,
        UseSaveLWF=False,
        OccupationFunction='FD',
        OccupationMPOrder=1,
        ):
        kwargs = locals()
        kwargs.pop('self')
        SolutionMethod.__init__(self, **kwargs)

class XC(LockedParameters):
    def __init__(self,
        functional='LDA',
        authors='PZ',
        ):
        LockedParameters.__init__(
                self,
                functional=functional,
                authors=authors,
                )
    def write_fdf(self, f):
        f.write(format_fdf('XC_functional', self['functional']))
        f.write(format_fdf('XC_authors', self['authors']))


class Specie(LockedParameters):
    def __init__(self,
                 symbol,
                 basis_set=DZP,
                 pseudopotential=None,
                 tag=None,
                 ghost=False,
                 ):
        kwargs = locals()
        kwargs.pop('self')
        LockedParameters.__init__(self, **kwargs)

class FDFArguments(LockedParameters):
    def __init__(self,
                 DM_Tolerance = 1e-4,
                 ):
        kwargs = locals()
        kwargs.pop('self')
        LockedParameters.__init__(self, **kwargs)

    def write_fdf(self, f):
        for key, value in self.iteritems():
            f.write(format_fdf(key, value))

class Siesta(FileIOCalculator):
    """  """
    implemented_properties = tuple(['energy', 'forces'])

    default_parameters = LockedParameters(
                 mesh_cutoff=(200, 'Ry'),
                 energy_shift=(100, 'meV'),
                 kpts=(1,1,1),
                 label='siesta',
                 atoms=None,
                 xc='LDA.PZ',
                 species=tuple(),
                 basis_set=DZP,
                 spin=UNPOLARIZED,
                 solution_method=Diag(),
                 pseudo_qualifier=None,
                 pseudo_path=None,
                 n_nodes=1,
                 restart=None,
                 ignore_bad_restart_file=False,
                 fdf_arguments=FDFArguments(),
                 )

    def __init__(self, **kwargs):
        parameters = self.get_default_parameters()
        functional, authors = parameters['xc'].split('.')
        parameters['xc'] = XC(functional=functional, authors=authors)
        parameters.update(kwargs)

        siesta = os.environ.get('SIESTA')
        label = parameters['label']
        if parameters['n_nodes']>1:
            command = 'mpirun -np %d %s < ./%s.fdf > ./%s.out'%(n_nodes, siesta, label, label)
        else:
            command = '%s < ./%s.fdf > ./%s.out'%(siesta, label, label)

        FileIOCalculator.__init__(
                self,
                command=command,
                **parameters
                )

    def __getitem__(self, key):
        return self.parameters[key]

    def species(self, atoms):
        symbols = np.array(atoms.get_chemical_symbols())
        tags = atoms.get_tags()
        species = list(self['species'])
        default_species = [specie for specie in species if specie['tag'] is None]
        default_symbols = [specie['symbol'] for specie in default_species]
        for symbol in symbols:
            if not symbol in default_symbols:
                specie = Specie(
                        symbol=symbol,
                        basis_set=self['basis_set'],
                        tag=None,
                        )
                default_species.append(specie)
                default_symbols.append(symbol)

        assert len(default_species) == len(np.unique(symbols))
        non_default_species = [specie for specie in species if not specie['tag'] is None]
        species_numbers = np.zeros(len(atoms), int)
        i = 1
        for specie in default_species:
            mask = symbols==specie['symbol']
            species_numbers[mask] = i
            i += 1

        for specie in non_default_species:
            mask1 = tags==specie['tag']
            mask2 = symbols==specie['symbol']
            mask = np.logical_and(mask1, mask2)
            if sum(mask) > 0:
                species_numbers[mask] = i
                i += 1

        all_species = default_species + non_default_species
        return all_species, species_numbers

    def calculate(self, atoms=None, properties=['energy'],
                  system_changes=all_changes):
        FileIOCalculator.calculate(self, atoms=atoms, properties=properties,
                  system_changes=system_changes)

    def write_input(self, atoms, properties=None, system_changes=None):
        FileIOCalculator.write_input(
                self,
                atoms=atoms,
                properties=properties,
                system_changes=system_changes,
                )
        if system_changes is None:
            return
        filename = self.label + '.fdf'
        fdf_dict = {}

        restart_fdf = ''
        if tuple(system_changes) == ('positions',):
            restart_fdf = 'DM.UseSaveDM  T\n'

        with open(filename, 'w') as f:
            f.write(restart_fdf)
            f.write(format_fdf('SystemName', self.label))
            f.write(format_fdf('SystemName', self.label))
            f.write(format_fdf('SystemLabel', self.label))
            self['solution_method'].write_fdf(f)
            self.__write_basis(f, atoms)
            self.__write_kpts(f)
            self.__write_structure(f, atoms)
            self['fdf_arguments'].write_fdf(f)

    def __write_structure(self, f, atoms):
        "Write STRUCT.fdf file"
        xyz=atoms.get_positions()
        species, species_numbers = self.species(atoms)
        basis_sizes = []
        for specie, number in zip(species, species_numbers):
            if specie['basis_set']['size'] != self['basis_set']['size']:
                basis_sizes.append((number, specie['basis_set']['size']))
        if len(basis_sizes) > 0:
            f.write(format_fdf('PAO.BasisSizes', basis_sizes))
        mask = np.zeros(len(atoms), bool)
        spec = np.zeros(0, bool)

        unit_cell = atoms.get_cell()

        f.write('\n')
        f.write(format_fdf('NumberOfAtoms', len(xyz)))
        f.write(format_fdf('NumberOfSpecies', len(species)))
        f.write(format_fdf('LatticeConstant', (1.0, 'Ang')))
        f.write('%block LatticeVectors\n')
        for i in range(3):
            for j in range(3):
                f.write(string.rjust('%.15f'%unit_cell[i,j],16)+ ' ')
            f.write('\n')
        f.write('%endblock LatticeVectors\n')

        self.__write_atomic_coordinates(f, atoms)

    def __write_atomic_coordinates(self, f, atoms):
        species, species_numbers = self.species(atoms)
        f.write('AtomicCoordinatesFormat  Ang\n')
        f.write('%block AtomicCoordinatesAndAtomicSpecies\n')
        for atom, number in zip(atoms, species_numbers):
            xyz = atom.position
            line=string.rjust('%.9f'%xyz[0],16)+' '
            line+=string.rjust('%.9f'%xyz[1],16)+' '
            line+=string.rjust('%.9f'%xyz[2],16)+' '
            line+=str(number)+'\n'
            f.write(line)
        f.write('%endblock AtomicCoordinatesAndAtomicSpecies\n')

        origin = tuple(-atoms.get_celldisp().flatten())
        f.write('%block AtomicCoordinatesOrigin\n')
        f.write('%.4f  %.4f  %.4f\n'%origin)
        f.write('%endblock AtomicCoordinatesOrigin\n')

    def __write_kpts(self, f):
        kpts = np.array(self['kpts'])
        f.write('\n')
        f.write('#KPoint grid\n')
        f.write('%block kgrid_Monkhorst_Pack\n')

        for i in range(3):
            s=''
            if i<len(kpts):
                number=kpts[i]
                displace=0.0
            else:
                number=1
                displace=0
            for j in range(3):
                if j==i:
                    write_this=number
                else:
                    write_this=0
                s+='%d  '%write_this
            s+='%1.1f\n'%displace
            f.write(s)
        f.write('%endblock kgrid_Monkhorst_Pack\n')
        f.write('\n')

    def __write_basis(self, f, atoms,):
        f.write(format_fdf('PAO.BasisSize', self['basis_set']['size']))
        f.write(format_fdf('PAO_EnergyShift', self['energy_shift']))

        species, species_numbers = self.species(atoms)
        f.write(format_fdf('MeshCutoff', self['mesh_cutoff']))
        self['spin'].write_fdf(f)
        self['xc'].write_fdf(f)

        if not self['pseudo_path'] is None:
            pseudo_path = self['pseudo_path']
        elif 'SIESTA_PP_PATH' in os.environ:
            pseudo_path = os.environ['SIESTA_PP_PATH']
        else:
            raise Exception

        f.write('%block ChemicalSpeciesLabel\n')
        for species_number, specie in enumerate(species):
            species_number += 1
            symbol = specie['symbol']
            atomic_number = atomic_numbers[symbol]

            if specie['pseudopotential'] is None:
                label = '.'.join([symbol, self.pseudo_qualifier()])
                pseudopotential = label + '.psf'
            else:
                pseudopotential = specie['pseudopotential']
                label = os.path.basename(pseudopotential)
                label = '.'.join(label.split('.')[:-1])

            if not os.path.isabs(pseudopotential):
                pseudopotential = join(pseudo_path, pseudopotential)

            if not os.path.exists(pseudopotential):
                raise RuntimeError('No pseudopotential for %s!' % symbol)

            name = os.path.basename(pseudopotential)
            name = name.split('.')
            name.insert(-1, str(species_number))
            if specie['ghost']:
                name.insert(-1, 'ghost')
                atomic_number= -atomic_number
            name = '.'.join(name)

            if join(os.getcwd(), name) != pseudopotential:
                if islink(name) or isfile(name):
                    os.remove(name)
                os.symlink(pseudopotential, name)

            label = '.'.join(np.array(name.split('.'))[:-1])

            f.write('%d %d %s\n'%(species_number, atomic_number, label))
        f.write('%endblock ChemicalSpeciesLabel\n')
        f.write('\n')

    def pseudo_qualifier(self):
        if self['pseudo_qualifier'] is None:
            return self['xc']['functional'].lower()
        else:
            return self['pseudo_qualifier']

    def read_results(self):
        self.read_energy()
        self.read_forces_stress()

    def read_energy(self):
        """Read results from SIESTA's text-output file."""
        with open(self.label + '.out', 'r') as f:
            text = f.read().lower()

        assert 'error' not in text
        lines = iter(text.split('\n'))

        # Get the number of grid points used:
        for line in lines:
            if line.startswith('initmesh: mesh ='):
                self.results['n_grid_point'] = [int(word) for word in line.split()[3:8:2]]
                break

        for line in lines:
            if line.startswith('siesta: etot    ='):
                self.results['energy'] = float(line.split()[-1])
                break
        else:
            raise RuntimeError

    def read_forces_stress(self):
        with open('FORCE_STRESS', 'r') as f:
            lines = f.readlines()

        stress_lines = lines[1:4]
        stress = np.empty((3, 3))
        for i in range(3):
            line = [s for s in stress_lines[i].strip().split(' ') if len(s)>0]
            stress[i] = map(float, line)

        self.results['stress'] = np.array(
                [stress[0, 0], stress[1, 1], stress[2, 2],
                 stress[1, 2], stress[0, 2], stress[0, 1]])

        start = 5
        self.results['forces'] = np.zeros((len(lines)-start, 3), float)
        for i in range(start, len(lines)):
            line = [s for s in lines[i].strip().split(' ') if len(s)>0]
            self.results['forces'][i-start] = map(float, line[2:5])


def getrecord(fileobj, dtype):
    """Used to read in binary files.
    """
    typetosize = {'l':4, 'f':4, 'd':8}# XXX np.int, np.float32, np.float64
    assert dtype in typetosize # XXX
    size = typetosize[dtype]
    record = array.array('l')
    trunk = array.array(dtype)
    record.fromfile(fileobj, 1)
    nofelements = int(record[-1]) / size
    trunk.fromfile(fileobj, nofelements)
    record.fromfile(fileobj, 1)
    data = np.array(trunk, dtype=dtype)
    if len(data) == 1:
        data = data[0]
    return data

def truncate_along_axis(h, s, direction, centers_ic, cutoff):
    """Truncate h and s such along a cartesian axis.

    Parameters:

    h: (N, N) ndarray
        Hamiltonian matrix.
    s: (N, N) ndarray
        Overlap matrix.
    direction: {'x', 'y', 'z'}
        Truncate allong a cartesian axis.
    centers_ic: (N, 3) ndarray
        Centers of the basis functions.
    cutoff: float
        The (direction-axis projected) cutoff distance.
    """
    dtype = h.dtype
    ni = len(centers_ic)
    d = 'xyz'.index(direction)
    pos_i = centers_ic[:, d]
    for i in range(ni):
        dpos_i = abs(pos_i - pos_i[i])
        mask_i = (dpos_i < cutoff).astype(dtype)
        h[i,:] *= mask_i
        h[:, i] *= mask_i
        s[i,:] *= mask_i
        s[:, i] *= mask_i

def get_nao(symbol, basis):
    """Number of basis functions.

    Parameters
    ==========
    symbol: str
        The chemical symbol.
    basis: str
        Basis function type.
    """
    ls = valence_config[symbol]
    nao = 0
    zeta = {'s':1, 'd':2, 't':3, 'q':4}
    nzeta = zeta[basis[0]]
    is_pol = 'p' in basis
    for l in ls:
        nao += (2 * l + 1) * nzeta
    if is_pol:
        l_pol = None
        l = -1
        while l_pol is None:
            l += 1
            if not l in ls:
                l_pol = l
        nao += 2 * l_pol + 1
    return nao

def get_bf_centers(symbols, positions, basis):
    """Centers of basis functions.

    Parameters
    ==========
    symbols: str, (N, ) array_like
        chemical symbol for each atom.
    positions: float, (N, 3) array_like
        Positions of the atoms.
    basis: {str,  dict}
        Basis set specification as either a string or a dictionary

    Examples
    ========
    >>> symbols = ['O', 'H']
    >>> positions = [(0, 0, 0), (0, 0, 1)]
    >>> basis = 'sz'
    >>> print get_bf_centers(symbols, positions, basis)
    [[0 0 0]
     [0 0 0]
     [0 0 0]
     [0 0 0]
     [0 0 1]]
    >>> basis = {'H':'dz', None:'sz'}
    >>> print get_bf_centers(symbols, positions, basis)
    [[0 0 0]
     [0 0 0]
     [0 0 0]
     [0 0 0]
     [0 0 1]
     [0 0 1]]

    """
    centers_ic = []
    dict_basis = False
    if isinstance(basis, dict):
        dict_basis = True
    for symbol, pos in zip(symbols, positions):
        if dict_basis:
            if symbol not in basis:
                bas = basis[None]
            else:
                bas = basis[symbol]
        else:
            bas = basis
        for i in range(get_nao(symbol, bas)):
            centers_ic.append(pos)
    return np.asarray(centers_ic)

def format_fdf(key, value):
    key = format_key(key)

    block = False
    if isinstance(value, list):
        block = True
    value = format_value(value)

    if block:
        return '%' + key + '\n' + value + '\n' + '%endblock ' + key + '\n'
    else:
        return '%s  %s\n'%(key, value)

def format_value(value):
    if isinstance(value, tuple):
        value = '%s %s'%value
    if isinstance(value, list):
        sub_values = map(format_value, value)
        value = '\n'.join(sub_values)

    return value

def format_key(key):
    return key.replace('_', '.')

valence_config = {
    'H': (0,),
    'C': (0, 1),
    'N': (0, 1),
    'O': (0, 1),
    'S': (0, 1),
    'Li': (0,),
    'Na': (0,),
    'Ni': (0, 2),
    'Cu': (0, 2),
    'Pd': (0, 2),
    'Ag': (0, 2),
    'Pt': (0, 2),
    'Au': (0, 2)}

def get_species(atoms, pseudopotentials, ghost_filter=None):
    if ghost_filter is None:
        ghost_filter = np.zeros(len(atoms), np.bool)

    species = []
    all_numbers = atoms.get_atomic_numbers()
    non_doped_numbers=all_numbers

    identities = set(zip(non_doped_numbers, ghost_filter))
    identities = sorted(list(identities))

    for number, ghost_state in identities:
        el = chemical_symbols[number]
        key = (el, number, ghost_state)
        species.append(key)

    return species

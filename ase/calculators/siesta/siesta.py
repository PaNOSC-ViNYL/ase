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
from collections import OrderedDict

from ase.data import chemical_symbols, atomic_numbers
from ase.units import Rydberg, fs, Bohr
from ase.io.siesta import read_rho, read_fdf, read_struct_out
from ase.io.cube import read_cube_data
from ase.calculators.siesta.basis_set import BasisSet, DZP
from ase.calculators.calculator import FileIOCalculator, all_changes
from ase.calculators.calculator import LockedParameters

class WithUnit:
    allowed = tuple()
    def __init__(self, value, unit):
        self.__value = value
        self.setUnit(unit)

    def setUnit(self, unit):
        assert unit in self.allowed
        self.__unit = unit

    def script(self):
        return '%s %s'%(self.__value, self.__unit)

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

class Specie(LockedParameters):
    def __init__(self,
                 symbol,
                 basis_set='DZP',
                 pseudopotential=None,
                 tag=None,
                 ghost=False,
                 ):
        kwargs = locals()
        kwargs.pop('self')
        LockedParameters.__init__(self, **kwargs)

class FDFArguments(LockedParameters):
    def __init__(self,
                 DM_Tolerance=1e-4,
                 ):
        kwargs = locals()
        kwargs.pop('self')
        LockedParameters.__init__(self, **kwargs)

    def write_fdf(self, f):
        for key, value in self.iteritems():
            f.write(format_fdf(key, value))

class Siesta(FileIOCalculator):
    """  """
    implemented_properties = tuple([
        'energy',
        'forces',
        'stress',
        'dipole',
        'eigenvalues',
        'density',
        'fermi_energy',
        ])

    default_parameters = LockedParameters(
                 mesh_cutoff=Energy(200, 'Ry'),
                 energy_shift=Energy(100, 'meV'),
                 kpts=(1,1,1),
                 label='siesta',
                 atoms=None,
                 xc='LDA.PZ',
                 species=tuple(),
                 basis_set='DZP',
                 spin='FULL',
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

    def get_potential_energy(self, atoms, force_consistent=False):
        if not 'energy' in self.results.keys():
            self.calculate(atoms=atoms, properties=['energy'])
        energy = self.results['energy']
        energy_free = self.results['energy free']
        if force_consistent:
            return energy
        else:
            return  (energy + energy_free) / 2

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

    def set(self, **kwargs):
        allowed_names = ['SZ', 'SZP', 'DZ', 'DZP']
        basis_set = kwargs.get('basis_set')
        if not basis_set is None and (not basis_set in allowed_names):
            raise Exception("Basis must be either %s, got %s"%(allowed_names, basis_set) )
        allowed_spins = ['UNPOLARIZED', 'COLLINEAR', 'FULL']
        spin = kwargs.get('spin')
        if not spin is None and (not spin in allowed_spins):
            raise Exception("Spin must be %s, got %s"%(allowed_spins, spin) )

        allowed_functionals = ['LDA', 'GGA', 'VDW']
        allowed_authors = {
                'LDA': ['PZ', 'CA', 'PW92'],
                'GGA': ['PW91', 'PBE', 'revPBE', 'RPBE',
                        'WC', 'AM05', 'PBEsol', 'PBEJsJrLO',
                        'PBEGcGxLO', 'PBEGcGxHEG', 'BLYP',
                        ],
                'VDW': ['DRSLL', 'LMKLL', 'KBM', 'C09', 'BH', 'VV'],
                }
        xc = kwargs.get('xc')
        if not xc is None:
            split = xc.split('.')
            if len(split) == 2:
                functional, authors = split
            elif len(split) == 1:
                functional = split[0]
                authors = None
            else:
                raise Exception("The xc argument must be of the format 'functional.authors', got %s"%xc)

            if not functional in allowed_functionals:
                raise Exception("Functional must be %s, got %s"%(allowed_functionals, functional) )
            if not authors is None:
                if not authors in allowed_authors[functional]:
                    raise Exception("Functional authors must be %s, got %s"%(allowed_authors[functional], authors) )

        d_parameters = self.get_default_parameters()
        for key, value in d_parameters.iteritems():
            if isinstance(value, WithUnit) and not isinstance(kwargs[key], WithUnit):
                new_value, unit = kwargs[key]
                kwargs[key] = value.__class__(value=new_value, unit=unit)

        FileIOCalculator.set(self, **kwargs)

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
        if len(system_changes) > 0:
            self.removeAnalysis()

        with open(filename, 'w') as f:
            if not 'numbers' in system_changes and \
                not 'initial_magmoms' in system_changes and \
                not 'initial_charges' in system_changes \
                :

                restart_fdf = 'DM.UseSaveDM  T\n'
                f.write(restart_fdf)
            if 'density' in properties:
                f.write(format_fdf('SaveRho', True))
            if 'hamiltonian' in properties or 'overlap' in properties:
                f.write(format_fdf('SaveHS', True))

            f.write(format_fdf('SystemName', self.label))
            f.write(format_fdf('SystemName', self.label))
            f.write(format_fdf('SystemLabel', self.label))
            self['solution_method'].write_fdf(f)
            self.__write_basis(f, atoms)
            self.__write_kpts(f)
            self.__write_structure(f, atoms)
            self['fdf_arguments'].write_fdf(f)

    def removeAnalysis(self):
        filename = self.label + '.RHO'
        if os.path.exists(filename):
            os.remove(filename)

    def __write_structure(self, f, atoms):
        "Write STRUCT.fdf file"
        xyz=atoms.get_positions()
        species, species_numbers = self.species(atoms)
        basis_sizes = []
        for specie, number in zip(species, species_numbers):
            if specie['basis_set'] != self['basis_set']:
                basis_sizes.append((number, specie['basis_set']))
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

        magmoms = atoms.get_initial_magnetic_moments()
        f.write('%block DM.InitSpin\n')
        for n, M in enumerate(magmoms):
            if M != 0:
                f.write('%d %.14f\n' % (n + 1, M))
        f.write('%endblock DM.InitSpin\n')

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
        f.write(format_fdf('PAO.BasisSize', self['basis_set']))
        f.write(format_fdf('PAO_EnergyShift', self['energy_shift']))

        species, species_numbers = self.species(atoms)
        f.write(format_fdf('MeshCutoff', self['mesh_cutoff']))
        if self['spin'] == 'UNPOLARIZED':
            f.write(format_fdf('SpinPolarized', False))
        elif self['spin'] == 'COLLINEAR':
            f.write(format_fdf('SpinPolarized', True))
        elif self['spin'] == 'FULL':
            f.write(format_fdf('SpinPolarized', True))
            f.write(format_fdf('NonCollinearSpin', True))

        split = self['xc'].split('.')
        if len(split) == 2:
            functional, authors = split
        else:
            functional = split[0]
            authors = None
        f.write(format_fdf('XC_functional', functional))
        if not authors is None:
            f.write(format_fdf('XC_authors', authors))

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
            return self['xc'].split('.')[0].lower()
        else:
            return self['pseudo_qualifier']

    def read_results(self):
        self.read_energy()
        self.read_forces_stress()
        self.read_eigenvalues()
        self.read_dipole()
        self.read_pseudo_density()

    def read_pseudo_density(self):
        filename = self.label + '.RHO'
        if isfile(filename):
            self.results['density'] = read_rho(filename)

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
                line = lines.next()
                self.results['energy free'] = float(line.split()[-1])
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

    def read_eigenvalues(self):
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
        with open(self.label+'.EIG', 'r') as f:
            text = f.read()
        lines = text.split('\n')
        e_fermi = float(lines[0].split()[0])
        tmp = lines[1].split()
        self.n_bands = int(tmp[0])
        n_spin_bands = int(tmp[1])
        self.spin_pol = n_spin_bands == 2
        lines = lines[2:-1]
        lines_per_kpt = (self.n_bands * n_spin_bands / 10 +
                         int((self.n_bands * n_spin_bands) % 10 != 0))
        eig = OrderedDict()
        for i in range(len(self.weights)):
            tmp = lines[i * lines_per_kpt:(i + 1) * lines_per_kpt]
            v = [float(v) for v in tmp[0].split()[1:]]
            for l in tmp[1:]:
                v.extend([float(t) for t in l.split()])
            if self.spin_pol:
                eig[(i, 0)] = np.array(v[0:self.n_bands])
                eig[(i, 1)] = np.array(v[self.n_bands:])
            else:
                eig[(i, 0)] = np.array(v)

        self.results['fermi_energy'] = e_fermi
        self.results['eigenvalues'] = eig

    def read_dipole(self):
        dipolemoment = np.zeros([1, 3])
        with open(self.label + '.out', 'r') as f:
            for line in f:
                if line.rfind('Electric dipole (Debye)') > -1:
                    dipolemoment = np.array([float(f) for f in line.split()[5:8]])

        # debye to e*Ang
        self.results['dipole'] = dipolemoment * 0.2081943482534

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
    if isinstance(value, WithUnit):
        value = value.script()
    elif isinstance(value, tuple):
        value = '%s %s'%value
    elif isinstance(value, list):
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

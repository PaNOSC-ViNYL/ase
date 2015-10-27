from __future__ import print_function
"""This module defines an ASE interface to SIESTA.

http://www.uam.es/departamentos/ciencias/fismateriac/siesta
"""
import os
import shutil
import sys
from ase.atoms import Atoms
import gzip
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
from ase.calculators.calculator import FileIOCalculator, all_changes, ReadError
from ase.calculators.calculator import LockedParameters

ALLOWED_BASIS_NAMES = ['SZ', 'SZP', 'DZ', 'DZP']
ALLOWED_SPINS = ['UNPOLARIZED', 'COLLINEAR', 'FULL']
ALLOWED_AUTHORS = {
                'LDA': ['PZ', 'CA', 'PW92'],
                'GGA': ['PW91', 'PBE', 'revPBE', 'RPBE',
                        'WC', 'AM05', 'PBEsol', 'PBEJsJrLO',
                        'PBEGcGxLO', 'PBEGcGxHEG', 'BLYP',
                        ],
                'VDW': ['DRSLL', 'LMKLL', 'KBM', 'C09', 'BH', 'VV'],
                }

Bohr2Ang = 0.529177
def handleAtomicNumbers(atomic_numbers):
    atomic_numbers = np.array(atomic_numbers)
    take = np.array(atomic_numbers) > 0
    atomic_numbers[take] = atomic_numbers[take]%200

    return atomic_numbers

def ReadXVFile(filename,InUnits='Bohr',OutUnits='Ang',ReadVelocity=False):
    "Returns tuple (vectors,speciesnumber,atomnumber,xyz,[v,]) from an XV-file"
    if (InUnits=='Bohr') and (OutUnits=='Ang'): convFactor = Bohr2Ang
    elif (InUnits=='Ang') and (OutUnits=='Bohr'): convFactor = Ang2Bohr
    elif (((InUnits=='Ang') and (OutUnits=='Ang')) \
       or ((InUnits=='Bohr') and (OutUnits=='Bohr'))): convFactor = 1
    else:pass
    try:
        file = open(filename,'r')
    except:
        file = gzip.open(filename+'.gz','r')

    # Read cell vectors (lines 1-3)
    vectors = []
    for i in range(3):
        data = string.split(file.readline())
        vectors.append([string.atof(data[j])*convFactor for j in range(3)])
    # Read number of atoms (line 4)
    numberOfAtoms = string.atoi(string.split(file.readline())[0])
    # Read remaining lines
    speciesnumber, atomnumber, xyz, V = [], [], [], []
    for line in file.readlines():
        if len(line)>5: # Ignore blank lines
            data = string.split(line)
            speciesnumber.append(string.atoi(data[0]))
            atomnumber.append(string.atoi(data[1]))
            xyz.append([string.atof(data[2+j])*convFactor for j in range(3)])
            V.append([string.atof(data[5+j])*convFactor for j in range(3)])
    file.close()
    if ReadVelocity:
        return vectors,speciesnumber,atomnumber,xyz, V
    else:
        return vectors,speciesnumber,atomnumber,xyz

def convertXVtoASE(filename):
    result=ReadXVFile(filename=filename,InUnits='Bohr',OutUnits='Ang',ReadVelocity=False)
    unit_cell=np.array(result[0])
    numbers = np.array(result[2])
    numbers = handleAtomicNumbers(numbers)
    positions=np.array(result[3])
    atoms = Atoms(numbers=numbers, positions=positions, cell=unit_cell)

    return atoms

class WithUnit:
    """ Abstract class to represent quantities with units."""
    # This is the units allowed by this type of quantity.
    allowed = tuple()

    def __init__(self, value, unit):
        """
        Set the value and unit of this quantity.
        """
        self.__value = value
        self.setUnit(unit)

    def setUnit(self, unit):
        if not unit in self.allowed:
            raise Exception("Received unit '%s', but only %s are allowed."%(unit, self.allowed) )
        self.__unit = unit

    def script(self):
        """
        Write the fdf script for this quantity.
        """
        return '%s %s'%(self.__value, self.__unit)

# All types of quantities with units in Siesta, with all allowed input units.
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
             'Hz', 'THz', 'cm-1', 'cm**-1', 'cm^-1', 'Bohr']

class Pressure(WithUnit):
    allowed=['Pa', 'MPa', 'GPa', 'atm', 'bar', 'Kbar', 'Mbar', 'Ry/Bohr**3', 'eV/Ang**3']
class Torque(WithUnit): allowed=['eV/deg', 'eV/rad', 'Ry/deg', 'Ry/rad', 'meV/deg',
                                 'meV/rad', 'mRy/deg', 'mRy/rad']


class SiestaParameters(LockedParameters):
    """
    Parameter class which can write out its own fdf-script.
    """
    def write_fdf(self, f):
        for key, value in self.iteritems():
            key = self.prefix() + '.' + key
            f.write(format_fdf(key, value))

class SolutionMethod(SiestaParameters):
    """
    Collection of parameters related to a specific solution method.
    """
    def identitier(self):
        """
        The string which begins all fdf-keywords in the group.
        """
        raise NotImplementedError

    def write_fdf(self, f):
        """
        Write the SolutionMethod keyword to the fdf-script as well as all
        parameters in this group.
        """
        f.write(format_fdf('SolutionMethod', self.identifier()))
        SiestaParameters.write_fdf(self, f)

class Diag(SolutionMethod):
    """
    Parameters related to the diagonalization solution method.
    """
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
    """
    Parameters related to the OrderN solution method.
    """
    def prefix(self):
        return 'ON'

    def identifier(self):
        return 'ON'

    def __init__(self,
        functional='Kim',
        MaxNumIter=1000,
        etol=1e-8,
        eta=Energy(0.0, 'eV'),
        eta_alpha=Energy(0.0, 'eV'),
        eta_beta=Energy(0.0, 'eV'),
        RcLWF=Energy(9.5, 'Bohr'),
        ChemicalPotential=False,
        ChemicalPotentialUse=False,
        ChemicalPotentialRc=Energy(9.5, 'Bohr'),
        ChemicalPotentialTemperature=Energy(0.05,'Ry'),
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
    """
    Parameters for specifying the behaviour for a single species in the
    calculation. If the tag argument is set to and integer then atoms with
    the specified element and tag will be a seperate species.

    Pseudopotential and basis set can be specified. Additionally the species
    can be set be a ghost species, meaning that they will not be considered
    atoms, but the corresponding basis set will be used.
    """
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

class Siesta(FileIOCalculator):
    """
    Calculator interface to the SIESTA code.
    """
    implemented_properties = tuple([
        'energy',
        'forces',
        'stress',
        'dipole',
        'eigenvalues',
        'density',
        'fermi_energy',
        ])

    # Dictionary of valid input vaiables.
    # Any additional variables will be written directly as to fdf-files.
    # Replace '.' with '_' in arguments.
    default_parameters = LockedParameters(
                 label='siesta',
                 mesh_cutoff=Energy(200, 'Ry'),
                 energy_shift=Energy(100, 'meV'),
                 kpts=(1,1,1),
                 atoms=None,
                 xc='LDA.PZ',
                 species=tuple(),
                 basis_set='DZP',
                 spin='COLLINEAR',
                 solution_method=Diag(),
                 pseudo_qualifier=None,
                 pseudo_path=None,
                 n_nodes=1,
                 restart=None,
                 ignore_bad_restart_file=False,
                 )

    def __init__(self, **kwargs):
        """
        ASE interface to the SIESTA code.

        Parameters:
            -label        : The base head of all created files.
            -mesh_cutoff  : tuple of (value, energy_unit)
                            The mesh cutoff energy for determining number of
                            grid points.
            -xc
            -energy_shift : tuple of (value, energy_unit)
                            The confining energy of the basis sets.
            -kpts         : Tuple of 3 integers, the k-points in different
                            directions.
            -atoms        : The Atoms object.
            -species      : None|list of Specie objects. The species objects
                            can be used to to specify the basis set,
                            pseudopotential and whether the species is ghost.
                            The tag on the atoms object and the element is used
                            together to identify the species.
            -basis_set    : "SZ"|"SZP"|"DZ"|"DZP", strings which specify the
                            type of functions basis set.
            -spin         : "UNPOLARIZED"|"COLLINEAR"|"FULL". The level of spin
                            description to be used.
            -solution_method : Diag object| OrderN object. This is the
                            internal method used by the siesta calculator.
            -pseudo_path  : None|path. This path is where
                            pseudopotentials are taken from.
                            If None is given, then then the path given
                            in $SIESTA_PP_PATH will be used.
            -pseudo_qualifier: None|string. This string will be added to the
                            pseudopotential path that will be retrieved.
                            For hydrogen with qualifier "abc" the
                            pseudopotential "H.abc.psf" will be retrieved.
            -n_nodes      : The number of nodes to use.
            -restart      : str.  Prefix for restart file.
                            May contain a directory.
                            Default is  None, don't restart.
            -ignore_bad_restart_file: bool.
                            Ignore broken or missing restart file.
                            By default, it is an error if the restart
                            file is missing or broken.
                            """
        parameters = self.get_default_parameters()
        parameters.update(kwargs)

        # Setup the siesta command based on number of nodes.
        siesta = os.environ.get('SIESTA')
        label = parameters['label']
        self.label = label
        if parameters['n_nodes']>1:
            command = 'mpirun -np %d %s < ./%s.fdf > ./%s.out'%(n_nodes, siesta, label, label)
        else:
            command = '%s < ./%s.fdf > ./%s.out'%(siesta, label, label)

        # Call the base class.
        FileIOCalculator.__init__(
                self,
                command=command,
                **parameters
                )

    def __getitem__(self, key):
        return self.parameters[key]

    def species(self, atoms):
        """
        Find all relevant species depending on the atoms object and
        species input.

        Parameters ::
            - atoms : An Atoms object.
        """
        # For each element use default specie from the species input, or set
        # up a default species  from the general default parameters.
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

        # Set default species as the first species.
        species_numbers = np.zeros(len(atoms), int)
        i = 1
        for specie in default_species:
            mask = symbols==specie['symbol']
            species_numbers[mask] = i
            i += 1

        # Set up the non-default species.
        non_default_species = [specie for specie in species if not specie['tag'] is None]
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
        """
        Set all parameters.
        """
        # Check the basis set input.
        basis_set = kwargs.get('basis_set')
        if not basis_set is None and (not basis_set in ALLOWED_BASIS_NAMES):
            raise Exception("Basis must be either %s, got %s"%(ALLOWED_BASIS_NAMES, basis_set) )

        # Check the spin input.
        spin = kwargs.get('spin')
        if not spin is None and (not spin in ALLOWED_SPINS):
            raise Exception("Spin must be %s, got %s"%(ALLOWED_SPINS, spin) )

        # Check the functional input.
        xc = kwargs.get('xc')
        if not xc is None:
            split = xc.split('.')
            if len(split) == 2:
                functional, authors = split
            elif len(split) == 1:
                functional = split[0]
                authors = ALLOWED_AUTHORS[functional]
            else:
                raise Exception("The xc argument must be of the format 'functional.authors', got %s"%xc)

            if not functional in ALLOWED_AUTHORS.keys():
                raise Exception("Functional must be %s, got %s"%(ALLOWED_AUTHORS, functional) )
            if not authors is None:
                if not authors in ALLOWED_AUTHORS[functional]:
                    raise Exception("Functional authors must be %s, got %s"%(ALLOWED_AUTHORS[functional], authors) )

        # Check that quantities with units have the correct units.
        d_parameters = self.get_default_parameters()
        for key, value in d_parameters.iteritems():
            if isinstance(value, WithUnit) and not isinstance(kwargs[key], WithUnit):
                new_value, unit = kwargs[key]
                kwargs[key] = value.__class__(value=new_value, unit=unit)

        FileIOCalculator.set(self, **kwargs)

    def write_input(self, atoms, properties=None, system_changes=None):
        """
        Write input (fdf)-file.
        """
        # Call base calculator.
        FileIOCalculator.write_input(
                self,
                atoms=atoms,
                properties=properties,
                system_changes=system_changes,
                )
        # Early exist.
        if system_changes is None and properties is None:
            return

        filename = self.label + '.fdf'
        fdf_dict = {}

        # On any changes, remove all analysis files.
        if len(system_changes) > 0:
            self.removeAnalysis()

        # Start writing the file.
        with open(filename, 'w') as f:
            # Use the saved density matrix if only 'cell' and 'positions'
            # haved changes.
            if not 'numbers' in system_changes and \
                not 'initial_magmoms' in system_changes and \
                not 'initial_charges' in system_changes \
                :
                restart_fdf = 'DM.UseSaveDM  T\n'
                f.write(restart_fdf)

            # Save density.
            if 'density' in properties:
                f.write(format_fdf('SaveRho', True))

            # Write system name and label.
            f.write(format_fdf('SystemName', self.label))
            f.write(format_fdf('SystemLabel', self.label))

            # Write solution method and all parameters.
            self['solution_method'].write_fdf(f)

            # Write the rest.
            self.__write_species(f, atoms)
            self.__write_kpts(f)
            self.__write_structure(f, atoms)
            self.__write_fdf_arguments(f)

    def read(self, restart):
        if not os.path.exists(restart):
            raise ReadError("The restart file '%s' does not exist"%restart)
        self.atoms = convertXVtoASE(restart)
        self.read_results()

    def __write_fdf_arguments(self, f):
        """
        Write all arguments not given as default directly as fdf-format.
        """
        d_parameters = self.get_default_parameters()
        for key, value in self.parameters.iteritems():
            if not key in d_parameters.keys():
                f.write(format_fdf(key, value))

    def removeAnalysis(self):
        """ Remove all analysis files"""
        filename = self.label + '.RHO'
        if os.path.exists(filename):
            os.remove(filename)

    def __write_structure(self, f, atoms):
        """
        Translate the Atoms object to fdf-format.

        Parameters:
            - f:     An open file object.
            - atoms: An atoms object.
        """
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

        xyz=atoms.get_positions()
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

        # Write magnetic moments.
        magmoms = atoms.get_initial_magnetic_moments()
        f.write('%block DM.InitSpin\n')
        for n, M in enumerate(magmoms):
            if M != 0:
                f.write('%d %.14f\n' % (n + 1, M))
        f.write('%endblock DM.InitSpin\n')

    def __write_atomic_coordinates(self, f, atoms):
        """
        Write atomic coordinates.

        Parameters:
            - f:     An open file object.
            - atoms: An atoms object.
        """
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
        """
        Write kpts.

        Parameters:
            - f : Open filename.
        """
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

    def __write_species(self, f, atoms):
        """
        Write input related the different species.

        Parameters:
            - f:     An open file object.
            - atoms: An atoms object.
        """
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
            raise Exception("Please set the environment variable 'SIESTA_PP_PATH'")

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
        """
        Get the extra string used in the middle of the pseudopotential.
        The retrieved pseudopotential for a specific element will be
        'H.xxx.psf' for the element 'H' with qualifier 'xxx'. If qualifier
        is set to None then the qualifier is set to functional name.
        """
        if self['pseudo_qualifier'] is None:
            return self['xc'].split('.')[0].lower()
        else:
            return self['pseudo_qualifier']

    def read_results(self):
        """
        Read the results.
        """
        self.read_energy()
        self.read_forces_stress()
        self.read_eigenvalues()
        self.read_dipole()
        self.read_pseudo_density()

    def read_pseudo_density(self):
        """
        Read the density if it is there.
        """
        filename = self.label + '.RHO'
        if isfile(filename):
            self.results['density'] = read_rho(filename)

    def read_energy(self):
        """
        Read energy from SIESTA's text-output file.
        """
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
                self.results['free_energy'] = float(line.split()[-1])
                break
        else:
            raise RuntimeError

    def read_forces_stress(self):
        """
        Read the forces and stress from the FORCE_STRESS file.
        """
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
        """
        Read eigenvalues from the '.EIG' file.
        This is done pr. kpoint.
        """
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
        """
        Read dipole moment.
        """
        dipolemoment = np.zeros([1, 3])
        with open(self.label + '.out', 'r') as f:
            for line in f:
                if line.rfind('Electric dipole (Debye)') > -1:
                    dipolemoment = np.array([float(f) for f in line.split()[5:8]])

        # debye to e*Ang
        self.results['dipole'] = dipolemoment * 0.2081943482534

def format_fdf(key, value):
    """
    Write an fdf key-word value pair.

    Parameters:
        - key   : The fdf-key
        - value : The fdf value.
    """
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
    """
    Format python values to fdf-format.

    Parameters:
        - value : The value to format.
    """
    if isinstance(value, WithUnit):
        value = value.script()
    elif isinstance(value, tuple):
        value = '%s %s'%value
    elif isinstance(value, list):
        sub_values = map(format_value, value)
        value = '\n'.join(sub_values)
    else:
        value = str(value)

    return value

def format_key(key):
    """ Fix the fdf-key replacing '_' with '.' """
    return key.replace('_', '.')

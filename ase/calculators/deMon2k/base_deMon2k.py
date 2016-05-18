from __future__ import print_function
"""This module defines an ASE interface to deMon2k.

http://www.demon-software.com 
"""
import os
from os.path import join, isfile, islink
import string
import numpy as np
from ase.units import Ry, eV, Bohr, Hartree
from ase.data import atomic_numbers
#from ase.calculators.siesta.import_functions import read_rho, xv_to_atoms
from ase.calculators.calculator import FileIOCalculator, ReadError
from ase.calculators.calculator import Parameters, all_changes
from ase.calculators.calculator import equal
#from ase.calculators.siesta.parameters import PAOBasisBlock, Specie
#from ase.calculators.siesta.parameters import format_fdf

meV = 0.001 * eV


class deMon2k_Parameters(Parameters):
    """Parameters class for the calculator.
    Documented in BaseSiesta.__init__

    """
    def __init__(
            self,
            label='rundir/deMon',
            atoms=None,
            restart=None,
            ignore_bad_restart_file=False,
            title="demon2k input file",
            scftype='RKS',
            multiplicity=None,
            charge=0,
            xc='VWN',
            guess='TB',
            print_out="",
            basis ={},
            ecps ={},
            mcps ={},
            auxis={},
            augment={},
            input_arguments=None):
        kwargs = locals()
        kwargs.pop('self')
        Parameters.__init__(self, **kwargs)


class Base_deMon2k(FileIOCalculator):
    """Calculator interface to the deMon2k code.
    """
    #allowed_xc = {}
    allowed_keywords = {}

    implemented_properties = (
        'energy',
        'forces',
        'dipole',
        'eigenvalues',
        'density',
        'fermi_energy')

    # Dictionary of valid input vaiables.
    #default_parameters = deMon2k_Parameters()

    def __init__(self, **kwargs):
        """ASE interface to the deMon2k code.

        Parameters:
            -xc           : The exchange-correlation potential. Can be set to
                            any allowed value for either the Siesta
                            XC.funtional or XC.authors keyword. Default "LDA"
            -basis_set    : "SZ"|"SZP"|"DZ"|"DZP", strings which specify the
                            type of functions basis set.
            -basis_path  : None|path. This path is where
                            basis sets, ecps and auxiliary basis sets are taken from.
                            If None is given, then then the path given
                            in $DEMON2K_BASIS_PATH will be used.
            -atoms        : The Atoms object.
            -restart      : str.  Prefix for restart file.
                            May contain a directory.
                            Default is  None, don't restart.
            -ignore_bad_restart_file: bool.
                            Ignore broken or missing restart file.
                            By default, it is an error if the restart
                            file is missing or broken.
            -input_arguments: Explicitly given input arguments. Dictonary using
                            deMon2k keywords as given in the manual. List values
                            are written as fdf blocks with each element on a
                            separate line, while tuples will write each element
                            in a single line.  ASE units are assumed in the
                            input.
        """
        # Put in the default arguments.
        #parameters = self.default_parameters.__class__(**kwargs)
        parameters = deMon2k_Parameters(**kwargs)
        
        # Setup the run command based on number of nodes.
        command = os.environ.get('DEMON2K_COMMAND')
        if command is None:
            mess = "The 'DEMON2K_COMMAND' environment is not defined."
            raise ValueError(mess)

        label = parameters['label']
        runfile = label + '.inp'
        outfile = label + '.out'

        try:
            command = command % (runfile, outfile)
        except TypeError:
            raise ValueError(
                "The 'DEMON2K_COMMAND' environment must " +
                "be a format string" +
                " with two string arguments.\n" +
                "Example : 'deMon < ./%s > ./%s'.\n" +
                "Got '%s'" % command)

        # Call the base class.
        FileIOCalculator.__init__(
            self,
            command=command,
            **parameters)

        #self.parameters =parameters
        
    def __getitem__(self, key):
        """Convenience method to retrieve a parameter as
        calculator[key] rather than calculator.parameters[key]

            Parameters:
                -key       : str, the name of the parameters to get.
        """
        return self.parameters[key]


    def set(self, **kwargs):
        """Set all parameters.

            Parameters:
                -kwargs  : Dictionary containing the keywords defined in
                           SiestaParameters.
        """
        # Put in the default arguments.
        kwargs = self.default_parameters.__class__(**kwargs)

        if 'parameters' in kwargs:
            filename = kwargs.pop('parameters')
            parameters = Parameters.read(filename)
            parameters.update(kwargs)
            kwargs = parameters

        changed_parameters = {}

        for key, value in kwargs.items():
            oldvalue = self.parameters.get(key)
            if key not in self.parameters or not equal(value, oldvalue):
                #if isinstance(oldvalue, dict):
                #    # Special treatment for dictionary parameters:
                #    for name in value:
                #        if name not in oldvalue:
                #            raise KeyError(
                #                'Unknown subparameter "%s" in '
                #                'dictionary parameter "%s"' % (name, key))
                #    oldvalue.update(value)
                #    value = oldvalue
                changed_parameters[key] = value
                self.parameters[key] = value

        return changed_parameters


        #FileIOCalculator.set(self, **kwargs)
        
        ## Check energy inputs.
        #for arg in ['mesh_cutoff', 'energy_shift']:
        #    value = kwargs.get(arg)
        #    if not (isinstance(value, (float, int)) and value > 0):
        #        mess = "'%s' must be a positive number(in eV), \
        #            got '%s'" % (arg, value)
        #        raise ValueError(mess)

#        # Check the basis set input.
#        basis_set = kwargs.get('basis_set')
#        allowed = self.allowed_basis_names
#        if not (isinstance(basis_set, PAOBasisBlock) or basis_set in allowed):
#            mess = "Basis must be either %s, got %s" % (allowed, basis_set)
#            raise Exception(mess)
#
#        # Check the spin input.
#        spin = kwargs.get('spin')
#        if spin is not None and (spin not in self.allowed_spins):
#            mess = "Spin must be %s, got %s" % (self.allowed_spins, spin)
#            raise Exception(mess)
#
#        # Check the functional input.
#        xc = kwargs.get('xc')
#        if isinstance(xc, (tuple, list)) and len(xc) == 2:
#            functional, authors = xc
#            if not functional in self.allowed_xc:
#                mess = "Unrecognized functional keyword: '%s'" % functional
#                raise ValueError(mess)
#            if not authors in self.allowed_xc[functional]:
#                mess = "Unrecognized authors keyword for %s: '%s'"
#                raise ValueError(mess % (functional, authors))
#
#        elif xc in self.allowed_xc:
#            functional = xc
#            authors = self.allowed_xc[xc][0]
#        else:
#            found = False
#            for key, value in self.allowed_xc.iteritems():
#                if xc in value:
#                    found = True
#                    functional = key
#                    authors = xc
#                    break
#
#            if not found:
#                raise ValueError("Unrecognized 'xc' keyword: '%s'" % xc)
#        kwargs['xc'] = (functional, authors)
#
#        # Check input_arguments.
#        input_arguments = kwargs['input_arguments']
#        if input_arguments is not None:
#            # Type checking.
#            if not isinstance(input_arguments, dict):
#                raise TypeError("input_arguments must be a dictionary.")
#
#            # Check if keywords are allowed.
#            input_keys = set(input_arguments.keys())
#            allowed_keys = set(self.allowed_input_keywords)
#            if not input_keys.issubset(allowed_keys):
#                offending_keys = input_keys.difference(allowed_keys)
#                raise ValueError("The 'input_arguments' dictionary " +
#                                 "argument does not allow " +
#                                 "the keywords: %s" % str(offending_keys))
#


    def calculate(self,
                  atoms=None,
                  properties=['energy'],
                  system_changes=all_changes):
        """Capture the RuntimeError from FileIOCalculator.calculate
        and add a little debug information from the Siesta output.

        See base FileIocalculator for documentation.
        """

        try:
            FileIOCalculator.calculate(
                self,
                atoms=atoms,
                properties=properties,
                system_changes=system_changes)

# this must be modified
        except RuntimeError as e:
            try:
                with open(self.label + '.out', 'r') as f:
                    lines = f.readlines()
                debug_lines = 10
                print('##### %d last lines of the Siesta output' % debug_lines)
                for line in lines[-20:]:
                    print(line.strip())
                print('##### end of siesta output')
                raise e
            except:
                raise e

    def write_input(self, atoms, properties=None, system_changes=None):
        """Write input (in)-file.
        See calculator.py for further details.

        Parameters:
            - atoms        : The Atoms object to write.
            - properties   : The properties which should be calculated.
            - system_changes : List of properties changed since last run.
        """
        # Call base calculator.
        FileIOCalculator.write_input(
            self,
            atoms=atoms,
            properties=properties,
            system_changes=system_changes)

        if system_changes is None and properties is None:
            return

        filename = self.label + '.inp'

        # On any changes, remove all analysis files.

#        # this must be modified
#        if system_changes is not None:
#            self.remove_analysis()

        # Start writing the file.
        with open(filename, 'w') as f:

            # write keyword argument keywords 
            value = self.parameters['title']
            self._write_argument('TITLE', value, f)

            f.write('#\n')
            
            value = self.parameters['scftype']
            self._write_argument('SCFTYPE', value, f)

            value = self.parameters['multiplicity']
            if value is not None:
                self._write_argument('MULTIPLICITY', value, f)
            
            value = self.parameters['charge']
            self._write_argument('CHARGE', value, f)
            
            value = self.parameters['xc']
            self._write_argument('VXCTYPE', value, f)

            value = self.parameters['guess']
            self._write_argument('GUESS', value, f)

            value = self.parameters['print_out']
            if not len(value) == 0:
                self._write_argument('PRINT', value, f)
                f.write('#\n')            

            # write general input arguments
            self._write_input_arguments(f)
            
            f.write('#\n')

            # write basis set etc 
            basis = self.parameters['basis']
            if not 'all' in basis:
                basis['all'] = 'DZVP'
            #if not len(basis) == 0:
            self._write_basis(f, atoms, basis, string='BASIS')

            ecps = self.parameters['ecps']
            if not len(ecps) == 0:
                self._write_basis(f, atoms, ecps, string='ECPS')

            mcps = self.parameters['mcps']
            if not len(mcps) == 0:
                self._write_basis(f, atoms, mcps, string='MCPS')

            auxis = self.parameters['auxis']
            if not len(auxis) == 0:
                self._write_basis(f, atoms, auxis, string='AUXIS')

            augment = self.parameters['augment']
            if not len(augment) == 0:
                self._write_basis(f, atoms, augment, string='AUGMENT')

            # write geometry
            self._write_atomic_coordinates(f, atoms)
                
                
    # this must be modified
    def read(self, filename):
        """Read parameters from file."""
        if not os.path.exists(filename):
            raise ReadError("The restart file '%s' does not exist" % filename)
        self.atoms = xv_to_atoms(filename)
        self.read_results()

    def _write_input_arguments(self, f):
        """Write directly given input-arguments.
        """
        input_arguments = self.parameters['input_arguments']

        # Early return
        if input_arguments is None:
            return

        for key, value in input_arguments.iteritems():
            #elif key in self.allowed_input_keywords:

            self._write_argument(key, value, f)

            
    def _write_argument(self, key, value, f):
        """ write an argument to file. 
        - key :  a string coresponding to the input keyword
        - value : the arguemnts, can be a string, a number or a list
        - f :  and open file
        """
        
        # for only one argument, write on same line
        if not isinstance(value, (tuple, list)):
            line = key
            line += "    " + str(value)
            f.write(line)
            f.write('\n')
        # for a list, write first argument on the first line, then the rest on new lines
        else:
            line = key
            if not isinstance(value[0], (tuple, list)):
                for i in range(len(value)):
                    line += "  " + str(value[i])
                f.write(line)
                f.write('\n')
            else:
                for i in range(len(value)):
                    for j in range(len(value[i])):
                        line += "  " + str(value[i][j])
                    f.write(line)
                    f.write('\n')
                    line = ""
                        
                    
                    
            #else:
            #    raise ValueError("%s not in allowed keywords." % key)

## this must be modified
#    def remove_analysis(self):
#        """ Remove all analysis files"""
#        filename = self.label + '.RHO'
#        if os.path.exists(filename):
#            os.remove(filename)

    # this has been modified
    def _write_atomic_coordinates(self, f, atoms):
        """Write atomic coordinates.
        
        Parameters:
        - f:     An open file object.
        - atoms: An atoms object.
        """

        #f.write('\n')
        f.write('#\n')
        f.write('# Atomic coordinates\n')
        f.write('#\n')
        f.write('GEOMETRY CARTESIAN ANGSTROM\n')

        for i in range(len(atoms)):
            xyz = atoms.get_positions()[i]
            chem_symbol = atoms.get_chemical_symbols()[i]
            chem_symbol += str(i+1)
            
            # if tag is set to 1 then we have a ghost atom, set nuclear charge to 0
            if(atoms.get_tags()[i] == 1):
                nuc_charge = str(0) 
            else:
                nuc_charge = str(atoms.get_atomic_numbers()[i])
            
            mass = atoms.get_masses()[i]
                
            line = '{0:6s}'.format(chem_symbol).rjust(10) + ' ' 
            line += '{0:.5f}'.format(xyz[0]).rjust(10) + ' ' 
            line += '{0:.5f}'.format(xyz[1]).rjust(10) + ' ' 
            line += '{0:.5f}'.format(xyz[2]).rjust(10) + ' ' 
            line += '{0:5s}'.format(nuc_charge).rjust(10) + ' ' 
            line += '{0:.5f}'.format(mass).rjust(10) + ' ' 
            
            f.write(line)
            f.write('\n')
        #f.write('#\n')

    # routine to write basis set inormation, including ecps and auxis
    def _write_basis(self, f, atoms, basis={}, string='BASIS'): 
        """Write basis set, ECPs, AUXIS, or AUGMENT basis
        
        Parameters:
        - f:     An open file object.
        - atoms: An atoms object.
        - basis: A dictionary specifying the basis set
        - string: 'BASIS', 'ECP','AUXIS' or 'AUGMENT'
        """

        # basis for all atoms
        line = '{0}'.format(string).ljust(10)

        if 'all' in basis:
            default_basis = basis['all']
            line += '({0})'.format(default_basis).rjust(16)
        
        f.write(line)
        f.write('\n')        

        # basis for all atomic species
        chemical_symbols = atoms.get_chemical_symbols()
        chemical_symbols_set = set(chemical_symbols)

        for i in range(chemical_symbols_set.__len__()):
            symbol = chemical_symbols_set.pop()

            if symbol in basis:
                line = '{0}'.format(symbol).ljust(10)
                line += '({0})'.format(basis[symbol]).rjust(16)
                f.write(line)
                f.write('\n')        

        # basis for individual atoms
        for i in range(len(atoms)):
            
            if i in basis:
                symbol = str(chemical_symbols[i])
                symbol += str(i+1)

                line = '{0}'.format(symbol).ljust(10)
                line += '({0})'.format(basis[i]).rjust(16)
                f.write(line)
                f.write('\n')        


        
                
# main method for analyisis, must be modified
    def read_results(self):
        """Read the results.
        """
        self.read_energy()
        #self.read_forces_stress()
        #self.read_eigenvalues()
        #self.read_dipole()
        #self.read_pseudo_density()

# might be useless
    def read_pseudo_density(self):
        """Read the density if it is there.
        """
        filename = self.label + '.RHO'
        if isfile(filename):
            self.results['density'] = read_rho(filename)

    def read_energy(self):
        """Read energy from deMon2k's text-output file.
        """
        print(self.label + '.out')
        with open(self.label + '.out', 'r') as f:
            text = f.read().upper()

        #assert 'error' not in text
        lines = iter(text.split('\n'))

        ## Get the number of grid points used:
        #for line in lines:
        #    if line.startswith('initmesh: mesh ='):
        #        n_points = [int(word) for word in line.split()[3:8:2]]
        #        self.results['n_grid_point'] = n_points
        #        break

        for line in lines:
            if line.startswith(' TOTAL ENERGY                ='):
                self.results['energy'] = float(line.split()[-1]) * Hartree
                #print(line)
                print(self.results['energy'])
                #line = lines.next()
                #self.results['free_energy'] = float(line.split()[-1])
                break
        else:
            raise RuntimeError

# this must be modified
    def read_forces_stress(self):
        """Read the forces and stress from the FORCE_STRESS file.
        """
        with open('FORCE_STRESS', 'r') as f:
            lines = f.readlines()

        stress_lines = lines[1:4]
        stress = np.empty((3, 3))
        for i in range(3):
            line = stress_lines[i].strip().split(' ')
            line = [s for s in line if len(s) > 0]
            stress[i] = map(float, line)

        self.results['stress'] = np.array(
            [stress[0, 0], stress[1, 1], stress[2, 2],
             stress[1, 2], stress[0, 2], stress[0, 1]])

        self.results['stress'] *= Ry / Bohr**3

        start = 5
        self.results['forces'] = np.zeros((len(lines) - start, 3), float)
        for i in range(start, len(lines)):
            line = [s for s in lines[i].strip().split(' ') if len(s) > 0]
            self.results['forces'][i - start] = map(float, line[2:5])

        self.results['forces'] *= Ry / Bohr

# this must be modified
    def read_eigenvalues(self):
        """Read eigenvalues from the '.EIG' file.
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
        with open(self.label + '.EIG', 'r') as f:
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
        eig = dict()
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

# this must be modified
    def read_dipole(self):
        """Read dipole moment.
        """
        dipole = np.zeros([1, 3])
        with open(self.label + '.out', 'r') as f:
            for line in f:
                if line.rfind('Electric dipole (Debye)') > -1:
                    dipole = np.array([float(f) for f in line.split()[5:8]])

        # debye to e*Ang
        self.results['dipole'] = dipole * 0.2081943482534

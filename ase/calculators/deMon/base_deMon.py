from __future__ import print_function
"""This module defines an ASE interface to deMon.

http://www.demon-software.com 
"""
import os
from os.path import join, isfile, islink
import string
import numpy as np
from ase.units import Ry, eV, Bohr, Hartree
from ase.data import atomic_numbers
import ase.data
from ase.calculators.calculator import FileIOCalculator, ReadError
from ase.calculators.calculator import Parameters, all_changes
from ase.calculators.calculator import equal
import ase.io
import subprocess
import pickle
import shutil

class Parameters_deMon(Parameters):
    """Parameters class for the calculator.
    Documented in BaseSiesta.__init__

    The options here are the most important ones that the user needs to be aware of. 
    Further options accepted by deMon can be set in the dictionary input_arguments.

    """
    def __init__(
            self,
            label='rundir', # relative path to the rundir
            atoms=None,
            restart=None,  # relative path to restart to dir for parameters and atoms object, not deMon restart
            basis_path=None, 
            ignore_bad_restart_file=False, # for the restart option
            deMon_restart_path='.',  #  relative path to the deMon restart dir it looks for deMon.mem and copies it to deMon.rst
            title="deMon input file", # title in input file. not so useful
            scftype='RKS', # 'RKS', 'UKS', 'ROKS' 
            forces=False,  # if True we will run a BOMD trajectory with one step in order to extract the forces 
            xc='VWN',
            guess='TB', # This controls the restart if 'RESTART' 
            print_out='MOE', # print out the orbital eigenvalues by default so that they can be parsed
            basis ={},
            ecps ={},
            mcps ={},
            auxis={},
            augment={},
            input_arguments=None):
        kwargs = locals()
        kwargs.pop('self')
        Parameters.__init__(self, **kwargs)


class Base_deMon(FileIOCalculator):
    """Calculator interface to the deMon code.
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
    #default_parameters = demon_Parameters()

    def __init__(self, **kwargs):
        """ASE interface to the deMon code.

        Parameters:
            -xc           : The exchange-correlation potential. Can be set to
                            any allowed value for either the Siesta
                            XC.funtional or XC.authors keyword. Default "LDA"
            -basis_set    : "SZ"|"SZP"|"DZ"|"DZP", strings which specify the
                            type of functions basis set.
            -basis_path  : None|path. This path is where
                            basis sets, ecps and auxiliary basis sets are taken from.
                            If None is given, then then the path given
                            in $DEMON_BASIS_PATH will be used.
            -atoms        : The Atoms object.
            -restart      : str.  Prefix for restart file.
                            May contain a directory.
                            Default is  None, don't restart.
            -ignore_bad_restart_file: bool.
                            Ignore broken or missing restart file.
                            By default, it is an error if the restart
                            file is missing or broken.
            -input_arguments: Explicitly given input arguments. Dictonary using
                            demon keywords as given in the manual. List values
                            are written as fdf blocks with each element on a
                            separate line, while tuples will write each element
                            in a single line.  ASE units are assumed in the
                            input.
        """
        # Put in the default arguments.
        #parameters = self.default_parameters.__class__(**kwargs)
        parameters = Parameters_deMon(**kwargs)
        
        # Setup the run command based on number of nodes.
        command = os.environ.get('DEMON_COMMAND')
        if command is None:
            mess = "The 'DEMON_COMMAND' environment is not defined."
            raise ValueError(mess)

        label = parameters['label']
        runfile = label + '.inp'
        outfile = label + '.out'

        try:
            command = command % (runfile, outfile)
        except TypeError:
            raise ValueError(
                "The 'DEMON_COMMAND' environment must " +
                "be a format string" +
                " with two string arguments.\n" +
                "Example : 'deMon < ./%s > ./%s'.\n" +
                "Got '%s'" % command)

        # set up basis_path


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

        if atoms is not None:
            self.atoms = atoms.copy()

        # I can't seem to manage to call the super super class...
        #self.FileIOCalculator.Calculator.calculate(self, atoms, properties, system_changes)
        #self.__class__.__bases__[0].__class__.__bases__[0].calculate(self, atoms, properties, system_changes)

        self.write_input(self.atoms, properties, system_changes)
        if self.command is None:
            raise RuntimeError('Please set $%s environment variable ' %
                               ('ASE_' + self.name.upper() + '_COMMAND') +
                               'or supply the command keyword')
        command = self.command #.replace('PREFIX', self.prefix)
        olddir = os.getcwd()

        # basis path
        basis_path = self.parameters['basis_path']
        if basis_path is None:
            basis_path = os.environ.get('DEMON_BASIS_PATH')

        if basis_path is None:
            raise RuntimeError('Please set basis_path keyword,' +
                               ' or the DEMON_BASIS_PATH' +
                               ' environment variable')

        try:
            # link restart file
            value = self.parameters['guess']
            if value.upper() == 'RESTART':
                value2 = self.parameters['deMon_restart_path']
                if os.path.exists(self.directory + '/deMon.rst')\
                        or  os.path.islink(self.directory + '/deMon.rst'):
                    os.remove(self.directory + '/deMon.rst')
                abspath = os.path.abspath(value2)
                
                if os.path.exists(abspath + '/deMon.mem') \
                        or  os.path.islink(abspath + '/deMon.mem'):
                    #os.symlink(abspath +'/deMon.mem', 
                    #           self.directory + '/deMon.rst')
                    shutil.copy(abspath +'/deMon.mem', self.directory + '/deMon.rst')
                else:
                    raise RuntimeError(
                        "{} doesn't exist".format(abspath + '/deMon.rst') )

            # link basis
            abspath = os.path.abspath(basis_path)
            
            if os.path.exists(self.directory + '/BASIS')\
                    or  os.path.islink(self.directory + '/BASIS'):
                os.remove(self.directory + '/BASIS')
                
            if os.path.exists(abspath + '/BASIS')\
                    or  os.path.islink(abspath + '/BASIS'):
                os.symlink(abspath +'/BASIS', 
                           self.directory + '/BASIS')
            else:
                raise RuntimeError(
                    "{} doesn't exist".format(abspath + '/BASIS') )

            # link auxis
            if os.path.exists(self.directory + '/AUXIS')\
                    or  os.path.islink(self.directory + '/AUXIS'):
                os.remove(self.directory + '/AUXIS')

            if os.path.exists(abspath + '/AUXIS')\
                    or  os.path.islink(abspath + '/AUXIS'):
                os.symlink(abspath +'/AUXIS', 
                       self.directory + '/AUXIS')
            else:
                raise RuntimeError(
                    "{} doesn't exist".format(abspath + '/AUXIS') )

            # link ecps
            if os.path.exists(self.directory + '/ECPS')\
                    or  os.path.islink(self.directory + '/ECPS'):
                os.remove(self.directory + '/ECPS')

            if os.path.exists(abspath + '/ECPS')\
                    or  os.path.islink(abspath + '/ECPS'):
                os.symlink(abspath +'/ECPS', 
                           self.directory + '/ECPS')
            else:
                raise RuntimeError(
                    "{} doesn't exist".format(abspath + '/ECPS') )

            # link mcps
            if os.path.exists(self.directory + '/MCPS')\
                    or  os.path.islink(self.directory + '/MCPS'):
                os.remove(self.directory + '/MCPS')

            if os.path.exists(abspath + '/MCPS')\
                    or  os.path.islink(abspath + '/ECPS'):
                os.symlink(abspath +'/MCPS', 
                           self.directory + '/MCPS')
            else:
                raise RuntimeError(
                    "{} doesn't exist".format(abspath + '/MCPS') )

            # go to directory and run calculation
            os.chdir(self.directory)
            errorcode = subprocess.call(command, shell=True)        
        finally:
            os.chdir(olddir)

        if errorcode:
            raise RuntimeError('%s returned an error: %d' %
                               (self.name, errorcode))

        if True: #try:
            self.read_results()
        else: #except:
            with open(self.directory + '/deMon.out', 'r') as f:
                lines = f.readlines()
            debug_lines = 10
            print('##### %d last lines of the deMon.out' % debug_lines)
            for line in lines[-20:]:
                print(line.strip())
            print('##### end of deMon.out')
            raise RuntimeError
            

    def set_label(self, label):
        """Set label directory """

        self.label = label

        self.directory = self.label # why have both?
        if self.directory == '':
            self.directory = os.curdir

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

        filename = self.label + '/deMon.inp'

# 2DO
# On any changes, remove all analysis files.
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

            value = self.parameters['xc']
            self._write_argument('VXCTYPE', value, f)

            value = self.parameters['guess']
            self._write_argument('GUESS', value, f)

            value = self.parameters['print_out']
            if not len(value) == 0:
                self._write_argument('PRINT', value, f)
                f.write('#\n')            

            # obtain forces through a single BOMD step
            value = self.parameters['forces']
            if value:
                self._write_argument('DYNAMICS', ['INT=1', 'MAX=1', 'STEP=0'], f)
                self._write_argument('TRAJECTORY', 'FORCES', f)
                self._write_argument('VELOCITIES', 'ZERO', f)

            # write general input arguments
            self._write_input_arguments(f)
            
            f.write('#\n')

            # write basis set, ecps, mcps, auxis, augment  
            basis = self.parameters['basis']
            if not 'all' in basis:
                basis['all'] = 'DZVP'
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

            # write pickle of Parameters 
            pickle.dump(self.parameters, open(self.directory +"/deMon_parameters.pckl", 'w'))

            # write xyz file for good measure. 
            ase.io.write(self.directory + "/deMon_atoms.xyz", self.atoms)
                
    def read(self, restart_path):
        """Read parameters from directory restart_path."""

        self.set_label(restart_path)

        if not os.path.exists(restart_path + '/deMon.inp'):
            raise ReadError("The restart_path file {} does not exist".format(restart_path))
        
        parameters = pickle.load(open(restart_path + '/deMon_parameters.pckl','r'))
        self.parameters = parameters

        self.atoms = self.deMon_inp_to_atoms(restart_path + '/deMon.inp')
        
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
            line = key.upper()
            line += "    " + str(value).upper()
            f.write(line)
            f.write('\n')

        # for a list, write first argument on the first line, then the rest on new lines
        else:
            line = key
            if not isinstance(value[0], (tuple, list)):
                for i in range(len(value)):
                    line += "  " + str(value[i].upper())
                f.write(line)
                f.write('\n')
            else:
                for i in range(len(value)):
                    for j in range(len(value[i])):
                        line += "  " + str(value[i][j]).upper()
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

    #
    # Analysis routines
    #

    def read_results(self):
        """Read the results from output files.
        """
        self.read_energy()
        self.read_forces(len(self.atoms))
        self.read_eigenvalues()
        #self.read_dipole()

    def read_energy(self):
        """Read energy from deMon's text-output file.
        """
        with open(self.label + '/deMon.out', 'r') as f:
            text = f.read().upper()

        lines = iter(text.split('\n'))

        for line in lines:
            if line.startswith(' TOTAL ENERGY                ='):
                self.results['energy'] = float(line.split()[-1]) * Hartree
                break
        else:
            raise RuntimeError

    def read_forces(self, natoms):
        """Read the forces from the deMon.trj file.
        """
        filename=self.label +'/deMon.trj'
        if isfile(filename):
            with open(filename, 'r') as f:
                lines = f.readlines()

            self.results['forces'] = np.zeros((natoms, 3), float)
        
            start = 8 + 2*natoms 
            for i in range(natoms):
                line = [s for s in lines[i+start].strip().split(' ') if len(s) > 0]
                self.results['forces'][i] = map(float, line[0:3])

            self.results['forces'] *= Hartree / Bohr
        
    def read_eigenvalues(self):
        """Read eigenvalues from the 'deMon.out' file.
        """
        assert os.access(self.label + '/deMon.out', os.F_OK)

        # Read eigenvalues 
        with open(self.label + '/deMon.out', 'r') as f:
            lines = f.readlines()

        # try  PRINT MOE
        eig_alpha, occ_alpha = self.read_eigenvalues_one_spin(lines, 'ALPHA MO ENERGIES',6)
        eig_beta, occ_beta = self.read_eigenvalues_one_spin(lines, 'BETA MO ENERGIES',6)

        # otherwise try PRINT MOS
        if len(eig_alpha) ==0 and len(eig_beta) == 0 :
            eig_alpha, occ_alpha = self.read_eigenvalues_one_spin(lines, 'ALPHA MO COEFFICIENTS',5)
            eig_beta, occ_beta = self.read_eigenvalues_one_spin(lines, 'BETA MO COEFFICIENTS',5)

        self.results['eigenvalues'] = np.array([eig_alpha, eig_beta]) * Hartree
        self.results['occupations'] = np.array([occ_alpha, occ_beta]) * Hartree


 
    def read_eigenvalues_one_spin(self,lines, string, neigs_per_line):
        """  utility method for retreiving eigenvalues after the string "string"
        with neigs_per_line eigenvlaues written per line
        """
        eig=[]
        occ=[]

        skip_line = False
        more_eigs = False

        # find line where the orbitals start
        for i in range(len(lines)):
            if lines[i].rfind(string) > -1:
                ii = i
                more_eigs =True
                break

        while more_eigs:
            # search for two empty lines in a row preceeding a line with numbers
            for i in range(ii+1, len(lines)):
                if len(lines[i].split()) == 0 and \
                        len(lines[i+1].split()) == 0 and \
                        len(lines[i+2].split()) > 0 :
                    ii = i+2
                    break

            # read eigenvalues, occupations
            line = lines[ii].split()
            if len(line) < neigs_per_line:
                #last row
                more_eigs = False
            if line[0] != str(len(eig) +1):
                more_eigs = False
                skip_line = True

            if not skip_line:
                line = lines[ii+1].split()
                for l in line:
                    eig.append(float(l))
                line = lines[ii+3].split()
                for l in line:
                    occ.append(float(l))
                ii = ii +3

        return eig, occ

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



    def deMon_inp_to_atoms(self, filename):
        """  routine to read deMon.inp and convert it to an atoms object
        
        """
        with open(filename, 'r') as f:
            lines = f.readlines()

        # find line where geometry starts
        for i in range(len(lines)):
            if lines[i].rfind('GEOMETRY') > -1:
                if lines[i].rfind('ANGSTROM'):
                    #print('Read coordinates in Ang units')
                    coord_units = 'Ang'
                elif lines.rfind('Bohr'):
                    #print('Read coordinates in Bohr units')
                    coord_units = 'Bohr'
                ii = i
                break

        chemical_symbols=[]
        xyz=[]
        atomic_numbers=[]
        masses=[]

        for i in range(ii+1, len(lines)):
            try:
                line = string.split(lines[i])
                
                for symbol in ase.data.chemical_symbols:
                    found = None
                    if line[0].upper().rfind(symbol.upper()) > -1:
                        found = symbol
                        break

                if found is not None:
                    chemical_symbols.append(found)
                else:
                    break

                xyz.append([float(line[1]), float(line[2]), float(line[3])] )
                
                if len(line) > 4:
                    atomic_numbers.append(int(line[4])) 
                
                if len(line) > 5:
                    masses.append(float(line[5])) 

            except:
                raise RuntimeError

        if coord_units == 'Bohr':
            xyz = xyz * Bohr

        natoms= len(chemical_symbols) 

        # set atoms object
        atoms = ase.Atoms(symbols=chemical_symbols, positions=xyz)

        # if atomic numbers were read in, set them 
        if(len(atomic_numbers) == natoms):
            atoms.set_atomic_numbers(atomic_numbers)
                
        # if masses were read in, set them
        if(len(masses) == natoms):
            atoms.set_masses(masses)
            
        return atoms

""" 
This module defines an ASE interface to gromacs
http://www.gromacs.org

It is VERY SLOW compared to standard Gromacs 
(due to slow formatted io required here).

Mainly stolen from aims.py lj.py  and vasp.py
Markus.Kaukonen@iki.fi
"""

from general import Calculator
import os
import glob
from ase import units
import numpy as np
from ase.io.gromacs import read_gromacs
from ase.io.gromos import write_gromos

string_keys = [
    'define',
    'integrator',
    'nsteps',
    'nstfout',
    'nstlog',
    'nstenergy',
    'energygrps',
    'nstlist',
    'ns_type',
    'pbc',
    'rlist',
    'coulombtype',
    'rcoulomb',
    'vdwtype',
    'rvdw',
    'rvdw_switch',
    'DispCorr',
]

class Gromacs(Calculator):
    """ A calculator using gromacs.org .
    Initializing variables
    and preparing gromacs topology and run input file
    """
    def __init__(self, init_structure_file='gromacs_init.gro', \
                     force_field='oplsaa', water_model='tip3p', \
                     **kwargs):
        self.name = 'Gromacs'
        self.init_structure_file = init_structure_file
        self.force_field = force_field 
        self.water_model = water_model 
        self.string_params = {}
        self.string_params_doc = {}
        for key in string_keys:
            self.string_params[key] = None
            self.string_params_doc[key] = None
        #add comments for gromacs input file
        self.string_params_doc['define'] = \
            'flexible/ rigid water'
        self.string_params_doc['integrator'] = \
            'md: molecular dynamics(Leapfrog), \n' + \
            '; md-vv: molecular dynamics(Velocity Verlet), \n' + \
            '; steep: steepest descent minimization, \n' + \
            '; cg: conjugate cradient minimization \n'
        self.positions = None
        self.atoms = None
        # storage for energy and forces
        self.energy = None
        self.forces = None
        
        self.set(**kwargs)
        #delete possible old files before a new calculation
        files = ['gromacs_confin.gro', 'gromacs_confout.gro', 
                 'gromacs.tpr', 'gromacs.top', 'gromacs.mdp',
                 'tmp_force.del', 'tmp_ene.del', 'energy.xvg',
                 'gromacs.log', 'gromacsEnergy.xvg','gromacsForce.xvg']
        for f in files:
            try:
                os.remove(f)
            except OSError:
                pass
        self.write_parameters()
        self.initialize()

    def set(self, **kwargs):
        """ Setting values for the parameters of the gromacs calculator """
        for key in kwargs:
            if self.string_params.has_key(key):
                self.string_params[key] = kwargs[key]
            else:
                raise TypeError('Parameter not defined: ' + key)

    def run(self):
        """ runs a 0-step gromacs-mdrun with the 
        current atom-configuration """
        delnames = glob.glob('#*')
        try:
            for delname in delnames:
                os.remove(delname)
        except:
            pass
        os.system('mdrun -s gromacs.tpr -o gromacs.trr ' + \
                      ' -e gromacs.edr -g gromacs.log -ffout ' + \
                      ' -rerun gromacs.g96 '+ ' > /dev/null 2>&1')

    def initialize(self):
        """ from coordinates (gromacs_init.gro)
            and gromacs run input file (gromacs.mdp)
            generate topology gromacs.top 
            (with a given force field and water model)
            and run input file gromacs.tpr 
        """
        #generate structure and topology files 
        os.system('pdb2gmx -f ' + self.init_structure_file + \
                      ' -o gromacs.gro ' + \
                      ' -p gromacs.top' + \
                      ' -ff ' + self.force_field + \
                      ' -water ' + self.water_model + ' > /dev/null 2>&1')
        #generate gromacs run input file
        os.system('grompp -f gromacs.mdp ' + \
                      ' -c gromacs.gro ' + \
                      ' -p gromacs.top' + \
                      ' -o gromacs.tpr' + ' > /dev/null 2>&1')
        #get the initial coordinates
        self.atoms = read_gromacs('gromacs.gro') 
        filename = 'inputGenergy.txt'
        output = open(filename,'w')
        output.write('Potential  \n')
        output.write('   \n')
        output.write('   \n')
        output.close()

        filename = 'inputGtraj.txt'
        output = open(filename, 'w')
        output.write('System  \n')
        output.write('   \n')
        output.write('   \n')
        output.close()

    def calculation_required(self, atoms):
        if ((self.atoms != atoms) or
        not(os.path.isfile('gromacs.log'))):
            return True
        else:
            return False

    def update(self, atoms):
        if self.calculation_required(atoms):
            # performs an update of the atoms 
            self.atoms = atoms.copy()
            #must be g96 format for accuracy, alternatively binary formats
            write_gromos('gromacs.g96', atoms)
            # does run to get forces and energies
            self.calculate()

    def get_potential_energy(self, atoms):
        self.update(atoms)
        return self.energy

    def get_forces(self, atoms):
        self.update(atoms)
        return self.forces

    def get_stress(self, atoms):
        return np.zeros(6)

    def set_atoms(self, atoms):
        self.atoms = atoms.copy()
    
    def calculate(self):
        """ runs one step gromacs-mdrun and 
        gets energy and forces
        """
        self.run()        
        self.energy = 0.0
        delnames = glob.glob('#*')
        try:
            for delname in delnames:
                os.remove(delname)
        except:
            pass
        # get energy
        try:
            os.remove('tmp_ene.del')
        except:
            pass

        os.system("g_energy_d -f gromacs.edr -dp"+\
                      " -o gromacsEnergy.xvg < inputGenergy.txt"+\
                      " > /dev/null 2>&1")
        os.system("tail -n 1 gromacsEnergy.xvg > tmp_ene.del")
        line = open('tmp_ene.del', 'r').readline()
        energy = float(line.split()[1])
        #We go for ASE units !
        self.energy = energy * units.kJ / units.mol 
        #self.energy = energy 
        # energies are about 100 times bigger in Gromacs units 
        # when compared to ase units

        #get forces
        try:
            os.remove('tmp_force.del')
        except:
            pass
        #os.system('gmxdump_d -f gromacs.trr > tmp_force.del 2>/dev/null')
        os.system('g_traj -f gromacs.trr -s gromacs.tpr -of '+\
                      ' -fp gromacsForce.xvg < inputGtraj.txt '+\
                      ' > /dev/null 2>&1')
        lines = open('gromacsForce.xvg', 'r').readlines()
        forces = []
        for line in lines:
            if (('#' in line) or ('@' in line)):
                pass
            else:
                #print line
                forces.append(np.array\
                                  ([float(f) for f in line.split()[1:]]))
        #We go for ASE units !
        self.forces = np.array(forces)/ units.nm * units.kJ / units.mol
        #self.forces = np.array(forces)

    def set_own(self, key, value, docstring=""):
        """Set own gromacs input file parameter."""
        self.string_params[key] = value
        self.string_params_doc[key] = docstring

    def write_parameters(self):
        """ Writes run-input file for gromacs (mdrun) """
        prefix = ';'
        filename = 'gromacs.mdp'
        output = open(filename,'w')
        output.write(prefix+\
            '=======================================================\n')
        output.write(prefix+'Gromacs input file \n')
        output.write(prefix+ \
            'Created using the Atomic Simulation Environment (ASE) \n')
        output.write(prefix+\
            '=======================================================\n')
        for key, val in self.string_params.items():
            if val is not None:
                if (self.string_params_doc[key] == None):
                    docstring = ''
                else:
                    docstring = self.string_params_doc[key]
                output.write('%-35s = %s ; %s\n' \
                                 % (key, val, docstring))
        output.close()


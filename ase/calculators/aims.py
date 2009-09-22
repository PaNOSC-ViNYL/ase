"""This module defines an ASE interface to FHI-aims.

Felix Hanke hanke@liverpool.ac.uk
Jonas Bjork j.bjork@liverpool.ac.uk
"""

from general import Calculator
import os
import sys
from os.path import join, isfile, islink

import numpy as np

import ase

float_keys = [
    'charge',
    'charge_mix_param',
    'ini_linear_mix_param',
    'ini_spin_mix_parma',
    'initial_moment',
    'MD_time_step',
    'prec_mix_param',
    'spin_mix_param',
]

exp_keys = [
    'sc_accuracy_eev',
    'sc_accuracy_etot',
    'sc_accuracy_forces',
    'sc_accuracy_rho',
]

string_keys = [
    'KS_method',
    'mixer',
    'output_level',
    'restart',
    'restart_read_only',
    'restart_write_only',
    'spin',
    'xc',
]

int_keys = [
    'empty_states',
    'ini_linear_mixing',
    'max_relaxation_steps',    
    'n_max_pulay',    
    'sc_iter_limit',
]

bool_keys = [
    'compute_forces',
    'compute_kinetic',
    'final_forces_cleaned',
    'MD_clean_rotations',
    'restart_relaxations',
    'use_dipole_correction',
    'vdw_correction_hirshfeld',
]

list_keys = [
    'hartree_convergence_parameter',
    'k_grid',
    'MD_run',
    'mixer_threshold',
    'occupation_type',
    'output',
    'preconditioner',
    'relativistic',
    'relax_geometry',
]

input_keys = [
    'run_command',
    'run_dir',
    'species_dir',
] 

class Aims(Calculator):
    def __init__(self, output_template = 'aims', track_output = False, **kwargs):
        self.name = 'Aims'
        self.float_params = {}
        self.exp_params = {}
        self.string_params = {}
        self.int_params = {}
        self.bool_params = {}
        self.s_bool_params = {}
        self.list_params = {}
        self.input_parameters = {}
        for key in float_keys:
            self.float_params[key] = None
        for key in exp_keys:
            self.exp_params[key] = None
        for key in string_keys:
            self.string_params[key] = None
        for key in int_keys:
            self.int_params[key] = None
        for key in bool_keys:
            self.bool_params[key] = None
        for key in list_keys:
            self.list_params[key] = None
        for key in input_keys:
            self.input_parameters[key] = None
        self.positions = None
        self.atoms = None
        self.run_counts = 0
        self.set(**kwargs)
        self.output_template = output_template
        self.track_output = track_output

    def set(self, **kwargs):
        for key in kwargs:
            if self.float_params.has_key(key):
                self.float_params[key] = kwargs[key]
            elif self.exp_params.has_key(key):
                self.exp_params[key] = kwargs[key]
            elif self.string_params.has_key(key):
                self.string_params[key] = kwargs[key]
            elif self.int_params.has_key(key):
                self.int_params[key] = kwargs[key]
            elif self.bool_params.has_key(key):
                self.bool_params[key] = kwargs[key]
            elif self.list_params.has_key(key):
                self.list_params[key] = kwargs[key]
            elif self.input_parameters.has_key(key):
                self.input_parameters[key] = kwargs[key]
            else:
                raise TypeError('Parameter not defined: ' + key)

    def update(self, atoms):
        if (self.positions is None or
            (self.atoms != atoms) or
            (self.float_params != self.old_float_params) or
            (self.string_params != self.old_string_params) or
            (self.int_params != self.old_int_params) or
            (self.input_parameters != self.old_input_parameters)
            ):
            self.calculate(atoms)

    def calculate(self, atoms):
        """Generate necessary files in the working directory.
        
        If the directory does not exist it will be created.

        """
        positions = atoms.get_positions()
        pbc = atoms.get_pbc()
        have_lattice_vectors = pbc[0] and pbc[1] and pbc[2]
        have_k_grid = not self.list_params['k_grid'] == None

        if have_lattice_vectors and not have_k_grid:
            raise RuntimeError("Found lattice vectors but no k-grid!")
        if not have_lattice_vectors and have_k_grid:
            raise RuntimeError("Found k-grid but no lattice vectors!")
        from ase.io.aims import write_aims
        write_aims('geometry.in', atoms)
        self.write_control()
        self.write_species()
        self.run()
        self.read(atoms)
        self.old_float_params = self.float_params.copy()
        self.old_string_params = self.string_params.copy()
        self.old_int_params = self.int_params.copy()
        self.old_input_parameters = self.input_parameters.copy()
        self.converged = self.read_convergence()

    def run(self):
        if (self.track_output):
            self.out = self.output_template+str(self.run_counts)+'.out'
            self.run_counts += 1
        else:
            self.out = self.output_template+'.out'            
        if self.input_parameters['run_command']:
            aims_command = self.input_parameters['run_command'] 
        elif os.environ.has_key('AIMS_COMMAND'):
            aims_command = os.environ['AIMS_COMMAND']
        else:
            raise RuntimeError("No specification for running FHI-aims. Aborting!")
        aims_command = aims_command + ' > ' 
        if self.input_parameters['run_dir']:
            # file exists: run_dir ??? 
            aims_command = aims_command + self.input_parameters['run_dir'] + '/'
        aims_command = aims_command + self.out
        exitcode = os.system(aims_command)
        if exitcode != 0:
            raise RuntimeError('FHI-aims exited with exit code: %d.  ' % exitcode)

    def write_control(self):
        """Writes the control.in file."""
        control = open('control.in', 'w')
        for key, val in self.float_params.items():
            if val is not None:
                control.write('%-30s%5.6f\n' % (key, val))
        for key, val in self.exp_params.items():
            if val is not None:
                control.write('%-30s%5.2e\n' % (key, val))
        for key, val in self.string_params.items():
            if val is not None:
                control.write('%-30s%s\n' % (key, val))
        for key, val in self.int_params.items():
            if val is not None:
                contol.write('%-30s%d\n' % (key, val))
        for key, val in self.bool_params.items():
            if val is not None:
                if key == 'vdw_correction_hirshfeld':
                    control.write('%-30s\n' % (key))
                elif val:
                    control.write('%-30s.true.\n' % (key))
                else:
                    control.write('%-30s.false.\n' % (key))
        for key, val in self.list_params.items():
            if val is not None:
                control.write('%-30s' % key)
                for ival in val:
                    control.write(str(ival)+' ')
                control.write('\n')
        for key, val in self.input_parameters.items():
            if val is not None:
                continue
        control.write('\n')
        control.close()

    def write_species(self):
        from ase.data import atomic_numbers

        if not self.input_parameters['species_dir']:
            raise RuntimeError('Missing species directory, THIS MUST BE SPECIFIED!')

        control = open('control.in', 'a')
        species_path = self.input_parameters['species_dir']
        symbols = self.atoms.get_chemical_symbols()
        symbols2 = []
        for n, symbol in enumerate(symbols):
            if symbol not in symbols2:
                symbols2.append(symbol)
        for symbol in symbols2:
            file = join(species_path, '%02i_%s_default' % (atomic_numbers[symbol], symbol))
            for line in open(file, 'r'):
                control.write(line)
        control.close()

    def read_energy(self, all=None):
        for line in open(self.out, 'r'):
            if line.rfind('Total energy corrected') > -1:
                E0 = float(line.split()[-2])
            elif line.rfind('Total energy uncorrected') > -1:
                F = float(line.split()[-2])
        energy_free, energy_zero = F, E0
        return [energy_free, energy_zero]

    def read_forces(self, atoms, filename='aims.out', all=False):
        """Method that reads forces from the output file.

        If 'all' is switched on, the forces for all ionic steps
        in the output file will be returned, in other case only the
        forces for the last ionic configuration are returned."""
        lines = open(self.out, 'r').readlines()
        forces = np.zeros([len(atoms), 3])
        for n, line in enumerate(lines):
            if line.rfind('Total atomic forces') > -1:
                for iatom in range(len(atoms)):
                    data = lines[n+iatom+1].split()
                    for iforce in range(3):
                        forces[iatom, iforce] = float(data[2+iforce])
        return forces

    def get_stress(self, atoms):
        raise NotImplementedError('Stresses are not currently available in FHI-aims, sorry. ')


# methods that should be quickly implemented some time, haven't had time yet:
    def read_fermi(self):
        """Method that reads Fermi energy from output file"""
        return

    def read_dipole(self):
        dipolemoment = np.array([0, 0, 0])
        return dipolemoment

    def read_magnetic_moment(self):
        return

    def read_convergence(self):
        return

    def read_eigenvalues(self, kpt=0, spin=0):
        return 


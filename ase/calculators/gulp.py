"""This module defines an ASE interface to GULP.

Written by:

Andy Cuko <andi.cuko@upmc.fr>
Antoni Macia <tonimacia@gmail.com>

EXPORT ASE_GULP_COMMAND="/path/to/gulp < PREFIX.gin > PREFIX.got"

Keywords
Options

"""
import os
import numpy as np
from ase.units import eV, Ang
#from ase.io.gulp_got import read_gulp
from ase.geometry import cellpar_to_cell, cell_to_cellpar
from ase.calculators.calculator import FileIOCalculator, ReadError

from ase.data import chemical_symbols
import re


class GULP(FileIOCalculator):
    implemented_properties = ['energy', 'forces']
    command = 'gulp < PREFIX.gin > PREFIX.got'
    default_parameters = dict(
        keywords='opti conp comp rfo',
        options=[],
        shel=[],
        library="ffsioh.lib",
        conditions=None
        )

#conditions=[['O', 'default', 'O1'], ['O', 'O2', 'H', '<', '1.6']]

    def __init__(self, restart=None, ignore_bad_restart_file=False,
                 label='gulp', atoms=None, optimized=None, Gnorm=1000.0, Steps=1000, conditions=None, **kwargs):
        """Construct GULP-calculator object."""
        FileIOCalculator.__init__(self, restart, ignore_bad_restart_file,
                                  label, atoms, **kwargs)
        self.optimized  = optimized
        self.Gnorm      = Gnorm
        self.Steps      = Steps
        self.conditions = conditions
        self.library_check()
        self.atom_types = []

    def set(self, **kwargs):
        changed_parameters = FileIOCalculator.set(self, **kwargs)
        if changed_parameters:
            self.reset()

    def write_input(self, atoms, properties=None, system_changes=None):
        FileIOCalculator.write_input(self, atoms, properties, system_changes)
        p = self.parameters

        # Build string to hold .gin input file:
        s = p.keywords
        s += '\ntitle\nASE calculation\nend\n\ncart\n'

#       IMPLEMENT HERE PERIODIC SYSTMES
#        for v, p in zip(atoms.cell, atoms.pbc):
#            if p:
#                s += ' {0} {1} {2}\n'.format(*v)
        # Write coordinates:
        if self.conditions is not None:
            c = self.conditions
            c.set_atoms(self.atoms)
            c.apply_rules()
            labels = c.get_atoms_labels()
            self.atom_types = c.get_atom_types()
        else:
            labels = self.atoms.get_chemical_symbols()
        for xyz, symbol in zip(atoms.positions, labels):
            s += ' {0:2} core {1}  {2}  {3}\n'.format(symbol, *xyz)
            if symbol in p.shel:
                s += ' {0:2} shel {1}  {2}  {3}\n'.format(symbol, *xyz)

        s += '\nlibrary {0}\n'.format(p.library)
        if p.options:
            for t in p.options:
                s += '%s\n' % t
        with open(self.prefix + '.gin', 'w') as f:
            f.write(s)

    def read_results(self):
        FileIOCalculator.read(self, self.label)
        if not os.path.isfile(self.label + '.got'):
            raise ReadError

        with open(self.label + '.got') as f:
            lines = f.readlines()

        cycles = -1
        self.optimized = None
        for i, line in enumerate(lines):
            if line.find('Final energy =') != -1:
                self.results['energy'] = float(line.split()[-2]) * eV

            elif line.find('Optimisation achieved') != -1:
                self.optimized = True

            elif line.find('Final Gnorm') != -1:
                self.Gnorm = float(line.split()[-1])

            elif line.find('Cycle:') != -1:
                cycles += 1

            elif line.find('Final Cartesian derivatives') != -1:
                s = i + 5
                forces = []
                while(True):
                    s = s + 1
                    if lines[s].find("------------") != -1:
                        break
                    if lines[s].find(" s ") != -1:
                        continue
                    g = lines[s].split()[3:6]
                    G = [-float(x) * eV/Ang for x in g]
                    forces.append(G)
                forces = np.array(forces)
                self.results['forces'] = forces

            elif line.find('Final cartesian coordinates of atoms') != -1:
                s = i + 5
                positions = []
                while(True):
                    s = s + 1
                    if lines[s].find("------------") != -1:
                        break
                    if lines[s].find(" s ") != -1:
                        continue
                    xyz = lines[s].split()[3:6]
                    XYZ = [float(x) * Ang for x in xyz]
                    positions.append(XYZ)
                positions = np.array(positions)
                self.atoms.set_positions(positions)

        self.Steps=cycles

    def get_opt_state(self):
        return self.optimized

    def get_opt_steps(self):
        return self.Steps

    def get_Gnorm(self):
        return self.Gnorm

    def library_check(self):
        if self.parameters['library'] is not None:
            try:
                path = os.environ['GULP_LIB'] + '/'
            except:
                raise RuntimeError("Be sure to have set correctly $GULP_LIB or to have the force field library.")

class Conditions:
    """ This class is made to return return an array similar to atoms.get_chemical_symbols()
    via get_atoms_labels() funtion but with atomic labels in stead of atomic symbols.
    This is usefull when you need to use calculators like gulp or lammps that uses
    force fileds. Some force fields can have different atom type for the same element.
    In this class you can create a set_rule() function that assigns labels according to
    structural criteria."""

    def __init__(self, atoms=None, rule=None):
        self.atoms = atoms
        self.atoms_symbols = atoms.get_chemical_symbols()
        self.atoms_labels  = atoms.get_chemical_symbols()
        self.rules = rule
        self.atom_types = []

    def min_distance_rule(self, sym1, sym2, ifyeslabel1 , ifyeslabel2, ifnolabel1):
        """This function is a rule that allows to define atom labels (like O1, O2, O_H etc..)
        starting from element symbols of an Atoms object that a force field can use and according
        to distance parameters.

        Example:
        atoms = read('some_xyz_format.xyz')
        a = Conditions(atoms)
        a.set_min_distance_rule("O", "H", "O2", "H", "O1")
        new_atoms_labels = a.get_atom_labels()

        In the example oxygens O are going to be labeled as O2 if they are close to a
        hydrogen atom othewise are labeled O1."""
        #self.atom_types is a list of element types  used instead of element
        #symbols in orger to track the changes made. Take care of this because
        # is very important.. gulp_read function that parse the output
        # has to know which atom_type it has to associate with which atom_symbol
        # Example: [['O','O1','O2'],['H', 'H_C', 'H_O']]
        # this beacuse Atoms oject accept only atoms symbols
        self.atom_types.append([sym1, ifyeslabel1 , ifnolabel1])
        self.atom_types.append([sym2, ifyeslabel2])

        dist_mat = self.atoms.get_all_distances()
        index_assiged_sym1=[]
        index_assiged_sym2=[]

        for i in range(len(self.atoms_symbols)):
            if self.atoms_symbols[i] is sym2:
                dist_12 = 1000
                index_assiged_sym2.append(i)
                for t in range(len(self.atoms_symbols)):
                    if self.atoms_symbols[t] is sym1 and dist_mat[i, t] < dist_12 and t not in index_assiged_sym1:
                        dist_12 = dist_mat[i, t]
                        closest_sym1_index = t
                index_assiged_sym1.append(closest_sym1_index)

        for s in range(len(self.atoms_symbols)):
            if s in index_assiged_sym1:
                self.atoms_labels[s] = ifyeslabel1
            elif s not in index_assiged_sym1 and self.atoms_symbols[s] is sym1:
                self.atoms_labels[s] = ifnolabel1
            elif s in index_assiged_sym2:
                self.atoms_labels[s] = ifyeslabel2
            else:
                pass

    def set_atoms(self, atoms):
        self.atoms = atoms

    def apply_rules(self):
        for rule in self.rules:
            getattr(self, rule)(*self.rules[rule])

    def get_atom_types(self):
        return self.atom_types

    def set_rule(self, rule):
        self.rules = rule

    def get_atoms_labels(self):
        labels = np.array(self.atoms_labels)
        return labels

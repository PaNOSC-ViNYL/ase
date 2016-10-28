"""This module defines an ASE interface to MOPAC.

Set $ASE_MOPAC_COMMAND to something like::
    
    LD_LIBRARY_PATH=/path/to/lib/ \
    MOPAC_LICENSE=/path/to/license \
    /path/to/MOPAC2012.exe PREFIX.mop 2> /dev/null

"""
import os
import re

import numpy as np

from ase.calculators.calculator import FileIOCalculator, ReadError, Parameters
from ase.units import kcal, mol


class MOPAC(FileIOCalculator):
    implemented_properties = ['energy', 'forces']
    command = 'mopac PREFIX.mop 2> /dev/null'

    default_parameters = dict(
        method='PM7',
        task='1SCF GRADIENTS',
        relscf=0.0001)

    methods = ['AM1', 'MNDO', 'MNDOD', 'PM3', 'PM6', 'PM6-D3', 'PM6-DH+',
               'PM6-DH2', 'PM6-DH2X', 'PM6-D3H4', 'PM6-D3H4X', 'PMEP', 'PM7',
               'PM7-TS', 'RM1']

    tasks = ['1SCF', 'RHF', 'UHF', 'SINGLET', 'DOUBLET', 'TRIPLET', ]

    def __init__(self, restart=None, ignore_bad_restart_file=False,
                 label='mopac', atoms=None, **kwargs):
        """Construct MOPAC-calculator object.

        Parameters
        ==========
        label: str
            Prefix for filenames (label.mop, label.out, ...)

        Examples
        ========
        Use default values to do a single SCF calculation and print
        the forces (task='1SCF GRADIENTS')

        >>> from ase.build import molecule
        >>> from ase.calculators.mopac import MOPAC
        >>> atoms = molecule('O2')
        >>> atoms.calc = MOPAC('label'='O2')
        >>> atoms.get_potential_energy()
        >>> eigs = calc.get_eigenvalues()
        >>> somos = calc.get_somo_levels()
        >>> homo, lumo = calc.get_homo_lumo_levels()

        Use the internal geometry optimization of Mopac
        >>> calc = MOPAC(label='H2', task='GRADIENTS')
        >>> atoms.set_calculator(calc)
        >>> atoms.get_potential_energy()

        """
        FileIOCalculator.__init__(self, restart, ignore_bad_restart_file,
                                  label, atoms, **kwargs)

    def set(self, **kwargs):
        changed_parameters = FileIOCalculator.set(self, **kwargs)
        if changed_parameters:
            self.reset()

    def write_input(self, atoms, properties=None, system_changes=None):
        FileIOCalculator.write_input(self, atoms, properties, system_changes)
        p = self.parameters

        # Build string to hold .mop input file:
        s = p.method + ' ' + p.task + ' '
                
        if p.relscf:
            s += 'RELSCF={0} '.format(p.relscf)
            
        # Write charge:
        charge = atoms.get_initial_charges().sum()
        if charge != 0:
            s += 'CHARGE={0} '.format(int(round(charge)))
        
        magmom = int(round(abs(atoms.get_initial_magnetic_moments().sum())))
        if magmom:
            s += (['DOUBLET', 'TRIPLET', 'QUARTET', 'QUINTET'][magmom - 1] +
                  ' UHF ')
            
        s += '\nTitle: ASE calculation\n\n'

        # Write coordinates:
        for xyz, symbol in zip(atoms.positions, atoms.get_chemical_symbols()):
            s += ' {0:2} {1} 1 {2} 1 {3} 1\n'.format(symbol, *xyz)

        for v, p in zip(atoms.cell, atoms.pbc):
            if p:
                s += 'Tv {0} {1} {2}\n'.format(*v)

        with open(self.label + '.mop', 'w') as f:
            f.write(s)

    def read(self, label):
        FileIOCalculator.read(self, label)
        if not os.path.isfile(self.label + '.out'):
            raise ReadError

        with open(self.label + '.out') as f:
            lines = f.readlines()

        self.parameters = Parameters(task='', method='')
        p = self.parameters
        parm_line = self.read_parameters(lines)
        for keyword in parm_line.split():
            if 'RELSCF' in keyword:
                p.relscf = float(keyword.split('=')[-1])
            elif keyword in self.methods:
                p.method = keyword
            else:
                p.task += keyword + ' '

        p.task.rstrip()
        # atoms
#        self.atoms =
#        self.parameters = {'task': , 'method': }
#        self.read_results()

    def get_index(self, lines, pattern):
        for i, line in enumerate(lines):
            if line.find(pattern) != -1:
                return i

    def read_atoms_from_file(self, lines):
        # first try to read from final point (last image)
        read_it = False
        i = self.get_index(lines, 'FINAL  POINT  AND  DERIVATIVES')
        if i is None:
            assert 0, 'Not implemented'

        i = self.get_index(lines[i:], 'CARTESIAN COORDINATES')
        j = i + 2
        symbols = []
        positions = []
        while line[j].strip():
            l = line.split()
            symbols += l[1]
            positions += [float(c) for c in l[2:]]

        atoms = Atoms(symbols=symbols, positions=positions)

    def read_parameters(self, lines):
        """Read the single line that defines a Mopac calculation
           lines:
        """
        for i, line in enumerate(lines):
            if line.find('CALCULATION DONE:') != -1:
                istart = i
                break

        lines1 = lines[istart:]
        for i, line in enumerate(lines1):
            if line.find('****') != -1:
                return lines1[i + 1]

    def read_results(self):
        FileIOCalculator.read(self, self.label)
        if not os.path.isfile(self.label + '.out'):
            raise ReadError

        with open(self.label + '.out') as f:
            lines = f.readlines()

        for i, line in enumerate(lines):
            if line.find('TOTAL ENERGY') != -1:
                self.results['energy'] = float(line.split()[3])
            elif line.find('FINAL HEAT OF FORMATION') != -1:
                self.final_hof = float(line.split()[5]) * kcal / mol
            elif line.find('NO. OF FILLED LEVELS') != -1:
                self.no_occ_levels = int(line.split()[-1])
            elif line.find('NO. OF ALPHA ELECTRON') != -1:
                self.no_alpha_electrons = int(line.split()[-1])
                self.no_beta_electrons = int(lines[i+1].split()[-1])
            elif line.find('FINAL  POINT  AND  DERIVATIVES') != -1:
                forces = [-float(line.split()[6])
                          for line in lines[i + 3:i + 3 + 3 * len(self.atoms)]]
                forces = np.array(forces).reshape((-1, 3)) * kcal / mol
                self.results['forces'] = forces
            elif line.find('EIGENVALUES') != -1:
                if line.find('ALPHA') != -1:
                    self.nspins = 2
                    j = i + 1
                    eigs_alpha = []
                    while lines[j].strip():
                        eigs_alpha += [float(e) for e in lines[j].split()]
                        j += 1
                elif line.find('BETA') != -1:
                    j = i + 1
                    eigs_beta = []
                    while lines[j].strip():
                        eigs_beta += [float(e) for e in lines[j].split()]
                        j += 1
                    eigs = np.array([eigs_alpha, eigs_beta]).reshape(2, 1, -1)
                    self.eigenvalues = eigs
                else:
                    self.nspins = 1
                    eigs = []
                    j = i + 1
                    while lines[j].strip():
                        eigs += [float(e) for e in lines[j].split()]
                        j += 1
                    self.eigenvalues = np.array(eigs).reshape(1, 1, -1)

    def get_eigenvalues(self, kpt=0, spin=0):
        print (self.eigenvalues.shape)
        return self.eigenvalues[spin, kpt]

    def get_homo_lumo_levels(self):
        if self.nspins == 1:
            n = self.no_occ_levels
            eigs = self.eigenvalues
            return np.array([eigs[0, 0, n - 1], eigs[0, 0, n]])
        else:
            na = self.no_alpha_electrons
            nb = self.no_beta_electrons
            eigs = self.eigenvalues
            eah, eal = eigs[0, 0, na - 1: na + 1]
            ebh, ebl = eigs[1, 0, nb - 1: nb + 1]
            return np.array([max(eah, ebh), min(eal, ebl)])

    def get_somo_levels(self):
        assert self.nspins == 2
        eigs = self.eigenvalues
        na, nb = self.no_alpha_electrons, self.no_beta_electrons
        return np.array([eigs[0, 0, na - 1], eigs[1, 0, nb - 1]])

    def get_final_heat_of_formation(self):
        return self.final_hof

"""This module defines an ASE interface to MOPAC.

Set $ASE_MOPAC_COMMAND to something like::
    
    LD_LIBRARY_PATH=/path/to/lib/ \
    MOPAC_LICENSE=/path/to/license \
    /path/to/MOPAC2012.exe PREFIX.mop 2> /dev/null

"""
import os
import re

import numpy as np

from ase.calculators.calculator import FileIOCalculator, ReadError
from ase.units import kcal, mol


class MOPAC(FileIOCalculator):
    implemented_properties = ['energy', 'forces']
    command = 'mopac PREFIX.mop 2> /dev/null'

    default_parameters = dict(
        method='PM7',
        task='1SCF GRADIENTS',
        relscf=0.0001)

    def __init__(self, restart=None, ignore_bad_restart_file=False,
                 label='mopac', atoms=None, **kwargs):
        """Construct MOPAC-calculator object."""
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

    def read_results(self):
        FileIOCalculator.read(self, self.label)
        if not os.path.isfile(self.label + '.out'):
            raise ReadError

        with open(self.label + '.out') as f:
            lines = f.readlines()

        p = re.compile(r"Empirical Formula:\s*.*=\s*(?P<na>[0-9]+)\s*atoms")
        for line in lines:
            m = re.search(p, line)
            if m:
                natoms = int(m.group('na'))
                break

        for i, line in enumerate(lines):
            if line.find('TOTAL ENERGY') != -1:
                self.results['energy'] = float(line.split()[3])
            elif line.find('FINAL HEAT OF FORMATION') != -1:
                self.final_hof = float(line.split()[5]) * kcal /  mol
            elif line.find('NO. OF FILLED LEVELS') != -1:
                self.no_occ_levels = int(line.split()[-1])
            elif line.find('NO. OF ALPHA ELECTRON') != -1:
                self.no_alpha_electrons = int(line.split()[-1])
                self.no_beta_electrons = int(lines[i+1].split()[-1])
            elif line.find('FINAL  POINT  AND  DERIVATIVES') != -1:
                forces = [-float(line.split()[6])
                          for line in lines[i + 3:i + 3 + 3 * natoms]]
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

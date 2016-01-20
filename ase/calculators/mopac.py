"""This module defines an ASE interface to MOPAC.

Set $ASE_MOPAC_COMMAND to something like::
    
    LD_LIBRARY_PATH=/path/to/lib/ \
    MOPAC_LICENSE=/path/to/license \
    /path/to/MOPAC2012.exe PREFIX.mop 2> /dev/null

"""
import os

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
            
        for i, line in enumerate(lines):
            if line.find('TOTAL ENERGY') != -1:
                self.results['energy'] = float(line.split()[3])
            elif line.find('FINAL  POINT  AND  DERIVATIVES') != -1:
                forces = [-float(line.split()[6])
                          for line in lines[i + 3:i + 3 + 3 * len(self.atoms)]]
                forces = np.array(forces).reshape((-1, 3)) * kcal / mol
                self.results['forces'] = forces
                break

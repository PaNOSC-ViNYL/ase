"""This module defines an ASE interface to NWchem

http://www.nwchem-sw.org/
"""
import os
import sys

import numpy as np

from ase.units import Hartree, Bohr
from ase.io.nwchem import write_nwchem, read_nwchem
from ase.calculators.general import Calculator

class NWchem(Calculator):
    def __init__(self, label='nwchem',
                 calculate_energy='nwchem', 
                 xc='LDA',
                 convergence = {'energy'  : 1.e-6,
                                'density' : 1.e-5,
                                'gradient': 1.e-5},
                 basis='3-21G',
                 charge=None,
                 ):
        self.label = label
        self.converged = False
        
        # set calculators for energy and forces
        self.calculate_energy = calculate_energy

        # turbomole has no stress
        self.stress = np.empty((3, 3))
        
        self.charge = charge
        self.xc = xc
        self.basis = basis
        self.convergence = convergence

        # atoms must be set
        self.atoms = None
        
    def execute(self, command):
        from subprocess import Popen, PIPE
        try:
            # the sub process gets started here
            proc = Popen([command], shell=True, stderr=PIPE)
            error = proc.communicate()[1]
            # check the error output
            if 'abnormally' in error:
                raise OSError(error)
#            print 'TM command: ', command, 'successfully executed'
        except OSError, e:
            print >> sys.stderr, 'Execution failed:', e
            sys.exit(1)

    def get_forces(self, atoms):
        self.get_potential_energy(atoms)
        return self.forces

    def get_potential_energy(self, atoms):
        # update atoms
        self.set_atoms(atoms)
        # if update of energy is neccessary
        if self.energy is None or self.forces is None:
            # write input file
            f = open(self.label + '.nw', 'w')
            if self.charge:
                f.write('charge ' + str(self.charge) + '\n')
            write_nwchem(f, atoms)
            basis = '\nbasis\n'
            basis += '  * library ' + self.basis + '\n'
            basis += 'end\n'
            f.write(basis)

            if self.xc == 'RHF':
                task = 'scf'
            else:
                task = 'dft'
                nwchem_xc_map = {
                    'LDA' : 'pw91lda',
                    'PBE' : 'xpbe96 cpbe96',
                    }
                if self.xc in nwchem_xc_map:
                    xc = nwchem_xc_map[self.xc]
                else:
                    xc = self.xc
                f.write('\ndft\n')
                f.write('  xc ' + xc + '\n')
                for key in self.convergence:
                    f.write('  convergence ' + key + ' ' +
                            str(self.convergence[key]) + '\n')
                f.write('end\n')

            f.write('\ntask ' + task + ' gradient\n')
            f.close()

            # calculate energy
            self.output = self.label + '.out'
            self.execute(self.calculate_energy + ' ' +
                         self.label + '.nw > ' + self.output)
            # read output
            self.atoms = read_nwchem(self.output)
            self.read_energy()
            self.read_forces()
        else:
            print 'taking old values (E)'

        return self.energy

    def read_energy(self):
        """Read Energy from nwchem output file."""
        text = open(self.output, 'r').read()
        lines = iter(text.split('\n'))

        # Energy:
        for line in lines:
            estring = 'Total '
            if self.xc == 'RHF':
                estring += 'SCF'
            else:
                estring += 'DFT'
            estring += ' energy'
            if line.find(estring) >=0:
                energy = float(line.split()[4])
        self.energy = energy

    def read_forces(self):
        """Read Forces from nwchem output file."""
        file = open(self.output, 'r')
        lines = file.readlines()
        file.close()

        for i, line in enumerate(lines):
            if line.find('ENERGY GRADIENTS') >=0:
                gradients = []
                for j in range(i + 4, i + 4 + len(self.atoms)):
                    word = lines[j].split()
                    gradients.append([float(word[k]) for k in range(5,8)])
        self.forces =  - np.array(gradients) * Hartree / Bohr

    def set_atoms(self, atoms):
        if self.atoms == atoms:
            return

        self.atoms = atoms
        self.energy = None
        self.forces = None

    def update(self, atoms):
        self.set_atoms(atoms)

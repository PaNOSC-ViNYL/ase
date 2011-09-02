"""This module defines an ASE interface to NWchem

http://www.nwchem-sw.org/
"""
import os
import sys

import numpy as np

from ase.units import Hartree, Bohr
from ase.io.nwchem import write_nwchem, read_nwchem
from ase.calculators.general import Calculator

class KPoint:
    def __init__(self, s):
        self.s = s
        self.eps_n = []
        self.f_n = []

class NWchem(Calculator):
    def __init__(self, label='nwchem',
                 calculate_energy='nwchem', 
                 xc='LDA',
                 convergence = {'energy'  : 1.e-6,
                                'density' : 1.e-5,
                                'gradient': 1.e-5},
                 maxiter = 120,
                 basis='3-21G',
                 ecp=None,
                 so=None,
                 charge=None,
                 multiplicity = 1,
                 spinorbit=False,
                 ):
        self.label = label
        self.converged = False
        
        # set calculators for energy and forces
        self.calculate_energy = calculate_energy

        # does nwchem have stress ???
        self.stress = np.empty((3, 3))
        
        self.charge = charge
        self.xc = xc
        self.basis = basis
        self.ecp = ecp
        self.so = so
        self.convergence = convergence
        self.maxiter = maxiter
        self.multiplicity = multiplicity
        self.spinorbit = spinorbit

        # atoms must be set
        self.atoms = None
        
    def execute(self, command):
        from subprocess import Popen, PIPE
        try:
            # the sub process gets started here
            proc = Popen([command], shell=True, stderr=PIPE)
            error = proc.communicate()[1]
            if error:
                raise OSError(error + '\ncheck ' + self.output)
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
            if self.charge is not None:
                f.write('charge ' + str(self.charge) + '\n')
            write_nwchem(f, atoms)

            def format_basis_set(string, tag='basis'):
                formatted = tag + '\n'
                lines = string.split('\n')
                if len(lines) > 1:
                    formatted += string
                else:
                    formatted += '  * library '  + string + '\n'
                return formatted + 'end\n'
            basis = format_basis_set(self.basis)
            if self.ecp is not None:
                basis += format_basis_set(self.ecp, 'ecp')
            if self.so is not None:
                basis += format_basis_set(self.so, 'so')
            f.write(basis)

            if self.xc == 'RHF':
                task = 'scf'
            else:
                if self.spinorbit:
                    task = 'sodft'
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
                f.write('  mult ' + str(self.multiplicity) + '\n')
                f.write('  xc ' + xc + '\n')
                f.write('  iterations ' + str(self.maxiter) + '\n')
                for key in self.convergence:
                    f.write('  convergence ' + key + ' ' +
                            str(self.convergence[key]) + '\n')
                f.write('end\n')

#            f.write('\ntask ' + task + ' gradient\n')
            f.write('\ntask ' + task + ' energy\n')
            f.close()

            # calculate energy
            self.output = self.label + '.out'
            self.execute(self.calculate_energy + ' ' +
                         self.label + '.nw > ' + self.output)
            # read output
            self.atoms = read_nwchem(self.output)
            self.read_energy()
#            self.read_forces()
        else:
            print 'taking old values (E)'

        return self.energy * Hartree

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
                break
        self.energy = energy

        # Eigenstates
        spin = -1
        kpts = []
        for line in lines:
            if line.find('Molecular Orbital Analysis') >= 0:
                spin += 1
                kpts.append(KPoint(spin))
            if spin >= 0:
                if line.find('Vector') >= 0:
                    line = line.lower().replace('d', 'e')
                    line = line.replace('=', ' ')
                    word = line.split()
                    kpts[spin].f_n.append(float(word[3]))
                    kpts[spin].eps_n.append(float(word[5]))
        self.kpts = kpts
        
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

    def get_eigenvalues(self, kpt=0, spin=0):
        """Return eigenvalue array."""
        return np.array(self.kpts[spin].eps_n) * Hartree
    
    def get_occupation_numbers(self, kpt=0, spin=0):
        """Return occupation number array."""
        return self.kpts[spin].f_n

    def get_number_of_spins(self):
        """Return the number of spins in the calculation.

        Spin-paired calculations: 1, spin-polarized calculation: 2."""
        return len(self.kpts)

    def get_spin_polarized(self):
        """Is it a spin-polarized calculation?"""
        return len(self.kpts) == 2
        
    def set_atoms(self, atoms):
        if self.atoms == atoms:
            return

        self.atoms = atoms
        self.energy = None
        self.forces = None

    def update(self, atoms):
        self.set_atoms(atoms)

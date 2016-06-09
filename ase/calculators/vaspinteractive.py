from __future__ import print_function
from subprocess import Popen, PIPE
from ase.calculators.general import Calculator
from ase.io import write
from os import environ
import atexit

import numpy as np


class VaspInteractive(Calculator):
    name = "VaspInteractive"
    implemented_properties = ['energy', 'forces', 'stress']

    def __init__(self, txt="interactive.log", print_log=True, process=None,
                 command=None, exitcleanly=False, path="./"):
        self.process = process
        self.path = path
        if txt is not None:
            self.txt = open(self.path + txt, "a")
        else:
            self.txt = None
        self.print_log = print_log
        self.positions = None
        self.converged = False
        self.command = command
        self.atoms = None
        self.exitcleanly = exitcleanly

    def set_atoms(self, atoms):
        if (atoms != self.atoms):
            self.converged = None
        self.atoms = atoms.copy()

    def get_atoms(self):
        atoms = self.atoms.copy()
        atoms.set_calculator(self)
        return atoms

    def get_property(self, name, atoms=None, allow_calculation=True):
        if name not in self.implemented_properties:
            raise NotImplementedError
        if atoms is None:
            atoms = self.atoms
        if self.calculation_required(atoms):
            if allow_calculation:
                self.update(atoms)
        if name == 'energy':
            return self.energy_free
        if name == 'forces':
            return self.forces
        if name == 'stress':
            return self.stress

    def get_executable(self):
        if self.command is not None:
            return "cd {0} && {1} ".format(self.path, self.command)
        elif 'VASP_COMMAND' in environ:
            return "cd {0} && {1}".format(self.path,
                                          environ['VASP_COMMAND'])
        else:
            raise RuntimeError('Please set either command in calculator'
                               ' or VASP_COMMAND environment variable')

    def _stdin(self, text, ending="\n"):
        if self.txt is not None:
            self.txt.write(text + ending)
        if self.print_log:
            print(text, end=ending)
        self.process.stdin.write(text + ending)

    def _stdout(self, text):
        if self.txt is not None:
            self.txt.write(text)
        if self.print_log:
            print(text, end="")

    def first_step(self, atoms):
        self._stdout("Writing Initial POSCAR\n")
        write(self.path + "POSCAR", atoms)
        self._stdout("Starting VASP for initial step...\n")
        self.process = Popen([self.get_executable()], stdout=PIPE,
                             stdin=PIPE, stderr=PIPE, shell=True)
        if self.exitcleanly:
            atexit.register(self.close)
        while self.process.poll() is None:
            text = self.process.stdout.readline()
            self._stdout(text)
            if "POSITIONS: reading from stdin" in text:
                self.converged = True
                break

    def continue_stepping(self, atoms):
        self._input_positions(atoms)
        while self.process.poll() is None:
            text = self.process.stdout.readline()
            self._stdout(text)
            if "POSITIONS: reading from stdin" in text:
                self.converged = True
                break

    def _input_positions(self, atoms):
        self._stdout("Inputting positions...\n")
        for atom in atoms.get_scaled_positions():
            self._stdin('{0} {1} {2}'.format(*atom))

    def close(self):
        self._stdout('Attemping to close VASP cleanly\n')
        stopcar = open(self.path + "STOPCAR", "w")
        stopcar.write("LABORT = .TRUE.")
        stopcar.close()
        self._input_positions(self.atoms)
        while self.process.poll() is None:
            text = self.process.stdout.readline()
            self._stdout(text)
            if "POSITIONS: reading from stdin" in text:
                self.converged = True
                break
        self._input_positions(self.atoms)
        while self.process.poll() is None:
            text = self.process.stdout.readline()
            self._stdout(text)
            if "deleting file STOPCAR" in text:
                self.converged = True
                break
        self._stdout("VASP has been closed\n")

    def calculation_required(self, atoms, quantities=None):
        return ((self.positions == atoms.positions).all() and self.converged)

    def update(self, atoms):
        if (self.positions is None):
            self.first_step(atoms)
        elif (((self.positions != atoms.positions).any()) or
              (not self.converged)):
            self.continue_stepping(atoms)
        self.read(atoms)
        self.results = {'energy': self.energy_free,
                        'forces': self.forces,
                        'stress': self.stress}

    def read_energy(self):
        [energy_free, energy_zero] = [0, 0]
        for line in open(self.path + 'vasprun.xml', 'r'):
            if '<i name="e_fr_energy">' in line:
                energy_free = float(line.split()[2])
            if '<i name="e_wo_entrp">' in line:
                energy_zero = float(line.split()[2])
        return [energy_free, energy_zero]

    def read_forces(self, atoms):
        start_tag = '<varray name="forces" >'
        end_tag = '</varray>'
        forces = []
        counter = None
        for line in open(self.path + 'vasprun.xml'):
            if start_tag in line:
                forces = [None] * len(atoms)
                counter = 0
            elif counter is None:
                pass
            elif counter == len(atoms):
                counter = None
                assert (end_tag in line)
            elif counter < len(atoms):
                forces[counter] = [float(x) for x in line.split()[1:4]]
                counter += 1
        return np.array(forces)

    def read_stress(self, atoms):
        start_tag = '<varray name="stress" >'
        end_tag = '</varray>'
        stress = [0, 0, 0, 0, 0, 0]
        counter = None
        for line in open(self.path + 'vasprun.xml'):
            if start_tag in line:
                counter = 0
            elif counter is None:
                pass
            elif counter == 3:
                counter = None
                assert (end_tag in line)
            elif counter == 0:
                values = [float(x) for x in line.split()[1:4]]
                stress[0] = values[0]  # XX
                counter += 1
            elif counter == 1:
                values = [float(x) for x in line.split()[1:4]]
                stress[1] = values[1]  # YY
                stress[5] = values[0]  # XY
                counter += 1
            elif counter == 2:
                values = [float(x) for x in line.split()[1:4]]
                stress[2] = values[2]  # ZZ
                stress[3] = values[1]  # YZ
                stress[4] = values[0]  # XZ
                counter += 1
        return np.array(stress)

    def read(self, atoms):
        self.positions = atoms.get_positions()
        self.energy_free, self.energy_zero = self.read_energy()
        self.forces = self.read_forces(atoms)
        self.stress = self.read_stress(atoms)
        self.atoms = atoms.copy()

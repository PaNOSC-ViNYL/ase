from __future__ import print_function
from subprocess import Popen, PIPE
from ase.calculators.calculator import Calculator
from ase.io import read, write
import os
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
        self.command = command
        self.atoms = None
        self.exitcleanly = exitcleanly

    def get_executable(self):
        if self.command is not None:
            return "cd {0} && {1} ".format(self.path, self.command)
        elif 'VASP_COMMAND' in os.environ:
            return "cd {0} && {1}".format(self.path,
                                          os.environ['VASP_COMMAND'])
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
        if os.path.isfile('STOPCAR'):
            os.remove('STOPCAR')
        self._stdout("Writing Initial POSCAR\n")
        write(os.path.join(self.path, "POSCAR"), atoms)
        self._stdout("Starting VASP for initial step...\n")
        self.process = Popen([self.get_executable()], stdout=PIPE,
                             stdin=PIPE, stderr=PIPE, shell=True)
        if self.exitcleanly:
            atexit.register(self.close)
        while self.process.poll() is None:
            text = self.process.stdout.readline()
            self._stdout(text)
            if "POSITIONS: reading from stdin" in text:
                break

    def continue_stepping(self, atoms):
        self._input_positions(atoms)
        while self.process.poll() is None:
            text = self.process.stdout.readline()
            self._stdout(text)
            if "POSITIONS: reading from stdin" in text:
                break

    def _input_positions(self, atoms):
        self._stdout("Inputting positions...\n")
        for atom in atoms.get_scaled_positions():
            self._stdin('{0} {1} {2}'.format(*atom))

    def close(self):
        if self.process is None:
            return

        self._stdout('Attemping to close VASP cleanly\n')
        stopcar = open(self.path + "STOPCAR", "w")
        stopcar.write("LABORT = .TRUE.")
        stopcar.close()
        self._input_positions(self.atoms)
        while self.process.poll() is None:
            text = self.process.stdout.readline()
            self._stdout(text)
            if "POSITIONS: reading from stdin" in text:
                break
        self._input_positions(self.atoms)
        while self.process.poll() is None:
            text = self.process.stdout.readline()
            self._stdout(text)
            if "deleting file STOPCAR" in text:
                break
        self._stdout("VASP has been closed\n")
        self.process = None

    def calculate(self, atoms=None, properties=['energy'],
                  system_changes=['positions', 'numbers', 'cell']):
        Calculator.calculate(self, atoms, properties, system_changes)

        if 'numbers' in system_changes:
            self.close()

        if self.process is None:
            self.first_step(atoms)
        else:
            if not system_changes:
                return
            self.continue_stepping(atoms)

        new = read('vasprun.xml', index=-1)

        self.results = {'free_energy': new.get_potential_energy(force_consistent=True),
                        'energy': new.get_potential_energy(),
                        'forces': new.get_forces(),
                        'stress': new.get_stress()}

    def __del__(self):
        self.close()

from __future__ import print_function
from subprocess import Popen, PIPE
from ase.calculators.calculator import Calculator
from ase.io import read, write
import os


class VaspInteractive(Calculator):
    name = "VaspInteractive"
    implemented_properties = ['energy', 'forces', 'stress']

    def __init__(self, txt="interactive.log", print_log=True, process=None,
                 command=None, path="./"):
        self.process = process
        self.path = path
        if txt is not None:
            self.txt = open(txt, "a")
        else:
            self.txt = None
        self.print_log = print_log

        if command is not None:
            self.command = command
        elif 'VASP_COMMAND' in os.environ:
            self.command = os.environ['VASP_COMMAND']
        elif 'VASP_SCRIPT' in os.environ:
            self.command = os.environ['VASP_SCRIPT']
        else:
            raise RuntimeError('Please set either command in calculator'
                               ' or VASP_COMMAND environment variable')

        if isinstance(self.command, str):
            self.command = self.command.split()

        self.atoms = None

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

    def _run_vasp(self, atoms):
        if self.process is None:
            if os.path.isfile('STOPCAR'):
                os.remove('STOPCAR')
            self._stdout("Writing Initial POSCAR\n")
            write(os.path.join(self.path, "POSCAR"), atoms)
            self._stdout("Starting VASP for initial step...\n")
            self.process = Popen(self.command, stdout=PIPE,
                                 stdin=PIPE, stderr=PIPE, cwd=self.path)
        else:
            self._stdout("Inputting positions...\n")
            for atom in atoms.get_scaled_positions():
                self._stdin(' '.join(map('{:19.16f}'.format, atom)))

        while self.process.poll() is None:
            text = self.process.stdout.readline()
            self._stdout(text)
            if "POSITIONS: reading from stdin" in text:
                break

    def close(self):
        if self.process is None:
            return

        self._stdout('Attemping to close VASP cleanly\n')
        with open(os.path.join(self.path, 'STOPCAR', 'w')) as stopcar:
            stopcar.write('LABORT = .TRUE.')

        self._run_vasp(self.atoms)
        self._run_vasp(self.atoms)
        while self.process.poll() is not None:
            time.sleep(1)
        self._stdout("VASP has been closed\n")
        self.process = None

    def calculate(self, atoms=None, properties=['energy'],
                  system_changes=['positions', 'numbers', 'cell']):
        Calculator.calculate(self, atoms, properties, system_changes)

        if not system_changes:
            return

        if 'numbers' in system_changes:
            self.close()

        self._run_vasp(atoms)

        new = read(os.path.join(self.path, 'vasprun.xml'), index=-1)

        self.results = {'free_energy': new.get_potential_energy(force_consistent=True),
                        'energy': new.get_potential_energy(),
                        'forces': new.get_forces(),
                        'stress': new.get_stress()}

    def __del__(self):
        self.close()

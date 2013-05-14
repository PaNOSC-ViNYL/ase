import sys
import os.path
import argparse
    
import numpy as np

from ase.io import read
from ase.db import connect
from ase.lattice import bulk
from ase.parallel import world
from ase.structure import molecule
from ase.utils import devnull, prnt
from ase.cli.plugin import PluginCommand
from ase.atoms import Atoms, string2symbols
from ase.data import ground_state_magnetic_moments
from ase.data import chemical_symbols, atomic_numbers, covalent_radii


def expand(names):
    """Expand ranges like H-Li to H, He and Li."""
    i = 0
    while i < len(names):
        name = names[i]
        if name.count('-') == 1:
            s1, s2 = name.split('-')
            Z1 = atomic_numbers.get(s1)
            Z2 = atomic_numbers.get(s2)
            if Z1 is not None and Z2 is not None:
                names[i:i + 1] = chemical_symbols[Z1:Z2 + 1]
                i += Z2 - Z1
        i += 1


class MyArgumentParser(argparse.ArgumentParser):
    def convert_arg_line_to_args(self, arg_line):
        for arg in arg_line.split():
            yield arg


class CLI:
    def __init__(self, args=None):
        if world.rank == 0:
            self.logfile = sys.stdout
        else:
            self.logfile = devnull

        self.command = None
        self.build_function = None
        self.args = args
        self.collection = None

    def log(self, *args, **kwargs):
        prnt(file=self.logfile, *args, **kwargs)

    def run(self):
        args = self.args

        command = self.get_command_object(args.command)

        if args.collection is not None:
            self.read_collection()
            if len(args.names) == 0:
                args.names = self.collection.keys()

        if args.plugin:
            f = open(args.plugin)
            script = f.read()
            f.close()
            namespace = {}
            exec script in namespace
            if 'names' in namespace and len(args.names) == 0:
                args.names = namespace['names']
            self.build_function = namespace.get('build')
            if 'calculate' in namespace:
                command = PluginCommand(namespace.get('calculate'))
                command.hook = self.hook

        command.logfile = self.logfile
        command.args = args

        expand(args.names)

        atoms = None
        for name in args.names:
            if atoms is not None:
                del atoms.calc  # release resources from last calculation
            atoms = self.build(name)
            command.run(atoms, name)

        command.finalize()

        return atoms

    def parse(self, args):
        # create the top-level parser
        parser = MyArgumentParser(fromfile_prefix_chars='@')
        add = parser.add_argument
        add('names', nargs='*')
        add('-q', '--quiet', action='store_const', const=0, default=1,
            dest='verbosity')
        add('-V', '--verbose', action='store_const', const=2, default=1,
            dest='verbosity')
        add('-t', '--tag',
            help='String tag added to filenames.')
        add('-M', '--magnetic-moment',
            metavar='M1,M2,...',
            help='Magnetic moment(s).  ' +
            'Use "-M 1" or "-M 2.3,-2.3".')
        add('--modify', metavar='...',
            help='Modify atoms with Python statement.  ' +
            'Example: --modify="atoms.positions[-1,2]+=0.1".')
        add('-v', '--vacuum', type=float, default=3.0,
            help='Amount of vacuum to add around isolated atoms '
            '(in Angstrom).')
        add('--unit-cell',
            help='Unit cell.  Examples: "10.0" or "9,10,11" ' +
            '(in Angstrom).')
        add('--bond-length', type=float,
            help='Bond length of dimer in Angstrom.')
        add('-c', '--collection')
        add('-x', '--crystal-structure',
            help='Crystal structure.',
            choices=['sc', 'fcc', 'bcc', 'hcp', 'diamond',
                     'zincblende', 'rocksalt', 'cesiumchloride',
                     'fluorite'])
        add('-a', '--lattice-constant', default='',
            help='Lattice constant(s) in Angstrom.')
        add('--orthorhombic', action='store_true',
            help='Use orthorhombic unit cell.')
        add('--cubic', action='store_true',
            help='Use cubic unit cell.')
        add('-r', '--repeat',
            help='Repeat unit cell.  Use "-r 2" or "-r 2,3,1".')
        add('--plugin')

        subparsers = parser.add_subparsers(dest='command',
                                           help='sub-command help')


        commands = ['run', 'optimize', 'eos', 'write',
                    'reaction', 'results', 'view', 'python']
        commands += self.hook.get('commands', [])
        for command in commands:
            cmd = self.get_command_object(command)
            cmd.add_parser(subparsers)

        self.args = parser.parse_args(args)
        if self.args.verbosity == 2:
            print(self.args)

    def get_command_object(self, name):
        classname = name.title().replace('-', '') + 'Command'
        if name in self.hook.get('commands', []):
            module = self.hook['command_module']
        else:
            module = 'ase.cli'
        name = name.replace('-', '_')
        module = __import__(module + '.' + name, {}, None, [classname])
        cmd = getattr(module, classname)()
        cmd.hook = self.hook
        return cmd

    def build(self, name):
        args = self.args
        if '.' in name:
            # Read from file:
            atoms = read(name)
        elif self.build_function:
            # Use build() function from plugin:
            atoms = self.build_function(name, self.args)
        elif self.args.crystal_structure:
            atoms = self.build_bulk(name)
        elif self.args.collection:
            atoms = self.collection[name]
        else:
            atoms = self.build_molecule(name)

        if args.magnetic_moment:
            magmoms = np.array(
                [float(m) for m in args.magnetic_moment.split(',')])
            atoms.set_initial_magnetic_moments(
                np.tile(magmoms, len(atoms) // len(magmoms)))

        if args.modify:
            exec args.modify in {'atoms': atoms}

        if args.repeat is not None:
            r = args.repeat.split(',')
            if len(r) == 1:
                r = 3 * r
            atoms = atoms.repeat([int(c) for c in r])

        return atoms

    def build_molecule(self, name):
        args = self.args
        try:
            # Known molecule or atom?
            atoms = molecule(name)
        except NotImplementedError:
            symbols = string2symbols(name)
            if len(symbols) == 1:
                Z = atomic_numbers[symbols[0]]
                magmom = ground_state_magnetic_moments[Z]
                atoms = Atoms(name, magmoms=[magmom])
            elif len(symbols) == 2:
                # Dimer
                if args.bond_length is None:
                    b = (covalent_radii[atomic_numbers[symbols[0]]] +
                         covalent_radii[atomic_numbers[symbols[1]]])
                else:
                    b = args.bond_length
                atoms = Atoms(name, positions=[(0, 0, 0),
                                               (b, 0, 0)])
            else:
                raise ValueError('Unknown molecule: ' + name)
        else:
            if len(atoms) == 2 and args.bond_length is not None:
                atoms.set_distance(0, 1, args.bond_length)

        if args.unit_cell is None:
            atoms.center(vacuum=args.vacuum)
        else:
            a = [float(x) for x in args.unit_cell.split(',')]
            if len(a) == 1:
                cell = [a[0], a[0], a[0]]
            elif len(a) == 3:
                cell = a
            else:
                a, b, c, alpha, beta, gamma = a
                degree = np.pi / 180.0
                cosa = np.cos(alpha * degree)
                cosb = np.cos(beta * degree)
                sinb = np.sin(beta * degree)
                cosg = np.cos(gamma * degree)
                sing = np.sin(gamma * degree)
                cell = [[a, 0, 0],
                        [b * cosg, b * sing, 0],
                        [c * cosb, c * (cosa - cosb * cosg) / sing,
                         c * np.sqrt(
                            sinb**2 - ((cosa - cosb * cosg) / sing)**2)]]
            atoms.cell = cell
            atoms.center()

        return atoms

    def build_bulk(self, name):
        args = self.args
        L = args.lattice_constant.replace(',', ' ').split()
        d = dict([(key, float(x)) for key, x in zip('ac', L)])
        atoms = bulk(name, crystalstructure=args.crystal_structure,
                     a=d.get('a'), c=d.get('c'),
                     orthorhombic=args.orthorhombic, cubic=args.cubic)

        M, X = {'Fe': (2.3, 'bcc'),
                'Co': (1.2, 'hcp'),
                'Ni': (0.6, 'fcc')}.get(name, (None, None))
        if M is not None and args.crystal_structure == X:
            atoms.set_initial_magnetic_moments([M] * len(atoms))

        return atoms

    def read_collection(self):
        colname = self.args.collection
        if os.path.isfile(colname):
            self.collection = connect(colname)
        else:
            module, name = colname.rsplit('.', 1)
            module = __import__(module, {}, None, [name])
            self.collection = getattr(module, name)
            if not isinstance(self.collection, dict):
                self.collection = self.collection()


def run(command=None, hook=None, **kwargs):
    if command is None:
        command = sys.argv[1:]
    elif ' ' in command:
        command = command.split()
    runner = CLI()
    runner.hook = hook or {}
    if isinstance(command, str):
        kwargs['command'] = command
        runner.args = FakeArguments(kwargs)
    else:
        runner.parse(command)
    atoms = runner.run()
    return atoms


class FakeArguments:
    def __init__(self, kwargs):
        self.args = kwargs

    def __getattr__(self, name):
        return self.args.get(name)

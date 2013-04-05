import os
import time
import argparse
import traceback

from ase.cli.command import Command
from ase.calculators.calculator import get_calculator, \
    names as calculator_names
import ase.db as db


def str2dict(s, namespace={}, sep='='):
    """Convert comma-separated key=value string to dictionary.

    Examples:

    >>> str2dict('xc=PBE,nbands=200,parallel={band:4}')
    {'xc': 'PBE', 'nbands': 200, 'parallel': {'band': 4}}
    >>> str2dict('a=1.2,b=True,c=ab,d=1,2,3,e={f:42,g:cd}')
    {'a': 1.2, 'c': 'ab', 'b': True, 'e': {'g': 'cd', 'f': 42}, 'd': (1, 2, 3)}
    """
    
    dct = {}
    s = (s + ',').split(sep)
    for i in range(len(s) - 1):
        key = s[i]
        m = s[i + 1].rfind(',')
        value = s[i + 1][:m]
        if value[0] == '{':
            assert value[-1] == '}'
            value = str2dict(value[1:-1], namespace, ':')
        else:
            try:
                value = eval(value, namespace)
            except (NameError, SyntaxError):
                pass
        dct[key] = value
        s[i + 1] = s[i + 1][m + 1:]
    return dct


class RunCommand(Command):
    db = None

    def add_parser(self, subparser):
        parser = subparser.add_parser('run', help='run ...')
        self.add_arguments(parser)

    def add_arguments(self, parser):
        add = parser.add_argument
        add('--after')
        calculator = self.default_calculator['name']
        if calculator == 'emt':
            help = ('Name of calculator to use: ' +
                    ', '.join(calculator_names) +
                    '.  Default is emt.')
        else:
            help=argparse.SUPPRESS
        add('-c', '--calculator', default=calculator, help=help)
        add('-p', '--parameters', default='',
            metavar='key=value,...',
            help='Comma-separated key=value pairs of ' +
            'calculator specific parameters.')
        add('-d', '--database',
            help='Use sqlite for a sqlite3 database and json for a simple ' +
            'json database.  Default is no database')
        add('-l', '--use-lock-file', action='store_true',
            help='Skip calculations where the json ' +
            'lock-file or result file already exists.')
        add('--properties', default='efsdMm',
            help='Default value is "efsdMm" meaning calculate energy, ' +
            'forces, stress, dipole moment, total magnetic moment and ' +
            'atomic magnetic moments.')
    
    def run(self, atoms, name):
        args = self.args

        if self.db is None:
            # Create database object:
            self.db = db.database(self.get_filename(ext=args.database),
                                  args.database,
                                  use_lock_file=args.use_lock_file)

        skip = False
        try:
            self.db.write(name, None)
        except db.KeyCollisionError:
            skip = True
        
        if not skip:
            self.set_calculator(atoms, name)

            tstart = time.time()
            try:
                data = self.calculate(atoms, name)
            except KeyboardInterrupt:
                raise
            except Exception:
                self.log(name, 'FAILED')
                traceback.print_exc(file=self.logfile)
            else:
                tstop = time.time()
                data['time'] = tstop - tstart
                self.db.write(name, atoms, data, overwrite=True)

    def set_calculator(self, atoms, name):
        args = self.args
        cls = get_calculator(args.calculator)
        namespace = self.default_calculator.get('namespace', {})
        parameters = str2dict(args.parameters, namespace)
        if getattr(cls, 'nolabel', False):
            atoms.calc = cls(**parameters)
        else:
            atoms.calc = cls(label=self.get_filename(name), **parameters)

    def calculate(self, atoms, name):
        args = self.args

        for p in args.properties:
            property, method = {'e': ('energy', 'get_potential_energy'),
                                'f': ('forces', 'get_forces'),
                                's': ('stress', 'get_stress'),
                                'd': ('dipole', 'get_dipole_moment'),
                                'M': ('magmom', 'get_magnetic_moment'),
                                'm': ('magmoms', 'get_magnetic_moments')}[p]
            try:
                x = getattr(atoms, method)()
            except NotImplementedError:
                pass

        data = {}
        if args.after:
            exec args.after in {'atoms': atoms, 'data': data}
        
        return data

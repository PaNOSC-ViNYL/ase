import os
import sys
import tempfile
import argparse

from ase.cli.run import RunCommand


class PythonCommand(RunCommand):
    def add_parser(self, subparser):
        parser = subparser.add_parser('python',
                                      help='Interactive Python session')
        self.add_arguments(parser)

    def add_arguments(self, parser):
        RunCommand.add_arguments(self, parser)
        parser.add_argument(
            '--interactive-python-session', action='store_true',
            help=argparse.SUPPRESS)
        
    def set_calculator(self, atoms, name):
        if self.args.interactive_python_session:
            RunCommand.set_calculator(self, atoms, name)

    def calculate(self, atoms, name):
        if self.args.interactive_python_session:
            return {}

        file = tempfile.NamedTemporaryFile()
        file.write('import os\n')
        file.write('if "PYTHONSTARTUP" in os.environ:\n')
        file.write('    execfile(os.environ["PYTHONSTARTUP"])\n')
        file.write('from ase.cli import run\n')
        file.write('atoms = run(%r)\n' % 
                   (' '.join(sys.argv[1:]) + ' --interactive-python-session'))
        file.flush()
        os.system('python -i %s' % file.name)
        return {}

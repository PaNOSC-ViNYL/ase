import numpy as np

from ase.cli.command import Command


class ResultsCommand(Command):
    def add_parser(self, subparser):
        parser = subparser.add_parser('results', help='results ...')

    def finalize(self):
        for name in self.args.names:
            if name in self.data:
                e = self.data[name].get('energy', 42)
            else:
                e = 117
            print '%2s %10.3f' % (name, e)

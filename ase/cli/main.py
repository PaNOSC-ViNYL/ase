import argparse

from ase import __version__
from ase.utils import import_module


commands = ['info']


def add_arguments(parser):
    parser.add_argument('-P', type=int, default=1)
    parser.add_argument('--version', action='version',
                        version='%(prog)s-{}'.format(__version__))
    #-vq--version

def main():
    parser = argparse.ArgumentParser(
        prog='ase',
        description='ASE command line tool')
    add_arguments(parser)
    subparsers = parser.add_subparsers(title='Subcommands',
                                       dest='command')#, help='H')

    for command in commands:
        module = import_module('ase.cli.' + command)
        subparser = subparsers.add_parser(command,
                                          description=module.description,
                                          help='hhhh')
        module.add_arguments(subparser)
        subparser.set_defaults(func=module.main)

    args = parser.parse_args()
    print(args)
    if 'func' in args:
        args.func(args)

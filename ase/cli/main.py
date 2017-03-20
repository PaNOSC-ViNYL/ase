from __future__ import print_function

import argparse
import sys

from ase import __version__
from ase.utils import import_module


commands = [
    ('info', 'ase.cli.info'),
    ('test', 'ase.test'),
    ('gui', 'ase.gui.ag'),
    ('run', 'ase.cli.run'),
    ('build', 'ase.cli.build'),
    ('db', 'ase.db.cli'),
    ('nomad-upload', 'ase.cli.nomad'),
    ('install-completion-script', 'ase.cli.complete')]


def add_arguments(parser):
    parser.add_argument('-v', '--verbose', action='store_true')
    parser.add_argument('-q', '--quiet', action='store_true')


def main():
    parser = argparse.ArgumentParser(
        prog='ase',
        description='ASE command line tool')
    #add_arguments(parser)
    parser.add_argument('--version', action='version',
                        version='%(prog)s-{}'.format(__version__))
    subparsers = parser.add_subparsers(title='Subcommands',
                                       dest='command')

    subparser = subparsers.add_parser('help',
                                      description='Help',
                                      help='hhhh')
    subparser.description = 'Help2'
    subparser.add_argument('helpcommand', nargs='?')

    functions = {}
    parsers = {}
    for command, module_name in commands:
        module = import_module(module_name)
        subparser = subparsers.add_parser(command,
                                          description=module.description,
                                          help='hhhh')
        add_arguments(subparser)
        module.add_arguments(subparser)
        functions[command] = module.main
        parsers[command] = subparser

    args = parser.parse_args()

    if args.command == 'help':
        if args.helpcommand is None:
            parser.print_help()
        else:
            parsers[args.helpcommand].print_help()
    elif args.command is None:
        parser.print_usage()
    else:
        try:
            functions[args.command](args)
        except KeyboardInterrupt:
            pass
        except Exception as x:
            if args.verbose:
                raise
            else:
                print('{}: {}'.format(x.__class__.__name__, x),
                      file=sys.stderr)
                print('To get a full traceback, use: ase --verbose',
                      file=sys.stderr)


def old():
    sys.argv[:1] = sys.argv[0].rsplit('-', 1)
    main()

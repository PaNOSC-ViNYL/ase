import argparse

from ase.utils import import_module


commands = ['info']


def add_arguments(parser):
    parser.add_argument('-P', type=int, default=1)


def main():
    parser = argparse.ArgumentParser()
    add_arguments(parser)
    subparsers = parser.add_subparsers()

    for command in commands:
        subparser = subparsers.add_parser(command)
        module = import_module('ase.cli.' + command)
        module.add_arguments(subparser)
        subparser.set_defaults(func=module.main)

    args = parser.parse_args()

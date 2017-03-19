#!/usr/bin/env python3
"""Bash completion for ase-db, ase-run, ase-build, ase-info and ase-gui.

Put this in your .bashrc::

    complete -o default -C /path/to/ase/cli/complete.py ase
"""

import os
import sys
from glob import glob


def match(word, *suffixes):
    return [w for w in glob(word + '*')
            if any(w.endswith(suffix) for suffix in suffixes)]


def generate():
    ...


commands = {'test': ['-c', '--calculators']}


def complete(word, previous, line, point):
    for w in line[:point - len(word)].strip().split()[1:]:
        if w[0].isalpha():
            if w in commands:
                command = w
                break
            return commands
    else:
        if word[:1] == '-':
            return ['-h', '--help', '-q', '--quiet', '-v', '--verbose',
                    '--version']
    if word[:1] == '-':
        return commands[comand]

    words = []

    if command == 'db':
        if previous == 'db':
            words = match(word, '.db', '.json')
    elif command == 'run':
        if previous == 'run':
            from ase.calculators.calculator import names as words
    elif command == 'build':
        if previous in ['-x', '--crystal-structure']:
            words = ['sc', 'fcc', 'bcc', 'hcp', 'diamond', 'zincblende',
                     'rocksalt', 'cesiumchloride', 'fluorite', 'wurtzite']
    return words


def add_arguments(parser):
    description = 'Add bash-completion script to ~/.bashrc.'
    return description


def main(args):
    print(__path__)
    cmd = 'complete -o default -C /path/to/ase/cli/complete.py ase'
    with open(os.path.expanduser('~/.bashrc')) as fd:
        if cmd in fd.readlines():
            print('Completion script already installed!')
            return
    with open(os.path.expanduser('~/.bashrc'), 'a') as fd:
        print(cmd, file=fd)


if __name__ == '__main__':
    word, previous = sys.argv[2:]
    line = os.environ['COMP_LINE']
    point = int(os.environ['COMP_POINT'])
    words = complete(word, previous, line, point)
    for w in words:
        if w.startswith(word):
            print(w)

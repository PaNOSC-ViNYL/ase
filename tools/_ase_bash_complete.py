#!/usr/bin/env python
"""Bash completion for ase-db, ase-run, ase-build, ase-info and ase-gui.

Put this in your .bashrc::
    
    complete -o default -C _ase_bash_complete.py ase-db ase-run \
    ase-build ase-info ase-gui
"""

import os
import sys
from glob import glob

command, word, previous = sys.argv[1:]
line = os.environ['COMP_LINE']
point = int(os.environ['COMP_POINT'])


def options(short, long):
    return ['-' + s for s in short] + ['--' + l for l in long.split()]


def match(word, *suffixes):
    return [w for w in glob(word + '*')
            if any(w.endswith(suffix) for suffix in suffixes)]
    
    
words = []

if command == 'ase-db':
    if word[:1] == '-':
        words = options(
            'hvqnliakycspwLj',
            'help verbose quiet count long insert-into add-from-file '
            'add-key-value-pairs limit offset delete '
            'delete-keys yes columns sort cut python csv '
            'open-web-browser json unique analyse')
    elif previous == 'ase-db':
        words = match(word, '.db', '.json')
elif command == 'ase-run':
    if word[:1] == '-':
        words = options(
            'htpdSfsEic',
            'help tag parameter database skip properties maximum-force '
            'constrain-tags maximum-stress equation-of-state modify after')
    elif previous == 'ase-run':
        from ase.calculators.calculator import names as words
elif command == 'ase-build':
    if previous in ['-x', '--crystal-structure']:
        words = ['sc', 'fcc', 'bcc', 'hcp', 'diamond', 'zincblende',
                 'rocksalt', 'cesiumchloride', 'fluorite', 'wurtzite']
    elif word[:1] == '-':
        words = options(
            'hMvxarg',
            'help magnetic-moment modify vacuum unit-cell bond-length '
            'crystal-structure lattice-constant orthorhombic cubic repeat gui')
elif command == 'ase-info':
    if word[:1] == '-':
        words = options('h', 'help')
    else:
        words = match(word, '.traj')
else:  # ase-gui
    if word[:1] == '-':
        words = options(
            'hnurRogtbs', 'help image-number show-unit-cell repeat verbose '
            'rotations output graph terminal aneb interpolate bonds scale')
        
for w in words:
    if w.startswith(word):
        print(w)

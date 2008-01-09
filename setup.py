#!/usr/bin/env python

# Copyright (C) 2007  CAMP
# Please see the accompanying LICENSE file for further information.

from distutils.core import setup
from glob import glob
from os.path import join

import os
import sys

long_description = """\
ASE is a python package providing an open source Atomic Simulation
Environment in the python scripting language."""


if sys.version_info < (2, 3, 0, 'final', 0):
    raise SystemExit, 'Python 2.3 or later is required!'

packages = ['ase',
            'ase.io',
            'ase.md',
            'ase.dft',
            'ase.gui',
            'ase.test',
            'ase.examples',
            'ase.optimize',
            'ase.calculators',
            'ase.gui.languages']

# set the version number 
version = '3.0.0'

setup(name = 'python-ase',
      version=version,
      description='Atomic Simulation Environment',
      url='http://www.fysik.dtu.dk/Campos/ase',
      maintainer='CAMd',
      maintainer_email='camd@fysik.dtu.dk',
      license='GPL',
      platforms=['linux'],
      packages=packages,
      scripts=['tools/ag'],
      long_description=long_description)



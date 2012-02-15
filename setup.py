#!/usr/bin/env python

# Copyright (C) 2007  CAMP
# Please see the accompanying LICENSE file for further information.

from distutils.core import setup
from glob import glob
from os.path import join

import os
import sys

from os import path
try:
    from subprocess import Popen, PIPE
except ImportError:
    from os import popen3
else:
    def popen3(cmd):
        p = Popen(cmd, shell=True, close_fds=True,
                  stdin=PIPE, stdout=PIPE, stderr=PIPE)
        return p.stdin, p.stdout, p.stderr

long_description = """\
ASE is a python package providing an open source Atomic Simulation
Environment in the python scripting language."""


if sys.version_info < (2, 3, 0, 'final', 0):
    raise SystemExit, 'Python 2.3 or later is required!'

packages = ['ase',
            'ase.cluster',
            'ase.cluster.data',
            'ase.io',
            'ase.md',
            'ase.dft',
            'ase.gui',
            'ase.gui.languages',
            'ase.data',
            'ase.test',
            'ase.tasks',
            'ase.utils',
            'ase.lattice',
            'ase.lattice.spacegroup',
            'ase.examples',
            'ase.optimize',
            'ase.optimize.test',
            'ase.visualize',
            'ase.visualize.vtk',
            'ase.transport',
            'ase.calculators',
            'ase.calculators.jacapo']

package_dir={'ase': 'ase'}

package_data={'ase': ['lattice/spacegroup/spacegroup.dat',
                      'gui/po/ag.pot',
                      'gui/po/makefile',
                      'gui/po/??_??/LC_MESSAGES/ag.po']}

# Compile makes sense only when building
if 'build' in sys.argv or 'build_ext' in sys.argv or 'install' in sys.argv:
    msgfmt = 'msgfmt'
    # Compile translation files (requires gettext)
    cmd = popen3(msgfmt + ' -V')[1]
    output = cmd.read().strip()
    cmd.close()
    if output:
        for pofile in glob('ase/gui/po/??_??/LC_MESSAGES/ag.po'):
            mofile = os.path.join(os.path.split(pofile)[0], 'ag.mo')
            os.system(msgfmt + ' -cv %s --output-file=%s' % (pofile, mofile))
        package_data.update({'ase': package_data['ase'] + \
                             ['gui/po/??_??/LC_MESSAGES/ag.mo']})

# Get the current version number:
execfile('ase/svnversion_io.py')  # write ase/svnversion.py and get svnversion
execfile('ase/version.py')        # get version_base
if svnversion:
    version = version_base + '.' + svnversion
else:
    version = version_base

setup(name = 'python-ase',
      version=version,
      description='Atomic Simulation Environment',
      url='https://wiki.fysik.dtu.dk/ase',
      maintainer='CAMd',
      maintainer_email='camd@fysik.dtu.dk',
      license='LGPLv2.1+',
      platforms=['linux'],
      packages=packages,
      package_dir=package_dir,
      package_data=package_data,
      scripts=['tools/ag', 'tools/ase', 'tools/ASE2ase.py', 'tools/testase.py'],
      long_description=long_description)

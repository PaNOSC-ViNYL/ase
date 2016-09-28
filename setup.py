#!/usr/bin/env python

# Copyright (C) 2007  CAMP
# Please see the accompanying LICENSE file for further information.

from __future__ import print_function
import os
import re
import sys
from distutils.core import setup
from distutils.command.build_py import build_py as _build_py
from distutils.command.sdist import sdist as _sdist
from glob import glob
from os.path import join


long_description = """\
ASE is a python package providing an open source Atomic Simulation
Environment in the Python language."""


if sys.version_info < (2, 6, 0, 'final', 0):
    raise SystemExit('Python 2.6 or later is required!')

    
class sdist(_sdist):
    """Fix distutils.
    
    Distutils insists that there should be a README or README.txt,
    but GitLab.com needs README.rst in order to parse it as reStructureText."""
    
    def warn(self, msg):
        if msg.startswith('standard file not found: should have one of'):
            self.filelist.append('README.rst')
        else:
            _sdist.warn(self, msg)
            
    
packages = []
for dirname, dirnames, filenames in os.walk('ase'):
    if '__init__.py' in filenames:
        packages.append(dirname.replace('/', '.'))

package_dir = {'ase': 'ase'}

package_data = {'ase': ['spacegroup/spacegroup.dat',
                        'collections/*.json',
                        'db/templates/*',
                        'db/static/*']}


class build_py(_build_py):
    """Custom distutils command to build translations."""
    def __init__(self, *args, **kwargs):
        _build_py.__init__(self, *args, **kwargs)
        # Keep list of files to appease bdist_rpm.  We have to keep track of
        # all the installed files for no particular reason.
        self.mofiles = []

    def run(self):
        """Compile translation files (requires gettext)."""
        _build_py.run(self)
        msgfmt = 'msgfmt'
        status = os.system(msgfmt + ' -V')
        if status == 0:
            for pofile in glob('ase/gui/po/*/LC_MESSAGES/ag.po'):
                dirname = join(self.build_lib, os.path.dirname(pofile))
                if not os.path.isdir(dirname):
                    os.makedirs(dirname)
                mofile = join(dirname, 'ag.mo')
                status = os.system('%s -cv %s --output-file=%s 2>&1' %
                                   (msgfmt, pofile, mofile))
                assert status == 0, 'msgfmt failed!'
                self.mofiles.append(mofile)

    def get_outputs(self, *args, **kwargs):
        return _build_py.get_outputs(self, *args, **kwargs) + self.mofiles

# Get the current version number:
with open('ase/__init__.py') as fd:
    version = re.search("__version__ = '(.*)'", fd.read()).group(1)
    
name = 'ase'  # PyPI name

# Linux-distributions may want to change the name:
if 0:
    name = 'python-ase'
    
scripts = ['tools/ase-gui', 'tools/ase-db', 'tools/ase-info',
           'tools/ase-build', 'tools/ase-run']
# provide bat executables in the tarball and always for Win
if 'sdist' in sys.argv or os.name in ['ce', 'nt']:
    for s in scripts[:]:
        scripts.append(s + '.bat')

setup(name=name,
      version=version,
      description='Atomic Simulation Environment',
      url='https://wiki.fysik.dtu.dk/ase',
      maintainer='ASE-community',
      maintainer_email='ase-developers@listserv.fysik.dtu.dk',
      license='LGPLv2.1+',
      platforms=['unix'],
      packages=packages,
      package_dir=package_dir,
      package_data=package_data,
      scripts=scripts,
      long_description=long_description,
      cmdclass={'build_py': build_py,
                'sdist': sdist},
      classifiers=[
          'Development Status :: 6 - Mature',
          'License :: OSI Approved :: '
          'GNU Lesser General Public License v2 or later (LGPLv2+)',
          'Operating System :: OS Independent',
          'Programming Language :: Python :: 2',
          'Programming Language :: Python :: 2.6',
          'Programming Language :: Python :: 2.7',
          'Programming Language :: Python :: 3',
          'Programming Language :: Python :: 3.4',
          'Programming Language :: Python :: 3.5',
          'Topic :: Scientific/Engineering :: Physics'])

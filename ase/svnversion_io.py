from __future__ import print_function
# Copyright (C) 2003  CAMP
# Please see the accompanying LICENSE file for further information.

import sys
import subprocess
from os import path

ON_POSIX = 'posix' in sys.builtin_module_names


def write_svnversion(svnversion, dir):
    svnversionfile = path.join(dir, 'svnversion.py')
    f = open(svnversionfile,'w')
    f.write('svnversion = "%s"\n' % svnversion)
    f.close()
    print('svnversion = ' +svnversion+' written to '+svnversionfile)
    # assert svn:ignore property if the installation is under svn control
    # because svnversion.py has to be ignored by svn!
    subprocess.call('svn propset svn:ignore svnversion.py ' + dir, shell=True)


def get_svnversion_from_svn(dir):
    # try to get the last svn version number from svnversion
    # assert we are in the project dir
    output = subprocess.Popen('svnversion -n ' + dir, shell=True,
                              stdout=subprocess.PIPE).stdout.read().decode()
    if not (output + ' ')[0].isdigit():
        output = None
    return output

svnversion = get_svnversion_from_svn(dir='ase')
if svnversion:
    write_svnversion(svnversion, dir='ase')

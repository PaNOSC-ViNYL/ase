# Copyright (C) 2003  CAMP
# Please see the accompanying LICENSE file for further information.

version = '3.3.1'

def get_ase_svnversion_from_import():
    try:
        # try to import the last svn version number from ase/svnversion.py
        from ase.svnversion import svnversion
    except:
        svnversion = None
    ##
    return svnversion

svnversion = get_ase_svnversion_from_import()
if svnversion:
    version = version+'.'+svnversion

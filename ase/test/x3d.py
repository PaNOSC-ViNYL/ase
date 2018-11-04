from __future__ import print_function

import warnings
from unittest.case import SkipTest
import sys

from ase import Atoms

if sys.platform == 'darwin':
    # The test crashes when part of the test suite, because fork()
    # is then called in a multithreaded environment.  That is a big
    # no-no on all platforms, but cause a quick crash on macOS.
    raise SkipTest('Test not compatible with macOS.')

try:
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', DeprecationWarning)
        from IPython.display import HTML
    from ase.visualize import x3d
except ImportError:
    raise SkipTest('cannot import HTML from IPython.displacy')
else:
    print('Testing x3d...')
    a = 3.6
    b = a / 2
    atoms = Atoms('Cu4',
                    positions=[(0, 0, 0),
                               (0, b, b),
                               (b, 0, b),
                               (b, b, 0)],
                    cell=(a, a, a),
                    pbc=True)
    my_obj = x3d.view_x3d(atoms)
    assert isinstance(my_obj, HTML)

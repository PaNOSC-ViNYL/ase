from __future__ import print_function
from ase import Atoms
from ase.calculators.emt import EMT
from ase.optimize import QuasiNewton
from ase.vibrations import Vibrations
from ase.thermochemistry import IdealGasThermo

import numpy as np

from ase.vibrations.franck_condon import FCOverlap
from math import factorial

def equal(x, y, tolerance=0, fail=True, msg=''):
    """Compare x and y."""

    if not np.isfinite(x - y).any() or (np.abs(x - y) > tolerance).any():
        msg = (msg + '%s != %s (error: |%s| > %.9g)' %
               (x, y, x - y, tolerance))
        if fail:
            raise AssertionError(msg)
        else:
            sys.stderr.write('WARNING: %s\n' % msg)


# -------------------
# FCOverlap

fco = FCOverlap()

# check factorial
assert(fco.factorial(8) == factorial(8))
# the second test is useful according to the implementation
assert(fco.factorial(5) == factorial(5))

# check T=0 and n=0 equality
S = np.array([1, 2.1, 34])
m = 5
assert(((fco.directT0(m, S) - fco.direct(0, m, S)) / fco.directT0(m, S) 
        < 1e-15).all())
#assert((fco.directT0(m, S) == fco.direct(0, m, S)).all()) 

# check symmetry
S = 2
n = 3
assert(fco.direct(n, m, S) == fco.direct(m, n, S))

# ---------------------------
#

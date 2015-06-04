# Copyright 2008, 2009 CAMd
# (see accompanying license files for details).

"""Atomic Simulation Environment."""

__all__ = ['Atoms', 'Atom']

from ase.atom import Atom
from ase.atoms import Atoms

import numpy as np
if np.version.version < '1.9':
    # Make isinstance(x, numbers.Integral) work also for np.intxx:
    import numbers
    numbers.Integral.register(np.integer)
    del numbers
del np

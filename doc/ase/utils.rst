.. module:: ase.utils

==============================
Utillity functions and classes
==============================

This module contains utility functions and classes.

.. toctree::

    xrdebye

.. autofunction:: ase.utils.opencew
.. autofunction:: ase.utils.gcd
.. autofunction:: ase.utils.seterr
.. autofunction:: ase.utils.plural
.. autofunction:: ase.utils.formula_hill
.. autofunction:: ase.utils.formula_metal
.. autofunction:: ase.utils.convert_string_to_fd
.. autofunction:: ase.utils.natural_cutoffs
.. autoclass:: ase.utils.timing.Timer
.. autoclass:: ase.utils.timing.timer


Symmetry equivalence checker
============================

This module compares two atomic structures to see if they are symmetrically equivalent. It is based on the recipe used in `XtalComp`__

__ https://doi.org/10.1016/j.cpc.2011.11.007

.. autoclass:: ase.utils.structure_comparator.SymmetryEquivalenceCheck
   :members:
               
Symmetry analysis
=================

https://atztogo.github.io/spglib/python-spglib.html


Phonons
=======

http://phonopy.sourceforge.net/

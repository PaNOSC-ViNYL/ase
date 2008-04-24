

-----------------------
The ``ase.data`` module
-----------------------

This module defines the following variables: ``atomic_masses``, ``atomic_names``, ``chemical_symbols``, ``covalent_radii``, ``cpk_colors`` and  ``reference_states``.  All of these are lists that should be indexed with an atomic number:

>>> from ase.data import *
>>> atomic_names[92]
'Uranium'
>>> atomic_masses[2]
4.0026000000000002

If you don't know the atomic number of some element, then you can look it up in the ``atomic_numbers`` dictionary:

>>> atomic_numbers['Cu']
29
>>> covalent_radii[29]
1.1699999999999999

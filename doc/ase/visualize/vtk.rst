.. module:: vtk
   :synopsis: Use of the Visualization Toolkit.

.. index:: vtk

===========
ASE-VTK
===========

For ASE, the :mod:`~visualize.vtk` interface consists of Python modules for automatic
visualization of positions, bonds, forces and volume data (e.g. wave functions)
from an :class:`~atoms.Atoms` object, provided such data is made available
by the calculator.

XXX more

Representing atoms
==================
.. autoclass:: ase.visualize.vtk.atoms.vtkAtoms
   :inherited-members:

Atom-centered data
==================
.. autoclass:: ase.visualize.vtk.module.vtkModule
   :inherited-members:

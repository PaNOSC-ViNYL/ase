.. module:: vtk
   :synopsis: Use of the Visualization Toolkit.

.. index:: vtk

===========
ASE-VTK
===========

For ASE, the :mod:`~visualize.vtk` interface consists of Python modules for 
automatic visualization of positions, bonds, forces and volume data (e.g. wave
functions) from an :class:`~atoms.Atoms` object, provided such data is made
available by the calculator.

Representing atoms
==================
.. autoclass:: ase.visualize.vtk.atoms.vtkAtoms
   :members:
   :show-inheritance:

Atom-centered data
------------------

The superclass :class:`~ase.visualize.vtk.grid.vtkAtomicPositions` implements
the basic concepts for representing atomic-centered data in VTK.

.. autoclass:: ase.visualize.vtk.grid.vtkAtomicPositions
   :members:
..   :show-inheritance:

Predefined shapes
------------------

The class :class:`~ase.visualize.vtk.module.vtkGlyphModule` implements
the lower-level objects for representing predefined shapes (glyphs) in VTK.

.. autoclass:: ase.visualize.vtk.module.vtkGlyphModule
   :members:
   :inherited-members:
..   :show-inheritance:


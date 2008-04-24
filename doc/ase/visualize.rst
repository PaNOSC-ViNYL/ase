Visualization
=============

.. automodule:: ase.visualize
.. autofunction:: view

This provides an interface to various visualization tools, such as RasMol, VMD, VTK...
If no viewer is specified, a simple built-in viewer is opened. You can plot your ``atoms``, which 
is an :class:`Atoms` object, simply by doing

.. highlight:: python

::

  >>> view(atoms)

RasMol
------

You can use the `RasMol`_ program to visualize your system by specifying

  >>> view(atoms, viewer='rasmol')

XXX

.. _RasMol: http://openrasmol.org/

VMD
---

You can use RasMol to visualize your system by specifying

  >>> view(atoms, viewer='vmd')

XXX


VTK
---

The Visualization ToolKit (VTK) is an open source, freely available software system for 3D 
visualization. For more information on VTK, go to the `vtk homepage`_.

XXX

.. _vtk homepage: http://public.kitware.com/VTK/


gOpenMol
--------

You can use RasMol to visualize your system by specifying

  >>> view(atoms, viewer='gopenmol')

XXX


PovrayPlotter
-------------

XXX



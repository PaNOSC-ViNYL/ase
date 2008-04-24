Visualization
=============

.. automodule:: ase.visualize
.. autofunction:: view

This provides an interface to various visualization tools, such as RasMol, VMD, VTK...
If no viewer is specified, a simple built-in viewer is opened. You can plot your ``atoms``, which 
are an instance of the :class:`Atoms` class, simply by doing

.. highlight:: python

::

  >>> view(atoms)

RasMol
------

  >>> view(atoms, viewer='rasmol')

XXX


VMD
---

XXX


VTK
---

The Visualization ToolKit (VTK) is an open source, freely available software system for 3D 
visualization. For more information on VTK, go to the `vtk homepage`_.

XXX

.. _vtk homepage: http://public.kitware.com/VTK/

PovrayPlotter
-------------

XXX



.. module:: visualize

Visualization
=============

.. function:: view(atoms, data=None, viewer=None, repeat=None)

This provides an interface to various visualization tools, such as
:mod:`ase.gui <gui>`, RasMol_, VMD_, VTK_, gOpenMol_, or
Avogadro_. The default viewer is the ase.gui, described in the
:mod:`gui` module. The simplest invocation is::

  >>> from ase import view
  >>> view(atoms)

where ``atoms`` is any :class:`Atoms` object.  Alternative viewers can
be used by specifying the optional keyword ``viewer=...`` - use one of
'ase.gui', 'gopenmol', 'vmd', or 'rasmol'.  The VMD and Avogadro
viewers can take an optional ``data`` argument to show 3D data, such
as charge density::

  >>> view(atoms, viewer='VMD', data=array)

If you do not wish to open an interactive gui, but rather visualize
your structure by dumping directly to a graphics file; you can use the
``write`` command of the :mod:`io` module, which can write 'eps',
'png', and 'pov' files directly, like this::

  >>> write('image.png', atoms)

.. _RasMol: http://openrasmol.org/
.. _VMD: http://www.ks.uiuc.edu/Research/vmd/
.. _VTK: http://public.kitware.com/VTK/
.. _gOpenMol: http://www.csc.fi/gopenmol/
.. _Avogadro: http://avogadro.openmolecules.net/


PrimiPlotter
------------

The PrimiPlotter is intended to do on-the-fly plotting of the
positions of the atoms during long molecular dynamics simulations.
The module :mod:ase.visualize.primiplotter contains the PrimiPlotter
and the various output modules, see below.


.. autoclass:: ase.visualize.primiplotter.PrimiPlotter
   :inherited-members:



FieldPlotter
------------

The FieldPlotter is intended to plot fields defined on the atoms in
large-scale simulations.  The fields could be e.g. pressure, stress or
temperature (kinetic energy), i.e. any quantity that in a given
simulation is best defined on a per-atom basis, but is best
interpreted as a continuum field.

The current version of FieldPlotter only works if the number of atoms
is at least 5-10 times larger than the number of pixels in the plot.

.. autoclass:: ase.visualize.fieldplotter.FieldPlotter
   :inherited-members:





.. module:: visualize

Visualization
=============

.. function:: view(atoms, data=None, viewer=None, repeat=None)

This provides an interface to various visualization tools, such as
:mod:`ase.gui <gui>`, RasMol_, VMD_, VTK_, or gOpenMol_. The default viewer is
the ase.gui, described in the :mod:`gui` module. The simplest
invocation is::

  >>> from ase import view
  >>> view(atoms)

where ``atoms`` is any :class:`Atoms` object.  Alternative viewers can
be used by specifying the optional keyword ``viewer=...`` - use one of
'ase.gui', 'gopenmol', 'vmd', or 'rasmol'.  The VMD viewer can take an
optional ``data`` argument to show 3D data::

  >>> view(atoms, viewer='VMD', data=array)

If you do not wish to open an interactive gui, but rather visualize
your structure by dumping directly to a grphics file; you can use the
``write`` command of the :mod:`io` module, which can write 'eps',
'png', and 'pov' files directly, like this::

  >>> write('image.png', atoms)

.. _RasMol: http://openrasmol.org/
.. _VMD: http://www.ks.uiuc.edu/Research/vmd/
.. _VTK: http://public.kitware.com/VTK/
.. _gOpenMol: http://www.csc.fi/gopenmol/


PovrayPlotter
-------------

XXX


PrimiPlotter
-------------

XXX



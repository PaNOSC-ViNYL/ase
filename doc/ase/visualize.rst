.. module:: visualize

Visualization
=============

.. automodule:: ase.visualize
.. autofunction:: view

This provides an interface to various visualization tools, such as
`ase.gui`_, RasMol_, VMD_, VTK_, or gOpenMol_. The default viewer is
the ase.gui, described in the :mod:`gui` module. The simplest
invocation is::

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

.. _ase.gui: gui.html
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



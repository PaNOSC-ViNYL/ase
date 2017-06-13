.. module:: ase.visualize

Visualization
=============

.. function:: view(atoms, data=None, viewer=None, repeat=None)

This provides an interface to various visualization tools, such as
:mod:`ase.gui`, RasMol_, VMD_, gOpenMol_, Avogadro_, or ParaView_. The default viewer is
the ase.gui, described in the :mod:`ase.gui` module. The simplest invocation
is:

.. testsetup::

    from ase import Atoms
    atoms = Atoms('Cu')

>>> from ase.visualize import view
>>> view(atoms)

where ``atoms`` is any :class:`~ase.Atoms` object.  Alternative viewers
can be used by specifying the optional keyword ``viewer=...`` - use one of
'ase.gui', 'gopenmol', 'vmd', 'rasmol', or 'paraview'.  The VMD and Avogadro viewers can
take an optional ``data`` argument to show 3D data, such as charge density:

>>> view(atoms, viewer='VMD', data=...)

If you do not wish to open an interactive gui, but rather visualize
your structure by dumping directly to a graphics file; you can use the
``write`` command of the :mod:`ase.io` module, which can write 'eps',
'png', and 'pov' files directly, like this:

>>> from ase.io import write
>>> write('image.png', atoms)

It is also possible to plot directly to a Matplotlib subplot object, which allows for
a large degree of customisability. :ref:`More information <matplotlib_plotting>`.

.. _RasMol: http://openrasmol.org/
.. _VMD: http://www.ks.uiuc.edu/Research/vmd/
.. _gOpenMol: http://www.csc.fi/gopenmol/
.. _Avogadro: http://avogadro.openmolecules.net/
.. _ParaView: http://www.paraview.org/


.. module:: ase.visualize.mlab
.. _iso surface:

Plotting iso-surfaces with Mayavi
---------------------------------

The :func:`ase.visualize.mlab.plot` function can be used from the
command-line::

    $ python -m ase.visualize.mlab abc.cube

to plot data from a cube-file or alternatively a wave function or an electron
density from a calculator restart file::

    $ python -m ase.visualize.mlab -C gpaw abc.gpw

Options:

.. include:: mlab_options.txt
    :start-after: Options:

.. autofunction:: ase.visualize.mlab.plot


PrimiPlotter
------------

The PrimiPlotter is intended to do on-the-fly plotting of the
positions of the atoms during long molecular dynamics simulations.
The module :mod:`ase.visualize.primiplotter` contains the
PrimiPlotter and the various output modules, see below.

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


.. _matplotlib_plotting:

Matplotlib
----------

>>> import matplotlib.pyplot as plt
>>> from ase.visualize.plot import plot_atoms
>>> from ase.lattice.cubic import FaceCenteredCubic
>>> slab = FaceCenteredCubic('Au', size=(2, 2, 2))
>>> fig, ax = plt.subplots()
>>> plot_atoms(slab, ax, radii=0.3, rotation=('90x,45y,0z'))
>>> fig.savefig("ase_slab.png")

.. image:: matplotlib_plot_atoms1.png

The data is plotted directly to a Matplotlib subplot object, giving a large
degree of customizability.


>>> import matplotlib.pyplot as plt
>>> from ase.visualize.plot import plot_atoms
>>> from ase.lattice.cubic import FaceCenteredCubic
>>> slab = FaceCenteredCubic('Au', size=(2, 2, 2))
>>> fig, axarr = plt.subplots(1, 4, figsize=(15, 5))
>>> plot_atoms(slab, axarr[0], radii=0.3, rotation=('0x,0y,0z'))
>>> plot_atoms(slab, axarr[1], scale=0.7, offset=(3, 4), radii=0.3, rotation=('0x,0y,0z'))
>>> plot_atoms(slab, axarr[2], radii=0.3, rotation=('45x,45y,0z'))
>>> plot_atoms(slab, axarr[3], radii=0.3, rotation=('0x,0y,0z'))
>>> axarr[0].set_title("No rotation")
>>> axarr[1].set_xlabel("X-axis, [$\mathrm{\AA}$]")
>>> axarr[1].set_ylabel("Y-axis, [$\mathrm{\AA}$]")
>>> axarr[2].set_axis_off()
>>> axarr[3].set_xlim(2, 6)
>>> axarr[3].set_ylim(2, 6)
>>> fig.savefig("ase_slab_multiple.png")

.. image:: matplotlib_plot_atoms2.png

>>> stem_image = mpimg.imread("stem_image.jpg")
>>> atom_pos = [(0.0, 0.0, 0.0), (0.5, 0.5, 0.5), (0.5, 0.5, 0.0)]
>>> srtio3 = crystal(['Sr','Ti','O'], atom_pos, spacegroup=221, cellpar=3.905, size=(3, 3, 3))
>>> fig, ax = plt.subplots()
>>> ax.imshow(stem_image, cmap='gray')
>>> plot_atoms(srtio3, ax, radii=0.3, scale=6.3, offset=(47, 54), rotation=('90x,45y,56z'))
>>> ax.set_xlim(0, stem_image.shape[0])
>>> ax.set_ylim(0, stem_image.shape[1])
>>> ax.set_axis_off()
>>> fig.tight_layout()
>>> fig.savefig("iomatplotlib3.png")

.. image:: matplotlib_plot_atoms3.png

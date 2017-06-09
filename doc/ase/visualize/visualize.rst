.. module:: ase.visualize

Visualization
=============

.. function:: view(atoms, data=None, viewer=None, repeat=None)

This provides an interface to various visualization tools, such as
:mod:`ase.gui`, RasMol_, VMD_, gOpenMol_, Avogadro_, ParaView_ or NGLView_. The default viewer is
the ase.gui, described in the :mod:`ase.gui` module. The simplest invocation
is:

.. testsetup::

    from ase import Atoms
    atoms = Atoms('Cu')

>>> from ase.visualize import view
>>> view(atoms)

where ``atoms`` is any :class:`~ase.Atoms` object.  Alternative viewers
can be used by specifying the optional keyword ``viewer=...`` - use one of
'ase.gui', 'gopenmol', 'vmd', 'rasmol', 'paraview', 'ngl'.  The VMD and Avogadro viewers can
take an optional ``data`` argument to show 3D data, such as charge density:

>>> view(atoms, viewer='VMD', data=...)

The nglview viewer additionally supports any indexible sequence of :class:`~ase.Atoms`
objects, e.g. lists of structures and :class:`~ase.io.Trajectory` objects.

If you do not wish to open an interactive gui, but rather visualize
your structure by dumping directly to a graphics file; you can use the
``write`` command of the :mod:`ase.io` module, which can write 'eps',
'png', and 'pov' files directly, like this:

>>> from ase.io import write
>>> write('image.png', atoms)

.. _RasMol: http://openrasmol.org/
.. _VMD: http://www.ks.uiuc.edu/Research/vmd/
.. _gOpenMol: http://www.csc.fi/gopenmol/
.. _Avogadro: http://avogadro.openmolecules.net/
.. _ParaView: http://www.paraview.org/
.. _NGLView: https://github.com/arose/nglview


.. module:: ase.visualize.nglview

Viewer for Jupyter notebooks
----------------------------

The nglview viewer is a dedicated viewer for the Jupyter_ notebook interface.
It uses an embeddable NGL_ WebGL molecular viewer. The viewer works only in
the web browser environment and embeds the live javascript output into the
notebook. To utilize this functionality you need to have NGLView_ and
ipywidgets_ packages installed in addition to the Jupyter_ notebook.

The basic usage provided by the :func:`ase.visualize.view` function exposes
only small fraction of the NGL_ widget capabilities. The simplest form:

>>> view(atoms, viewer='ngl')

creates interactive ngl viewer widget with the few additional control widgets
added on the side. The object returned by the above call is a reference to
the `.gui` member of the :class:`ase.visualize.nglview.NGLDisplay` containing
actual viewer (`.view` member), a reference to control widgets box
(`.control_box` member) and
:func:`ase.visualize.view.nglview.NGLDisplay.custom_colors` method. The
notebook interface is not blocked by the above call and the returned object
may be further manipulated by the following code in the separate cell (the
`_` variable contains output from the previous cell):

>>> v=_
>>> v.custom_colors({'Mn':'green','As':'blue'})
>>> v.view._remote_call("setSize", target="Widget", args=["400px", "400px"])
>>> v.view.center_view()
>>> v.view.background='#ffc'
>>> v.view.parameters=dict(clipDist=-200)

The `.view` member exposes full API of the NGLView_ widget. The
`.control_box` member is a :class:`ipywidgets.HBox` containing
:class:`nglview.widget.NGLWidget` and :class:`ipywidgets.VBox` with control
widgets. For the full documentation of these objects consult the NGLView_,
NGL_ and ipywidgets_ websites.

.. _Jupyter: https://www.jupyter.org/
.. _NGL: https://github.com/arose/ngl
.. _ipywidgets: https://github.com/jupyter-widgets/ipywidgets

.. autoclass:: ase.visualize.nglview.NGLDisplay
   :inherited-members:

.. autofunction:: ase.visualize.nglview.view_ngl
.. automethod:: ase.visualize.nglview.NGLDisplay.custom_colors


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

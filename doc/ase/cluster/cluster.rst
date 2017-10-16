.. module:: ase.cluster

==========================
Nanoparticles and clusters
==========================

There are modules for creating nanoparticles (clusters) with a given crystal
structure by specifying either the number of layers in different directions,
or by making a Wulff construction.

Examples
========

Layer specification
-------------------

This example sets up a nanoparticle of copper in the FCC crystal structure,
by specifying 6 layers in the (100) directions, 9 in the (110) directions and
5 in the (111) directions::

  import ase
  from ase.cluster.cubic import FaceCenteredCubic

  surfaces = [(1, 0, 0), (1, 1, 0), (1, 1, 1)]
  layers = [6, 9, 5]
  lc = 3.61000
  atoms = FaceCenteredCubic('Cu', surfaces, layers, latticeconstant=lc)

|culayer|

.. |culayer| image:: culayer.png


Wulff construction
------------------

To set up a Wulff construction, the surface energies should be
specified, in units of energy per area (*not* energy per atom).  The
actual unit used does not matter, as only the ratio between surface
energies is important.  In addition, the approximate size of the
nanoparticle should be given.  As the Wulff construction is build from
whole layers, it is not possible to hit the desired particles size
exactly::

  from ase.cluster import wulff_construction

  surfaces = [(1, 0, 0), (1, 1, 0), (1, 1, 1)]
  esurf = [1.0, 1.1, 0.9]   # Surface energies.
  lc = 3.61000
  size = 1000  # Number of atoms
  atoms = wulff_construction('Cu', surfaces, esurf,
                             size, 'fcc',
                             rounding='above', latticeconstant=lc)

Note that the Wulff construction currently only work with cubic
lattices.


Creating a nanoparticle
=======================

The :mod:`ase.cluster` module contains a number of sub-modules for
defining clusters, one for each crystal structure.  They are all
called the same way, by specifying the element, the number of layers
in different directions, and optionally the lattice constant.

The layer specification is the only part that may not be intuitive.
It is given as two arrays, one specifying the Miller indices of the
surfaces, and one specifying the number of layers from the center of
the cluster to the respective surfaces.

The surface specification allows for one or more surfaces of a given
family of surfaces to be different from the other surfaces.  This can
be used e.g. to create a cluster where one part has been truncated by
a substrate.  This is done by *first* specifying the number of layers
for the family of surfaces, and *later* specifying the number of
layers for a given surface.  Consider the surface specification

::

   surfaces = [(1, 0, 0), (1, 1, 1), (1, -1, 1)]
   layers = [6, 5, -1]
   atoms = FaceCenteredCubic('Cu', surfaces, layers)

Here we first ask for 6 layers for the {100} surface family, i.e. the
directions (100), (010), (001), (-1,0,0), etc.  Then we ask for 5
layers for the {111} family of surfaces.  Finally, we change the
number of layers for the (-1,1,1) surface.  This is interpreted as a
single surface, since it is part of a family that has already been
specified.  Asking for a negative number of layers is allowed, this
cause the particle to be truncated *before* its center point.  The
result is seen below.

|truncated|

.. |truncated| image:: truncated.png



The functions for creating nanoparticles take the following
arguments:

  ``symbols``: A string specifying the element (or a tuple of strings
  for compounds).

  ``surfaces``: A list of surfaces, as explained above.

  ``layers``:  A corresponding list of the number of layers to be
  included.

  ``vacuum=0.0``:  The amount of vacuum to include around the particle.
  Defaults to 0.

  ``latticeconstant=None``:  The lattice constant of the lattice.  If
  not specified, the experimental value from :mod:`ase.data` is used.


Possible crystal structures
---------------------------

You select the crystal structure by selecting the right function for
creating the nanoparticle.  Currently, these modules only work for the
three cubic crystal structures: FaceCenteredCubic, BodyCenteredCubic,
and SimpleCubic.  Other structures are implemented, but do currently
not work correctly.


Wulff constructions
===================

As an alternative to specifying the number of layers, a Wulff
construction can be used to create a nanoparticle (with cubic
symmetry).  The function can be imported as::

  from ase.cluster import wulff_construction

.. autofunction:: ase.cluster.wulff_construction

=====
Tools
=====

Graphs
------

Allows to graph different quantities for a given trajectory. A 'save' button
also gives the opportunity to save the data to file.

This example plots the maximal force for each image i and could help in
investigating the convergence properties for relaxations:

::

  i, e-min(E), fmax

These are the symbols that can be used:

================ ==================================================
 Symbol          Interpretation
================ ==================================================
e                total energy
epot             potential energy
ekin             kinetic energy
fmax             maximum force
fave             average force
d(n1,n2)         distance between two atoms
R[n,0-2]         position of atom number n
i                current image number
E[i]             energy of image number i
F[n,0-2]         force on atom number n
M[n]             magnetic moment of atom number n
A[0-2,0-2]       unit-cell basis vectors
s                path length
a(n1,n2,n3)      tangle between atoms n1, n2 and n3, centered on n2
dih(n1,n2,n3,n4) dihedral angle between n1, n2, n3, and n4
T                temperature (requires velocity)
================ ==================================================


Movie
-----

Allows to play the current trajectory as a movie using a number of
different settings. Default duration is 5 s.


Constraints
-----------

Allows to set (or remove) constraints based on the currently selected atoms.


Render scene
------------

(Currently disabled)

Graphical interface to the ASE povray interface, ideally it requires
that povray is installed on your computer to function, but it also can
be used just to export the complete set of povray files.

The texture of each atom is adjustable: The default texture is applied
to all atoms, but then additional textures can be defined based on
selections (``Create new texture from current selection``). These can
be obtained either from selecting atoms by hand or by defining a
selection with a boolean expression, for example ``Z==6 and x>5 and
y<0`` will select all carbons with coordinates x>5 and y<0. The
available commands are listed in the ``Help on textures`` window.

A movie-making mode (``render all N frames``) is also available. After
rendering, the frames can be stitched together using the ``convert``
unix program e.g.

::

    localhost:doc hanke$ convert -delay 4.17 temp.*.png temp.gif

For this particular application it might be a good idea to use a white
background instead of the default transparent option.


Move atoms
----------

Allows selected atoms to be moved using the arrow keys. The direction
is always parallel to the plane of the screen. Two possible movements
are available: Just pressing the arrow keys will move by 0.1
Angstrom, ``shift`` + arrow keys will move by 0.01 Angstrom.

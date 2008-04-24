Bader Analysis
--------------

Henkelman et al. have implemented a fast and robust algorithm for
calculating the electronic charges on individual atoms in molecules or
crystals using the Bader scheme [#bader]_, [#improved_bader]_. In that
method the electron density is analyzed and so-called zero-flux
surfaces are used to divide a system into atoms. The program can be
downloaded from http://theory.cm.utexas.edu/bader/ .

This algorithm is very well suited for large solid state physical
systems as well as large biomolecular systems. The computational time
depends only on the size of the 3D grid used in the calculations and
typically it takes less than a minute to do the analysis. To get more
accurate results in the charge analysis it is recommended to use a
finer grid than the default values, however, usually it is enough to
calculate the electronic structure with a finer grid after the
geometry optimization has been done with the default values.

The program takes cube input files. The program does *not* support
units, and assumes atomic units.

All ase calculators should have a ``get_pseudo_density`` method which
can be used to get the density.

A simple python script for making a cube file, ready for the Bader
program, could be:

>>> from ase import *
>>> density = calc.get_pseudo_density() * Bohr**3
>>> write('filename.cube', atoms, data=density)

Some calculators also have a method called
``get_all_electron_density``, in which case this is preferable to
``get_pseudo_density``.

If you use pseudo densities, you have to be carefull. There can be a
problem with e.g. O-H bonds (and possibly C-H bonds too).


.. [#bader] G. Henkelman, A. Arnaldsson, and H. Jónsson, A fast and
            robust algorithm for Bader decomposition of charge
            density, *Comput. Mater. Sci.* **36** 254-360 (2006)

.. [#improved_bader] E. Sanville, S. D. Kenny, R. Smith, and
                     H. Henkelman, An improved grid-baed algorithm for
                     Bader charge allocation, *J. Comp. Chem.* **28**
                     899-908 (2007)

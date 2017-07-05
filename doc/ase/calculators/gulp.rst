.. module:: ase.calculators.gulp

====
GULP
====

GULP_, the General Utility Lattice Program, is a forcefield code.


Instructions
============

Make sure that the following variables are set correctly:

 * ``$GULP_LIB`` : path to the folder containing the potential files
 * ``$ASE_GULP_COMMAND="/path/to/gulp < PREFIX.gin > PREFIX.got"``

The latter defaults to ``gulp < PREFIX.gin > PREFIX.got``, assuming that
the gulp command is in ``$PATH``.

GULP uses several variables to write the input file.  The most
important ones are:

 * ``keywords``: string of space-separated keywords, e.g.: ``'angle compare comp conp'``
 * ``options``: list of strings with options, e.g.: ``['xtol opt 5.0']``
 * ``library``:  filename of the potential library file (e.g. ``'reaxff.lib'``).  Recall that the search path for this library is given by ``$GULP_LIB``.

If the potential uses different atom types, one must define an instance of conditions
with the atom object to label them and add rules to rename atoms. The only rule suported
right now is ``min_distance_rule`` which selects atoms of one type close to atoms of another type and
renames them. For example::

  c = Conditions(atoms)
  c.min_distance_rule('O', 'H', ifcloselabel1='O2', ifcloselabel2='H',
                      elselabel1='O1')
  calc = GULP(conditions=c)

This will assign, for each H in the system, a corresponding O atom as O2,
selecting those O which are the closest to H. The rest of oxygens will
be labeled O1.

Finally, if the potential file requires the use of shells, the
variable ``shel`` must be used::

  GULP(..., shel=['O'])

Example
=======

Here is an example of how to use the GULP calculator.

.. literalinclude:: gulp_example.py


The script performs a single-point calculation and then an optimization
of the H8Si8O20 double ring
structure within ASE using the ffsioh potential.  The
method ``get_optimizer()`` returns an object which works like an ASE optimizer,
but actually triggers an optimization using GULP's internal optimizer.
This works on GULP calculators that have the ``'opti'`` keyword.


_GULP: http://gulp.curtin.edu.au/gulp/

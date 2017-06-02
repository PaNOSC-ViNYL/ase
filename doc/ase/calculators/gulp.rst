.. module:: ase.calculators.gulp

====
GULP
====

GULP_, the General Utility Lattice Program, is a forcefield code.
Brief description.

Instructions
============

It is mandatory to have set-up the following variables:

 * $GULP_LIB : path to the folder containing the potential files 
 * ASE_GULP_COMMAND="/path/to/gulp < PREFIX.gin > PREFIX.got"

GULP uses several variables to write the input file, the most important ones are:

keywords = string of keywords, I.e: 'angle compare comp conp'
options = [] list of strings with options, I.e: ['xtol opt 5.0']
library = String with the name of the potential on the library file

If the potential uses different atom types, one must define an instance of conditions
with the atom object to label them and add rules to rename atoms. The only rule suported 
right now is "min_distance_rule" which selects the Elements1 closest to Elements2 and 
renames them by Tag1 while keeping the rest as Tag2. I.e: 
c = Conditions(atoms)
c.min_distance_rule('O', 'H', ifcloselabel1='O2', ifcloselabel2='H', elselabel1='O1')
calc = GULP(conditions = c)
Would assign as many O2 types as the number of H in the system selecting those O which
are the closest to H. The rest of oxygens will be labeled O1.

Finally, if the potential file requires the use of shells, the variable shel must be used:

shel = ['O']

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

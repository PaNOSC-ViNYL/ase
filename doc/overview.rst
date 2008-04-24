================
A quick overview
================

This section gives a quick overview of what ASE can do.  For more details:

* Look at the documentation for the individual :mod:`modules <ase>`.
* See the automatically generated documentation: :epydoc:`ase`.
* Browse the `source code`_ online.


.. _source code: http://trac.fysik.dtu.dk/projects/ase/browser/trunk


-----
Atoms
-----

The ``Atoms`` object is a collection of atoms.  Here is how to define
a CO molecule::

  from ase import Atom, Atoms
  d = 1.1
  co = Atoms('CO', positions=[(0, 0, 0), (0, 0, d)])

Here, the first argument specifies the type of the atoms and we used
the ``positions`` keywords to specify their positions.  Other
possible keywords are: ``numbers``, ``tags``, ``momenta``, ``masses``,
``magmoms`` and ``charges``.

Here is how you could make an infinite gold wire with a bond length of
2.9 Å::

  from ase import *
  d = 2.9
  L = 10.0
  wire = Atoms('Au',
               positions=[(L / 2, L / 2, 0)],
               cell=(L, L, d),
               pbc=(0, 0, 1))

Here, two more optional keyword arguments were used:

``cell``: Unit cell size
  This can be a sequence of three numbers for
  an orthorhombic unit cell or three by three numbers for a general
  unit cell (a sequence of three sequences of three numbers).  The
  default value is ``[1.0, 1.0, 1.0]``.

``pbc``: Boundary conditions
  The default value is ``False`` - a value of ``True`` would give
  pbc boundary conditions along all three axes.  It is possible
  to give a sequence of three booleans to specify periodicity along
  specific axes.

You can also use the following methods to work with the unit cell and the
boundary conditions:

``get_cell()``:
  Returns a three by three array.

``set_cell(cell, fix=False)``:
  Change the size of the unit cell.  If the optional argument ``fix``
  is ``True`` (defaults to ``False``), then the positons of the atoms
  are fixed, otherwise the atoms are moved so that their positions
  relative to the unit cell are kept.

``get_pbc()``:
  Return periodic boundary condition flags as three booleans.

``set_pbc()``:
  Set the periodic boundary condition flag.

Here is how you would do bulk ruthenium (hcp)::

  from math import sqrt
  a = 2.70
  b = a * sqrt(3) / 2
  c = 1.59 * a
  bulk = Atoms('Ru2',
               positions=[(0,     0,     0    ),
                          (a / 2, b / 3, c / 2)],
               pbc=True,
               cell=[(a,     0, 0),
                     (a / 2, b, 0),
                     (0,     0, c)])

In addition, an ``Atoms`` instance has the following methods:
``append``, ``center``, ``copy``, ``distance``, ``extend``,
``get_center_of_mass``, ``pop``, ``rattle``, ``repeat``, ``rotate``
and ``translate``.  See `Atom Manipulations` for more details.


-----
Units
-----

The units used for length, energy and mass are Å, eV and atomic mass
units.  To convert to/from other units, use the constants:  ``nm``,
``Bohr``, ``Hartree``, ``Rydberg``, ``kJ``, ``kcal``, ``mol``, ``fs``,
``kB``.

>>> 2 * Bohr
1.0583545150138329
>>> 25 * Rydberg
340.14244569396635
>>> 100 * kJ/mol
1.0364272141304978
>>> 300 * kB
0.025852157076770025
>>> 0.1 * fs
0.0010180506973856182



------------
Input-output
------------

Writing the atomic positions to a file is done with the ``write``
function::

  write('CO.xyz', co)

This will write a file in the xyz-format.  Other possible formats:

========  ===========================
format    description
========  ===========================
``xyz``   Simple xyz-format
``cube``  Gaussian cube file
``pdb``   Protein data bank file
``traj``  ASE's own trajectory format
``py``    Python script
========  ===========================


The ``write`` function will choose the format from the given filename
or from the optional ``format`` argument::

  write('CO.1', co, format='cube')

Reading from file is done like this::

  co = read('CO.xyz')

Some file-formats can take several configurations (trajectories), and
the default behavior of the ``read`` function is to return the last
configuration::

  read('A.traj', 0)   # first configuration
  read('A.traj', -1)  # last configuration
  read('A.traj')      # same as above

In addition to the ``xyz``, ``cube``, and ``traj`` formats, ASE can read and understand the following types of files:

=======  =================================
format   description
=======  =================================
``nc``   Old ASE-2 NetCDF trajectory files
``gpw``  GPAW restart files
``txt``  GPAW text output
``nc``   Dacapo NetCDF output
``out``  Dacapo text output
=======  =================================



Gaussian Cube file format
-------------------------

The Gaussian Cube file format describes volumetric data as well as
atom positions, it originates from the Gaussian software package.  The
volume data should be a 3 dimensional ndarray describing the
volumetric data for the unit cell, given in the Atoms object::

  write('x.cube', co, data=a)

Here ``a`` is the ndarray.  If the array has complex numbers, then the
absolute vale is written.  Use::

  from numpy import angle
  write('xp.cube', co, data=angle(a))

to write the phases.

Reading back in the data from a cube file is done like this::

  from ase.io.cube import read_cube_data
  co, a = read_cube_data('x.cube')

As can be seen, the ``read_cube_data`` function returns both the atoms
object and the ndarray.




Trajectories
------------

Molecular trajectories are useful for storing result from molecular
dynamics runs, structure optimization runs or the configurations along
a minimum energy path from reactant to product.

::

  traj = PickleTrajectory('CO.traj', 'w', co)
  for i in range(10):
      # do something to the CO molecule
      traj.write()

This will write 10 configurations to the 'CO.traj' file.  Read it like this::

  traj = PickleTrajectory('CO.traj', 'r')
  co = traj[-1]
  co = traj[9]  # same as above




-----------
Calculators
----------- 

In order to calculate froces and energies, you need to attach a calculator object to your atoms object:

>>> co.set_calculator(EMT())
>>> co.get_potential_energy()
-0.1577706320763923
>>> co.get_forces()
array([[  0.        ,   0.        , -16.76090913],
       [  0.        ,   0.        ,  16.76090913]])
>>> co.positions
array([[ 0. ,  0. ,  0. ],
       [ 0. ,  0. ,  1.1]])
>>> co.positions[1, 2] = 1.2
>>> co.get_forces()
array([[ 0.        ,  0.        ,  1.38699718],
       [ 0.        ,  0.        , -1.38699718]])

Here, we used an effective medium theory calculator to calculate
energies and forces.  There are currently five :mod:`calculators` that
can be used with ASE: :class:`EMT`, Asap_, GPAW_, Dacapo_,
:class:`Siesta` (and MMTK - work in progress).
  
.. _Asap: http://wiki.fysik.dtu.dk/Asap
.. _Dacapo: http://wiki.fysik.dtu.dk/dacapo
.. _GPAW: http://wiki.fysik.dtu.dk/gpaw




-----------------------
The ``ase.data`` module
-----------------------

This module defines the following variables: ``atomic_masses``, ``atomic_names``, ``chemical_symbols``, ``covalent_radii``, ``cpk_colors`` and  ``reference_states``.  All of these are lists that should be indexed with an atomic number:

>>> from ase.data import *
>>> atomic_names[92]
'Uranium'
>>> atomic_masses[2]
4.0026000000000002

If you don't know the atomic number of some element, then you can look it up in the ``atomic_numbers`` dictionary:

>>> atomic_numbers['Cu']
29
>>> covalent_radii[29]
1.1699999999999999



----------------------
Structure optimization
----------------------





QuasiNewton
-----------

The ``QuasiNewton`` object is one of the minimizers in the ASE
package.  Let's try to use it to optimize the structure of a water
molecule.  We start with the experimental geometry::

  from ase import *
  d = 0.9575
  t = pi / 180 * 104.51
  water = Atoms('H2O',
                positions=[(d, 0, 0),
                           (d * cos(t), d * sin(t), 0),
                           (0, 0, 0)],
                calculator=EMT())
  dyn = QuasiNewton(water)
  dyn.run(fmax=0.05)
  QuasiNewton:   0        6.445801      51.6847
  QuasiNewton:   1        2.418583      27.2946
  QuasiNewton:   2        0.551767      12.1607
  QuasiNewton:   3       -0.039301       4.0520
  QuasiNewton:   4       -0.128045       0.8479
  QuasiNewton:   5       -0.132312       0.0397

When doing structure optimization, it is useful to write the
trajectory to a file, so that the progress of the optimization run can
be followed during or after the run::

  traj = PickleTrajectory('H2O.traj', 'w', water)
  dyn = QuasiNewton(water)
  dyn.attatch(traj.write)
  dyn.run(fmax=0.05)
  
Use the command ``ag H2O.traj`` to see what is going on (more here: ase.gui_).

The ``attach`` method takes an optional argument ``interval=n`` that can
be used to tell the structure optimizer object to write the
configuration to the trajectory file only every ``n`` steps.


Hessian ...

Restart ...


LBFGS
-----

...

FIRE
----

...

MDMin
-----

The MDmin algorithm is a modification of the usual velocity-Verlet
molecular dynamics algorithm.  Newtons second law is solved
numerically, but after each time step the dot product between the
forces and the momenta is checked.  If it is zero, the system has just
passed through a (local) minimum in the potential energy, the kinetic
energy is large and about to decrease again.  At this point, the
momentum is set to zero.  Unlike a "real" molecular dynamics, the
masses of the atoms are not used, instead all masses are set to one.

The MDmin algorithm exists in two flavors, one where each atom is
tested and stopped individually (QuickMinAtomByAtom in the old ASE),
and one where all coordinates are treated as one long vector, and all
momenta are set to zero if the dotproduct between the momentum vector
and force vector (both of length 3N) is zero (QuickMinAllCoordinates
in the old ASE).  This module implements the latter version.

Although the algorithm is primitive, it performs very well because it
takes advantage of the physics of the problem.  Once the system is so
near the minimum that the potential energy surface is approximately
quadratic it becomes advantageous to switch to a minimization method
with quadratic convergence, such as `Conjugate Gradient` or `Quasi
Newton`.



------------------
Molecular dynamics
------------------

Let us look at an example: Molecular dynamics with a water molecule:

>>> from ase import *
>>> d = 1.1
>>> n2 = Atoms('N2', positions=[(0, 0, 0), (d, 0, 0)],
...            calculator=EMT())
>>> dyn = VelocityVerlet(n2)
>>> for i in range(10):
...     pot = n2.get_potential_energy()
...     kin = n2.get_kinetic_energy()
...     print '%2d: %.5f eV, %.5f eV, %.5f eV' % (i, pot + kin, pot, kin)
...     dyn.run(dt=1.0 * fs, steps=20)
... 
 0: 0.04217 eV, 0.04217 eV, 0.00000 eV
 1: 0.04216 eV, 0.02159 eV, 0.02057 eV
 2: 0.04216 eV, 0.00009 eV, 0.04206 eV
 3: 0.04216 eV, 0.01637 eV, 0.02580 eV
 4: 0.04217 eV, 0.04045 eV, 0.00171 eV
 5: 0.04217 eV, 0.03297 eV, 0.00920 eV
 6: 0.04216 eV, 0.00585 eV, 0.03631 eV
 7: 0.04216 eV, 0.00497 eV, 0.03718 eV
 8: 0.04217 eV, 0.03392 eV, 0.00825 eV
 9: 0.04217 eV, 0.03804 eV, 0.00413 eV

The ``dyn`` object has a method called ``run(dt, steps)`` that takes
two arguments:  A time step for the
integration of Newtons equation and the number of steps to take.
Here, we take 20 steps of length 1 fs.  Since we didn't set any
initial velocities for the nitrogen molecule, the kinetic energy
starts at 0.0 eV.  Notice also that the total energy is conserved to
within 0.1 meV at this time step.


Initial velocities
------------------

...


More advanced MD algorithms
---------------------------

...




-----------
Constraints
-----------

Applying constraints to the atomic positions can be useful in many
cases.  Let's look at a simple example:  We want to relax the bond
length of a nitrogen molecule with the first atom fixed at the
position (0, 0, 0):

>>> from ase import *
>>> from ase import *
>>> d = 1.1
>>> n2 = Atoms('N2', positions=[(0, 0, 0), (d, 0, 0)],
...            calculator=EMT(),
...      constraint=FixAtoms(indices=[0]))
>>> QuasiNewton(n2).run(fmax=0.01)
QuasiNewton:   0        0.042171       2.9357
QuasiNewton:   1        0.001205       0.4725
QuasiNewton:   2        0.000009       0.0396
QuasiNewton:   3        0.000000       0.0006
>>> print n2.get_positions()
[[ 0.          0.          0.        ]
 [ 1.12958567  0.          0.        ]]

It's also possible to attach a constraint to an ``Atoms`` object using
the ``set_constraint()`` method.  These three are euivalent::

  n2.set_constraint(FixAtoms(indices=[0]))
  n2.set_constraint(FixAtoms(mask=[True, False]))
  n2.set_constraint(FixAtoms(mask=[1, 0]))


Fix atoms
---------

We have just seen how to use the ``FixAtoms`` constraint.  Use
``mask=[...]`` where the list contains one boolean flag for each atom
indicating wheter this atom should be fixed or free (use 1 or ``True``
to fix the atoms and 0 or ``False`` for atoms free to move).
Alternatively, use ``indices=[...]``, where the list contains the
indices ofthe atoms that should be fixed - as always in Python code,
the first atoms has index zero.

Fix a bond length
-----------------

The ``FixBondLength`` constraint can fix a distance between two atoms.
You construct the constraint  like this::

  constraint = FixBondLength(5, 6)
  molecule.set_constraint(constraint)

This will fix the distance between atoms number 5 and 6.






-------------
Visualization
-------------

Use the ``view`` function to visualize the atoms::

  view(atoms)

This will pop up an `ase.gui`_ window.  Alternative viewers can be used
by specifying the optional keyword ``virwer=...`` - use one of
'ase.gui', 'gopenmol', 'vmd', or 'rasmol'.  The VMD viewer can take an
optional ``data`` argument to show 3D data::

  view(atoms, viewer='VMD', data=array)


.. _ase.gui: :mod:`ase.gui`


VTK
---

...


PNG and EPS files
-----------------

...

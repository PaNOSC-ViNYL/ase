. _overview:

================
A quick overview
================

This section gives a quick overview of what ASE can do.  For more details:

* Look at the documentation for the individual :mod:`modules <ase>`.
* See the automatically generated documentation: :epydoc:`ase`.
* Browse the `source code`_ online.
* Read the :ref:`tutorials`.


.. _source code: http://trac.fysik.dtu.dk/projects/ase/browser/trunk


-----
Atoms
-----

The ``Atoms`` object is a collection of atoms.  Here is how to define
a N2 molecule by directly specifying the position of two nitrogen
atoms::

  from ase import Atoms
  d = 1.10
  molecule = Atoms('2N', positions=[(0., 0., 0.), (0., 0., d)])

You can also build crystals using, for example, the lattice module
which returns ``Atoms`` objects corresponding to common crystal
structures. Let us make a Cu (111) surface::

  from ase.lattice.surface import *
  slab = fcc111('Cu', size=(4,4,2), vacuum=10.0)



-----------
Calculators
----------- 

There are currently five :mod:`calculators` that can be used with ASE:
:class:`EMT`, Asap_, GPAW_, Dacapo_, :class:`Siesta` (Abinit and MMTK
- work in progress).
  
.. _Asap: http://wiki.fysik.dtu.dk/Asap
.. _Dacapo: http://wiki.fysik.dtu.dk/dacapo
.. _GPAW: http://wiki.fysik.dtu.dk/gpaw

In this overview we use an effective medium theory (EMT) calculator,
as it is very fast and hence useful for getting started.

We can attach a calculator to the previously created ``Atoms`` objects::

  from ase import EMT
  slab.set_calculator(EMT())
  molecule.set_calculator(EMT()) 

and use it to calculate the total energies for the systems by using
the ``get_potential_energy`` method from the ``Atoms`` class::

  e_slab = slab.get_potential_energy()
  e_N2 = molecule.get_potential_energy()


--------------------
Structure relaxation
--------------------

Let's use the ``QuasiNewton`` minimizer to optimize the structure of
the N2 molecule adsorbed on the Cu surface. First add the adsorbate to
the Cu slab, for example in the on-top position::
  
  h = 1.85
  add_adsorbate(slab, molecule, h, 'ontop')

In order to speed up the relaxation, let us keep the Cu atoms fixed in
the slab by using the ``FixAtoms`` constraint. Only the N2 molecule is
then allowed to relax to the equilibrium structure::

  constraint = FixAtoms(mask=[a.symbol != 'N' for a in slab])
  slab.set_constraint(constraint)

Now attach the ``QuasiNewton`` minimizer to the system and save the
trajectory file. Run the minimizer with the convergence criteria that
the force on all atoms should be less than some ``fmax``::

  dyn = QuasiNewton(slab, trajectory='ontop.traj')
  dyn.run(fmax=0.05)


------------
Input-output
------------

Writing the atomic positions to a file is done with the ``write``
function::

  write('slab.xyz', slab)

This will write a file in the xyz-format.  Possible formats are:

========  ===========================
format    description
========  ===========================
``xyz``   Simple xyz-format
``cube``  Gaussian cube file
``pdb``   Protein data bank file
``traj``  ASE's own trajectory format
``py``    Python script
========  ===========================

Reading from a file is done like this::

  slab_from_file = read('slab.xyz')

If the file contains several configurations, the default behavior of
the ``read`` function is to return the last configuration. However, we
can load a specific configuration by doing::

  read('slab.traj')      # last configuration
  read('slab.traj', -1)  # same as above
  read('slab.traj', 0)   # first configuration


-------------
Visualization
-------------

The simplest way to visualize the atoms is the ``view`` function::

  view(slab)

This will pop up a :mod:`gui` window.  Alternative viewers can be used
by specifying the optional keyword ``viewer=...`` - use one of
'ase.gui', 'gopenmol', 'vmd', or 'rasmol'.  The VMD viewer can take an
optional ``data`` argument to show 3D data::

  view(slab, viewer='VMD', data=array)


------------------
Molecular dynamics
------------------

Let us look at the nitrogen molecule as an example of molecular
dynamics with the ``VelocityVerlet`` algorithm. We first create the
:class:`VelocityVerlet` object giving it the molecule and the time
step for the integration of Newton's law. We then perform the dynamics
by calling its :meth:`run` methodand giving it the number of steps to
take::

  dyn = VelocityVerlet(molecule, dt=1.0 * fs)
  for i in range(10):
     pot = molecule.get_potential_energy()
     kin = molecule.get_kinetic_energy()
     print '%2d: %.5f eV, %.5f eV, %.5f eV' % (i, pot + kin, pot, kin)
     dyn.run(steps=20)


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
0.009822693531550318



-----------------------
The ``ase.data`` module
-----------------------

This module defines the following variables: ``atomic_masses``,
``atomic_names``, ``chemical_symbols``, ``covalent_radii``,
``cpk_colors`` and ``reference_states``.  All of these are lists that
should be indexed with an atomic number:

>>> atomic_names[92]
'Uranium'
>>> atomic_masses[2]
4.0026000000000002

If you don't know the atomic number of some element, then you can look
it up in the ``atomic_numbers`` dictionary:

>>> atomic_numbers['Cu']
29
>>> covalent_radii[29]
1.1699999999999999











Gaussian Cube file format
-------------------------

The Gaussian Cube file format describes volumetric data as well as
atom positions, it originates from the Gaussian software package.  The
volume data should be a 3 dimensional :term:`ndarray` describing the
volumetric data for the unit cell, given in the Atoms object::

  write('x.cube', co, data=a)

Here *a* is the ndarray.  If the array has complex numbers, then the
absolute vale is written.  Use::

  write('xp.cube', co, data=angle(a))

to write the phases.

Reading back in the data from a cube file is done like this::

  from ase.io.cube import read_cube_data
  co, a = read_cube_data('x.cube')

As can be seen, the ``read_cube_data`` function returns both the atoms
object and the ndarray.


-----------
Constraints
-----------

Applying constraints to the atomic positions can be useful in many
cases.  Let's look at a simple example:  We want to relax the bond
length of a nitrogen molecule with the first atom fixed at the
position (0, 0, 0):

>>> d = 1.1
>>> n2 = Atoms('N2', positions=[(0, 0, 0), (d, 0, 0)],
...            calculator=EMT(),
...            constraint=FixAtoms(indices=[0]))
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



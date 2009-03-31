.. _overview:

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

The :class:`~ase.atoms.Atoms` object is a collection of atoms.  Here
is how to define a N2 molecule by directly specifying the position of
two nitrogen atoms::

  from ase import Atoms
  d = 1.10
  molecule = Atoms('2N', positions=[(0., 0., 0.), (0., 0., d)])

You can also build crystals using, for example, the lattice module
which returns :class:`~ase.atoms.Atoms` objects corresponding to
common crystal structures. Let us make a Cu (111) surface::

  from ase.lattice.surface import *
  slab = fcc111('Cu', size=(4,4,2), vacuum=10.0)



-----------
Calculators
----------- 

There are currently five :mod:`calculators` that can be used with ASE:
:mod:`~calculators.emt`, Asap_, GPAW_, Dacapo_,
:mod:`~calculators.siesta` (Abinit and MMTK - work in
progress).
  
.. _Asap: http://wiki.fysik.dtu.dk/Asap
.. _Dacapo: http://wiki.fysik.dtu.dk/dacapo
.. _GPAW: http://wiki.fysik.dtu.dk/gpaw

In this overview we use the effective medium theory (EMT) calculator,
as it is very fast and hence useful for getting started.

We can attach a calculator to the previously created
:class:`~ase.atoms.Atoms` objects::

  from ase import EMT
  slab.set_calculator(EMT())
  molecule.set_calculator(EMT()) 

and use it to calculate the total energies for the systems by using
the :meth:`~ase.atoms.Atoms.get_potential_energy` method from the
:class:`~ase.atoms.Atoms` class::

  e_slab = slab.get_potential_energy()
  e_N2 = molecule.get_potential_energy()


--------------------
Structure relaxation
--------------------

Let's use the :mod:`QuasiNewton <optimize.qn>` minimizer to optimize the
structure of the N2 molecule adsorbed on the Cu surface. First add the
adsorbate to the Cu slab, for example in the on-top position::
  
  h = 1.85
  add_adsorbate(slab, molecule, h, 'ontop')

In order to speed up the relaxation, let us keep the Cu atoms fixed in
the slab by using :class:`~constraints.FixAtoms` from the
:mod:`~ase.constraints` module. Only the N2 molecule is then allowed
to relax to the equilibrium structure::

  constraint = FixAtoms(mask=[a.symbol != 'N' for a in slab])
  slab.set_constraint(constraint)

Now attach the :mod:`QuasiNewton <optimize.qn>` minimizer to the
system and save the trajectory file. Run the minimizer with the
convergence criteria that the force on all atoms should be less than
some ``fmax``::

  dyn = QuasiNewton(slab, trajectory='ontop.traj')
  dyn.run(fmax=0.05)


------------
Input-output
------------

Writing the atomic positions to a file is done with the
:func:`~ase.io.write` function::

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
the :func:`~ase.io.write` function is to return the last
configuration. However, we can load a specific configuration by
doing::

  read('slab.traj')      # last configuration
  read('slab.traj', -1)  # same as above
  read('slab.traj', 0)   # first configuration


-------------
Visualization
-------------

The simplest way to visualize the atoms is the :func:`~visualize.view`
function::

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
dynamics with the :class:`VelocityVerlet <md.verlet.VelocityVerlet>`
algorithm. We first create the :class:`VelocityVerlet
<md.verlet.VelocityVerlet>` object giving it the molecule and the time
step for the integration of Newton's law. We then perform the dynamics
by calling its :meth:`run` method and giving it the number of steps to
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










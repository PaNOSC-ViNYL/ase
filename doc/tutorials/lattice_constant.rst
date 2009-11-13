.. _lattice_constant:

=========================
Finding lattice constants
=========================

.. seealso::

   :ref:`eos`.


HCP
===

Let's try to find the `a` and `c` lattice constants for HCP nickel
using the :mod:`EMT <emt>` potential.  First, we import the things we need
like the :func:`~ase.structure.bulk` function::

  from ase import *
  from ase.structure import bulk

Then we make a good intial guess for `a` and `c`::

  a0 = 2.5
  c0 = sqrt(8 / 3.0) * a0

We want to try three values for `a` and three for `c`::

  eps = 0.01
  strains = np.array([1 - eps, 1, 1 + eps])

and put the results in a trajectory::

  traj = PickleTrajectory('Ni.traj', 'w')

Finally, we do the 9 calculations::

  for a in a0 * strains:
      for c in c0 * strains:
          ni = bulk('Ni', 'hcp', a=a, covera=c / a)
          ni.set_calculator(EMT())
          ni.get_potential_energy()
          traj.write(ni)


Analysis
--------

We fit the energy to this expression:

.. math:: c_0 + c_1 a + c_2 c + c_3 a^2 + c_4 ac + c_5 c^2,

using the function:

.. autofunction:: ase.optimize.fitpoly2

Now, we need to extract the data from the trajectory.  Try this:

>>> from ase.structure import bulk
>>> ni = bulk('Ni', 'hcp', a=2.5, covera=4.0 / 2.5)
>>> ni.cell
array([[ 2.5       ,  0.        ,  0.        ],
       [ 1.25      ,  2.16506351,  0.        ],
       [ 0.        ,  0.        ,  4.        ]])

So, we can get `a` and `c` from ``ni.cell[0, 0]`` and ``ni.cell[2,
2]``::

  from ase import *
  from ase.optimize import fitpoly2
  energies = []
  aandc = []
  for atoms in read('Ni.traj@:'):
      energies.append(atoms.get_potential_energy())
      a = atoms.cell[0, 0]
      c = atoms.cell[2, 2]
      aandc.append((a, c))

  print fitpoly2(aandc, energies)

The result is `a=2.468` Å and `c=4.026` Å.  Using those as initial
guess, you get `a=2.470` Å and `c=4.008` Å which are the correct
numbers.

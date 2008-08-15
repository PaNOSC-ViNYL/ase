.. _siesta2:

----------------------------
Siesta 2: Molecular Dynamics
----------------------------

This exercise is intended to illustrate a Molecular Dynamics run with
the SIESTA calculator. The system is a Si(001) surface, in the 2x1
reconstruction with asymmetric dimers. The simulation cell contains
two dimers. An H2 molecule approaches the surface, above one of the
dimers, and dissociates, ending up with a H atom bonded to each of the
Si atoms in the dimer (and thus leading to a symmetric dimer). You can
get the ``xyz`` file with the initial geometry :svn:`here
<doc/exercises/siesta2/geom.xyz?format=raw>`.

.. literalinclude:: siesta2.py

Note that both H atoms are given an initial velocity towards the
surface through the lines::

  p = atoms.get_momenta()
  p[0,2]= -1.5 
  p[1,2]= -1.5 
  atoms.set_momenta(p)

Run the program, and check the results. You can visualize the dynamics
using the trajectory file with the help of the ASE :mod:`gui`. For
example you can visualize the behaviour of the potential and total
energies during the dynamics by doing::

  $ ag si001+h2.traj -g i,e,e+ekin

Check that the total energy is a conserved quantity in this
microcanonical simulation.In order to get a better picture you can
open VMD from the ``View`` menu and repeat the unit cell along the x
and y directions. Also try to visualize the atoms in the
ball-and-sticks fashion.

.. note::

  You can repeat your structure directly from the command line by
  typing :command:`ag si001+h2.traj -r 2,2,1`. This will repeat it 2
  times along x and y and just once along the z direction.

You can also try to repeat the simulation with a different
:mod:`dynamics <md>`. For instance you can model a canonical ensemble
with the Langevin dynamics, which couples the system to a heat bath::

  dyn = Langevin(atoms, 1 * fs, kB * 300, 0.002, trajectory=traj)







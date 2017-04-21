.. _md_tutorial:

==================
Molecular dynamics
==================

.. note::

  These examples *can* be used without Asap installed, then
  the ase.EMT calculator (implemented in Python) is used, but nearly
  superhuman patience is required.

Here we demonstrate now simple molecular dynamics is performed.  A
crystal is set up, the atoms are given momenta corresponding to a
temperature of 300K, then Newtons second law is integrated numerically
with a time step of 5 fs (a good choice for copper).

.. literalinclude:: moldyn1.py

Note how the total energy is conserved, but the kinetic energy quickly
drops to half the expected value.  Why?


Instead of printing within a loop, it is possible to use an "observer"
to observe the atoms and do the printing (or more sophisticated
analysis).

.. literalinclude:: moldyn2.py

Constant temperature MD
=======================

Often, you want to control the temperature of an MD simulation.  This
can be done with the Langevin dynamics module.  In the previous
examples, replace the line ``dyn = VelocityVerlet(...)`` with::

  dyn = Langevin(atoms, 5*units.fs, T*units.kB, 0.002)

where T is the desired temperature in Kelvin.  You also need to import
Langevin, see the class below.

The Langevin dynamics will then slowly adjust the total energy of the
system so the temperature approaches the desired one.

As a slightly less boring example, let us use this to melt a chunk of
copper by starting the simulation without any momentum of the atoms
(no kinetic energy), and with a desired temperature above the melting
point.  We will also save information about the atoms in a trajectory
file called moldyn3.traj.

.. literalinclude:: moldyn3.py

After running the simulation, you can study the result with the
command

::

  ase-gui moldyn3.traj

Try plotting the kinetic energy.  You will *not* see a well-defined
melting point due to finite size effects (including surface melting),
but you will probably see an almost flat region where the inside of
the system melts.  The outermost layers melt at a lower temperature.

.. note::

  The Langevin dynamics will by default keep the position and momentum
  of the center of mass unperturbed. This is another improvement over
  just setting momenta corresponding to a temperature, as we did before.


Isolated particle MD
====================

When simulating isolated particles with MD, it is sometimes preferable
to set random momenta corresponding to a specific temperature and let the
system evolve freely. With a relatively high temperature, the is however
a risk that the collection of atoms will drift out of the simulation box
because the randomized momenta gave the center of mass a small but
non-zero velocity too.

Let us see what happens when we propagate a nanoparticle for a long time:

.. literalinclude:: moldyn4.py

After running the simulation, use :ref:`ase-gui` to compare the results
with how it looks if you comment out either the line that says `Stationary(atoms)`, `ZeroRotation(atoms)` or both.

::

  ase-gui moldyn4.traj

Try playing the movie with a high frame rate and set frame skipping to a
low number. Can you spot the subtle difference?


.. _TIPnP Water Box Equillibration:

Equilibrating A TIPnP Water Box
===============================

This tutorial shows how to use the TIP3P and TIP4P force fields in
ASE.


Since the TIPnP type water interpotentials are for rigid
molecules, there are no intramolecular force terms, and we need to
constrain all internal degrees of freedom. For this, we're
using the RATTLE-type constraints of the :ref:`FixBondLengths` class to
constrain all internal atomic distances (O-H1, O-H2, and H1-H2) for
each molecule.

The box is equillibrated with the Langevin thermostat.


For efficiency, we first equillibrate a smaller box, and then repeat that
once more for the final equillibration. However, the potentials are not
parallelized, and are mainly included for testing and for use with QM/MM
tasks, so expect to let it run for some time.


The following is for TIP3P:

.. literalinclude:: tip3p_equil.py

.. note::

  The temperature calculated by ASE is assuming all degrees of freedom
  are available to the system. Since the constraints have removed the 3
  vibrational modes from each water, the shown temperature will be 2/3
  of the actual value.

The procedure for the TIP4P force field is the same, with the following
exceptions:

- the atomic sequence **must** be OHH, OHH, ...
- charges are set automatically

So to perform the same task using TIP4P, you simply have to import
that calculator instead:

::

    from ase.calculators.tip4p import TIP4P, rOH, thetaHOH


And remove the following line from the above script:

::

  set_tip3p_charges(atoms)


More info about the TIP4P potential: :mod:`ase.calculators.tip4p`

==================
Molecular dynamics
==================

.. module:: md

.. contents::

Typical computer simulations involve moving the atoms around, either
to optimize a structure (energy minimization) or to do molecular
dynamics.    This chapter discusses molecular dynamics,
energy minimization algorithms will be discussed in the next chapter XXX.


The dynamics is an object that operate on the atoms by moving them
according to their forces - it integrates Newton's second law
numerically.  A typical molecular dynamics simulation will use the
`Velocity Verlet dynamics`_.  You create the
:class:`VelocityVerlet` object, giving it the atoms and a time step, and then
you perform dynamics by calling its ``Run(n)`` method, ``n`` is the
number of time steps you want performed::

  from ase.md.verlet import VelocityVerlet
  dyn = VelocityVerlet(atoms, 5.0*femtosecond)
  dyn.run(1000)

A number of different algorithms can be used to perform molecular
dynamics, with slightly different results.  

Choosing the time step
======================

Al the dynamics objects need a time step.  Choosing it too small will
waste computer time, choosing it too large will make the dynamics
unstable, typically the energy increases dramatically (the system
"blows up").  If the time step is only a little to large, the lack of
energy conservation is most obvious in `Velocity Verlet dynamics`_,
where energy should otherwise be conserved.

Experience has shown that 5 femtoseconds is a good choice for most metallic
systems.  Systems with light atoms (e.g. hydrogen) and/or with strong
bonds (carbon) will need a smaller time step.

All the dynamics objects documented here are sufficiently related to
have the same optimal time step.


Constant NVE simulations (the microcanonical ensemble)
======================================================

Newton's second law preserves the total energy of the system, and a
straightforward integration of Newton's second law therefore leads to
simulations preserving the total energy of the system (E), the number
of atoms (N) and the volume of the system (V).  The most appropriate
algorithm for doing this is velocity Verlet dynamics, since it gives
very good long-term stability of the total energy even with quite
large time steps.  Fancier algorithms such as Runge-Kutta may give
very good short-term energy preservation, but at the price of a slow
drift in energy over longer timescales, causing trouble for long
simulations.

In a typical NVE simulation, the temperature will remain approximately
constant, but if significant structural changes occurs they may result
in temperature changes.  If external work is done on the system, the
temperature is likely to rise significantly.

Velocity Verlet dynamics
------------------------

.. class:: VelocityVerlet(atoms, timestep)


``VelocityVerlet`` is the only dynamics implementing the NVE ensemble.
It requires two arguments, the atoms and the time step.  Choosing
a two large time step will immediately be obvious, as the energy will
increase with time, often very rapidly.

Example: XXXX



Constant NVT simulations (the canonical ensemble)
=================================================

Since Newton's second law conserves energy and not temperature,
simulations at constant temperature will somehow involve coupling the
system to a heat bath.  This cannot help being somewhat artificial.
Two different approaches are possible within the ASE.  In Langevin
dynamics, each atom is coupled to a heat bath through a fluctuating
force and a friction term.  In Nosé-Hoover dynamics, a term
representing the heat bath through a single degree of freedom is
introduced into the Hamiltonian.

Langevin dynamics
-----------------

The ``Langevin`` module implements Langevin dynamics, where a (small)
friction term and a fluctuating force are added to Newton's second law
which is then integrated numerically.  The temperature of the heat
bath and magnitude of the friction is specified by the user, the
amplitude of the fluctuating force is then calculated to give that
temperature.  This procedure has some physical justification: in a
real metal the atoms are (weakly) coupled to the electron gas, and the
electron gas therefore acts like a heat bath for the atoms.  If heat
is produced locally, the atoms locally get a temperature that is
higher than the temperature of the electrons, heat is transferred to
the electrons and then rapidly transported away by them.  A Langevin
equation is probably a reasonable model for this process.

A disadvantage of using Langevin dynamics is that if significant heat
is produced in the simulation, then the temperature will stabilize at
a value higher than the specified temperature of the heat bath, since
a temperature difference between the system and the heat bath is
necessary to get a finite heat flow.  Another disadvantage is that the
fluctuating force is stochastic in nature, so repeating the simulation
will not give exactly the same trajectory.

When the ``Langevin`` object is created, you must specify a time step,
a temperature (in energy units) and a friction.  Typical values for
the friction are 0.01-0.02 atomic units.

::

  from ase.units import kB
  from ase.md.langevin import Langevin
  # Room temperature simulation
  dyn = Langevin(atoms, 5*femtosecond, kB*300, 0.002)

Both the friction and the temperature can be replaced with arrays
giving per-atom values.  This is mostly useful for the friction, where
one can choose a rather high friction near the boundaries, and set it
to zero in the part of the system where the phenomenon being studied
is located.



Nosé-Hoover dynamics
--------------------

In Nosé-Hoover dynamics, an extra term is added to the Hamiltonian
representing the coupling to the heat bath.  From a pragmatic point of
view one can regard Nosé-Hoover dynamics as adding a friction term to
Newton's second law, but dynamically changing the friction coefficient
to move the system towards the desired temperature.  Typically the
"friction coefficient" will fluctuate around zero.

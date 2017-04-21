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
exception: the atomic sequence **must** be OHH, OHH, ... .

So to perform the same task using TIP4P, you simply have to import
that calculator instead:

::

    from ase.calculators.tip4p import TIP4P, rOH, angleHOH

More info about the TIP4P potential: :mod:`ase.calculators.tip4p`

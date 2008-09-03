==================
Molecular dynamics
==================

Note that these examples *can* be used without Asap installed, then
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

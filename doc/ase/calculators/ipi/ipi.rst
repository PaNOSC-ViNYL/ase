.. module:: ase.calculators.ipi

===========================================
Communication with calculators over sockets
===========================================

ASE can use sockets to communicate efficiently with certain external
codes using the protocol of `i-PI <http://ipi-code.org/>`_.  This may
significantly speed up geometry optimizations, dynamics and other
algorithms in which ASE moves the atoms while the external code
calculates energies, forces, and stress.  Note that ASE does not
require i-PI, but simply uses the same protocol.

The reference article for i-PI is `Ceriotti, More, Manolopoulos, Comp. Phys. Comm. 185, 1019-1026 (2014) <http://dx.doi.org/10.1016/j.cpc.2013.10.027>`_.


Introduction
------------

Normally, file-IO calculators in ASE launch a new process to calculate
every atomic step.  This is inefficient since the codes will need to
either start from scratch or perform significant IO between steps.

Some codes can run in "driver mode" where a server provides atomic
coordinates through a socket connection, and the code returns
energies, forces, and stress to the server.  That way the startup
overhead is eliminated, and the codes can reuse and extrapolate
wavefunctions and other quantities for increased efficiency.

Which codes can be used with ASE/i-PI?
--------------------------------------

Below is a list of codes that can run as i-PI clients, and whether ASE
provides a calculator that supports doing so.

================ =========================================
Program name     Supported by ASE calculator
================ =========================================
Quantum Espresso Yes
FHI-aims         Yes
Siesta           Yes, presumably (untested)
DFTB+            Yes, presumably (untested)
Yaff             No; there is no ASE calculator for Yaff
cp2k             No; ASE uses cp2k shell instead
Lammps           No; ASE uses lammpsrun/lammpslib instead
ASE              Yes - ASE implements a client as well
================ =========================================

The codes that are "not supported" by ASE can still be used as
clients, but you will need to generate the input files and launch the
client programs yourself.

Codes may require different commands, keywords, or compilation options
in order to run in driver mode.  See the code's documentation for
details.  The i-PI documentation may also be useful.

How to use ASE/i-PI
-------------------

Example using Quantum Espresso

.. literalinclude:: ipi/ase_ipi_qe.py

Use ASE as a client
-------------------

ASE can run as an i-PI client using the IPIClient class.

Module documentation
--------------------

.. autoclass:: ase.calculators.ipi.IPICalculator

.. autoclass:: ase.calculators.ipi.IPIClient

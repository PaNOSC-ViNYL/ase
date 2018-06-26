.. module:: ase.calculators.socketio

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

ASE provides such a server in the form of a calculator.

Which codes can be used with socket I/O calculators?
----------------------------------------------------

Below is a list of codes that can run as clients, and whether ASE
provides a calculator that supports doing so.

================ =========================================
Program name     Supported by ASE calculator
================ =========================================
Quantum Espresso Yes
FHI-aims         Yes
Siesta           Yes
DFTB+            Yes, presumably (untested)
Yaff             No; there is no ASE calculator for Yaff
cp2k             No; ASE uses cp2k shell instead
Lammps           No; ASE uses lammpsrun/lammpslib instead
ASE              Yes - ASE provides a client as well
GPAW             Yes, using the ASE client
================ =========================================

The codes that are "not supported" by ASE can still be used as
clients, but you will need to generate the input files and launch the
client programs yourself.

Codes may require different commands, keywords, or compilation options
in order to run in driver mode.  See the code's documentation for
details.  The i-PI documentation may also be useful.

How to use the ASE socket I/O interface
---------------------------------------

Example using Quantum Espresso

.. literalinclude:: example_espresso.py

.. note::

   It is wise to ensure smooth termination of the connection.  This
   can be done by calling ``calc.close()`` at the end or, more
   elegantly, by enclosing using the ``with`` statement as done in all
   examples here.

Example using FHI-aims

.. literalinclude:: example_aims.py

Example using Siesta

.. literalinclude:: example_siesta.py

For codes other than these, see the next section.

Run server and client manually
------------------------------

ASE can run as a client using the SocketClient class.  This may be
useful for controlling calculations remotely or using a serial process
to control a parallel one.

This example will launch a server without (necessarily) launching any client:

.. literalinclude:: example_server.py

Run it and then run the client:

.. literalinclude:: example_client_gpaw.py

This also demonstrates how to use the interface with GPAW.
Instead of running the client script, it is also possible
to run any other program that acts as a client.  This
includes the codes listed in the compatibility table above.

Module documentation
--------------------

.. autoclass:: ase.calculators.socketio.SocketIOCalculator

.. autoclass:: ase.calculators.socketio.SocketClient

The SocketServer allows launching a server without the need
to create a calculator:

.. autoclass:: ase.calculators.socketio.SocketServer

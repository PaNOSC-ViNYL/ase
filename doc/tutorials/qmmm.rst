.. _qmmm:

=========================
ASE for QM/MM Simulations
=========================

QM/MM Simulations couple two (or more) descriptions to get total energy
and forces for the entire system in an efficiant manner. The method paper 
on our implementation is currently being written. 
General background can be found `here <http://www.python.org/>`_, here, and here. Examples of what this
code has been used for can be seen here, and here. 


This section will show you how to setup up various QM/MM simulations.
We will be using GPAW_ for the QM part. Other QM calculators should
be straightforwardly compatible with the subtractive-scheme SimpleQMMM
calculator, but for the Excplicit Interaction EIQMMM calculator, you
would need to be able to put an electrostatic external potential into
the calculator for the QM subsystem. 

.. _GPAW: http://wiki.fysik.dtu.dk/gpaw

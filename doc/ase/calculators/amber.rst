.. module:: ase.calculators.amber

Amber
=======

Introduction
------------------------

Amber is a powerfull classical simulations package
(http://ambermd.org).  It is not free, academic license costs $500.
Ase-Amber has been tested only for amber16 (2016). It can bee usefull
as MM part of QM/MM calculations since amber supports fast netCDF
fileIO.

.. _Amber: http://ambermd.org


Example
-------------------------

You need input files
(instructions for md, structure and topology, respectively):
:download:`mm.in<./mm.in>` :download:`2h2o.pdb<./2h2o.pdb>` :download:`mm.top<./2h2o.top>`

The actual ase-amber script:
	  
.. literalinclude:: ./test_amber.py

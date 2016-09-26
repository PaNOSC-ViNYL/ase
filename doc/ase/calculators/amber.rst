.. module:: ase.calculators.amber

Amber
=====

Introduction
------------

Amber_ is a powerfull classical simulations package.  It is not free, academic
license costs $500. Ase-Amber has been tested only for amber16 (2016). It can
bee usefull as MM part of QM/MM calculations since amber supports fast netCDF
fileIO.

.. _Amber: http://ambermd.org


Example
-------

You need input files (instructions for md, structure and topology,
respectively): :download:`mm.in <../../../ase/test/amber/mm.in>`,
:download:`2h2o.pdb <../../../ase/test/amber/2h2o.pdb>`, :download:`mm.top
<../../../ase/test/amber/2h2o.top>`.

The actual ase-amber script:

.. literalinclude:: ../../../ase/test/amber/amber.py

MMTK
====

Introduction
---------------

The Molecular Modelling Toolkit (MMTK) by Konrad Hinsen is a powerful program library for 
molecular simulation.

MMTK's functionality includes construction of molecular systems, including support for proteins 
and nucleic acids, dynamics, normal mode analysis and much more.

The ASE calculator interface enables the use of ASE dynamics/filters etc, with the MMTK 
force-fields, in this way the ASE functionality, like the Nudged Elastic Band method and Langevin 
dynamics, can be added to the MMTK functionality.

For general information on MMTK see the Molecular Modelling Toolkit homepage.


MMTK Calculator
---------------

The MMTK ASE calculator can not be attached to a general ListOfAtoms, but must use a ListOfAtom 
constructed from the MMTK atoms. The reason is that in constructing the force-fields additional 
information for each atoms must be present; a Carbon atom is not just a Carbon atom. This 
information is usually supplied in a PDB file.

Setting up the force-field starting from PDB files, the MMTK universe must be handle by MMTK 
commands.

The examples in the next section shows examples of setting up a MMTK universe.

Once the MMTK universe is constructed a ASE :class:`Atoms` can be constructed and attached to a 
MMTK calculator:

.. highlight:: python

::

  >>> atoms = MMTKListOfAtoms(universe)
  >>> atoms.set_calculator(MMTK())


Examples
--------
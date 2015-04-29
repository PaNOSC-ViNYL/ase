==============
Bader Analysis
==============

An example of Bader analysis: `The water molecule`_.

.. _The water molecule:: https://wiki.fysik.dtu.dk/gpaw/tutorials/bader/bader.html

You can attach the output charges from the bader program to the atoms
for further processing::

  from ase import *
  from ase.io.bader import attach_charges

  # define the molecule as above
  atoms = molecule('H2O')
  atoms.set_cell([7.5, 9, 9])
  atoms.center()

  # the next two lines are equivalent (only one needed)
  attach_charges(atoms)
  attach_charges(atoms, 'ACF.dat')

  for atom in atoms:
      print('Atom', atom.symbol, 'Bader charge', atom.charge)

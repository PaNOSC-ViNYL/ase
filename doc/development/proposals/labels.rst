================
Labels for atoms
================

**WIP**

This proposal describes how to introduce labels for the atoms in an Atoms
object.


Why?
====

Some atoms are special:

* ghost atoms
* atoms with a core hole
* atoms that need special basis set
* ...

and some atoms are not atoms at all (methyl group, ...).


Proposal
========

Introduce labels and kinds.  A label is a string and a kind is the combination
of a chemical symbol and a label like:

* 'H:ghost'
* 'X:methyl' (same as 'methyl')
* 'Cu:' (same as 'Cu')

We add symbols, labels and kinds list-of-string like attributes to the Atoms
object as well as new labels and kinds keyword arguments to the Atoms
constructor.  Currently, the first argument to the Atoms constructor
('symbols') will accept a chemical formula as a string or a list of chemical
symbols or atomic numbers.  Extend this to also accept kinds.

.. note::

    Not really related to this label issue, but since we are anyway adding new
    attributes, we should also add masses, tags, initial_magmoms, momenta,
    initial_charges and mass.

Examples
========

>>> a1 = Atoms(['N', 'C', 'methyl'])
>>> a1.positions[:, 0] = [0, 1.2, 2.6]
>>> a1.masses[a1.labels == 'methyl'] = Atoms('CH3').mass
>>> a1.numbers
array([7, 6, 0])
>>> a1.symbols
Symbols(['N', 'C', 'X'])
>>> a1.labels
Labels(['', '', 'methyl'])
>>> a1.kinds
Kinds(['N', 'C', 'X:methyl'])
>>> a.masses[a.labels == 'deuterium'] = 2.0

A DFT code could use the kinds to select pseudo-potentials:

>>> n2 = molecule('N2')
>>> n2.labels[0] = 'core-hole'
>>> n2.calc = DFT(xc='LDA',
...               pp={'N:core-hole': '~/mypps/N-LDA.core-hole.pp'})


List-like objects
=================

...


Atom objects?
=============

...

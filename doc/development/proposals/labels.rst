================
Labels for atoms
================

**WORK IN PROGRESS**

This proposal describes how to introduce labels for the atoms in an
:class:`~ase.Atoms` object in a backwards compatible way.


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

Introduce *labels* and *kinds*.  A label is a string and a kind is the
combination of these three:

* A chemical symbol or an atomic number.  Default: ``None``
  (same as ``'X'`` or ``0``)
* A label.  Default ``None`` (same as ``''``)
* A mass.  Default ``None`` (use standard value)

A *kind* can be represented as a ``tuple`` or a ``str``:

* ``(Z, 'label', mass)``
* ``('symbol', 'label', mass)``
* ``'symbol'``
* ``'symbol:label'``
* ``'label'``

Examples:

* ``(1, 'ghost', 0.0)`` (same as ``'H:ghost'``)
* ``('H', '', None)`` (same as ``'H'``)
* ``('H', '', 2.0)``
* ``'methyl'`` (same as ``(None, 'methyl', None)``)

We add ``symbols``, ``labels`` and ``kinds`` list-of-string like attributes to
the :class:`~ase.Atoms` object as well as new ``labels`` and ``kinds`` keyword
arguments to the Atoms constructor.  Currently, the first argument to the
Atoms constructor (``symbols``) will accept a chemical formula as a string or a
list of chemical symbols or atomic numbers.  We extend this to also accept
kinds.

.. note::

    Not really related to this label issue, but since we are anyway adding new
    attributes, we should also add ``masses``, ``tags``, ``initial_magmoms``,
    ``momenta``, ``initial_charges`` and ``mass``.


Examples
========

>>> a = Atoms(['N', 'C', (None, 'methyl', 9.0)])
>>> a.number_of_species
3
>>> a.positions[:, 0] = [0, 1.2, 2.6]
>>> a.masses[a.labels == 'methyl'] = 10
>>> a.numbers
array([7, 6, 0])
>>> a.symbols  # special list-like object tied to a
Symbols(['N', 'C', 'X'])
>>> a.get_chemical_symbols()  # simple list
['N', 'C', 'X']
>>> a.labels
Labels(['', '', 'methyl'])
>>> a.kinds
Kinds(['N', 'C', ('X', 'methyl', 10.0)])

Here are 50 hydrogen molecules:

>>> h = Atoms('H100', positions=...)
>>> h.labels[::2] = 'deuterium'
>>> h.masses[h.labels == 'deuterium'] = 2.0

or equivalently:

>>> h.kinds[::2] = ('H', 'deuterium', 2.0)

A DFT code could use the kinds to select pseudo-potentials:

>>> n2 = molecule('N2')
>>> n2.labels[0] = 'core-hole'
>>> n2.calc = DFT(xc='LDA',
...               pp={'N:core-hole': '~/mypps/N-LDA.core-hole.upf'})


List-like objects
=================

New ``Labels``, ``Kinds`` and ``Symbols`` list-like objects will
be introduced that can handle all the indexing operations in a storage
efficient way.  A statement like ``a.symbols[0] = 'He'`` must somehow lead to
``a.numbers[0] == 2`` and other magic.


Atom objects?
=============

``Atom.label`` and ``Atom.kind``?


I/O
===

???


Questions
=========

Tags?

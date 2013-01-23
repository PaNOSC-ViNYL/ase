=============================
Calculator interface proposal
=============================

All ASE calculators should behave similarly if there is no good reason
for them not to.  This should make it simpler for both users and developers.

This proposal tries to define how a good ASE calculator should behave.
The goal is to have ASE calculators:

* that share more code
* are more uniform to use
* are better tested

Setting some standards is a good thing, but we should also be careful
not to set too strict rules that could limit each calculator to the
lowest common denominator.

So far, this proposal is mostly ideas and questions.  


Behavior
========

When a calculator calculates the energy, forces, stress tensor, atomic
magnetic moments or something else, it should store a copy of the
system (atomic numbers, atomic positions, unit cell and boundary
conditions).  When asked again, it should return the value already
calculated if the system hasn't been changed.

If calculational parameters such as plane wave cutoff or XC-functional
has been changed the calculator should throw away old calculated
values.


Standards parameters
====================

The standard keywords are: ``xc``, ``kpts``, ``smearing_method``,
``width``, ``nbands`` and ``ranks``.  Each calculator can have
additional parameters.  The units are eV and Å.

Magnetic moments and charges are read from the
:class:`~ase.atoms.Atoms` object.


kpts
----

* ``None``: Gamma-point

* ``(n1,n2,n3)``: Monkhorst-Pack grid

* ``(n1,n2,n3,'gamma')``: Shifted Monkhorst-Pack grid that includes `\Gamma`

* ``[(k11,k12,k13),(k21,k22,k23),...]``: Explicit list in units of the
  reciprocal lattice vectors

* Should there be a k-point density parameter?  Maybe ``kpts=3.5`` as
  in 3.5 `\vec k`-points per Å\ `^{-1}`?

* Should there be a way to get a reasonable number of k-points for the
  given cell?


xc
---

* Should all calculators agree on a default functional?  LDA?

* Which LDA?

* What about wave function methods where XC-functional doesn't make
  sense?  Hartree-Fock, MP2, ...


smearing_method
---------------

Smearing method can be one of

* ``'Fermi-Dirac'``
* ``'Gaussian'``
* ``'Methfessel-Paxton'`` of some order?

Default value?


width
-----

The width parameter used for the chosen smearing method (in eV).
Default value?  Is there a better way to specify the smearing method
than using ``smearing_method`` and ``width`` keywords?


nbands
------

Each band can be occupied by two electrons.  Should it be a
requirement that this parameters is optional?

Alternative names: ``states``, ``nband``, ...


ranks
-----

If a calculator can run in parallel on a subset of the available
CPU's, it should support a ``ranks`` keyword that can be set to a list
of ranks.  Useful for parallel NEB calculations.


Restart
=======

A calculator should be able to prefix all output files with a given
label or run the calculation in a directory.  There are three
possibilities for the first argument of the constructor of a
calculator object:

* name of a file containing all results of a calculation
* a prefix used for several files containing results
* name of a directory containing result files with fixed names

Each calculator can decide what the default value is: ``None`` for no
output, ``'-'`` for standard output or something else.  All others
parameters are given as keyword arguments.

Example:  Do a calculation with ABC calculator and write results to
:file:`si.abc`:

>>> atoms = ...
>>> atoms.calc = ABC('si.abc', xc='LDA')
>>> atoms.get_potential_energy()
-1.2

Start from previous calculation:

>>> atoms = ABC.read('si.abc')  # or read_atoms?
>>> atoms.calc
<ABC-calculator>
>>> atoms.get_potential_energy()
-1.2

The :meth:`read()` method is equivalent to:

>>> atoms = ABC('si.abc').get_atoms()

If we do:

>>> atoms = ABC.read('si.abc')
>>> atoms.rattle()            # change positions and/or
>>> atoms.calc.set(xc='PBE')  # change a calculator-parameter
>>> atoms.get_potential_energy()
-0.7

then the :file:`si.abc` will be overwritten or maybe appended to.


Implementation
==============

There will be a hierarchy of classes:

* Common base class for all calculators: ``Calculator``.  Takes care
  of file read/write logic, handles setting of parameters and checks
  for state changes.

* Specialized DFT class: ``ElectronicStructureCalculator`` or
  ``DFTCalculator``?  Special treatment of ``xc``, ``nbands``,
  ``smearing_method`` and ``width`` parameters.

* Class for handling `\vec k`-points: ``KPointCalculator``. Handles
  ``kpts`` logic.


Other stuff
===========

Support for ASE's :ref:`command line tool`.

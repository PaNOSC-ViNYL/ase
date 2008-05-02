============
Introduction
============

ASE is an Atomistic Simulation Environment written in the
Python_ programming language with the aim of setting up, stearing, and
analyzing atomistic simulations. The ASE has been constructed with a
number of "design goals" as:


:Simplicity in use:
  Setting up an atomistic total energy calculation or molecular
  dynamics simulation with ASE is simple and straightforward. The Python
  scripts are almost self-explanatory and look like well
  commented input files to an atomistic simulation program.

:Flexibility in use:
  Since ASE is based on the Python scripting language it is possible
  without any code modifications to perform very complicated simulation
  tasks. For example a sequence of calculations may be performed with
  the use of simple "for-loop" constructions or simulations of different
  types (:term:`DFT` and classical molecular mechanics potentials) may
  be coupled together.

:Simplicity in development:
  The ASE defines a set of interfaces for different objects, i.e. an
  :class:`~ase.atoms.Atoms` object is required to posses a method with the name
  :meth:`get_positions` which returns the coordinates of
  the atom. By following a few such standard interfaces it is easy for
  new users to get access to all of the functionality of ASE.

:Flexibility in development:
  The Python code in ASE is structured in different modules intended for
  different purposes. There are :mod:`calculators` for calculating
  energies, forces and stresses, :mod:`md` and :mod:`optimize` modules
  for controlling the motion of atoms, :mod:`constraint <constraints>`
  objects and filters for performing :mod:`nudged-elastic-band <neb>`
  calculations etc. The modularity of the code and the documented
  interfaces make it simple to contribute new functionality to ASE.

:Pythonic:
  It fits nicely into the rest of the Python world, and *it fist your brain*.

:Open to participation:
  All of the code in ASE is carrying the :term:`GNU General Public License`
  and people are invited to participate in using and :ref:`developing the
  code <devel>`.


.. _Python: http://www.python.org




.. _ml:

Mailing list
============

XXX  campos, ase, checkins, ???

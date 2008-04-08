Introduction
============

Campos ASE is an Atomistic Simulation Environment written in the
python programming language with the aim of setting up, stearing, and
analyzing atomistic simulations. The ASE has been constructed with a
number of "design goals" as:


Simplicity in use.
  Setting up an atomistic total energy calculation or molecular
  dynamics simulation with ASE is simple and straightforward. The python
  scripts are almost self-explanatory and look like well
  commented input files to an atomistic simulation program.

Flexibility in use.
  Since ASE is based on the python scripting language it is possible
  without any code modifications to perform very complicated
  simulation tasks. For example a sequence of calculations may be performed with
  the use of simple "for-loop" constructions or simulations of different
  types (DFT and classical molecular dynamics) may be coupled
  together.

Simplicity in development
  The ASE defines a set of interfaces for different objects, i.e. an
  "atom" object is required to posses a method with the name
  "GetCartesianPosition" which returns the Cartesian coordinates of
  the atom. By following a few such standard interfaces it is easy for
  new users to get access to all of the functionality of ASE.

Flexibility in development
  The python code in ASE is structured in different modules intended
  for different purposes. There are "DFT-calculators" for performing
  DFT calculations, "dynamics" for controlling the motion of atoms,
  "filters" for constraining motion or performing nudged-elastic-band
  calculations etc. The modularity of the code and the documented
  interfaces make it simple to contribute new functionality to ASE.

Open to participation
  All of the code in ASE is carrying the GNU General Public License
  and people are invited to participate in using and developing the
  code.

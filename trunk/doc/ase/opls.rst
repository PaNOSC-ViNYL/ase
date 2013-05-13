==========================================
Setting up an OPLS force field calculation
==========================================

.. module:: opls
   :synopsis: OPLS force field

In order to facilitate the definition of structures for the use
of OPLS force fields, there are some helper classes.

Modified xyz
============

Suppose, we define the ethanal molecule as an modified xyz file (``172_mod.xyz``)::

  7
  # ethanal from ChemSpider
  O       1.613900000000000     -0.762100000000000     -0.000000000000000
  CT      -0.327900000000000      0.522700000000000      0.000000000000000
  C       0.392400000000000     -0.722900000000000     -0.000000000000000
  HC      -0.960000000000000      0.580900000000000      0.887500000000000
  HC      -0.960000000000000      0.580900000000000     -0.887500000000000
  HC       0.346400000000000      1.382000000000000      0.000000000000000
  H1      -0.104900000000000     -1.581400000000000     -0.000000000000000

Then we can read and view the structure using::

  from ase.visualize import view
  from ase.io.opls import OPLSStructure

  s = OPLSStructure('172_mod.xyz') # 172_mod.xyz if the file name for the structure above
  view(s) # view with real elements
  elements = { 'CT' : 'Si', 'HC' : 'H', 'H1' : 'He' }
  view(s.colored(elements)) # view with fake elements

Defining the force field
========================

The definitions of the force field can be stored in an Amber like style (``172_defs.par``)::

  # the blocks are separated by empty lines
  # comments are allowed 
  #
  # one body - LJ-parameters and charge
  CT 0.0028619844 3.50  0.000
  C  0.0028619844 3.50  0.000
  O  0.0073717780 3.12 -0.683
  HC 0.0013009018 2.50  0.000
  H1 0 0 0

  # bonds
  C -CT  317.0    1.522       JCC,7,(1986),230; AA
  C -O   570.0    1.229       JCC,7,(1986),230; AA,CYT,GUA,THY,URA
  CT-HC  340.0    1.090       changed from 331 bsd on NMA nmodes; AA, SUGARS
  C -H1  340.0    1.090       

  # angles
  CT-C -O     80.0      120.40
  HC-CT-HC    35.0      109.50
  CT-C -H1    70.0      117.00

  # dihedrals - here empty

  # cutoffs
  C -CT 2.0
  C -O  1.8
  CT-HC 1.4
  C -H1 1.4
  C3-O1 1.8 # extra stuff, should not bother

We can write LAMMPS input using the information above::

  from ase.io.opls import OPLSff, OPLSStructure

  s = OPLSStructure('172_mod.xyz')
  opls = OPLSff('172_defs.par')
  opls.write_lammps(s)

which writes the LAMMPS input files ``lammps_atoms`` defining atoms, bonds, etc., and
``lammps_opls`` defining the corresponding OPLS force field. 

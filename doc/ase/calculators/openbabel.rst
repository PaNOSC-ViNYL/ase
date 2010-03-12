.. module:: openbabel

=========
OpenBabel
=========

Introduction
============

OpenBabel_ is a chemical toolbox designed to speak the many languages of
chemical data. It's an open, collaborative project allowing anyone to search,
convert, analyze, or store data from molecular modeling, chemistry,
solid-state materials, biochemistry, or related areas.  

The :class:`~ase.calculators.openbabel.OBForceField` calculator use the force
fields (currently UFF and ghemical) included in OpenBabel.

To use this calculator you need to have the OpenBabel python bindings
installed:

- Ubuntu/Fedora: python-openbabel

.. _OpenBabel: http://www.openbabel.org

OBForceField Calculator Class
=============================

.. class:: OBForceField(force_field='UFF', bonds=None)

Here is a detailed list of all the keywords for the calculator:

================ ========= ================  =================================================
keyword          type      default value     description
================ ========= ================  =================================================
``force_field``  ``str``   ``'UFF'``         Force field ('UFF' or 'ghemical')
``bonds``        ``list``  ``None``          List of bonds in molecule given as:
                                             [[begin atom idx, end atom idx, bond order], ...]
                                             If None OpenBabel will try to guess it. 
================ ========= ================  =================================================

Examples
========

Automatic bond detection 
------------------------

Here is an example of how to calculate the total energy CO::
        
  #!/usr/bin/env python
  from ase import molecule
  from ase.calculators.openbabel import OBForceField
  
  atoms = molecule('CO')

  calc = OBForceField()
  atoms.set_calculator(calc)
  e = atoms.get_potential_energy()

Adding bonds manually
---------------------

If we want to relax e.g. CO2 starting with a very large interatomic distances,
OpenBabel will not detect the bonds and we have to put them manually::

  #!/usr/bin/env python
  from ase import Atoms, QuasiNewton
  from ase.calculators.openbabel import OBForceField

  atoms = Atoms('OCO', [[0.0, 0.0, 0.0],
                        [3.0, 0.0, 0.0],
                        [6.0, 0.0, 0.0]])

  calc = OBForceField()
  atoms.set_calculator(calc)
  relax = QuasiNewton(atoms)
  relax.run(fmax=0.05)
  # This results in a even larger interatomic distances
  # Therefore we add two double bonds

  bonds = [[0, 1, 2],
           [1, 2, 2]]
  calc = OBForceField(bonds=bonds)
  atoms.set_calculator(calc)
  relax = QuasiNewton(atoms)
  relax.run(fmax=0.05)
  # This gives the expected results

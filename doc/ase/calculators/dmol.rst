.. module:: ase.calculators.dmol

=====
DMol3
=====

DMol3 is an atomic orbital DFT code.

Environment variables
=====================
DMOL_COMMAND should point to the RunDmol script

.. highlight:: bash
 
::

  $ export DMOL_COMMAND="./RunDmol.sh -np 16"

DMol3 Calculator
================
The DMol3 calculator is a FileIOCalculator. The default setting used by the
dmol interface is

.. class:: DMol3(functional='pbe', symmetry='on')

The dmol calculator supports the calculate gradient function in DMol3, meaning
the internal relaxation is not supported. Forces and potential energy are the
supported properties

.. code-block:: python

    implemented_properties = ['energy', 'forces']


.. note::

   DMol3 often reorients the atomic system. Therefore it's recommended to use
   the calculator with care. Forces are reoriented to match the atoms object,
   however properties like k-points and density files (.grd) may be misoriented
   when reading.

.. note::

   Only 3D periodic systems (pbc = [True, True, True]) and fully non-periodic
   systems are supported by the DMol3 calculator.

Example
=======

.. code-block:: python

   from ase.build import molecule
   from ase.calculators.dmol import DMol3

   atoms = molecule('H2O')
   calc = DMol3(symmetry='auto',
                spin_polarization='unrestricted',
                charge=0,
                basis='dnp',
                pseudopotential='none',
                functional='pbe',
                scf_density_convergence=1.0e-7)
   atoms.set_calculator(calc)
   atoms.get_potential_energy()

File formats
============

The supported dmol file formats (for write/read) are

* .car
* .incoor
* .arc

For molecules and systems without periodic boundary conditions, the .car
format is used, while for periodic systems the .incoor format, which allows
specification of the unit cell, is used. The .arc files are trajectory files
from internal relaxation runs in DMol3 (which is not supported by this
calculator)



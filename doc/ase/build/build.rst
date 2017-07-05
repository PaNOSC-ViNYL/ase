.. module:: ase.build

================
Building things
================

Quick links:

* Simple bulk crystals: :func:`~ase.build.bulk`

* Simple molecules: :func:`~ase.build.molecule`

* Special surfaces:

  * fcc: :func:`~ase.build.fcc100`, :func:`~ase.build.fcc110`,
    :func:`~ase.build.fcc111`, :func:`~ase.build.fcc211`,
    :func:`~ase.build.fcc111_root`

  * bcc: :func:`~ase.build.bcc100`, :func:`~ase.build.bcc110`,
    :func:`~ase.build.bcc111`
    * - :func:`~ase.build.bcc111_root`

  * hcp: :func:`~ase.build.hcp0001`, :func:`~ase.build.hcp10m10`,
    :func:`~ase.build.hcp0001_root`

  * diamond: :func:`~ase.build.diamond100`, :func:`~ase.build.diamond111`

* `MX_2` (2H or 1T): :func:`~ase.build.mx2`

* Other surface tools: :func:`~ase.build.surface`,
  :func:`~ase.build.add_adsorbate`, :func:`~ase.build.add_vacuum`,
  :func:`~ase.build.root_surface`

* 1D: :func:`~ase.build.nanotube`, :func:`~ase.build.graphene_nanoribbon`

* Other tools: :func:`~ase.build.cut`, :func:`~ase.build.stack`,
  :func:`~ase.build.sort`, :func:`~ase.build.minimize_tilt`,
  :func:`~ase.build.niggli_reduce`, :func:`~ase.build.rotate`,
  :func:`~ase.build.minimize_rotation_and_translation`,
  :func:`~ase.build.get_deviation_from_optimal_cell_shape`,
  :func:`~ase.build.find_optimal_cell_shape`,
  :func:`~ase.build.find_optimal_cell_shape_pure_python`,
  :func:`~ase.build.make_supercell`



.. toctree::
   :maxdepth: 2

   surface
   tools
   
.. seealso::

   * The :mod:`ase.lattice` module.  The module contains functions for
     creating most common crystal structures with arbitrary orientation.
     The user can specify the desired Miller index along the three axes
     of the simulation, and the smallest periodic structure fulfilling
     this specification is created.  Both bulk crystals and surfaces can
     be created.

   * The :mod:`ase.cluster` module.  Useful for creating nanoparticles
     and clusters.
     
   * The :mod:`ase.spacegroup` module

   * The :mod:`ase.geometry` module



Molecules
=========

The G2-database of common molecules is available:

.. autofunction:: molecule

Example::

>>> from ase.build import molecule
>>> atoms = molecule('H2O')

The list of available molecules is those from the :data:`ase.collections.g2`
database:

>>> from ase.collections import g2
>>> g2.names
['PH3', 'P2', 'CH3CHO', 'H2COH', 'CS', 'OCHCHO', 'C3H9C', 'CH3COF',
 'CH3CH2OCH3', 'HCOOH', 'HCCl3', 'HOCl', 'H2', 'SH2', 'C2H2',
 'C4H4NH', 'CH3SCH3', 'SiH2_s3B1d', 'CH3SH', 'CH3CO', 'CO', 'ClF3',
 'SiH4', 'C2H6CHOH', 'CH2NHCH2', 'isobutene', 'HCO', 'bicyclobutane',
 'LiF', 'Si', 'C2H6', 'CN', 'ClNO', 'S', 'SiF4', 'H3CNH2',
 'methylenecyclopropane', 'CH3CH2OH', 'F', 'NaCl', 'CH3Cl',
 'CH3SiH3', 'AlF3', 'C2H3', 'ClF', 'PF3', 'PH2', 'CH3CN',
 'cyclobutene', 'CH3ONO', 'SiH3', 'C3H6_D3h', 'CO2', 'NO',
 'trans-butane', 'H2CCHCl', 'LiH', 'NH2', 'CH', 'CH2OCH2',
 'C6H6', 'CH3CONH2', 'cyclobutane', 'H2CCHCN', 'butadiene', 'C',
 'H2CO', 'CH3COOH', 'HCF3', 'CH3S', 'CS2', 'SiH2_s1A1d', 'C4H4S',
 'N2H4', 'OH', 'CH3OCH3', 'C5H5N', 'H2O', 'HCl', 'CH2_s1A1d',
 'CH3CH2SH', 'CH3NO2', 'Cl', 'Be', 'BCl3', 'C4H4O', 'Al', 'CH3O',
 'CH3OH', 'C3H7Cl', 'isobutane', 'Na', 'CCl4', 'CH3CH2O', 'H2CCHF',
 'C3H7', 'CH3', 'O3', 'P', 'C2H4', 'NCCN', 'S2', 'AlCl3', 'SiCl4',
 'SiO', 'C3H4_D2d', 'H', 'COF2', '2-butyne', 'C2H5', 'BF3', 'N2O',
 'F2O', 'SO2', 'H2CCl2', 'CF3CN', 'HCN', 'C2H6NH', 'OCS', 'B', 'ClO',
 'C3H8', 'HF', 'O2', 'SO', 'NH', 'C2F4', 'NF3', 'CH2_s3B1d', 'CH3CH2Cl',
 'CH3COCl', 'NH3', 'C3H9N', 'CF4', 'C3H6_Cs', 'Si2H6', 'HCOOCH3', 'O',
 'CCH', 'N', 'Si2', 'C2H6SO', 'C5H8', 'H2CF2', 'Li2', 'CH2SCH2', 'C2Cl4',
 'C3H4_C3v', 'CH3COCH3', 'F2', 'CH4', 'SH', 'H2CCO', 'CH3CH2NH2', 'Li',
 'N2', 'Cl2', 'H2O2', 'Na2', 'BeH', 'C3H4_C2v', 'NO2']

plus ``Be2``, ``C7NH5``, ``BDA``, ``biphenyl`` and ``C60`` (for historical
reasons).


.. _bulk-crystal-section:

Common bulk crystals
====================

.. autofunction:: bulk

examples:

>>> from ase.build import bulk
>>> a1 = bulk('Cu', 'fcc', a=3.6)
>>> a2 = bulk('Cu', 'fcc', a=3.6, orthorhombic=True)
>>> a3 = bulk('Cu', 'fcc', a=3.6, cubic=True)
>>> a1.cell
array([[ 0. ,  1.8,  1.8],
       [ 1.8,  0. ,  1.8],
       [ 1.8,  1.8,  0. ]])
>>> a2.cell
array([[ 2.546,  0.   ,  0.   ],
       [ 0.   ,  2.546,  0.   ],
       [ 0.   ,  0.   ,  3.6  ]])
>>> a3.cell
array([[ 3.6,  0. ,  0. ],
       [ 0. ,  3.6,  0. ],
       [ 0. ,  0. ,  3.6]])

|a1| |a2| |a3|

.. |a1| image:: a1.png
.. |a2| image:: a2.png
.. |a3| image:: a3.png


.. _nanotubes-section:

Nanotubes
=========

.. autofunction:: nanotube

examples:

>>> from ase.build import nanotube
>>> cnt1 = nanotube(6, 0, length=4)
>>> cnt2 = nanotube(3, 3, length=6, bond=1.4, symbol='Si')

|cnt1| |cnt2|

.. |cnt1| image:: cnt1.png
.. |cnt2| image:: cnt2.png


.. _nanoribbons-section:

Graphene nanoribbons
====================

.. autofunction:: graphene_nanoribbon

examples:

>>> from ase.build import graphene_nanoribbon
>>> gnr1 = graphene_nanoribbon(3, 4, type='armchair', saturated=True)
>>> gnr2 = graphene_nanoribbon(2, 6, type='zigzag', saturated=True,
...                            C_H=1.1, C_C=1.4, vacuum=6.0,
...                            magnetic=True, initial_mag=1.12)

|gnr1| |gnr2|

.. |gnr1| image:: gnr1.png
.. |gnr2| image:: gnr2.png

ASE contains a number of modules for setting up atomic structures,
mainly molecules, bulk crystals and surfaces.  Some of these modules
have overlapping functionality, but strike a different balance between
flexibility and ease-of-use.




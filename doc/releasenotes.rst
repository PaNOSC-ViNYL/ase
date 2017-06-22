.. _releasenotes:

=============
Release notes
=============


Git master branch
=================

:git:`master <>`.


Version 3.14.0
==============

20 June 2017: :git:`3.14.0 <../3.14.0>`.

* Python 2.6 no longer supported.

* The command-line tools :program:`ase-???` have been replaced by a
  single :program:`ase` command with sub-commands (see :ref:`cli`).
  For help, type::

      $ ase --help
      $ ase sub-command --help

* The old :program:`ase-build` command which is now called
  :program:`ase build` will no longer add vacuum by default.  Use
  ``ase build -V 3.0`` to get the old behavior.

* All methods of the :class:`~ase.Atoms` object that deal with angles now
  have new API's that use degrees instead of radians as the unit of angle
  (:meth:`~ase.Atoms.get_angle`, :meth:`~ase.Atoms.set_angle`,
  :meth:`~ase.Atoms.get_dihedral`, :meth:`~ase.Atoms.set_dihedral`,
  :meth:`~ase.Atoms.rotate_dihedral`, :meth:`~ase.Atoms.rotate`,
  :meth:`~ase.Atoms.euler_rotate`).

  The old way of calling these methods works as always, but will give
  you a warning.  Example:

  >>> water.get_angle(0, 1, 2)  # new API
  104.52
  >>> water.get_angle([0, 1, 2])  # old API
  /home/jensj/ase/ase/atoms.py:1484: UserWarning: Please use new API (which will return the angle in degrees): atoms_obj.get_angle(a1,a2,a3)*pi/180 instead of atoms_obj.get_angle([a1,a2,a3])
  1.8242181341844732

  Here are the changes you need to make in order to get rid of warnings:

  Old API:

  >>> a1 = atoms.get_angle([0, 1, 2])
  >>> atoms.set_angle([0, 1, 2], pi / 2)
  >>> a2 = atoms.get_dihedral([0, 1, 2, 3])
  >>> atoms.set_dihedral([0, 1, 2, 3], pi / 6)
  >>> atoms.rotate_dihedral([0, 1, 2, 3], 10.5 * pi / 180)
  >>> atoms.rotate('z', pi / 4)
  >>> atoms.rotate_euler(phi=phi, theta=theta, psi=psi)

  New API:

  >>> a1 = atoms.get_angle(0, 1, 2) * pi / 180
  >>> atoms.set_angle(0, 1, 2, angle=90)
  >>> a2 = atoms.get_dihedral(0, 1, 2, 3) * pi / 180
  >>> atoms.set_dihedral(0, 1, 2, 3, angle=30)
  >>> atoms.rotate_dihedral(0, 1, 2, 3, angle=10.5)
  >>> atoms.rotate(45, 'z')
  >>> atoms.euler_rotate(phi=phi * 180 / pi,
  ...                    theta=theta * 180 / pi,
  ...                    psi=psi * 180 / pi)

* The web-interface to the :mod:`ase.db` module now uses Bootstrap and looks
  much nicer.  Querying the database is also much easier.  See
  https://cmrdb.fysik.dtu.dk for an example.

* The PostgreSQL backend for :mod:`ase.db` can now contain more than one ASE
  database.

* An ASE database can now have :ref:`metadata` describing the data.
  Metadata is a dict with any of the following keys: ``title``,
  ``key_descriptions``, ``default_columns``, ``special_keys`` and
  ``layout``.

* :data:`ase.data.atomic_masses` has been updated to IUPAC values from
  2016. Several elements will now have different weights which will affect
  dynamic calculations. The old values can be recovered like this:

  >>> from ase.data import atomic_masses_legacy
  >>> atoms.set_masses(atomic_masses_legacy[atoms.numbers])

* New :func:`ase.data.isotopes.download_isotope_data` function for getting
  individual isotope masses from NIST.

* New :func:`ase.eos.calculate_eos` helper function added.

* Added DeltaCodesDFT data: :data:`ase.collections.dcdft`.

* :mod:`ase.gui` can now load and display any sequence of :class:`~ase.Atoms`
  objects; it is no longer restricted to sequences with a constant number
  of atoms or same chemical composition.

* Trajectory files can now store any sequence of :class:`~ase.Atoms`
  objects.  Previously, atomic numbers, masses, and constraints were
  only saved for the first image, and had to apply for all subsequent ones.

* Added calculator interface for DMol\ :sup:`3`.

* Added calculator interface for GULP.

* Added file formats .car, .incoor, and .arc, related to DMol\ :sup:`3`.

* New function for interpolating from Monkhors-Pack sampled values in the BZ
  to arbitrary points in the BZ:
  :func:`ase.dft.kpoints.monkhorst_pack_interpolate`.

* New *band-structure* command for the :program:`ase` :ref:`cli`.

* Two new functions for producing chemical formulas:
  :func:`ase.utils.formula_hill` and :func:`ase.utils.formula_metal`.

* The :func:`ase.dft.bandgap.get_band_gap` function is now deprecated.  Use
  the new one called :func:`ase.dft.bandgap.bandgap` (it's more flexible and
  returns also band indices).

* New :mod:`Viewer for Jupyter notebooks <ase.visualize.nglview>`.


Version 3.13.0
==============

7 February 2017: :git:`3.13.0 <../3.13.0>`.

* The default unit-cell when you create an :class:`~ase.Atoms` object has
  been changed from ``[[1,0,0],[0,1,0],[0,0,1]]`` to
  ``[[0,0,0],[0,0,0],[0,0,0]]``.

* New :attr:`ase.Atoms.number_of_lattice_vectors` attribute equal to,
  big surprise, the number of non-zero lattice vectors.

* The :meth:`ase.Atoms.get_cell` method has a new keyword argument
  ``complete``.  Use ``atoms.get_cell(complete=True)`` to get a complete
  unit cell with missing lattice vectors added at right angles to the
  existing ones.  There is also a function :func:`ase.geometry.complete_cell`
  that will complete a unit cell.

* :func:`~ase.build.graphene_nanoribbon` no longer adds 2.5 Å of vacuum by
  default.

* All functions that create molecules, chains or surfaces
  (see the :mod:`ase.build` module) will no longer add "dummy" lattice
  vectors along the non-periodic directions.  As an example, the surface
  functions will generate unit cells of the type
  ``[[a1,a2,0],[b1,b2,0],[0,0,0]]``.  In order to define all three lattice
  vectors, use the ``vacuum`` keyword that all
  of the 0-d, 1-d and 2-d functions have or, equivalently, call the
  :meth:`~ase.Atoms.center` method.

* Many of the :ref:`surface generating functions <surfaces>` have changed
  their behavior when called with ``vacuum=None`` (the default).  Before, a
  vacuum layer equal to the interlayer spacing would be added on the upper
  surface of the slab. Now, the third axis perpendicular to the surface will be
  undefined (``[0, 0, 0]``).  Use ``vacuum=<half-the-interlater-distance>`` to
  get something similar to the old behavior.

* New :func:`ase.geometry.is_orthorhombic` and
  :func:`ase.geometry.orthorhombic` functions added.

* :mod:`ase.gui` now works on Python 3.

* NEB-tools class has been renamed to :class:`~ase.neb.NEBTools`.

* :mod:`Optimizers <ase.optimize>` now try force-consistent energies if
  possible (instead of energies extrapolated to 0.0 K).


Version 3.12.0
==============

24 October 2016: :git:`3.12.0 <../3.12.0>`.

* New :class:`ase.constraints.ExternalForce` constraint.

* Updated :mod:`ase.units` definition to CODATA 2014. Additionally, support
  for older versions of CODATA was added such that the respective units can
  be created by the user when needed (e.g. interfacing codes with different
  CODATA versions in use).

* New :mod:`ase.calculators.checkpoint` module.  Adds restart and rollback
  capabilities to ASE scripts.

* Two new flawors of :class:`~ase.neb.NEB` calculations have been added:
  ``method='eb'`` and ``method='improvedtangent'``.

* :func:`ase.io.write` can now write XSD files.

* Interfaces for deMon, Amber and ONETEP added.

* New :ref:`defects` tutorial and new super-cell functions:
  :func:`~ase.build.get_deviation_from_optimal_cell_shape`,
  :func:`~ase.build.find_optimal_cell_shape`,
  :func:`~ase.build.find_optimal_cell_shape_pure_python`,
  :func:`~ase.build.make_supercell`.

* New :class:`~ase.dft.band_structure.BandStructure` object.  Can identify
  special points and create nice plots.

* Calculators that inherit from :class:`ase.calculators.calculator.Calculator`
  will now have a :meth:`~ase.calculators.calculator.Calculator.band_structure`
  method that creates a :class:`~ase.dft.band_structure.BandStructure` object.

* Addition to :mod:`~ase.geometry` module:
  :func:`~ase.geometry.crystal_structure_from_cell`.

* New functions in :mod:`ase.dft.kpoints` module:
  :func:`~ase.dft.kpoints.parse_path_string`,
  :func:`~ase.dft.kpoints.labels_from_kpts` and
  :func:`~ase.dft.kpoints.bandpath`.

* Helper function for generation of Monkhorst-Pack samplings and BZ-paths:
  :func:`ase.calculators.calculator.kpts2ndarray`.

* Useful class for testing band-structure stuff:
  :class:`ase.calculators.test.FreeElectrons`.

* The ``cell`` attribute of an :class:`~ase.Atoms` object and the ``cell``
  keyword for the :class:`~ase.Atoms` constructor and the
  :meth:`~ase.Atoms.set_cell` method now accepts unit cells given ase
  ``[a, b, c, alpha, beta, gamma]``, where the three angles are in degrees.
  There is also a corresponding :meth:`~ase.Atoms.get_cell_lengths_and_angles`
  method.

* Galician translation of ASE's GUI.

* Two new preconditioned structure optimizers available.  See
  :mod:`ase.optimize.precon`.

* Trajectory files now contain information about the calculator and also
  information from an optimizer that wrote the trajectory.


Version 3.11.0
==============

10 May 2016: :git:`3.11.0 <../3.11.0>`.

* Special `\mathbf{k}`-points from the [Setyawana-Curtarolo]_ paper was added:
  :data:`ase.dft.kpoints.special_points`.

* New :mod:`ase.collections` module added.  Currently contains the G2 database
  of molecules and the S22 set of weakly interacting dimers and complexes.

* Moved modules:

  * ``ase.utils.eos`` moved to :mod:`ase.eos`
  * ``ase.calculators.neighborlist`` moved to :mod:`ase.neighborlist`
  * ``ase.lattice.spacegroup`` moved to :mod:`ase.spacegroup`

* The ``InfraRed`` that used to be in the ``ase.infrared`` or
  ``ase.vibrations.infrared`` modules is now called
  :class:`~ase.vibrations.Infrared` and should be imported from the
  :mod:`ase.vibrations` module.

* Deprecated modules: ``ase.structure``, ``ase.utils.geometry``,
  ``ase.utils.distance``, ``ase.lattice.surface``.  The functions from these
  modules that will create and manipulate :class:`~ase.Atoms` objects are now
  in the new :mod:`ase.build` module.  The remaining functions have been moved
  to the new :mod:`ase.geometry` module.

* The ``ase.lattice.bulk()`` function has been moved to :func:`ase.build.bulk`.

* Two new functions: :func:`~ase.geometry.cell_to_cellpar` and
  :func:`~ase.geometry.cellpar_to_cell`.

* We can now :func:`~ase.io.read` and :func:`~ase.io.write` magres files.

* :class:`~ase.neb.NEB` improvement:  calculations for molecules can now be
  told to minimize ratation and translation along the path.


Version 3.10.0
==============

17 Mar 2016: :git:`3.10.0 <../3.10.0>`.

* :ref:`old trajectory` files can no longer be used.  See :ref:`convert`.

* New iterator function :func:`ase.io.iread` for iteratively reading Atoms
  objects from a file.

* The :func:`ase.io.read` function and command-line tools can now read ``.gz``
  and ``.bz2`` compressed files.

* Two new decorators :func:`~ase.parallel.parallel_function` and
  :func:`~ase.parallel.parallel_generator` added.

* Source code moved to https://gitlab.com/ase/ase.

* Preliminary :mod:`ase.calculators.qmmm` module.

* Improved :mod:`~ase.calculators.tip3p.TIP3P` potential.

* Velocity Verlet will now work correctly with constraints.

* ASE's GUI no longer needs a special GTK-backend for matplotlib to work.
  This will make installation of ASE much simpler.

* We can now :func:`~ase.io.read` and :func:`~ase.io.write` JSV files.

* New :func:`ase.dft.kpoints.get_special_points` function.

* New :func:`ase.geometry.get_duplicate_atoms` function for finding and
  removing atoms on top of each other.

* New: A replacement :mod:`Siesta <ase.calculators.siesta>` calculator was
  implemented. It closely follows the
  :class:`ase.calculators.calculator.FileIOCalculator` class which should
  ease further development. Handling pseudopotentials, basis sets and ghost
  atoms have been made much more flexible in the new version.


Version 3.9.1
=============

21 July 2015: :git:`3.9.1 <../3.9.1>`.

* Added function for finding maximally-reduced Niggli unit cell:
  :func:`ase.build.niggli_reduce`.

* Octopus interface added (experimental).


Version 3.9.0
=============

28 May 2015: :git:`3.9.0 <../3.9.0>`.

* Genetic algorithm implemented; :mod:`ase.ga`. This can be used
  for the optimization of: atomic cluster structure, materials
  properties by use of template structures. Extension to other projects
  related to atomic simulations should be straightforward.

* The ``ase.lattice.bulk`` function can now build the Wurtzite structure.

* The :class:`ase.utils.timing.Timer` was moved from GPAW to ASE.

* New :mod:`ase.db` module.

* New functions: :func:`ase.build.fcc211` and
  :func:`ase.visualize.mlab.plot`.

* New :class:`~ase.Atoms` methods:
  :meth:`ase.Atoms.get_distances()` and
  :meth:`ase.Atoms.get_all_distances()`.

* :ref:`bash completion` can now be enabled.

* Preliminary support for Python 3.

* Wrapping: new :meth:`ase.Atoms.wrap` method and
  :func:`ase.geometry.wrap_positions` function.  Also
  added ``wrap=True`` keyword argument to
  :meth:`ase.Atoms.get_scaled_positions` that can be used to turn
  off wrapping.

* New improved method for initializing NEB calculations:
  :meth:`ase.neb.NEB.interpolate`.

* New pickle-free future-proof trajectory file format added:
  :ref:`new trajectory`.

* We can now do :ref:`phase diagrams`.

* New :func:`ase.build.mx2` function for 1T and 2H metal
  dichalcogenides and friends.

* New :func:`ase.dft.bandgap.get_band_gap` function

* :class:`~ase.calculators.cp2k.CP2K` interface.


Version 3.8.0
=============

22 October 2013: :git:`3.8.0 <../3.8.0>`.

* ASE's :mod:`gui <ase.gui>` renamed from ``ag`` to ``ase-gui``.
* New :ref:`STM <stm>` module.
* Python 2.6 is now a requirement.
* The old ``ase.build.bulk`` function is now deprecated.
  Use the new one instead (:func:`ase.lattice.bulk`).
* We're now using BuildBot for continuous integration:
  https://ase-buildbot.fysik.dtu.dk/waterfall
* New interface to the JDFTx code.


Version 3.7.0
=============

13 May 2013: :git:`3.7.0 <../3.7.0>`.

* ASE's GUI can now be configured to be more friendly to visually
  impaired users: :ref:`high contrast`.

* The :class:`ase.neb.NEB` object now accepts a list of spring constants.

* *Important backwards incompatible change*: The
  :func:`ase.build.surface` function now returns a
  right-handed unit cell.

* Mopac, NWChem and Gaussian interfaces and EAM potential added.

* New :meth:`~ase.Atoms.set_initial_charges` and
  :meth:`~ase.Atoms.get_initial_charges` methods.  The
  :meth:`~ase.Atoms.get_charges` method will now ask the
  calculator to calculate the atomic charges.

* The :ref:`aep1` has been implemented and 6 ASE calculators are now
  based on the new base classes.

* ASE now runs on Windows and Mac.

* :ref:`mhtutorial` added to ASE.


Version 3.6.0
=============

24 Feb 2012: :git:`3.6.0 <../3.6.0>`.

* ASE GUI translations added, available: da_DK, en_GB, es_ES.

* New function for making surfaces with arbitrary Miller indices with
  the smallest possible surface unit cell:
  ase.build.surface()

* New ase.lattice.bulk() function.  Will replace old
  ase.build.bulk() function.  The new one will produce a more
  natural hcp lattice and it will use experimental data for crystal
  structure and lattice constants if not provided explicitely.

* New values for ase.data.covalent_radii from Cordeo *et al.*.

* New command line tool: :ref:`cli` and tests based on it:
  abinit, elk, fleur, nwchem.

* New crystal builder for ase-gui

* Van der Waals radii in ase.data

* ASE's GUI (ase-gui) now supports velocities for both graphs and coloring

* Cleaned up some name-spaces:

  * ``ase`` now contains only :class:`~ase.Atoms` and
    :class:`~ase.atom.Atom`
  * ``ase.calculators`` is now empty


Version 3.5.1
=============

24 May 2011: :git:`3.5.1 <../3.5.1>`.

* Problem with parallel vibration calculations fixed:
  `Ticket #80 <https://trac.fysik.dtu.dk/projects/ase/ticket/80>`_.


Version 3.5.0
=============

13 April 2011: :git:`3.5.0 <../3.5.0>`.

* Improved EMT potential:  uses a
  :class:`~ase.neighborlist.NeighborList` object and is
  now ASAP_ compatible.

* :mod:`BFGSLineSearch <ase.optimize.bfgslinesearch>` is now the default
  (``QuasiNewton==BFGSLineSearch``).

* There is a new interface to the LAMMPS molecular dynamics code.

* New :mod:`ase.phonons` module.

* Van der Waals corrections for DFT, see GPAW_ usage.

* New :class:`~ase.io.bundletrajectory.BundleTrajectory` added.

* Updated GUI interface:

  * Stability and usability improvements.
  * Povray render facility.
  * Updated expert user mode.
  * Enabled customization of colours and atomic radii.
  * Enabled user default settings via :file:`~/.ase/gui.py`.

* :mod:`Database library <ase.data>` expanded to include:

  * The s22, s26 and s22x5 sets of van der Waals bonded dimers and
    complexes by the Hobza group.
  * The DBH24 set of gas-phase reaction barrier heights by the Truhlar
    group.

* Implementation of the Dimer method.


.. _ASAP: http://wiki.fysik.dtu.dk/asap
.. _GPAW: https://wiki.fysik.dtu.dk/gpaw/documentation/xc/vdwcorrection.html


Version 3.4.1
=============

11 August 2010: :git:`3.4.1 <../3.4.1>`.

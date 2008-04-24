Maximally localized Wannier functions
-------------------------------------

.. default-role math



Wannier functions are maximally localized orbitals
constructed from a set of extended eigenstates. This page describes
how to construct the Wannier orbitals using the class **Wannier**.

Band structure and orbital analysis, based on the
set of localized Wannier orbitals are discussed in the chapter
band structure analysis XXX

For the theory behind this method see the paper
Partly Occupied Wannier Functions XXX
Thygesen, Hansen and Jacobsen, Phys. Rev. Lett, Vol.94, 26405 (2005).


The class Wannier
`````````````````

A simple initialization of the :class:`Wannier` looks like this::


     >>> from ASE.Utilities.Wannier import Wannier
     >>> atoms = Calculator.read_atoms('ethylene.nc')
     >>> calc = atoms.get_calculator()
     >>> wannier = Wannier(numberofwannier=6,calculator=calc)


The first argument to the Wannier constructor must be the
number of Wannier functions to be constructed, in this case 6.
The next argument is a calculator containing all the information about
the electronic structure calculation, i.e. number of bands, k-points,
wavefunctions and so on.

For examples of how to use the **Wannier** class, see the `Wannier tutorial`_.

.. _Wannier tutorial: http://www.fysik.dtu.dk/campos/ASE/tut/wannier.html

XXX tip:

.. note::
   For a calculation using **k**-points the relevant unit cell for
   eg. visualization of the Wannier orbitals is not the original unit cell,
   but rather a larger unit cell defined by repeating the original
   unit cell by the number of **k**-points in each direction.
   We will refer to this unit cell as the large unit cell.
   Note that for a :math:`\Gamma`-point calculation the large unit cell
   coinsides with the original unit cell.
   The large unitcell defines also the periodicity of the Wannier
   orbitals.

.. note::
   For calculations using **k**-points, make sure that the
   gamma-point is included in the **k**-point grid. Moreover you must shift all
   **k**-points by a small amount (but not less than 2e-5 in) in e.g. the x direction, before performing
   the Dacapo calculation. If this is not done the symmetry program in
   Dacapo will use time-reversal symmetry to reduce the number of
   **k**-points by a factor 2. The shift can be performed like this::

                kpoints = calc.get_b_z_k_points()
                kpoints[:,0] += 2e-5
                calc.set_b_z_k_points(kpoints)


Below is the full list of keyword arguments:

``numberofwannier``: Number of Wannier orbitals.
  The number of Wannier orbitals to be constructed
  must be <= ``numberofbands``.

``calculator``: Calculator holding the eigenstates, the unit cell and
  providing the method GetWannierLocalizationMatrix.

``numberofbands``: Total number of bands defining the external
  localization space.
  This number is optional. It must be <= the
  total number of bands in the DFT calculation. Default is the total
  number of bands used in the DFT calculation.

``spin``: Select specific spin for a spin-polarized calculation (0 or 1).

``occupationenergy``: Eigenstates below this energy are included in
  the internal localization space. In practice this means that all
  eigenstates below this energy can be exactly reproduced in terms of
  the resulting Wannier functions.
  The energy is relative to the Fermi-level, and 0 is default (all
  occupied states are included in the localization space).

``numberoffixedstates``: Fix the number of states to be included in
  the internal localization space starting from the lowest eigenstate.
  This keyword provides another way of specifying how many
  states should be included, and overrides ``occupationenergy`` if
  this is also set. Default is None meaning that
  the number of fixed states is set by the ``occupationenergy``
  keyword.

``initialwannier``: Setup an initial set of Wannier orbitals.
  *initialwannier* can  set up a  starting guess for the Wannier functions.
  This is important to speed up convergence in particular for large systems
  For transition elements with **d** electrons you will always find 5 highly
  localized **d**-orbitals centered at the atom.
  Placing 5 **d**-like orbitals with a radius of
  0.4 Angstroms and center at atom no. 7, and 3 **p**-like orbitals with a
  radius of 0.4 Angstroms and center at atom no. 27 looks like this::

     initialwannier = [[[7],2,0.4],[[27],1,0.4]]

  Placing only the l=2, m=-2 and m=-1 orbitals at atom no. 7 looks like this::

     initialwannier = [[[7],2,-2,0.4],[[7],2,-1,0.4]]

  I.e. if you do not specify the m quantum number all allowed values are used.
  Instead of placing an orbital at an atom, you can place it at a specified
  position. For example the following::

     initialwannier = [[[0.5,0.5,0.5],0,0.5]]

  places an **s** orbital with radius 0.5 Angstroms at the position (0.5,0.5,0.5)
  in scaled coordinates of the unit cell.

.. note::
   For the ethylene molecule (12 valence electrons)  we assume that we
   have the Dacapo output file ethylene.nc:

     >>> atoms = Calculator.read_atoms('ethylene.nc')
     >>> calc = atoms.get_calculator()
     >>> wannier = Wannier(numberofwannier=6,calculator=calc)

   This will construct 6 Wannier orbitals using the 6 occupied
   valence states, corresponding to the 12 valence electrons.


Below is a list of the most important methods of the :class:`Wannier`:

``Localize(step=0.5,tolerance=1.0e-08)``: Perform the localization of the
  Wannier orbitals. This method will localize the Wannier orbitals, i.e will try to
  maximize the localization functional. If the functional value is not increasing
  decrease the *step* size from the default value 0.5.

``GetWannierFunctionDOS(n,energies,width)``: Get projected density of states for WF.
  Returns the projected density of states (PDOS) for Wannier function n. The calculation
  is performed over the energy grid specified in energies. The PDOS is produced as a sum
  of Gaussians centered at the points of the energy grid and with the specified width.

``WriteWannierFunctionDOSToNetCDFFile(filename,n,energies,width)``:
  Same as GetWannierFunctionDOS, but writes the output to a NetCDF file.

``GetElectronicState(wannierindex,repeat=None)``: Returns an ``ElectronicState`` instance
  corresponding to the Wannier orbital with index *wannierindex*. The keyword repeat can be
  a list of 3 integers [n1,n2,n3], specifying how many times the unit cell is repeated
  along the unit cell basis vectors.

``GetCentersAsAtoms``: Returns a Atoms object with the Wannier centers.
  The chemical element is set to 'X'.

``TranslateAllWannierFunctionsToCell(cell)``: Move all Wannier orbitals to a specific unit cell.
  There exists an arbitrariness  in the positions of the Wannier orbitals relative to the
  unit cell. This method can move all orbitals to the unit cell specified by *cell*.
  For a gamma-point calculation, this has no effect. For
  a **k**-point calculation the periodicity of the orbitals are given by the large unit cell
  defined by repeating the original unitcell by the number of **k**-points in each direction.
  In this case it is usefull to move the orbitals away from the boundaries of the large cell
  before plotting them. For a bulk calculation with, say 10x10x10 **k** points, one could move
  the orbitals to the cell [2,2,2].
  In this way the pbc boundary conditions will not be noticed.

``WriteCube(wannierindex,filename,repeat=(7,7,7),real=False)``: Write a Cube formatted file.
  A Cube formatted file is written for the given wannier index.
  *repeat* can be used to repeat the unitcell, this is only relevant for calculations using
  **k**-points. In this case ``repeat``, will default be
  the number of **k**-points in each directions, i.e for a 11x11x11
  **k**-point set, repeat will be (11x11x11). This cell size represents the
  periodicity of the Wannier orbitals.

  Localized Wannier functions can often be chosen to be real.
  If the keyword *real* is set to *True*, the complex Wannier function will be transformed
  into a real one by multiplication be a suitable phase factor.
  In VMD you can use this to add two *isosurfaces* using  +- isosurface value, to get an
  approximation for the sign of the Wannier function.

``Save/ReadZIBlochMatrix(filename)``: Save and read ZI bloch matrix.
  These methods save and restore the localization matrix generated from the initial
  set of bloch function. This can save time, since the ZI matrix must be provided each time
  a localization is performed. If a ZI matrix is not read from file, it will be calculated.

``Save/ReadRotation(filename)``:
  These methods can be used to save and restore the unitary matrices used to produce a set
  of Wannier functions. The method produces two files: filename_rot.pickle and
  filename_coeff.pickle.


.. note::

   You can save your localized Wannier orbital like this::

     >>> wannier = Wannier(..)
     >>> wannier.save_z_i_bloch_matrix('fe_bloch.pickle')
     >>> wannier.localize(tolerance=0.000001)
     >>> wannier.save_rotation('fe')

   and read them in again like this::

     >>> wannier = Wannier(..)
     >>> wannier.read_z_i_bloch_matrix('fe_bloch.pickle')
     >>> wannier.read_rotation('fe')
     >>> wannier.localize(tolerance=0.000001)

   Localize should now converge in one step.

.. module:: dft.wannier
   :synopsis: Maximally localized Wannier functions

=====================================
Maximally localized Wannier functions
=====================================

.. default-role:: math


This page describes how to construct the Wannier orbitals using the
class :class:`Wannier`. The page is organized as follows:

* `Introduction`_: A short summary of the basic theory.
* `The Wannier class`_ : A description of how the Wannier class is
  used, and the methods defined within.
* `Band structure and orbital analysis`_ : A description of how to
  analyse the band structure using localized Wannier orbitals.

Introduction
============

The point of Wannier functions is the transform the extended Bloch
eigenstates of a DFT calculation, into a smaller set of states
designed to facilitate the analysis of e.g. chemical bonding. This is
achived by designing the Wannier functions to be localized in real
space instead of energy (which would be the eigen states).

The standard Wannier transformation is a unitary rotation of the Bloch
states. This implies that the Wannier functions (WF) span the same
Hilbert space as the Bloch states, i.e. they have the same eigenvalue
spectrum, and the original Bloch states can all be exactly reproduced
from a linear combination of the WF. For maximally localized Wannier
functions (MLWF), the unitary transformation is chosen such that the
spread of the resulting WF is minimized.

The standard choice is to make a unitary transformation of the
occupied bands only, thus resulting in as many WF as there are
occupied bands. If you make a rotation using more bands, the
localization will be improved, but the number of wannier functions
increase, thus making orbital based analysis harder.

The class defined here allows for construction of *partly* occupied
MLWF. In this scheme the transformation is still a unitary rotation
for the lowest states (the *fixed space*), but it uses a dynamically
optimized linear combination of the remaining orbitals (the *active
space*) to improve localization. This implies that e.g. the
eigenvalues of the Bloch states contained in the fixed space can be
exactly reproduced by the resulting WF, whereas the largest
eigenvalues of the WF will not necessarily correspond to any "real"
eigenvalues (this is irrelevant, as the fixed space is usually chosen
large enough, i.e. high enough above the fermilevel, that the
remaining DFT eigenvalues are meaningless anyway).

For the theory behind this method see the paper "Partly Occupied
Wannier Functions" Thygesen, Hansen and Jacobsen, *Phys. Rev. Lett*,
Vol. **94**, 26405 (2005).


The Wannier class
=================

The class is accessed by::

  from ase.dft import Wannier

The basic initialization is::

  Wannier(nwannier,
          calc,
          nbands=None,
          occupationenergy=None,
          fixedstates=None,
          spin=0)

The required arguments are:

``nwannier``: The number of Wannier functions you wish to construct.
  This must be at least half the number of electrons in
  the system and at most equal to the number of bands in the
  calculation.

``calc``: A converged DFT calculator class.
  A calculator containing all the information about the electronic
  structure calculation, i.e. number of bands, k-points, wavefunctions
  and so on. The calculator you are using *must* provide the method
  ``get_wannier_localization_matrix``.

The optional arguments are:

``nbands``: Bands to include in localization.
  You can specify ``nbands`` as a smaller number than the
  number of bands in the calculator. This is useful if you do not
  believe the highest bands of your calculation are well
  converged. This value defaults to the number of bands in the
  calculation.

``spin``: The spin channel to be considered.
  The Wannier code treats each spin channel independently.

``occupationenergy``: Maximal energy of the fixed part of Hilbert space.
  Eigenstates below this energy are included in the *fixed space*. In
  practice this means that all eigenstates below this energy can be
  exactly reproduced in terms of the resulting Wannier functions.  The
  energy is relative to the Fermi-level.

``fixedstates``: Number of states in the fixed part of Hilbert space.
  Fix the number of states to be included in the fixed space, starting
  from the lowest eigenstate.  This keyword provides another way of
  specifying how many states should be fixed, and overrides
  ``occupationenergy`` if this is also set. Default is None meaning
  that the number of fixed states is equated to ``nwannier``.

Below is a list of the most important methods of the :class:`Wannier`:

.. function:: initialize(calc, initialwannier=None, seed=None, bloch=False, first=False)

  This method initializes the class, and should always be called, if
  the calculation is not restarted from file (using the method
  ``load``).
  
  The keywords are described below

  ``initialwannier``: Setup an initial set of Wannier orbitals.
    *initialwannier* can set up a starting guess for the Wannier
    functions.  This is important to speed up convergence in
    particular for large systems For transition elements with **d**
    electrons you will always find 5 highly localized **d**-orbitals
    centered at the atom.  Placing 5 **d**-like orbitals with a radius
    of 0.4 Angstroms and center at atom no. 7, and 3 **p**-like
    orbitals with a radius of 0.4 Angstroms and center at atom no. 27
    looks like this::

       initialwannier = [[[7],2,0.4],[[27],1,0.4]]

    Placing only the l=2, m=-2 and m=-1 orbitals at atom no. 7 looks
    like this::

       initialwannier = [[[7],2,-2,0.4],[[7],2,-1,0.4]]

    I.e. if you do not specify the m quantum number all allowed values
    are used.  Instead of placing an orbital at an atom, you can place
    it at a specified position. For example the following::

       initialwannier = [[[0.5,0.5,0.5],0,0.5]]

    places an **s** orbital with radius 0.5 Angstroms at the position
    (0.5,0.5,0.5) in scaled coordinates of the unit cell.

  ``seed``: The seed used for any randomly generated initial rotations.

  ``bloch``: Use Bloch states for initial guess
    If ``True``, sets the initial guess for the rotation matrix to be
    identity, i.e. the Bloch states are used.

  ``first``: First time initialize is called?
    Indicates if the localization matrix of the Bloch states should be
    re-evaluated. Setting this to False is useful if you wish to
    randomize the rotation matrices (or use initialwannier, or
    bloch=True) after having initialized once, or loaded the
    localization matrix from file.

.. function:: localize(step=0.25, tolerance=1.0e-08)

.. function:: dump(file)

   Save the rotation-, coefficient-, and wannier localization matrices
   to indicated ``file`` (string).

.. function:: load(file)

   Load the rotation-, coefficient-, and wannier localization matrices
   from indicated ``file`` (string).

.. function:: get_function(calc, index, repeat=None)

   Returns an array with the funcion values of the indicated Wannier
   function on a grid with the size of the *repeated* unit cell.
   
   For a calculation using **k**-points the relevant unit cell for
   eg. visualization of the Wannier orbitals is not the original unit
   cell, but rather a larger unit cell defined by repeating the
   original unit cell by the number of **k**-points in each direction.
   We will refer to this unit cell as the large unit cell.  Note that
   for a `\Gamma`-point calculation the large unit cell coinsides with
   the original unit cell.  The large unitcell defines also the
   periodicity of the Wannier orbitals.

.. function:: get_centers()

.. function:: get_radii()

.. function:: get_pdos(calc, w, energies, width)

  Returns the projected density of states (PDOS) for Wannier function
  ``w``. The calculation is performed over the energy grid specified
  in energies. The PDOS is produced as a sum of Gaussians centered at
  the points of the energy grid and with the specified width.

.. function:: translate_all_to_cell(cell=[0, 0, 0])

  Move all Wannier orbitals to a specific unit cell.  There exists an
  arbitrariness in the positions of the Wannier orbitals relative to
  the unit cell. This method can move all orbitals to the unit cell
  specified by ``cell``.  For a `\Gamma`-point calculation, this has no
  effect. For a **k**-point calculation the periodicity of the
  orbitals are given by the large unit cell defined by repeating the
  original unitcell by the number of **k**-points in each direction.
  In this case it is usefull to move the orbitals away from the
  boundaries of the large cell before plotting them. For a bulk
  calculation with, say 10x10x10 **k** points, one could move the
  orbitals to the cell [2,2,2].  In this way the pbc boundary
  conditions will not be noticed.

.. function:: write_cube(calc, index, fname, repeat=None, real=True)

.. function:: get_hopping(R, calc)

   Returns the matrix H(R)_nm = <0,n|H|R,m>, where R is the
   cell-distance (in units of the basis vectors of the small cell) and
   n,m are indices of the Wannier functions.

.. function:: get_hamiltonian(calc, k=0):
   
   Get Wannier function Hamiltonian at existing k-vector of index k

.. function:: get_hamiltonian_kpoint(kpt_c, calc):
   
   Get Wannier function Hamiltonian at some new arbitrary k-vector.

   XXX Not implemented yet!

For examples of how to use the **Wannier** class, see the `Wannier tutorial`_.

.. _Wannier tutorial: https://wiki.fysik.dtu.dk/ase/tutorials/wannier.html

.. note:: For calculations using **k**-points, make sure that the
   `\Gamma`-point is included in the **k**-point grid. Moreover you must
   shift all **k**-points by a small amount (but not less than 2e-5
   in) in e.g. the x direction, before performing the Dacapo
   calculation. If this is not done the symmetry program in Dacapo
   will use time-reversal symmetry to reduce the number of
   **k**-points by a factor 2. The shift can be performed like this::

                kpoints = calc.get_b_z_k_points()
                kpoints[:,0] += 2e-5
                calc.set_b_z_k_points(kpoints)

..  The cutoff distance truncates the Wannier orbitals at the
    specified distance. This distance should be smaller than half the
    length of large unitcell. The truncation is necessary because the
    Wannier functions will always be periodic (with a periodicity
    given by the large cell), and thus in order to describe completely
    localized orbitals the WFs must be truncated.


    ``GetWFHamiltonian()``: The Hamiltonian matrix in the basis of the
        Wannier orbitals are returned.  We will refer to this
        Hamiltonian as **H** in that follows. The Hamiltonian refers
        to the large unit cell, and its dimension is therefore
        (N_w*N_k)x(N_w*N_k), where N_w is the number of Wannier
        functions in a unit cell and N_k is the number of **k**
        points. Periodic boundary conditions are imposed on the
        boundaries of the large cell.

.. default-role::

.. _defects:

======================================================
Tools for defect calculations
======================================================

This section gives an (incomplete) overview of features in ASE that
help in the preparation and analysis of supercell calculations as most
commonly employed in the computation of defect properties.

Supercell creation
==============================

Background
---------------------

One is usually interested in defect properties in the so-called dilute
limite, i.e. under conditions, in which defect-defect interactions are
negligible. While alternative approaches in particular embedding
techniques \cite{} exist, the most common approach is to use
supercells. To this end, one creates a supercell by a *suitable* (see
below) repetition of the primitive unit cell, after which a defect,
e.g., a vacancy or an impurity atom, is inserted.

.. image:: supercell-approach.png

The calculation thus corresponds to a periodic arrangement of
defects. Accordingly, care must be taken to keep the interactions
between defects as small as possible, which generally calls for large
supercells. It is furthermore indicated to maximize the defect-defect
separation in *all* directions, which is in principle achieved if the
supercell used has a cubic (or close to cubic) shape.

.. image:: cubic-supercell.png

While this is readily achievable for cubic materials (diamond,
zincblende, face-centered cubic ...) by using simple repetitions of
the *conventional* unit cell, for countless materials of lower
symmetry the choice of a supercell is not necessarily straightforward.

Furthermore, electrostatic and strain interactions between periodic
images die of very slowly with distance. While various correction
schemes have been developed, the most reliable approach is still
finite-size extrapolation using supercells of different size.


Finding optimal supercell shapes
------------------------------------------


The above considerations illustrate the need for a more systematic
approach to supercell construction. A simple scheme to construct
"optimal" supercells is described in [Erhart]_. Optimality here
implies that one identifies the supercell that for a given size
(number of atoms) most closely approximates a cube. This approach
ensures that the defect separation is large and that the electrostatic
interactions exhibit a systematic scaling.

The ideal cubic cell metric for a given volume :math:`Omega` is
simply given by :math:`\mathbf{h}_\text{cub} = \Omega^{1/3}
\mathbf{I}`, which in general does not satisfy the crystallographic
boundary conditions. The :math:`l_2`-norm provides a convenient
measure of the deviation of any other cell metric from a cubic
shape. The "acubicity" can thus be defined as

.. math:: \Delta_c(\mathbf{h}) = ||\mathbf{h} - \mathbf{h}_\text{cub}||_2

Any cell metric that is compatible with the crystal symmetry can be
written in the form

.. math:: \mathbf{h} = \mathbf{P} \mathbf{h}_p

where

.. math:: \mathbf{P} \in \mathbb{Z}^{3\times3}.

The matrix :math:`\mathbf{P}`, which yields the optimal cell shape for
a given cell size can then be defined as

.. math:: \mathbf{P}_\text{opt} = \underset{\mathbf{P}}{\operatorname{argmin}} \left\{ \Delta_c\left(\mathbf{Ph}_p\right) | \det\mathbf{P} = N_{uc}\right\},

where :math:`N_{uc}` defines the size of the supercell in terms of the
number of primitive unit cells. This approach is general and can be
applied to generate suitable supercells for arbitrary primitive cell
metrics.

Implementation
------------------------------------------

The algorithm described above requires finding a matrix
:math:`\mathbf{P}` from a *discrete* set of integer matrices. This can
be achieved by a brute force loop over a subset of :math:`\mathbf{P}`
matrices

.. math:: \mathbf{P} \in [n_{min},n_{max}]^{3\times3}.

Since this involves looping over :math:`(n_{max}-n_{min}+1)^9`
matrices, which is a staggeringly large number already for a moderate
range of integers, a plain python implementation is rather
slow. Instead the algorithm has been implemented as inline C using
`scipy.weave
<http://docs.scipy.org/doc/scipy/reference/tutorial/weave.html>`_. For
this reason the respective code is currently not included as an ASE
module but is made available as ascript that can be executed from the
command line, :download:`find_optimal_cell_shape.py`. Using inline C
reduces the execution time by at least two orders of magnitude and
:program:`find_optimal_cell_shape.py` takes only a few seconds to run
on a standard desktop given a trial range of :math:`n_{max},n_{min}
\lessim 16`.

For illustration consider the following example. First we set up a
primitive face-centered cubic (fcc) unit cell::

    ase-build -x fcc -a 4.05 Al Al.xyz

after which we call :program:`find_optimal_cell_shape.py` to find the
cell metric :math:`\mathbf{h}` and transformation matrix :math:`P`
that yield the most cubic shape for a supercell comprising 32
primitive unit cells::
     
     find_optimal_cell_shape.py 32 Al.xyz

This will produce the following output::

     target size in number of unit cells: 32
     primitive cell read from file: Al.xyz
     primitive cell metric:
     [[ 0.     2.025  2.025]
      [ 2.025  0.     2.025]
      [ 2.025  2.025  0.   ]]
     searching over 387420489 matrices with coefficients in the range from -4 to 4
     best acubicity value: 4.37527e-24
     optimal P:
     [[-2.  2.  2.]
      [ 2. -2.  2.]
      [ 2.  2. -2.]]
     cell:
     [[ 8.1  0.   0. ]
      [ 0.   8.1  0. ]
      [ 0.   0.   8.1]]
     determinant of optimal P:  32.0

and thus

.. math:: \mathbf{P}_\text{opt} = \left(\begin{array} -2 & 2 & 2 \\ 2 & -2 & 2 \\ 2 & 2 & -2 \end{array}\right) \\
          \mathbf{h}_{opt} = \left(\begin{array} 8.1 & 0 & 0 \\ 0 & 8.1 & 0 \\ 0 & 0 & 8.1 \end{array}\right)\,\text{\AA},

which is the expected outcome as it corresponds to a
:math:`2\times2\times2` repetition of the *conventional* (4-atom) unit
cell. On the other hand repeating this exercise with::

      find_optimal_cell_shape.py 90 Al.xyz

yields a less obvious results, namely

.. math:: \mathbf{P}_\text{opt} = \left(\begin{array} -3 & 3 & 3 \\ 3 & -3 & 3 \\ 2 & 3 & -3 \end{array}\right) \\
          \mathbf{h}_{opt} = \left(\begin{array} 12.24 & 0 & 0 \\ 0 & 12.24 & 0 \\ 0 & -2.04 & 10.2 \end{array}\right)\,\text{\AA},

which indeed corresponds to a reasonably cubic cell shape.

Since this procedure requires only knowledge of the cell metric (and not the atomic positions) for standard metrics, e.g., fcc, bcc, and simple cubic one can generte series of shapes that are usable for *all* structures with the respective metric. For example the :math:`\mathbf{P}_\text{opt}` matrices that optimize the shape of a supercell build using a primitive FCC cell are directly applicable to diamond and zincblende lattices.

For convenience the :math:`\mathbf{P}_\text{opt}` matrices for the aforementioned lattices have already been generated for :math:`N_{uc}\leq800` and are provide here as dictionaries in python pickle format for download: face-centered cubic (:download:`Popt-fcc.pkl`), body-centerd cubic (:download:`Popt-bcc.pkl`), simple cubic (:download:`Popt-sc.pkl`).

.. [Erhart] P. Erhart, B. Sadigh, A. Schleife, and D. Åberg.
   First-principles study of codoping in lanthanum bromide,
   Phys. Rev. B, Vol **91**, pp. 165206 (2012); Appendix C

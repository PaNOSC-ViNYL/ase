.. _defects:

======================================================
Tools for defect calculations
======================================================

This section gives an (incomplete) overview of features in ASE that
help in the preparation and analysis of supercell calculations as most
commonly employed in the computation of defect properties.

.. contents::

Supercell creation
==============================

Background
---------------------

Defect properties are most commonly investigated in the so-called
dilute limit, i.e. under conditions, in which defect-defect
interactions are negligible. While alternative approaches in
particular embedding techniques exist, the most common approach is to
use supercells. To this end, one creates a supercell by a *suitable*
(see below) repetition of the primitive unit cell, after which a
defect, e.g., a vacancy or an impurity atom, is inserted. This
procedure can be schematically depicted as follows:

.. image:: supercell-1.svg
   :width: 30%
.. image:: supercell-2.svg
   :width: 30%
.. image:: supercell-3.svg
   :width: 30%

The calculation thus corresponds to a periodic arrangement of
defects. Accordingly, care must be taken to keep the interactions
between defects as small as possible, which generally calls for large
supercells. It is furthermore indicated to maximize the defect-defect
separation in *all* directions, which is in principle achieved if the
supercell used has a suitable shape. Consider for illustration the
following three 2D lattices with identical unit cell area but
different lattice symmetry:

.. image:: periodic-images-1.svg
   :width: 30%
.. image:: periodic-images-2.svg
   :width: 30%
.. image:: periodic-images-3.svg
   :width: 30%

In the case of the square lattice, each defect has :math:`Z_1=4`
nearest neighbors at a distance of :math:`r_1=a_0`, where
:math:`a_0=\sqrt{A}` with :math:`A` being the unit cell area. By
comparison in a rectangular lattice with an aspect ratio of 2:1, the
defects are much closer to each other with :math:`r_1 = 0.5 a_0` and
:math:`Z_1=2`. The largest defect-defect distance (at constant unit
cell area) is obtained for the hexagonal lattice, which also
correponds to the most closely packed 2D arrangement. Here, one
obtains :math:`r_1=\sqrt{2}/\sqrt[4]{3}=1.075 a_0` and
:math:`Z_1=6`. For defect calculation supercells corresponding to
hexagonal or square lattices have thus clear advantages. This argument
can be extended to 3D: Square lattices in 2D correspond to cubic
lattices (supercells) in 3D with :math:`r_1=a_0` and
:math:`Z_1=6`. The 3D analogue of the hexagonal 2D lattice are
hexagonal and cubic close packed structures, both of which yield
:math:`r_1 = \sqrt{3}/2 a_0` and :math:`Z_1=12`.

It is straightforward to construct cubic or face-centered cubic (fcc,
cubic closed packed) supercells for cubic materials (including e.g,
diamond and zincblende) by using simple repetitions of the
conventional or primitive unit cells. For countless materials of lower
symmetry the choice of a supercell is, however not necessarily so
simple. The algorithm below represents a general solution to this
issue.

In the case of semiconductors and insulators with small dielectric
constants, defect-defect interactions are particularly pronounced due
to the weak screening of long-ranged electrostatic interactions. While
various correction schemes have been proposed, the most reliable
approach is still finite-size extrapolation using supercells of
different size. In this case care must be taken to use a sequence of
self-similar supercells in order for the extrapolation to be
meaningful. To motivate this statement consider that the leading
(monopole-monopole) term :math:`E_{mp}`, which scales with :math:`1/r`
and is proportional to the (ionic) dielectric constant
:math:`\epsilon_0`. The :math:`E_{mp}` term is geometry dependent and
in the case of simple lattices the dependence is easily expressed by
the Madelung constant. The geometry dependence implies that different
(super)cell shapes fall on different lines when plotting e.g., the
formation energy as a function of :math:`N^{-1/3}` (equivalent to an
effective inverse cell size, :math:`1/L \propto N^{-1/3}`. For
extrapolation one should therefore only use geometrically equivalent
cells or at least cells that are as self-similar to each other as
possibly, see Fig.~10 in [Erhart]_ for a very clear example. In this
context there is therefore also a particular need for supercells of a
particular shape.



Algorithm for finding optimal supercell shapes
-----------------------------------------------

The above considerations illustrate the need for a more systematic
approach to supercell construction. A simple scheme to construct
"optimal" supercells is described in [Erhart]_. Optimality here
implies that one identifies the supercell that for a given size
(number of atoms) most closely approximates the desired shape, most
commonly a simple cubic or fcc metric. This approach ensures that the
defect separation is large and that the electrostatic interactions
exhibit a systematic scaling.

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
metrics. Specifically, in order to obtain supercells that resemble the
shape of a fcc primitive cell one simply needs to replace the above
definition of the "acubicity" with the following expression

.. math:: \Delta(\mathbf{h}) = ||\mathbf{h} - \mathbf{h}_\text{fcc}||_2

where

.. math:: \mathbf{h}_\text{fcc} = \Omega^{1/3} \frac{1}{2} \left(\begin{array}{rrr} 0 & 1 & 1 \\ 1 & 0 & 1 \\ 1 & 1 & 0 \end{array}\right).



Implementation of algorithm
------------------------------------------

The algorithm described above requires finding a matrix
:math:`\mathbf{P}` from a *discrete* set of integer matrices. This can
be achieved by a brute force loop over a subset of :math:`\mathbf{P}`
matrices [#algorithm-note]_

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
on a standard desktop given a trial range of :math:`n_{max}-n_{min}
\lesssim 16`.

For illustration consider the following example. First we set up a
primitive face-centered cubic (fcc) unit cell::

    ase-build -x fcc -a 4.05 Al Al.xyz

after which we call :program:`find_optimal_cell_shape.py` to find the
cell metric :math:`\mathbf{h}` and transformation matrix
:math:`\mathbf{P}` that yield the cell that closely resembles a simple
cubic (sc) cell metric for a supercell comprising 32 primitive unit
cells::
     
     find_optimal_cell_shape.py 32 Al.xyz sc

This will produce the following output::

     target size in number of unit cells: 32
     primitive cell read from file: Al.xyz
     primitive cell metric:
     [[ 0.     2.025  2.025]
      [ 2.025  0.     2.025]
      [ 2.025  2.025  0.   ]]
     effective target dimension: 8.1
     target metric:
     [[ 8.1  0.   0. ]
      [ 0.   8.1  0. ]
      [ 0.   0.   8.1]]
     searching over 10604499373 matrices with coefficients in the range from -6 to 6
     smallest |P h_p - h_target|_2 value: 0.000000
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

.. math:: \mathbf{P}_\text{opt} = \left(\begin{array}{rrr} -2 & 2 & 2 \\ 2 & -2 & 2 \\ 2 & 2 & -2 \end{array}\right) \quad
          \mathbf{h}_{opt} = \left(\begin{array}{ccc} 2 a_0 & 0 & 0 \\ 0 & 2 a_0 & 0 \\ 0 & 0 & 2 a_0 \end{array}\right),

where :math:`a_0` =4.05 Å is the lattice constant. This is indeed the
expected outcome as it corresponds to a :math:`2\times2\times2`
repetition of the *conventional* (4-atom) unit cell. On the other hand
repeating this exercise with::

      find_optimal_cell_shape.py 90 Al.xyz

yields a less obvious result, namely

.. math:: \mathbf{P}_\text{opt} = \left(\begin{array}{rrr} -3 & 3 & 3 \\ 3 & -3 & 3 \\ 2 & 3 & -3 \end{array}\right) \quad
          \mathbf{h}_{opt} = \left(\begin{array}{ccc} 3 a_0 & 0 & 0 \\ 0 & 3 a_0 & 0 \\ 0 & -0.5 a_0 & 2.5 a_0 \end{array}\right),

which indeed corresponds to a reasonably cubic cell shape.

Since this procedure requires only knowledge of the cell metric (and
not the atomic positions) for standard metrics, e.g., fcc, bcc, and
simple cubic one can generate series of shapes that are usable for
*all* structures with the respective metric. For example the
:math:`\mathbf{P}_\text{opt}` matrices that optimize the shape of a
supercell build using a primitive FCC cell are directly applicable to
diamond and zincblende lattices.

For convenience the :math:`\mathbf{P}_\text{opt}` matrices for the
aforementioned lattices have already been generated for
:math:`N_{uc}\leq800` and are provided here as dictionaries in `json
<https://en.wikipedia.org/wiki/JSON>`_ format.

 * Transformation of face-centered cubic metric to simple cubic-like shapes: :download:`Popt-fcc-sc.json`
 * Transformation of face-centered cubic metric to face-centered cubic-like shapes: :download:`Popt-fcc-fcc.json`
 * Transformation of body-centered cubic metric to simple cubic-like shapes: :download:`Popt-bcc-sc.json`
 * Transformation of body-centered cubic metric to face-centered cubic-like shapes: :download:`Popt-bcc-fcc.json`
 * Transformation of simple cubic metric to simple cubic-like shapes: :download:`Popt-sc-sc.json`
 * Transformation of simple cubic metric to face-centered cubic-like shapes: :download:`Popt-sc-fcc.json`




Generation of supercell
------------------------------------------

Once 


.. [#algorithm-note] One can easily conceive of more refined variants
   but here we stick to the most obvious/trivial version as it
   suffices for the task at hand.

.. [Erhart] P. Erhart, B. Sadigh, A. Schleife, and D. Åberg.
   First-principles study of codoping in lanthanum bromide,
   Phys. Rev. B, Vol **91**, pp. 165206 (2012); Appendix C

.. _defects:

=============================
Tools for defect calculations
=============================

This section gives an (incomplete) overview of features in ASE that
help in the preparation and analysis of supercell calculations as most
commonly employed in the computation of defect properties.

.. contents::

Supercell creation
==================

Background
----------

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
effective inverse cell size, :math:`L^{-1} \propto N^{-1/3}`. For
extrapolation one should therefore only use geometrically equivalent
cells or at least cells that are as self-similar to each other as
possibly, see Fig. 10 in [Erhart]_ for a very clear example. In this
context there is therefore also a particular need for supercells of a
particular shape.


Algorithm for finding optimal supercell shapes
----------------------------------------------

The above considerations illustrate the need for a more systematic
approach to supercell construction. A simple scheme to construct
"optimal" supercells is described in [Erhart]_. Optimality here
implies that one identifies the supercell that for a given size
(number of atoms) most closely approximates the desired shape, most
commonly a simple cubic or fcc metric (see above). This approach
ensures that the defect separation is large and that the electrostatic
interactions exhibit a systematic scaling.

The ideal cubic cell metric for a given volume :math:`\Omega` is simply
given by :math:`\Omega^{1/3} \mathbf{I}`, which in general does not
satisfy the crystallographic boundary conditions. The :math:`l_2`-norm
provides a convenient measure of the deviation of any other cell
metric from a cubic shape. The optimality measure can thus be defined
as

.. math:: \Delta_\text{sc}(\mathbf{h}) = ||\mathbf{h} - \Omega^{1/3} \mathbf{1}||_2,

Any cell metric that is compatible with the crystal symmetry can be
written in the form

.. math:: \mathbf{h} = \mathbf{P} \mathbf{h}_p

where :math:`\mathbf{P} \in \mathbb{Z}^{3\times3}` and
:math:`\mathbf{h}_p` is the primitive cell metric.  This approach can
be readily generalized to arbitrary target cell metrics. In order to
obtain a measure that is size-independent it is furthermore convenient
to introduce a normalization, which leads to the expression
implemented here, namely

.. math:: \bar{\Delta}(\mathbf{Ph}_p) = ||Q\mathbf{Ph}_p - \mathbf{h}_\text{target}||_2,

where :math:`Q = \left(\det\mathbf{h}_\text{target} \big/
\det\mathbf{h}_p\right)^{1/3}` is a normalization factor.  The
matrix :math:`\mathbf{P}_\text{opt}` that yields the optimal cell
shape for a given cell size can then be obtained by

.. math:: \mathbf{P}_\text{opt} = \underset{\mathbf{P}}{\operatorname{argmin}} \left\{ \bar\Delta\left(\mathbf{Ph}_p\right) | \det\mathbf{P} = N_{uc}\right\},

where :math:`N_{uc}` defines the size of the supercell in terms of the
number of primitive unit cells.


Implementation of algorithm
---------------------------

The algorithm described above has been implemented in two versions:

* :func:`~ase.build.find_optimal_cell_shape`: An
  implementation that invokes inline-C code via the `scipy.weave
  <http://docs.scipy.org/doc/scipy/reference/tutorial/weave.html>`_
  interface as is substantially (about one to two orders of magnitude)
  faster than the pure python version. `scipy
  <https://www.scipy.org/>`_ is available on many systems for
  scientific computing and should be preferable if available.

* :func:`~ase.build.find_optimal_cell_shape_pure_python`:
  A pure python implementation that is, however, quite slow.

For illustration consider the following example. First we set up a
primitive face-centered cubic (fcc) unit cell, after which we call
:func:`~ase.build.find_optimal_cell_shape` to obtain a
:math:`\mathbf{P}` matrix that will enable us to generate a supercell
with 32 atoms that is as close as possible to a simple cubic shape::
                                 
   from ase.build import bulk
   from ase.build import find_optimal_cell_shape, get_deviation_from_optimal_cell_shape
   import numpy as np
   conf = bulk('Au')
   P1 = find_optimal_cell_shape(conf.cell, 32, 'sc')
 
This yields

.. math:: \mathbf{P}_1 = \left(\begin{array}{rrr} -2 & 2 & 2 \\ 2 & -2 & 2 \\ 2 & 2 & -2 \end{array}\right) \quad
          \mathbf{h}_1 = \left(\begin{array}{ccc} 2 a_0 & 0 & 0 \\ 0 & 2 a_0 & 0 \\ 0 & 0 & 2 a_0 \end{array}\right),

where :math:`a_0` =4.05 Å is the lattice constant. This is indeed the
expected outcome as it corresponds to a :math:`2\times2\times2`
repetition of the *conventional* (4-atom) unit cell. On the other hand
repeating this exercise with::

   P2 = find_optimal_cell_shape(conf.cell, 495, 'sc')

yields a less obvious result, namely

.. math:: \mathbf{P}_2 = \left(\begin{array}{rrr} -5 & 5 & 5 \\ 5 & -4 & 5 \\ 5 & 5 & -4 \end{array}\right) \quad
          \mathbf{h}_2 = a_0 \left(\begin{array}{ccc}  5 & 0 & 0 \\ 0.5 & 5 & 0.5 \\ 0.5 & 0.5 & 5 \end{array}\right),

which indeed corresponds to a reasonably cubic cell shape. One can
also obtain the optimality measure :math:`\bar{\Delta}` by executing::

   dev1 = get_deviation_from_optimal_cell_shape(np.dot(P1, conf.cell)
   dev2 = get_deviation_from_optimal_cell_shape(np.dot(P2, conf.cell)

which yields :math:`\bar{\Delta}(\mathbf{P}_1)=0` and
:math:`\bar{\Delta}(\mathbf{P}_2)=0.201`.

Since this procedure requires only knowledge of the cell metric (and
not the atomic positions) for standard metrics, e.g., fcc, bcc, and
simple cubic one can generate series of shapes that are usable for
*all* structures with the respective metric. For example the
:math:`\mathbf{P}_\text{opt}` matrices that optimize the shape of a
supercell build using a primitive FCC cell are directly applicable to
diamond and zincblende lattices.

For convenience the :math:`\mathbf{P}_\text{opt}` matrices for the
aforementioned lattices have already been generated for
:math:`N_{uc}\leq2000` and are provided here as dictionaries in `json
<https://en.wikipedia.org/wiki/JSON>`_ format.

 * Transformation of face-centered cubic metric to simple cubic-like shapes: :download:`Popt-fcc2sc.json`
 * Transformation of face-centered cubic metric to face-centered cubic-like shapes: :download:`Popt-fcc2fcc.json`
 * Transformation of body-centered cubic metric to simple cubic-like shapes: :download:`Popt-bcc2sc.json`
 * Transformation of body-centered cubic metric to face-centered cubic-like shapes: :download:`Popt-bcc2fcc.json`
 * Transformation of simple cubic metric to simple cubic-like shapes: :download:`Popt-sc2sc.json`
 * Transformation of simple cubic metric to face-centered cubic-like shapes: :download:`Popt-sc2fcc.json`

The thus obtained :math:`\bar{\Delta}` values are shown as a function
of the number of unit cells :math:`N_{uc}` in the panel below, which
demonstrates that this approach provides access to a large number of
supercells with e.g., simple cubic or face-centered cubic shapes that
span the range between the "exact" solutions, for which
:math:`\bar{\Delta}=0`. The algorithm is, however, most useful for
non-cubic cell shapes, for which finding several reasonably sized cell
shapes is more challenging as illustrated for a hexagonal material
(LaBr\ :sub:`3`) in [Erhart]_.

.. image:: score-size-sc2sc.svg
   :width: 30%
.. image:: score-size-fcc2sc.svg
   :width: 30%
.. image:: score-size-bcc2sc.svg
   :width: 30%
.. image:: score-size-sc2fcc.svg
   :width: 30%
.. image:: score-size-fcc2fcc.svg
   :width: 30%
.. image:: score-size-bcc2fcc.svg
   :width: 30%

   
Generation of supercell
-----------------------

Once the transformation matrix :math:`\mathbf{P}` it is
straightforward to generate the actual supercell using e.g., the
:func:`~ase.build.cut` function. A convenient interface is provided by
the :func:`~ase.build.make_supercell` function, which is invoked as
follows::

   from ase.build import bulk
   from ase.build import find_optimal_cell_shape
   from ase.build import make_supercell
   conf = bulk('Au')
   P = find_optimal_cell_shape(conf.cell, 495, 'sc')
   supercell = make_supercell(conf, P)


   
.. [Erhart] P. Erhart, B. Sadigh, A. Schleife, and D. Åberg.
   First-principles study of codoping in lanthanum bromide,
   Phys. Rev. B, Vol **91**, 165206 (2012),
   `doi: 10.1103/PhysRevB.91.165206 <http://dx.doi.org/10.1103/PhysRevB.91.165206>`_; Appendix C

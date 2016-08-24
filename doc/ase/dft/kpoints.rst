.. module:: ase.dft.kpoints
   :synopsis: Brillouin zone sampling

=======================
Brillouin zone sampling
=======================

The **k**-points are always given relative to the basis vectors of the
reciprocal unit cell.


Monkhorst-Pack
--------------

.. autofunction:: monkhorst_pack

The k-points are given as [MonkhorstPack]_:

.. math::

    \sum_{i=1,2,3} \frac{2n_i -N_i - 1}{2N_i} \mathbf{b}_i,

where `n_i=1,2,...,N_i`, ``size`` = `(N_1, N_2, N_3)` and the
`\mathbf{b}_i`'s are reciprocal lattice vectors.

.. autofunction:: get_monkhorst_pack_size_and_offset

Example:

>>> from ase.dft.kpoints import *
>>> monkhorst_pack((4, 1, 1))
array([[-0.375,  0.   ,  0.   ],
       [-0.125,  0.   ,  0.   ],
       [ 0.125,  0.   ,  0.   ],
       [ 0.375,  0.   ,  0.   ]])
>>> get_monkhorst_pack_size_and_offset([[0, 0, 0]])
(array([1, 1, 1]), array([ 0.,  0.,  0.]))


.. [MonkhorstPack]
    Hendrik J. Monkhorst and James D. Pack:
    *Special points for Brillouin-zone integrations*,
    Phys. Rev. B 13, 5188–5192 (1976)


Special points in the Brillouin zone
------------------------------------

.. data:: special_points

Special points from [Setyawana-Curtarolo]_:

.. list-table::
    :widths: 10 30 35
    
    * - Cubic
      - `\Gamma`-`X`-`M`-`\Gamma`-`R`-`X`, `M`-`R`
      - .. image:: cubic.svg
            :width: 25 %
    * - FCC
      - `\Gamma`-`X`-`W`-`K`-`\Gamma`-`L`-`U`-`W`-`L`-`K`, `U`-`X`
      - .. image:: fcc.svg
            :width: 25 %
    * - BCC
      - `\Gamma`-`H`-`N`-`\Gamma`-`P`-`H`, `P`-`N`
      - .. image:: bcc.svg
            :width: 25 %
    * - Tetragonal
      - `\Gamma`-`X`-`M`-`\Gamma`-`Z`-`R`-`A`-`Z`, `X`-`R`, `M`-`A`
      - .. image:: tetragonal.svg
            :width: 25 %
    * - Orthorhombic
      - `\Gamma`-`X`-`S`-`Y`-`\Gamma`-`Z`-`U`-`R`-`T`-`Z`, `Y`-`T`, `U`-`X`, `S`-`R`
      - .. image:: orthorhombic.svg
            :width: 25 %
    * - Hexagonal
      - `\Gamma`-`M`-`K`-`\Gamma`-`A`-`L`-`H`-`A`, `L`-`M`, `K`-`H`
      - .. image:: hexagonal.svg
            :width: 25 %
    * - Monoclinic
      - `\Gamma`-`Y`-`H`-`C`-`E`-`M_1`-`A`-`X`-`H_1`, `M`-`D`-`Z`, `Y`-`D`
      - .. image:: monoclinic.svg
            :width: 25 %

.. [Setyawana-Curtarolo]
    High-throughput electronic band structure calculations:
    Challenges and tools

    Wahyu Setyawana, Stefano Curtarolo
    
    Computational Materials Science,
    Volume 49, Issue 2, August 2010, Pages 299–312
    
    http://dx.doi.org/10.1016/j.commatsci.2010.05.010

You can find the special points in the Brillouin zone:

>>> from ase.lattice import bulk
>>> from ase.dft.kpoints import get_special_points, get_bandpath
>>> si = bulk('Si', 'diamond', a=5.459)
>>> points = get_special_points('fcc', si.cell)
>>> GXW = [points[k] for k in 'GXW']
>>> kpts, x, X = get_bandpath(GXW, si.cell, 100)
>>> print(kpts.shape, len(x), len(X))
(100, 3) 100 3

.. autofunction:: get_special_points
.. autofunction:: get_bandpath


High symmetry paths
-------------------

.. data:: special_paths

The ``special_paths`` dictionary contains suggestions for high symmetry
paths in the BZ from the [Setyawana-Curtarolo]_ paper.

>>> from ase.dft.kpoints import special_paths
>>> paths = special_paths['bcc']
>>> paths
[['G', 'H', 'N', 'G', 'P', 'H'], ['P', 'N']],
>>> points = get_special_points('bcc', [[-1, 1, 1], [1, -1, 1], [1, 1, -1]])
>>> points
{'H': [0.5, -0.5, 0.5], 'N': [0, 0, 0.5], 'P': [0.25, 0.25, 0.25], 'G': [0, 0, 0]}
>>> kpts = [points[k] for k in paths[0]]  # G-H-N-G-P-H
>>> kpts
[[0, 0, 0], [0.5, -0.5, 0.5], [0, 0, 0.5], [0, 0, 0], [0.25, 0.25, 0.25], [0.5, -0.5, 0.5]]


Chadi-Cohen
-----------

Predefined sets of **k**-points:

.. data:: cc6_1x1
.. data:: cc12_2x3
.. data:: cc18_sq3xsq3
.. data:: cc18_1x1
.. data:: cc54_sq3xsq3
.. data:: cc54_1x1
.. data:: cc162_sq3xsq3
.. data:: cc162_1x1


Naming convention: ``cc18_sq3xsq3`` is 18 **k**-points for a
sq(3)xsq(3) cell.

Try this:

>>> import numpy as np
>>> import pylab as plt
>>> from ase.dft.kpoints import cc162_1x1
>>> B = [(1, 0, 0), (-0.5, 3**0.5 / 2, 0), (0, 0, 1)]
>>> k = np.dot(cc162_1x1, B)
>>> plt.plot(k[:, 0], k[:, 1], 'o')
[<matplotlib.lines.Line2D object at 0x9b61dcc>]
>>> p.show()

.. image:: cc.png

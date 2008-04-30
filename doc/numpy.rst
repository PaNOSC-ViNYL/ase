Numeric arrays in Python
========================

:term:`ndarray`.

XXX Links to scipy.org here


ASE makes heavy use of an extension to Python called :term:`NumPy`.  The ``Numeric`` module defines an ``array`` type that can
hold large arrays of uniform multidimensional numeric data.  An
``array`` is similar to a ``list`` or a ``tuple``, but it is a lot
more powerful and efficient.


>>> import Numeric as num
>>> a = num.zeros((3, 2), num.Float)
>>> a[:, 1] = 1
>>> a[1] = 2
>>> a
array([[ 0.,  1.],
       [ 2.,  2.],
       [ 0.,  1.]])
>>> a.shape
(3, 2)

ase import *: np, array,zeros, ase.core import *
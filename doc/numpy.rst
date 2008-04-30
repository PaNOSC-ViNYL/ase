.. _numpy:

Numeric arrays in Python
========================

:term:`ndarray`.

XXX Links to scipy.org here


ASE makes heavy use of an extension to Python called NumPy_.  The
NumPy module defines an ``ndarray`` type that can hold large arrays of
uniform multidimensional numeric data.  An ``array`` is similar to a
``list`` or a ``tuple``, but it is a lot more powerful and efficient.


>>> import numpy as np
>>> a = np.zeros((3, 2))
>>> a[:, 1] = 1
>>> a[1] = 2
>>> a
array([[ 0.,  1.],
       [ 2.,  2.],
       [ 0.,  1.]])
>>> a.shape
(3, 2)

XXX ``from ase import *``: np, array, zeros

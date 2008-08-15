.. _numpy:

Numeric arrays in Python
========================

Links to NumPy's webpage:

* `Numpy example list`_
* `Numpy functions by category`_
* http://mentat.za.net/numpy/refguide
* `Numpy guide book <http://www.tramy.us/>`_

.. _Numpy example list: http://www.scipy.org/Numpy_Example_List_With_Doc
.. _Numpy functions by category:
                        http://www.scipy.org/Numpy_Functions_by_Category


ASE makes heavy use of an extension to Python called NumPy.  The
NumPy module defines an :term:`ndarray` type that can hold large arrays of
uniform multidimensional numeric data.  An array is similar to a
``list`` or a ``tuple``, but it is a lot more powerful and efficient.

XXX More examples from everyday ASE-life here ...

>>> import numpy as np
>>> a = np.zeros((3, 2))
>>> a[:, 1] = 1.0
>>> a[1] = 2.0
>>> a
array([[ 0.,  1.],
       [ 2.,  2.],
       [ 0.,  1.]])
>>> a.shape
(3, 2)
>>> a.ndim
2


All functions from the ``numpy`` module's namespace are available in
the ASE namespace as well - thus, ``from ase import *`` will import
``array``, ``zeros`` and so on.  The numpy module will furthermore be
aliased as ``np``.

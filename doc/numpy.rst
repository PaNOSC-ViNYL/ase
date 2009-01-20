.. _numpy:

Numeric arrays in Python
========================

Links to NumPy's webpage:

* `Numpy and Scipy Documentation`_
* `Numpy guide book <http://www.tramy.us/numpybook.pdf>`_
* `Numpy example list`_
* `Numpy functions by category`_


.. _Numpy and Scipy Documentation: http://docs.scipy.org/doc
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


.. note::

  When you do ``from ase import *``, you will get the ``np`` alias
  automatically, so there is no need to do a ``import numpy as np``
  also.  Try this:

  >>> from ase import *
  >>> dir(np)

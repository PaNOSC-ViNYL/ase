------
Python
------

This section will give a very brief introduction to the Python
language.  For more information on the Python language, go to the
`Python homepage`_.

.. _Python homepage: http://www.python.org



Executing Python code
---------------------

You can execute Python code interactively by starting the interpreter
like this::

  [hugo@cpu1 hugo]$ python
  Python 2.2.1 (#1, Aug 30 2002, 12:15:30)
  [GCC 3.2 20020822 (Red Hat Linux Rawhide 3.2-4)] on linux2
  Type "help", "copyright", "credits" or "license" for more information.
  >>> print 'hello'
  hello

You can also put the ``print 'hello'`` line in a file (``hello.py``)
and execute it as a Python script::

  [hugo@cpu1 hugo]$ python hello.py
  hello

Or like this::

  [hugo@cpu1 hugo]$ python -i hello.py
  hello
  >>> print 'hi!'
  hi!

Finally, you can put ``#!/usr/bin/env python`` in the first line of
the ``hello.py`` file, make it executable (``chmod +x hello.py``) and
execute it like any other executable.





Types
-----

Python has the following predefined types:

===========  =====================  ==========================
type         description            example
===========  =====================  ==========================
``bool``     boolean                ``False``
``int``       integer                ``117``
``float``    floating point number  ``1.78``
``complex``  complex number         ``0.5 + 2.0j``
``str``      string                 ``'abc'``
``tuple``    tuple                  ``(1, 'hmm', 2.0)``
``list``     list                   ``[1, 'hmm', 2.0]``
``dict``     dictionary             ``{'a': 7.0, 23: True}``
``file``     file                   ``file('stuff.dat', 'w')``
===========  =====================  ==========================

A ``dict`` object is mapping from keys to values:

>>> d = {'s': 0, 'p': 1}
>>> d['d'] = 2
>>> d
{'p': 1, 's': 0, 'd': 2}
>>> d['p']
1

A ``list`` object is an ordered collection of arbitrary objects:

>>> l = [1, ('gg', 7), 'hmm']
>>> l[1]
('gg', 7)
>>> l.append(1.2)
>>> l[-2]
'hmm'

A ``tuple`` behaves like a ``list`` - except that it can't be modified
inplace.  Objects of types ``list`` and ``dict`` are *mutable* - all
the other types listed in the table are *immutable*, which means that
once an object has been created, it can not change.

.. note::
   List and dictionary objects *can* change.  Variables in
   Python are references to objects.  This is demonstrated here:

   >>> a = ['q', 'w']
   >>> b = a
   >>> a.append('e')
   >>> a
   ['q', 'w', 'e']
   >>> b
   ['q', 'w', 'e']



Numeric Python
--------------

ASE makes heavy use of an extension to Python called `Numeric
Python`_.  The ``Numeric`` module defines an ``array`` type that can
hold large arrays of uniform multidimensional numeric data.  An
``array`` is similar to a ``list`` or a ``tuple``, but it is a lot
more powerful and efficient.

.. _Numeric:
.. _Numeric Python: http://numpy.sf.net

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



Loops
-----

A loop in Python can be done like this:

>>> things = ['a', 7]
>>> for x in things:
...     print x
...
a
7

The ``things`` object could be any sequence.  Strings, tuples, lists,
dictionaries, Numeric arrays and files are sequences.






Functions and classes
---------------------

A function is defined like this:

>>> def f(x, y):
...     return x + 2 * x * y
...
>>> f(1, 2)
5

A class is defined like this:

>>> class C:
...     def __init__(self, x):
...         self.x = x
...     def M(self, y):
...         return f(self.x, y)
...

The ``__init__()`` function is called a *constructor*.  You can think
of a class as a template for creating user defined objects:

>>> o = C(1)
>>> o.M(2)
5

Here we just called the method ``M`` of the object ``o`` (``o`` is an
instance of the class ``C``).






Importing modules
-----------------

If you put the definitions of the function ``f`` and the class ``C``
in a file ``stuff.py``, then you can use that code from another piece
of code::

  from stuff import f, C
  print f(1, 2)
  print C(1).M(2)

or::

  import stuff
  print stuff.f(1, 2)
  print stuff.C(1).M(2)

or::

  import stuff as st
  print st.f(1, 2)
  print st.C(1).M(2)

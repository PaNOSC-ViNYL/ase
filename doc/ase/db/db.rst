.. module:: db

====================
A database for atoms
====================

ASE has its own database that can be used for storing and retrieving atoms and
associated data in a compact and convenient way.  
    
.. note::

    This is work in progress.  Use at your own risk!
    
There are currently three back-ends:

JSON_:
    Simple human-readable text file with a ``.json`` extension.
SQLite3_:
    Self-contained, server-less, zero-configuration database.  Lives in a file
    with a ``.db`` extension.
PostgreSQL_:
    Server based database.

The JSON and SQLite3 back-ends work "out of the box", whereas PostgreSQL
requires a server.

There is a `command-line tool`_ called :program:`ase-db` that can be
used to query and manipulate databases and also a `Python interface`_.

.. _JSON: http://www.json.org/
.. _SQLite3: http://www.sqlite.org/
.. _PostgreSQL: http://www.postgresql.org/

.. contents::
    

What's in the database?
=======================

Every row in the database contains:

* all the information stored in the :class:`~ase.atoms.Atoms` object
  (positions, atomic numbers, ...)
* calculator name and parameters (if a calculator is present)
* already calculated properties such as energy and forces
  (if a calculator is present)
* keywords and key-value pairs (for finding the calculation again)
* an integer ID (unique for each database) starting with 1 and always
  increasing for each new row
* a unique ID which is a 128 bit random number which should be globally
  unique (at least in the lifetime of our universe)


Command-line tool
=================

The :program:`ase-db` command can be used to query databases and for manipulating keywords and key-value pairs.  Try::
    
    $ ase-db --help
    
Example: Show all rows of SQLite database abc.db::
 
    $ ase-db abc.db 
    id|age|user |formula|calc|energy| fmax|pbc|keywords|keyvals      | mass
     1|6s |jensj|H2     |emt | 1.419|9.803|000|molecule|relaxed=False|2.016
     2|5s |jensj|H2     |emt | 1.071|0.000|000|molecule|relaxed=True |2.016
     3|5s |jensj|H      |emt | 3.210|0.000|000|        |             |1.008
    
Show all details for a single row::
    
    $ ase-db abc.db id=3 -l
    id: 3
    formula: H
    user: jensj
    age: 174s
    calculator: emt
    energy: 3.210 eV
    maximum atomic force: 0.000 eV/Ang
    magnetic moment: 0
    periodic boundary conditions: [False False False]
    unit cell [Ang]:
         1.000     0.000     0.000
         0.000     1.000     0.000
         0.000     0.000     1.000
    volume: 1.000 Ang^3
    mass: 1.008 au
                       
Here are some example query strings:
    
.. list-table::
    :widths: 25 75
    
    * - v3
      - has 'v3' keyword
    * - abc=H
      - has key 'abc' with value 'H'
    * - v3,abc=H
      - both of the above
    * - calculator=nwchem
      - calculations done with NWChem
    * - 2.2<bandgap<3.0
      - 'bandgap' key has value between 2.2 and 3.0
    * - natoms>=10
      - 10 or more atoms
    * - H<3
      - less than 3 hydrogen atoms
    * - id=2345
      - specific id
    * - age<1h
      - not older than 1 hour
    * - age>1y
      - older than 1 year

These names are special:

.. list-table::
    :widths: 25 75
    
    * - id
      - integer identifier
    * - natoms
      - number of atoms
    * - energy
      - potential energy
    * - charge
      - total charge
    * - magmom
      - total magnetic moment
    * - calculator
      - name of calculator
    * - user
      - who did it
    * - age
      - age of calculation (use s, m, h, d, w, M and y for second, minute,
        hour, day, week, month and year respectively)

        
Integration with other parts of ASE
===================================

ASE's :func:`ase.io.read` function can also read directly from databases:
    
>>> from ase.io import read
>>> a = read('abc.db@42')
>>> a = read('abc.db@id=42')  # same thing
>>> b = read('abc.db@v3,abc=H')

Also the :program:`ase-gui` program can read from databases using the
same syntax.
        

Python Interface
================

.. module:: ase.db.core

First, we :func:`~ase.db.connect` to the database:
    
>>> from ase.db import connect
>>> c = connect('abc.db')

or

>>> import ase.db
>>> c = ase.db.connect('abc.db')

Do a calculation for a hydrogen molecule:
    
>>> from ase import Atoms
>>> from ase.calculators.emt import EMT
>>> h2 = Atoms('H2', [(0, 0, 0), (0, 0, 0.7)])
>>> h2.calc = EMT()
>>> h2.get_forces()
array([[ 0.        ,  0.        , -9.80290573],
       [ 0.        ,  0.        ,  9.80290573]])

Write a row to the database with keyword 'molecule' and a key-value pair
('relaxed', ``False``):
    
>>> c.write(h2, ['molecule'], relaxed=False)
1

The :meth:`~Database.write` method returns an integer id.

Do one more calculation and write results:
    
>>> from ase.optimize import BFGS
>>> BFGS(h2).run(fmax=0.01)
BFGS:   0  12:49:25        1.419427       9.8029
BFGS:   1  12:49:25        1.070582       0.0853
BFGS:   2  12:49:25        1.070544       0.0236
BFGS:   3  12:49:25        1.070541       0.0001
>>> c.write(h2, ['molecule'], relaxed=True)
2

Loop over selected rows using the :meth:`~Database.select` method:
    
>>> for d in c.select('molecule'):
...     print d.forces[0, 2], d.relaxed
... 
-9.8029057329 False
-9.2526347333e-05 True

The :meth:`~Database.select` method will generate dictionaries that one can
loop over.  The dictionaries are special in the sense that keys can be
accessed as attributes also (``d.relaxed == d['relaxed']``).

Write the energy of an isolated hydrogen atom to the database:

>>> h = Atoms('H')
>>> h.calc = EMT()
>>> h.get_potential_energy()
3.21
>>> c.write(h)
3

Select a single row with the :meth:`~Database.get` method:
    
>>> d = c.get(relaxed=1, calculator='emt')
>>> for k, v in d.items():
...     print '%-25s: %s' % (k, v)
... 
user                     : jensj
key_value_pairs          : {u'relaxed': True}
timestamp                : 14.0195850909
energy                   : 1.07054126233
relaxed                  : True
calculator_parameters    : {}
cell                     : [[ 1.  0.  0.] [ 0.  1.  0.] [ 0.  0.  1.]]
numbers                  : [1 1]
forces                   : [[ ... ]]
positions                : [[ ... ]]
keywords                 : [u'molecule']
pbc                      : [False False False]
id                       : 2
unique_id                : 22e4a2d64800987dd8d398c6cc9fd085
calculator               : emt

Calculate the atomization energy and :meth:`~Database.update` a row in
the database:
    
>>> e2 = d.energy
>>> e1 = c.get(H=1).energy
>>> ae = 2 * e1 - e2
>>> print ae
5.34945873767
>>> id = c.get(relaxed=1).id
>>> c.update(id, atomization_energy=ae)
(0, 1)

Delete a single row:
    
>>> del c[c.get(relaxed=0).id]

or use the :meth:`~Database.delete` method to delete several rows.


More details
------------

.. autofunction:: ase.db.connect

.. autoclass:: ase.db.core.Database
    :members:
    :member-order: bysource

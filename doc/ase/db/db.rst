.. module:: ase.db

====================
A database for atoms
====================

ASE has its own database that can be used for storing and retrieving atoms and
associated data in a compact and convenient way.

There are currently three back-ends:

JSON_:
    Simple human-readable text file with a ``.json`` extension.
SQLite3_:
    Self-contained, server-less, zero-configuration database.  Lives in a file
    with a ``.db`` extension.
PostgreSQL_:
    Server based database.

The JSON and SQLite3 back-ends work "out of the box", whereas PostgreSQL
requires a :ref:`server`.

There is a command-line tool called :ref:`ase-db` that can be
used to query and manipulate databases and also a `Python interface`_.

.. _JSON: http://www.json.org/
.. _SQLite3: http://www.sqlite.org/
.. _PostgreSQL: http://www.postgresql.org/


.. contents::


What's in the database?
=======================

Every row in the database contains:

* all the information stored in the :class:`~ase.Atoms` object
  (positions, atomic numbers, ...)
* calculator name and parameters (if a calculator is present)
* already calculated properties such as energy and forces
  (if a calculator is present)
* key-value pairs (for finding the calculation again)
* an integer ID (unique for each database) starting with 1 and always
  increasing for each new row
* a unique ID which is a 128 bit random number which should be globally
  unique (at least in the lifetime of our universe)
* constraints (if present)
* user-name
* creation and modification time


.. _ase-db:

ase db
======

The :ref:`ase-db` command-line tool can be used to query databases and for
manipulating key-value pairs.  Try::

    $ ase db --help

Example: Show all rows of SQLite database abc.db:

.. literalinclude:: ase-db.txt

Show all details for a single row:

.. literalinclude:: ase-db-long.txt

.. seealso::

    * :ref:`cli`


Querying
--------

Here are some example query strings:

.. list-table::
    :widths: 25 75

    * - Cu
      - contains copper
    * - H<3
      - less than 3 hydrogen atoms
    * - Cu,H<3
      - contains copper and has less than 3 hydrogen atoms
    * - H2O
      - At least two hydrogens and at least one oxygen
    * - formula=H2O
      - Exactly two hydrogens and one oxygen
    * - v3
      - has 'v3' key
    * - abc=bla-bla
      - has key 'abc' with value 'bla-bla'
    * - v3,abc=bla-bla
      - both of the above
    * - calculator=nwchem
      - calculations done with NWChem
    * - 2.2<bandgap<3.0
      - 'bandgap' key has value between 2.2 and 3.0
    * - natoms>=10
      - 10 or more atoms
    * - id=2345
      - specific id
    * - age<1h
      - not older than 1 hour
    * - age>1y
      - older than 1 year
    * - pbc=TTT
      - Periodic boundary conditions along all three axes
    * - pbc=TTF
      - Periodic boundary conditions along the first two axes (F=False, T=True)

These names are special:

.. list-table::
    :widths: 25 75

    * - id
      - integer identifier
    * - natoms
      - number of atoms
    * - pbc
      - Periodic boundary conditions
    * - formula
      - formula
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

Also the :ref:`ase-gui` program can read from databases using the
same syntax.


.. _ase-db-web:

Browse database with your web-browser
=====================================

You can use your web-browser to look at and query databases like this::

    $ ase db abc.db -w
    $ firefox http://0.0.0.0:5000/

Click individual rows to see details.  See the CMR_ web-page for an example of
how this works.

.. _CMR: https://cmrdb.fysik.dtu.dk/


Python Interface
================

.. module:: ase.db.core

First, we :func:`connect` to the database:

>>> from ase.db import connect
>>> db = connect('abc.db')

or

>>> import ase.db
>>> db = ase.db.connect('abc.db')

Let's do a calculation for a hydrogen molecule and write some results to a
database:

>>> from ase import Atoms
>>> from ase.calculators.emt import EMT
>>> h2 = Atoms('H2', [(0, 0, 0), (0, 0, 0.7)])
>>> h2.calc = EMT()
>>> h2.get_forces()
array([[ 0.   ,  0.   , -9.803],
       [ 0.   ,  0.   ,  9.803]])

Write a row to the database with a key-value pair (``'relaxed'``, ``False``):

>>> db.write(h2, relaxed=False)
1

The :meth:`~Database.write` method returns an integer id.

Do one more calculation and write results:

>>> from ase.optimize import BFGS
>>> BFGS(h2).run(fmax=0.01)  # doctest: +SKIP
BFGS:   0  12:49:25        1.419427       9.8029
BFGS:   1  12:49:25        1.070582       0.0853
BFGS:   2  12:49:25        1.070544       0.0236
BFGS:   3  12:49:25        1.070541       0.0001
>>> db.write(h2, relaxed=True)
2

Loop over selected rows using the :meth:`~Database.select` method:

>>> for row in db.select(relaxed=True):
...     print(row.forces[0, 2], row.relaxed)  # doctest: +SKIP
-9.8029057329 False
-9.2526347333e-05 True

The :meth:`~Database.select` method will generate :ref:`row objects`
that one can loop over.

Write the energy of an isolated hydrogen atom to the database:

>>> h = Atoms('H')
>>> h.calc = EMT()
>>> h.get_potential_energy()
3.21
>>> db.write(h)
3

Select a single row with the :meth:`~Database.get` method:

>>> row = db.get(relaxed=1, calculator='emt')
>>> for key in row:
...    print('{0:22}: {1}'.format(key, row[key]))  # doctest: +SKIP
pbc                   : [False False False]
relaxed               : True
calculator_parameters : {}
user                  : jensj
mtime                 : 15.3439399027
calculator            : emt
ctime                 : 15.3439399027
positions             : [[ ... ]]
id                    : 2
cell                  : [[ 0.  0.  0.] [ 0.  0.  0.] [ 0.  0.  0.]]
forces                : [[ ... ]]
energy                : 1.07054126233
unique_id             : bce90ff3ea7661690b54f9794c1d7ef6
numbers               : [1 1]

Calculate the atomization energy and :meth:`~Database.update` a row in
the database:

>>> e2 = row.energy
>>> e1 = db.get(H=1).energy
>>> ae = 2 * e1 - e2
>>> print(ae)
5.34945873767
>>> id = db.get(relaxed=1).id
>>> db.update(id, atomization_energy=ae)
1

Delete a single row:

>>> del db[db.get(relaxed=0).id]

or use the :meth:`~Database.delete` method to delete several rows.


Description of a row
--------------------

The first 9 keys (from "id" to "positions") are always present --- the rest
may be there:

=====================  =================================  ============  ======
key                    description                        datatype      shape
=====================  =================================  ============  ======
id                     Local database id                  int
unique_id              Globally unique hexadecimal id     str
ctime                  Creation time                      float
mtime                  Modification time                  float
user                   User name                          str
numbers                Atomic numbers                     int           (N,)
pbc                    Periodic boundary condition flags  bool          (3,)
cell                   Unit cell                          float         (3, 3)
positions              Atomic positions                   float         (N, 3)
initial_magmoms        Initial atomic magnetic moments    float         (N,)
initial_charges        Initial atomic charges             float         (N,)
masses                 Atomic masses                      float         (N,)
tags                   Tags                               int           (N,)
momenta                Atomic momenta                     float         (N, 3)
constraints            Constraints                        list of dict
energy                 Total energy                       float
forces                 Atomic forces                      float         (N, 3)
stress                 Stress tensor                      float         (6,)
dipole                 Electrical dipole                  float         (3,)
charges                Atomic charges                     float         (N,)
magmom                 Magnetic moment                    float
magmoms                Atomic magnetic moments            float         (N,)
calculator             Calculator name                    str
calculator_parameters  Calculator parameters              dict
=====================  =================================  ============  ======


Extracting Atoms objects from the database
------------------------------------------

If you want an :class:`~ase.Atoms` object insted of an
:class:`~ase.db.row.AtomsRow` object, you should use the
:meth:`~Database.get_atoms` method:

>>> h2 = db.get_atoms(H=2)

or if you want the original EMT calculator attached:

>>> h2 = db.get_atoms(H=2, attach_calculator=True)


Add additional data
-------------------

When you write a row to a database using the :meth:`~Database.write` method,
you can add key-value pairs where the values can be
strings, floating point numbers, integers and booleans:

>>> db.write(atoms, functional='LDA', distance=7.2)

More complicated data can be written like this:

>>> db.write(atoms, ..., data={'parents': [7, 34, 14], 'stuff': ...})

and accessed like this:

>>> row = db.get(...)
>>> row.data.parents
[7, 34, 14]


.. _row objects:

Row objects
-----------

There are three ways to get at the columns of a row:

1) as attributes (``row.key``)

2) indexing (``row['key']``)

3) the :meth:`~ase.db.row.AtomsRow.get` method (``row.get('key')``)

The first two will fail if there is no ``key`` column whereas the last will
just return ``None`` in that case.  Use ``row.get('key', ...)`` to use
another default value.

.. autoclass:: ase.db.row.AtomsRow
    :members:
    :member-order: bysource


Writing and updating many rows efficiently
------------------------------------------

If you do this::

    db = connect('mols.db')
    for mol in molecules:
        db.write(mol, ...)

the database will make sure that each molecule is written to permanent
starage (typically a harddisk) before it moves on to the next molecule.  This
can be quite slow.  To speed this up, you can write all the molecules in a
single transaction like this::

    with connect('mols.db') as db:
        for mol in molecules:
            db.write(mol, ...)

When the for-loop is done, the database will commit (or roll back if there
was an error) the transaction.

Similarly, the :meth:`~Database.update` method will do up to
``block_size=1000`` rows in one transaction::

    # slow:
    for row in db.select(...):
        db.update(row.id, foo='bar')  # a single id
    # faster:
    ids = [row.id for row in db.select(...)]
    db.update(ids, foo='bar')  # list of id's


Writing rows in parallel
------------------------

Say you want to run a series of jobs and store the calculations in one
database::

    for name in many_molecules:
        mol = read(name)
        calculate_something(mol)
        db.write(mol, name=name)

With four extra lines (see the :meth:`~Database.reserve` method)::

    for name in many_molecules:
        id = db.reserve(name=name)
        if id is None:
            continue
        mol = read(name)
        calculate_something(mol)
        db.write(mol, name=name)
        del db[id]

you will be able to run several jobs in parallel without worrying about two
jobs trying to do the same calculation.  The :meth:`~Database.reserve` method
will write an empty row with the ``name`` key and return the ID of that row.
Other jobs trying to make the same reservation will fail.  While the jobs are
running, you can keep an eye on the ongoing (reserved) calculations by
identifying empty rows::

    $ ase db many_results.db natoms=0


More details
------------

Use this function for getting a connection to a database:

.. autofunction:: connect

Here is a description of the database object:

.. autoclass:: ase.db.core.Database
    :members:
    :member-order: bysource
    :exclude-members: write, reserve, update

    .. decorators hide these three from Sphinx, so we add them by hand:

    .. automethod:: write(atoms, key_value_pairs={}, data={}, **kwargs)
    .. automethod:: reserve(**key_value_pairs)
    .. automethod:: update(ids, delete_keys=[], block_size=1000, **add_key_value_pairs)

    .. attribute:: metadata

       Dictionary


.. _metadata:

Metadata
--------

It's can be useful to add a discription of the data that's in a database.
This description will be used when browsing the database with a
:ref:`web-browser <ase-db-web>`.
You can add the desciption using the :attr:`~Database.metadata` attribute:

>>> db.metadata = {
...     'title': 'Project 1',
...     'key_descriptions':
...         {'v0': ('Voltage', 'Longer description ...', 'V')},
...     'default_columns': ['id', 'formula', 'v0']}

ASE already knows all about the following keys:

.. csv-table::
    :file: known-keys.csv
    :header-rows: 1
    :widths: 2 3 4 2

You can also write/read to/from JSON using::

    $ ase db proj1.db --set-metadata metadata.json
    $ ase db proj1.db --show-metadata > metadata.json


.. _server:

Running a PostgreSQL server
===========================

.. highlight:: bash

With your PostgreSQL server up and running, you should run the following
command as the ``postgres`` user::

    $ sudo -u postgres psql

Then you create an 'ase' user and one database for each project you have::

    postgres=# create user ase login password 'pw';
    postgres=# create database project1;
    postgres=# create database project2;

Show databases and quit::

    postgres=# \l
    postgres=# \q

You should now be able to
query the database using a URI like
``postgresql://user:pw@host:port/dbname``::

    $ PROJ1=postgresql://ase:pw@localhost:5432/project1
    $ ase db $PROJ1

If you have some data in a ``data.db`` SQLite3 file, then you can insert that
into the PostgreSQL database like this::

    $ ase db data.db --insert-into $PROJ1

Now you can start the Flask_\ -app ``ase.db.app``.  You can use Flask's own
web-server or use any WSGI_ compatible server.  We will use
Twisted_ in the example below. Set the $ASE_DB_APP_CONFIG environment variable
to point to a Python configuration file containing something similar to this::

    ASE_DB_NAMES = ['postgresql://ase:pw@localhost:5432/project1',
                    'postgresql://ase:pw@localhost:5432/project2',
                    ...]
    ASE_DB_HOMEPAGE = '<a href="https://home.page.org">HOME</a>'

and then start the server with::

    $ ASE_DB_APP_CONFIG=... twistd web --wsgi=ase.db.app.app --port=tcp:8000

.. note::

    Please review the code carefully before exposing the ``ase.db.app`` to
    the internet or `bad things <https://xkcd.com/327/>`__ could happen.

.. _Flask: http://flask.pocoo.org/
.. _WSGI: https://www.python.org/dev/peps/pep-3333/
.. _Twisted: https://twistedmatrix.com/

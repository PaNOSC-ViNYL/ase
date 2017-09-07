.. _db tutorial:

======================
Using the ASE database
======================

In this tutorial we will adsorb C, N and O on 7 different FCC(111) surfaces
with 1, 2 and 3 layers and we will use database files to store the results.


Bulk
----

First, we equilibium bolk FCC lattice constants for the seven elements where
the :mod:`emt` potential works well:

.. literalinclude:: bulk.py

.. highlight:: bash

Run the :download:`bulk.py` script and look at the results::

    $ python3 bulk.py
    $ ase db bulk.db -c +bm  # show also the bulk-modulus column
    ...
    $ ase gui bulk.db

The `bulk.db` is a SQLite database in a single file::

    $ file bulk.db
    bulk.db: SQLite 3.x database

If you want to see what's inside you can convert the database file to a json
file and open that in your text editor::

    $ ase db bulk.db --insert-into bulk.json
    ??????????

or, you can look at a single row like this::

    $ ase db bulk.db Cu -j > cu.json

The json file format is human readable, but much less efficient to work with
compared to a SQLite file.


Adsorbates
----------

Now we do the adsorbtion calculations (run the :download:`ads.py` script).

.. literalinclude:: ads.py

We now have a new database file with 63 rows::

    $ ase db ads.db -n
    63 rows

These 63 calculation only take a few seconds with EMT.  Suppose you want to use
DFT and send the calculations to a supercomputer. In that case you may want to
run several calculations in different jobs on the supercomputer.  In addition,
some of the jobs could time out and not finish.  It's a good idea to modify
the script a bit for this scenario.  We replace the `run` and `write` lines
in the inner loop by this:

.. highlight:: python

::

    id = db2.reserve(layers=n, surf=symb, ads=ads)
    if id is not None:
        atoms = run(symb, a, n, ads)
        db2.write(atoms, layers=n, surf=symb, ads=ads)
        del db2[id]

The :meth:`~ase.db.Database.reserve` method will check if there is a row with
the keys `layers?n`, `surf=symb` and `ads=ads`.  If there is, then the
calculation will be skipped.  If there is not, then an empty row with those
keys-values will be written and the calculation will start.  When done, the
real row will be written and the empty one will be removed.  This modified
script can run several jobs all running in parallel and no calculation will
be done twice.

.. highlight:: bash

In case a calculation crashes, you will see empty rows in the database::

    $ ase db ads.db natoms=0

Delete them, fix the problem and run the script again::

    $ ase db ads.db natoms=0 --delete
    ???????
    $ python ads.py  # or sbatch ...
    $ ase db ads.db natoms=0


Reference energies
------------------

Let's also calculate the energy of the clean surfaces and the isolated
adsorbates:

.. literalinclude:: refs.py

::

    $ python refs.py
    $ ase db ads.db -n
    87 rows

Say we want those 24 reference energies (clean surfaces and isolated adsorbates) in a `refs.db` file instead of the big `ads.db` file.  We could change the
`refs.py` script and run the calculations again, but we can also manipulate the files using the `ase db` tool like this::

    $ ase db -L 24 --sort -id --insert-into refs.db
    $ ase db -L 24 --sort -id --delete --yes
    $ ase db ads.db -n
    63 rows
    $ ase db refs.db -n
    24 rows


Analysis
--------

Now we have what we need to calulate the adsorbiton energies and heights:

.. literalinclude:: ea.py

Here are the results for three layers of Pt::

    $ ase db ads.db Pt,layers=3 -c formula,ea,height

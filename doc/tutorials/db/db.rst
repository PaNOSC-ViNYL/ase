.. _db tutorial:

======================
Using the ASE database
======================

In this tutorial we will adsorb C, N and O on 7 different FCC(111) surfaces
with 1, 2 and 3 layers and we will use database files to store the results.

.. seealso::

    The :mod:`ase.db` module documentation.


|cu1o| |cu2o| |cu3o|

.. |cu1o| image:: cu1o.png
.. |cu2o| image:: cu2o.png
.. |cu3o| image:: cu3o.png


Bulk
----

First, we calculate the equilibrium bulk FCC lattice constants for the seven
elements where the :mod:`EMT <ase.calculators.emt>` potential works well:

.. literalinclude:: bulk.py

.. highlight:: bash

Run the :download:`bulk.py` script and look at the results::

    $ python3 bulk.py
    $ ase db bulk.db -c +bm  # show also the bulk-modulus column
    id|age|formula|calculator|energy| fmax|pbc|volume|charge|   mass|   bm
     1|10s|Al     |emt       |-0.005|0.000|TTT|15.932| 0.000| 26.982|0.249
     2| 9s|Ni     |emt       |-0.013|0.000|TTT|10.601| 0.000| 58.693|1.105
     3| 9s|Cu     |emt       |-0.007|0.000|TTT|11.565| 0.000| 63.546|0.839
     4| 9s|Pd     |emt       |-0.000|0.000|TTT|14.588| 0.000|106.420|1.118
     5| 9s|Ag     |emt       |-0.000|0.000|TTT|16.775| 0.000|107.868|0.625
     6| 9s|Pt     |emt       |-0.000|0.000|TTT|15.080| 0.000|195.084|1.736
     7| 9s|Au     |emt       |-0.000|0.000|TTT|16.684| 0.000|196.967|1.085
    Rows: 7
    Keys: bm
    $ ase gui bulk.db

The :file:`bulk.db` is an SQLite3_ database in a single file::

    $ file bulk.db
    bulk.db: SQLite 3.x database

.. _SQLite3: http://www.sqlite.org/

If you want to see what's inside you can convert the database file to a json
file and open that in your text editor::

    $ ase db bulk.db --insert-into bulk.json
    Added 0 key-value pairs (0 pairs updated)
    Inserted 7 rows

or, you can look at a single row like this::

    $ ase db bulk.db Cu -j
    {"1": {
     "calculator": "emt",
     "energy": -0.007036492048371201,
     "forces": [[0.0, 0.0, 0.0]],
     "key_value_pairs": {"bm": 0.8392875566787444},
     ...
     ...
    }

The json file format is human readable, but much less efficient to work with
compared to a SQLite3 file.


Adsorbates
----------

Now we do the adsorption calculations (run the :download:`ads.py` script).

.. literalinclude:: ads.py

We now have a new database file with 63 rows::

    $ ase db ads.db -n
    63 rows

These 63 calculations only take a few seconds with EMT.  Suppose you want to
use DFT and send the calculations to a supercomputer. In that case you may
want to run several calculations in different jobs on the computer.  In
addition, some of the jobs could time out and not finish.  It's a good idea
to modify the script a bit for this scenario.  We add a couple of lines to
the inner loop:

.. highlight:: python

::

    for row in db1.select():
        a = row.cell[0, 1] * 2
        symb = row.symbols[0]
        for n in [1, 2, 3]:
            for ads in 'CNO':
                id = db2.reserve(layers=n, surf=symb, ads=ads)
                if id is not None:
                    atoms = run(symb, a, n, ads)
                    db2.write(atoms, layers=n, surf=symb, ads=ads)
                    del db2[id]

The :meth:`~ase.db.core.Database.reserve` method will check if there is a row
with the keys ``layers=n``, ``surf=symb`` and ``ads=ads``.  If there is, then
the calculation will be skipped.  If there is not, then an empty row with
those keys-values will be written and the calculation will start.  When done,
the real row will be written and the empty one will be removed.  This
modified script can run in several jobs all running in parallel and no
calculation will be done twice.

.. highlight:: bash

In case a calculation crashes, you will see empty rows left in the database::

    $ ase db ads.db natoms=0 -c ++
    id|age|user |formula|pbc|charge| mass|ads|layers|surf
    17|31s|jensj|       |FFF| 0.000|0.000|  N|     1|  Cu
    Rows: 1
    Keys: ads, layers, surf

Delete them, fix the problem and run the script again::

    $ ase db ads.db natoms=0 --delete
    Delete 1 row? (yes/No): yes
    Deleted 1 row
    $ python ads.py  # or sbatch ...
    $ ase db ads.db natoms=0
    Rows: 0


Reference energies
------------------

Let's also calculate the energy of the clean surfaces and the isolated
adsorbates (:download:`refs.py`):

.. literalinclude:: refs.py

::

    $ python refs.py
    $ ase db ads.db -n
    87 rows

Say we want those 24 reference energies (clean surfaces and isolated
adsorbates) in a :file:`refs.db` file instead of the big :file:`ads.db` file.
We could change the :file:`refs.py` script and run the calculations again,
but we can also manipulate the files using the ``ase db`` tool.  First, we
move over the clean surfaces::

    $ ase db ads.db ads=clean --insert-into refs.db
    Added 0 key-value pairs (0 pairs updated)
    Inserted 21 rows
    $ ase db ads.db ads=clean --delete --yes
    Deleted 21 rows

and then the three atoms (``pbc=FFF``, no periodicity)::

    $ ase db ads.db pbc=FFF --insert-into refs.db
    Added 0 key-value pairs (0 pairs updated)
    Inserted 3 rows
    $ ase db ads.db pbc=FFF --delete --yes
    Deleted 3 rows
    $ ase db ads.db -n
    63 rows
    $ ase db refs.db -n
    24 rows


Analysis
--------

Now we have what we need to calculate the adsorption energies and heights
(:download:`ea.py`):

.. literalinclude:: ea.py

Here are the results for three layers of Pt::

    $ python3 ea.py
    $ ase db ads.db Pt,layers=3 -c formula,ea,height
    formula|    ea|height
    Pt3C   |-3.715| 1.504
    Pt3N   |-5.419| 1.534
    Pt3O   |-4.724| 1.706
    Rows: 3
    Keys: ads, ea, height, layers, surf

.. note::

    While the EMT description of Ni, Cu, Pd, Ag, Pt, Au and Al is OK, the
    parameters for C, N and O are not intended for real work!

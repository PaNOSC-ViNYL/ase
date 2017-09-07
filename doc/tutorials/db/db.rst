.. _db tutorial:

======================
Using the ASE database
======================

In this tutorial we will adsorb C, N and O on 7 different FCC(111) surfaces
with 1, 2 and 3 layers and we will use database files to store the results.

.. literalinclude:: bulk.py
.. literalinclude:: ads.py
.. literalinclude:: refs.py
.. literalinclude:: ea.py

    $ ase gui Ag.traj
ase gui ads.db@CuN -r2,2,1
ase gui ads.db -r2,2,1
ase db ads.db
ase db ads.db id=33 --delete -y
gpaw -T run ads.db@AlC,layers=1 -p kpts=4,4,1 -w alc.gpw --modify="atoms.center(vacuum=4, axis=2)"
ase band-structure -k GKMG alc.gpw
ase db ads.db ads=none -i refs.db
ase db ads.db ads=none --delete
ase db ads.db "age<5m" -i refs.db
ase db ads.db "age<5m" --delete -y
ase db ads.db -c formula,ea,h
ase db ads.db AlC -c formula,ea,h
ase build Al -xfcc|gpaw run -p kpts=11,11,11,mode=pw -w al.gpw
gpaw dos al.gpw -p -w0

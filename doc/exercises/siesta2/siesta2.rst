.. _siesta2:

----------------------------
Siesta 2: Molecular Dynamics
----------------------------

This exercise is intended to illustrate a Molecular Dynamics run with
the SIESTA calculator. The system is a Si(001) surface, in the 2x1
reconstruction with asymmetric dimers. The simulation cell contains
two dimers. An H2 molecule approaches the surface, above one of the
dimers, and dissociates, ending up with a H atom bonded to each of the
Si atoms in the dimer (and thus leading to a symmetric dimer). You can
get the ``xyz`` file with the initial geometry :svn:`here
<doc/exercises/siesta2/geom.xyz?format=raw>`.

.. literalinclude:: siesta2.py

Note that both H atoms are given an initial velocity towards the
surface through the lines::

  p = atoms.get_momenta()
  p[0,2]= -1.5 
  p[1,2]= -1.5 
  atoms.set_momenta(p)

Run the program, and check the results. You can visualize the dynamics
using the trajectory file with the help of the ASE :mod:`gui`.




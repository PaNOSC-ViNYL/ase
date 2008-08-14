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
get the ``xyz`` file with the initial geometry :svn:`geom.xyz`.

.. literalinclude:: siesta2.py

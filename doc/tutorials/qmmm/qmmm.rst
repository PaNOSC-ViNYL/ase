.. _qmmm:

=========================
ASE for QM/MM Simulations
=========================

QM/MM Simulations couple two (or more) descriptions to get total energy
and forces for the entire system in an efficiant manner. The method paper 
on our implementation is currently being written. 
General background can be found `here <https://link.springer.com/article/10.1007/s00214-006-0143-z/>`_, 
`here <http://onlinelibrary.wiley.com/doi/10.1002/anie.200802019/abstract>`_, and
`here <https://www.elsevier.com/books/combining-quantum-mechanics-and-molecular-mechanics-some-recent-progresses-in-qm-mm-methods/sabin/978-0-12-380898-1>`_ . 
Examples of what this code has been used for can be seen `here <http://pubs.acs.org/doi/abs/10.1021/jz500850s>`_, 
and `here <http://pubs.acs.org/doi/abs/10.1021/acs.inorgchem.6b01840>`_. 

This section will show you how to setup up various QM/MM simulations.
We will be using GPAW_ for the QM part. Other QM calculators should
be straightforwardly compatible with the subtractive-scheme SimpleQMMM
calculator, but for the Excplicit Interaction EIQMMM calculator, you
would need to be able to put an electrostatic external potential into
the calculator for the QM subsystem. 

.. _GPAW: http://wiki.fysik.dtu.dk/gpaw


You might also be interested in the solvent MM potentials included in ASE. 
The tutorial on :ref:`TIPnP Water Box Equillibration` could be relevant to have a look at. 


Electrostatic Embedding QM/MM
-----------------------------

Let's first do a very simple electrostatic embedding QM/MM single point energy calculation
on the water dimer. The necessary inputs are described in the class help text:

.. autoclass:: ase.calculators.qmmm.EIQMMM


The following script will calculate the QM/MM single point energy of the water dimer from
the :ref:`s22`, using LDA and TIP3P. 

.. literalinclude:: water_dimer.py

Here, we have just used the TIP3P LJ parameters for the QM part as well. If this is a good
idea or not isn't trivial. The LJInteractions module needs combined parameters for all 
possible permutations of atom types in your system, that have LJ parameters. A list of 
combination rules can be found `here <http://www.sklogwiki.org/SklogWiki/index.php/Combining_rules>`_.
Here's a code snippet of how to combine LJ parameters of atom types A and B via the Lorentz-Berthelot rules::

   import itertools as it

   parameters = {'A': (epsAA, sigAA),
                 'B': (epsBB, sigBB)}

   def lorenz_berthelot(p):
       combined = {}
       for comb in it.product(p.keys(), repeat=2):
          combined[comb] = ((p[comb[0]][0] * p[comb[1]][0])**0.5,
                           (p[comb[0]][1] + p[comb[1]][1])/2)
       return combined

   combined = lorenz_berthelot(parameters)
   interaction = LJInteractions(combined)

This will (somewhat redundantly) yield::

    >>>combined
    {('A', 'A'): (epsAA, sigAA),
     ('A', 'B'): (epsAB, sigAB),
     ('B', 'A'): (epsAB, sigAB),
     ('B', 'B'): (epsBB, sigBB)}
 

Current limitations
    - No QM/MM border over bonds
    - No QM PBCs


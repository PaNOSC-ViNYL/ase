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


You might also be interested in the tutorial on the solvent MM potentials included in ASE. 
The tutorial on :ref:`TIPnP Water Box Equillibration` could be relevant to have a look at. 

Current limitations
    - No QM/MM border over bonds
    - No QM PBCs


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
The total energy expression for the full QM/MM system is:

.. math::  E_\mathrm{TOT} = E_\mathrm{QM} + E_\mathrm{I} + E_\mathrm{MM}.

The MM region is modelled using point charge force fields, with charges 
:math:`q_i` and :math:`\tau_i` denoting their spatial coordinates, so the 
QM/MM coupling term :math:`E_\mathrm{I}` will be

.. math:: E_\mathrm{I} = \sum_{i=1}^C q_i \int \frac{n({\bf r})}{\mid\!{\bf r} -
                         \tau_i\!\mid}\mathrm{d}{\bf r} + 
                         \sum_{i=1}^C\sum_{\alpha=1}^A 
                         \frac{q_i Z_{\alpha}}{\mid\!{\bf R}_\alpha - \tau_i\!\mid} + E_\mathrm{RD}

where :math:`n({\bf r})` is the spatial electronic density of the quantum 
region, :math:`Z_\alpha` and :math:`{\bf R}_\alpha` are the charge and coordinates 
of the nuclei in the QM region, respectively, and :math:`E_\mathrm{RD}` is the term describing 
the remaining, non-Coulomb interactions between the two subsystems.

For the MM point-charge external potential in GPAW, we use the 
total pseudo-charge density :math:`\tilde{\rho}({\bf r})` for the coupling, and since
the Coloumb integral is evaluated numerically on the real space grid, 
thus the coupling term ends up like this:

.. math:: E_\mathrm{I} = \sum_{i=1}^C q_i \sum_{g} \frac{\tilde{\rho}({\bf r})}{\mid\!{\bf r}_g  - \tau_i\!\mid} v_g + E_\mathrm{RD}

Currently, the term for :math:`E_{\mathrm{RD}}` implemented is a Lennard-Jones-type
potential:

.. math:: E_\mathrm{RD} = \sum_i^C \sum_\alpha^A 
                          4\epsilon\left[ \left(\frac{\sigma}{\mid\!{\bf R}_\alpha 
                          - \tau_i\!\mid}\right)^{12} 
                          - \left(\frac{\sigma}{\mid\!{\bf R}_\alpha 
                          - \tau_i\!\mid}\right)^{6} \right]

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


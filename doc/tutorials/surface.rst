=====================
Nitrogen on ruthenium
=====================

In this tutorial we calculate the adsorption energy of a nitrogen
molecule on a ruthenium surface. This is done by calculating the total
energy for the isolated slab and for the isolated molecule. The
adsorbate is then added to the slab and relaxed, and the total energy
for this composite system is calculated. The adsorption energy is
obtained as the sum of the isolated energies minus the energy of the
composite system.

.. literalinclude:: N2Cu.py

Here is a picture of the system after the relaxation::

  view(slab)  

.. image:: surface.png

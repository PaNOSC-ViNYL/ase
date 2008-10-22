.. module:: dft.dos
   :synopsis: Density of states

=================
Density of states
=================

Example::

  calc = ...
  dos = DOS(calc, width=0.2)
  d = dos.get_dos()
  e = dos.get_energies()

You can plot the result like this::

  import pylab as plt
  plt.plot(e, d)
  plt.xlabel('energy [eV]')
  plt.ylabel('DOS')
  plt.show()

.. image:: dos.png


More details
------------

.. autoclass:: ase.dft.dos.DOS
   :members: get_energies, get_dos


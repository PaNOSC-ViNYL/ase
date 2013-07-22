.. _stm-tutorial:

==============================
Tutorial: STM images - Al(100)
==============================

The STM is a revolutionary experimental surface probe that has
provided direct local insight into the surface electronic
structure. Sometimes the interpretation of STM topographs are not
straightforward and therefore theoretically modeled STM images may
resolve conflicting possibilities and point to an underlying atomistic
model. ASE includes python modules for generating
Tersoff-Hamann STM topographs. Use of the STM code is illustrated here for a
Al(100) surface.

Let's make a 2 layer Al(100) fcc surface using the :mod:`lattice` module::

  from ase.lattice.surface import fcc100
  atoms = fcc100('Al', size=(1, 1, 2))
  atoms.center(vacuum=4.0, axis=2)

Now a calculator must be defined, in this tutorial we will make a STM
image from the GPAW calculator.

For the GPAW code the calculator for the Al(100) surface can be
defined like this::

  from gpaw import GPAW
  calc = GPAW(mode='pw',
              kpts=(4, 4, 1),
	      txt='Al100.out')
  atoms.set_calculator(calc)
  energy = atoms.get_potential_energy() 
  calc.write('Al100.gpw', 'all')


2-d scans
=========

In this section we will make simulated STM scans and contour plot
using matplotlib. First initialize the :class:`STM` object and get the
averaged current at `z=8.0` Å (for our setup, the top layer is a
`z=6.025` Å)::

  from ase.dft.stm import STM
  from gpaw import GPAW

  calc = GPAW('Al100.gpw')
  atoms = calc.get_atoms()

  stm = STM(atoms, symmetries=[0, 1, 2])
  z = 8.0
  bias = 1.0
  c = stm.get_averaged_current(bias, z)

From the current we make a scan to get a 2-d array of constant current
heights::

  h = stm.scan(bias, c)

Finally we make a contour plot::

  import pylab as p
  p.contourf(h, 40)
  p.hot()
  p.colorbar()
  p.show()	
  
If you want to plot more than a single unit cell, use this::

  import numpy as np
  h = np.tile(h, (3, 3))


Linescans
=========

Here is how to make a line-scan::

  x, y = stm.linescan(bias, c, [0, 0], [2.5, 2.5])
  p.plot(x,y)
  p.show()

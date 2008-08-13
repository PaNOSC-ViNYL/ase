.. _stm-tutorial:

==============================
Tutorial: STM images - Al(100)
==============================

The STM is a revolutionary experimental surface probe that has
provided direct local insight into the surface electronic
structure. Sometimes the interpretation of STM topographs are not
straightforward and therefore theoretically modeled STM images may
resolve conflicting possibilities and point to an underlying atomistic
model. The CAMPOS code includes python modules for generating
Tersoff-Hamann STM topographs. The STM code is illustrated here for a
Al(100) in a 2x2 unit cell.

Let's make the Al(100) fcc surface by using the :mod:`lattice` module::

  from ase.lattice.surface import *
  atoms = fcc100('Al', size=(2,2,2))

Now a calculator must be defined, in this tutorial we will make a STM
image from the GPAW calculator.

For the GPAW code the calculator for the Al(100) surface can be
defined like this::

  from gpaw import GPAW
  calc = GPAW(gpts=(20,20,48),nbands=28,
  	kpts=(4,4,1),out='Al100.out')
  atoms.set_calculator(calc)
  energy = atoms.get_potential_energy() 
  calc.write('Al100.gpw')


3D Visualization
==========================

XXX

Now a ElectronicStates object must be defined from the output files. This can be done like this for the GPAW code::

    from ASE.Utilities.ElectronicStates import ElectronicStates
    electronicstates = ElectronicStates(filename = 'Al100.nc')

and for the dacapo code::

    from Dacapo.ElectronicStates import ElectronicStates
    electronicstates = ElectronicStates(filename='Al100.nc')

dacapo must import its own version of the ElectronicStates class, while the GPAW code can use the generic ASE version.

Now generate the STMTool object from the electronicstates object::

    from  ASE.Utilities.STMTool import STMTool
    stm = GetSTMTool(electronicstates,contourvalue=1.0e-6, channel="s", normalaxis=2)

The initialization is done with three keyword arguments (keyword=actualvalue) - more options are possible, see section STMTool in the manual.

The most important options are:

contourvalue
    XXX
channel
    XXX
normalaxis
    XXX
smoothfactor
    XXX

Linescans and matplotlib
==========================

In this section we will make simulated STM linescans and contour plot using matplotlib.

The example below make a contour plot and a linescan plot. The position of the linescan is plotted on the contour plot. The simulated linescans are made using the STMLineScan class:

XXX

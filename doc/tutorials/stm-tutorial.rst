==============================
Tutorial: STM images - Al(100)
==============================

The STM is a revolutionary experimental surface probe that has provided direct local insight into the surface electronic structure. Sometimes the interpretation of STM topographs are not straightforward and therefore theoretically modeled STM images may resolve conflicting possibilities and point to an underlying atomistic model. The CAMPOS code includes python modules for generating Tersoff-Hamann STM topographs. Using this functionality requires VTK installed. The STM code is illustrated here for a Al(100) in a 2x2 unit cell.

The setup for the atoms in the Al(100) in a 2x2 unitcell is described below.

After the import statements, the section starting with a0 = ... defines some major variables in the script ; this is not at all required, but increases the readability and reusability of the script. The part starting at atoms=Atoms(..) illustrates a three-fold python loop that iteratively builds up the slab in a layerwise fashion::

    from ase import *
    from numpy import sqrt,array
    a0     = 4.05         # cubic fcc lattice constant
    N      = 2            # repetition along x
    M      = 2            # repetition along y
    layers = 2            # slab layers
    electronsperatom = 3
    vaclay = 5            # interlayer dist = a0/2
    
    atoms   = Atoms([],pbc=True)
    for n in range(layers):
        for i in range(N):
            for j in range(M):
                scaledpos = [(i+(n%2)/2.)/sqrt(2.),(j+(n%2)/2.)/sqrt(2.),-n/2.]
                atoms.append(Atom('Al', a0*array(scaledpos)))
     
    unitcell = [[N/sqrt(2.), 0.0,        0.0],
               [0.0,        M/sqrt(2.), 0.0],
               [0.0,        0.0,        (vaclay+layers)/2.]]
    
    atoms.set_cell(a0*array(unitcell),fix=True)

Now a calculator must be defined, in this tutorial we will make a STM image from both the GPAW real space grid code and from the dacapo planewave code.

For the GPAW code the calculator for the Al(100) surface can be defined like this::

    from gpaw import Calculator
    calc = Calculator(gpts=(20,20,48),nbands=28,
                      kpts=(4,4,1),out='Al100.out')
    atoms.set_calculator(calc)
    energy = atoms.get_potential_energy() 
    calc.write('Al100.gpw')

The calculator for the dacapo code is similar::

    from Dacapo import Dacapo
    calc = Dacapo(planewavecutoff = 150,nbands = 28,
                  kpts=(4,4,1),xc = 'LDA',
                  usesymm=True,
                  out = 'Al100.nc',txtout = 'Al100.txt')
    atoms.SetCalculator(calc)
    energy = atoms.GetPotentialEnergy()


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
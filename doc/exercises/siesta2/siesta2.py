#!/usr/bin/env python
from ase import *

# Read in the geometry from a xyz file, set the cell and center
atoms = read('geom.xyz')
atoms.set_cell([7.66348,7.66348,7.66348*2])
atoms.set_pbc((1,1,1))
atoms.center()

# Set the magnetic moments
p = atoms.get_momenta()
p[0,2]=-1.5
p[1,2]=-1.5
atoms.set_momenta(p)

# Keep some atoms fixed during the simulation
atoms.set_constraint(FixAtoms(indices=range(18,38)))

# Set the calculator and attach it to the system
calc = Siesta('si001_h2_siesta.txt',basis='SZ',xc='PBE',meshcutoff=50*Ry)
atoms.set_calculator(calc)

# Set the VelocityVerlet algorithm and run it
traj = PickleTrajectory('si001_h2_siesta.traj','w',atoms)
dyn = VelocityVerlet(atoms,dt=1.0 * fs,trajectory=traj)
dyn.run(steps=100)

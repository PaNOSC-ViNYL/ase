#!/usr/bin/env python
""" An example for using gromacs calculator in ase.
    Atom positions are relaxed.
    A sample call:
   ./gromacs_example_relax.py hishBOX.gro
"""

from ase import Atoms
from ase.visualize import view
from ase.calculators.gromacs import Gromacs
from ase.optimize import BFGS
from ase.optimize.bfgslinesearch import BFGSLineSearch
from ase.optimize import MDMin
from ase.optimize.sciopt import SciPyFminBFGS, SciPyFminCG 
from ase.optimize.sciopt import SciPyFmin
from ase.optimize.sciopt import SciPyOptimizer

import sys
from ase.io import read, write
from ase import units

infilename = sys.argv[1]
calc = Gromacs(
    init_structure_file=infilename, 
    force_field='oplsaa', 
    water_model='tip3p',
    define = '-DFLEXIBLE',
    integrator = 'md',
    nsteps = '0',
    nstfout = '1',
    nstlog = '1',
    nstenergy = '1',
    energygrps = 'System',
    nstlist = '1',
    ns_type = 'grid',
    pbc = 'xyz',
    rlist = '1.15',
    coulombtype = 'PME-Switch',
    rcoulomb = '0.8',
    vdwtype = 'shift',
    rvdw = '0.8',
    rvdw_switch = '0.75',
    DispCorr = 'Ener')

system = read('gromacs.gro')

system.set_calculator(calc)
e_system = system.get_potential_energy()
f_system = system.get_forces()

# SciPyFminCG == non-linear Polak-Ribiere, should be the same as gromacs 
dyn = SciPyFminCG(system)
dyn.run(fmax=0.01)

from ase.io import read,write
from ase import Atoms
from ase.io import crystal
from ase.optimize import BFGS
from ase.calculators.crystal import CRYSTAL
import time

a0 = 5.43
bulk = Atoms('Si2', [(0, 0, 0),
                     (0.25, 0.25, 0.25)],
             pbc=True)
b = a0 / 2
bulk.set_cell([(0, b, b),
               (b, 0, b),
               (b, b, 0)], scale_atoms=True)

bulk.set_calculator(CRYSTAL(label='Si2',
                         guess = 'True',
                         basis = 'custom',
                         xc = 'PBE',
                         kpts = (2,2,2),
                         otherkey = ['SCFDIR','ANDERSON',['MAXCYCLES','500'],['TOLDEE','6'],['TOLINTEG','7 7 7 7 14'],['FMIXING','90']],
                         ))

dyn = BFGS(bulk, trajectory='init.traj')
dyn.run(fmax=0.01)

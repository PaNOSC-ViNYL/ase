from ase.io import read,write
from ase.io import crystal
from ase.optimize import BFGS
from ase.calculators.crystal import CRYSTAL
import time

geom = read('tio2_slab.xyz')
geom.set_calculator(CRYSTAL(label='tio2_slab',
                         guess = 'True',
                         basis = 'sto-3g',
                         xc = 'PBE',
                         kpts = (2,2,1),
                         otherkey = ['SCFDIR','ANDERSON',['MAXCYCLES','500'],['TOLDEE','6'],['TOLINTEG','7 7 7 7 14'],['FMIXING','90']],
                         ))

dyn = BFGS(geom, trajectory='init.traj')
dyn.run(fmax=0.01)

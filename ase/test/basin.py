import numpy as np
from ase import *
from ase.calculators.lj import LennardJones
from ase.optimize.basin import BasinHopping

N = 13
R = N**(1./3.)
pos = np.random.uniform(-R, R, (N, 3))
s = Atoms('He' + str(N),
          positions = pos)
s.set_calculator(LennardJones())

bh = BasinHopping(s, temperature=100 * kB, dr=0.5)
bh.run(30)

E, smin = bh.get_minimum()
view(smin)

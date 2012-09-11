import numpy as np
from math import pi, sqrt
from ase import Atoms
from ase.calculators.lj import LennardJones
from ase.optimize.basin import BasinHopping
from ase.io import PickleTrajectory, read
from ase.units import kB

# Global minima from
# Wales and Doye, J. Phys. Chem. A, vol 101 (1997) 5111-5116
E_global = {
    4: -6.000000, 
    5: -9.103852, 
    6: -12.712062,
    7: -16.505384, 
}
N = 7
R = N**(1./3.)
pos = np.random.uniform(-R, R, (N, 3))
s = Atoms('He' + str(N),
          positions = pos)
s.set_calculator(LennardJones())

ftraj = 'lowest.traj'
traj = PickleTrajectory(ftraj, 'w', s)
bh = BasinHopping(s, 
                  temperature=100 * kB, dr=0.5, 
                  optimizer_logfile=None)
bh.attach(traj)
bh.run(10)

Emin, smin = bh.get_minimum()
print "N=", N, 'minimal energy found', Emin, 
print ' global minimum:', E_global[N]

# recalc energy
smin.set_calculator(LennardJones())
E = smin.get_potential_energy()
assert abs(E - Emin) < 1e-15
traj.close()
smim = read(ftraj)
E = smin.get_potential_energy()
assert abs(E - Emin) < 1e-15


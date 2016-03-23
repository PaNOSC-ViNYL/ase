from ase.io import read, write
from ase.calculators.lj import LennardJones
from ase.optimize import FIRE
from ase.neb import NEB, NEBtools
from ase import Atoms, Atom

nimages=3
fmax = 0.01
# Define coordinates for initial and final states
initial = Atoms('O4', 
           [( 1.94366484,  2.24788196,  2.32204726),
       ( 3.05353823,  2.08091038,  2.30712548),
       ( 2.63770601,  3.05694348,  2.67368242),
       ( 2.50579418,  2.12540646,  3.28585811)])
           
final = Atoms('O4',
           [( 1.9550137 ,  2.22270649,  2.33191017),
       ( 3.07439495,  2.13662682,  2.31948449),
       ( 2.4473055 ,  1.26930465,  2.65964947),
       ( 2.52788189,  2.1899024 ,  3.29728667)])
           
           
final.set_cell((5,5,5))
initial.set_cell((5,5,5))
final.set_calculator(LennardJones())
initial.set_calculator(LennardJones())

images = [initial]

#Set calculator
for i in range(nimages):
    image = initial.copy()
    image.set_calculator(LennardJones())
    images.append(image)

images.append(final)

# Define the NEB and make a linear interpolation
# with removing translational 
# and rotational degrees of freedom
neb = NEB (images, tr = True)
neb.interpolate(method = 'idpp')

qn = FIRE(neb, logfile = 'LJ4_NEB-TR.log',dt=0.005, maxmove=0.05, dtmax=0.1)
qn.run(steps=20)

# Switch to CI-NEB, still removing the external degrees of freedom
# Also spesify the linearly varying spring constants

neb = NEB (images, climb=True, tr = True)
qn = FIRE(neb, logfile = 'LJ4_CINEB-TR.log',dt=0.005, maxmove=0.05, dtmax=0.1)
qn.run(fmax=fmax)

images = neb.images

nebtools = NEBtools(images)
Ef_neb_rt, dE_neb_rt = nebtools.get_barrier(fit=False)
nsteps_neb_rt = qn.nsteps

##########################################

# normal NEB
# Define coordinates for initial and final states
initial = Atoms('O4', 
           [( 1.94366484,  2.24788196,  2.32204726),
       ( 3.05353823,  2.08091038,  2.30712548),
       ( 2.63770601,  3.05694348,  2.67368242),
       ( 2.50579418,  2.12540646,  3.28585811)])
           
final = Atoms('O4',
           [( 1.9550137 ,  2.22270649,  2.33191017),
       ( 3.07439495,  2.13662682,  2.31948449),
       ( 2.4473055 ,  1.26930465,  2.65964947),
       ( 2.52788189,  2.1899024 ,  3.29728667)])
           
           
final.set_cell((5,5,5))
initial.set_cell((5,5,5))
final.set_calculator(LennardJones())
initial.set_calculator(LennardJones())

images = [initial]

#Set calculator
for i in range(nimages):
    image = initial.copy()
    image.set_calculator(LennardJones())
    images.append(image)

images.append(final)

# Define the NEB and make a linear interpolation
# with removing translational 
# and rotational degrees of freedom
neb = NEB (images, tr = False)
neb.interpolate(method = 'idpp')

qn = FIRE(neb, logfile = 'LJ4_NEB.log',dt=0.005, maxmove=0.05, dtmax=0.1)
qn.run(steps=20)

# Switch to CI-NEB, still removing the external degrees of freedom
# Also spesify the linearly varying spring constants

neb = NEB (images, climb=True)
qn = FIRE(neb, logfile = 'LJ4_CINEB.log',dt=0.005, maxmove=0.05, dtmax=0.1)
qn.run(fmax=fmax)

images = neb.images

nebtools = NEBtools(images)
Ef_neb, dE_neb= nebtools.get_barrier(fit = False)
nsteps_neb = qn.nsteps


assert abs(Ef_neb-Ef_neb_rt) < 1e-2
assert nsteps_neb > nsteps_neb_rt

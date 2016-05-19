from ase.io import Trajectory, read
from ase.neb import NEB, NEBtools
from ase.calculators.morse import MorsePotential
from ase.optimize import BFGS


fmax = 0.05
nimages = 3

print([a.get_potential_energy() for a in Trajectory('H.traj')])
images = [Trajectory('H.traj')[-1]]
for i in range(nimages):
    images.append(images[0].copy())
images[-1].positions[6, 1] = 2 - images[0].positions[6, 1]
neb = NEB(images)
neb.interpolate()
if 0:  # verify that initial images make sense
    from ase.visualize import view
    view(neb.images)

for image in images:
    image.set_calculator(MorsePotential())

dyn = BFGS(neb, trajectory='mep.traj')  # , logfile='mep.log')

dyn.run(fmax=fmax)

for a in neb.images:
    print(a.positions[-1], a.get_potential_energy())

neb.climb = True
dyn.run(fmax=fmax)

# Check NEB tools.
nt_images = read('mep.traj@-4:')
nebtools = NEBtools(nt_images)
nt_fmax = nebtools.get_fmax(climb=True)
Ef, dE = nebtools.get_barrier()
print(Ef, dE, fmax, nt_fmax)
assert nt_fmax < fmax
assert abs(Ef - 1.389) < 0.001
